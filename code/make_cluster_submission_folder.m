clear all;
clc;

cluster_dir = sprintf('./cluster_sub');
if exist(cluster_dir,'dir') == 0
    mkdir(cluster_dir);
end;
 
cluster_job_type = 'batch'; % 'batch' or 'array'
total_num_sims = 500;
array_start = 1;
% array_len = 20;
feasible_flag = 0;
dgp_type_vec = (1:3);
id_type_vec = [0 1 2];
T_vec = [100 250 500 1000]; % [100 250 500 1000];

sim_file_name = 'simulate';
infile = fullfile(sprintf('%s.m', sim_file_name));

for dgp_type = dgp_type_vec
    
if feasible_flag == 0
    if (dgp_type == 1 || dgp_type == 2)
        Js_vec{100} = 500; Js_vec{250} = 250; Js_vec{500} = 100; Js_vec{1000} = 50;
        array_len{100} = 20; array_len{250} = 20; array_len{500} = 50; array_len{1000} = 50;
    elseif dgp_type == 3
        Js_vec{100} = 100; Js_vec{250} = 100; Js_vec{500} = 50; Js_vec{1000} = 50;
        array_len{100} = 50; array_len{250} = 50; array_len{500} = 50; array_len{1000} = 50;        
    end
    if (dgp_type == 1 || dgp_type == 3)
        innovation_type_vec_group{100} = 7; %(1:7);
        innovation_type_vec_group{250} = (1:7);
        innovation_type_vec_group{500} = (1:8);
        innovation_type_vec_group{1000} = (1:9);
    elseif dgp_type == 2
        innovation_type_vec_group{100} = (1:2);
        innovation_type_vec_group{250} = (1:2);
        innovation_type_vec_group{500} = (1:2);
        innovation_type_vec_group{1000} = (1:2);
    end
elseif feasible_flag == 1
    if (dgp_type == 1 || dgp_type == 2)
        Js_vec{100} = 100; Js_vec{250} = 100; Js_vec{500} = 50; Js_vec{1000} = 50;
        array_len{100} = 50; array_len{250} = 50; array_len{500} = 50; array_len{1000} = 50;
    elseif dgp_type == 3
        Js_vec{100} = 100; Js_vec{250} = 100; Js_vec{500} = 50; Js_vec{1000} = 50;
        array_len{100} = 50; array_len{250} = 50; array_len{500} = 50; array_len{1000} = 50;
    end
    if (dgp_type == 1 || dgp_type == 3)
        innovation_type_vec_group{100} = (1:7);
        innovation_type_vec_group{250} = [(1:6), 7];
        innovation_type_vec_group{500} = [(1:6), 8];
        innovation_type_vec_group{1000} = [(1:6), 9];
    elseif dgp_type == 2
        innovation_type_vec_group{100} = (1:2);
        innovation_type_vec_group{250} = (1:2);
        innovation_type_vec_group{500} = (1:2);
        innovation_type_vec_group{1000} = (1:2);
    end
end

% module load matlab/2016b; sbatch --job-name=sim_dgp1_id1_e1_T100 -o out.%j.txt -t 450 -n 12 --wrap="matlab -nodisplay -nojvm -nodesktop -nosplash -singleCompThread -r sim_dgp1_id1_e1_T100"

for indT = 1:length(T_vec)
    T = T_vec(indT);
    innovation_type_vec = innovation_type_vec_group{T};
    if strcmp(cluster_job_type, 'array') == 1
        J = floor(total_num_sims / array_len{T});
        array_end = array_start + array_len{T} - 1;
    elseif strcmp(cluster_job_type, 'batch') == 1
        J = Js_vec{T};
    end    
    
    outputname_ll = sprintf('ll_d%d_T%d', dgp_type, T);
    fname = sprintf('%s/%s.txt', cluster_dir, outputname_ll);
    FID = fopen(fname, 'w');
    % fprintf(FID, 'module load matlab/2016b; ');

        for innovation_type = innovation_type_vec
            for id_type = id_type_vec
                
                %allocate 25 min/sim (for longleaf) * J * T/100
                wall_time = 25*J*T/100*(1-feasible_flag) + 25*J*T/100*(11*11*feasible_flag);
                job_name = sprintf('MCd%di%de%df%dT%d', dgp_type, id_type, innovation_type, feasible_flag, T);
                if strcmp(cluster_job_type, 'array') == 1
                    m_file_name = sprintf('%s(%d,%d,%d,%d,%d,%d,$SLURM_ARRAY_TASK_ID)', sim_file_name, dgp_type, id_type, innovation_type, feasible_flag, J, T);
                    fprintf(FID, 'sbatch --job-name=%s -o out.%%j.txt -t %d --array=[%d-%d%%1] --wrap=\"matlab -nodisplay -nojvm -nodesktop -nosplash -singleCompThread -r ''%s''\"; ', job_name, wall_time, array_start, array_end, m_file_name);
                elseif strcmp(cluster_job_type, 'batch') == 1
                    m_file_name = sprintf('%s(%d,%d,%d,%d,%d,%d)', sim_file_name, dgp_type, id_type, innovation_type, feasible_flag, J, T);
                    fprintf(FID, 'sbatch --job-name=%s -o out.%%j.txt -t %d --wrap=\"matlab -nodisplay -nojvm -nodesktop -nosplash -singleCompThread -r ''%s''\"; ', job_name, wall_time, m_file_name);
                end
                
            end
        end

    fclose(FID);
end % T_vec loop

end % dgp_type loof

outfile = fullfile(cluster_dir, sprintf('%s.m', sim_file_name));
copyfile(infile, outfile);

tempfile = 'class_innovations.m';
infile = fullfile(tempfile);
outfile = fullfile(cluster_dir, tempfile);
copyfile(infile, outfile);

tempfile = 'class_dgp.m';
infile = fullfile(tempfile);
outfile = fullfile(cluster_dir, tempfile);
copyfile(infile, outfile);

tempfile = 'class_tests.m';
infile = fullfile(tempfile);
outfile = fullfile(cluster_dir, tempfile);
copyfile(infile, outfile);




