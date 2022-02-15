% data 2: output latex file with includes for all tables and figures.
clear all;
clc;


main_dir = 'fig';
table_dir = 'tables';
fig_dir = 'figs';
fig_scale = '.5';
output_name_vec{1} = 'RF_p';
output_name_vec{2} = 'RF_cv';

T_vec = [100 250 500 1000];
dgp_type_vec = 3;
id_type_vec = [0 1 2];
%innovation_type_vec = (1:9);
feasible_vec = 0; % [0 1];

for ind_on = 1:length(output_name_vec)
output_name = output_name_vec{ind_on};

%% Tables
data_dir = sprintf('./%s/%s', main_dir, table_dir);
if exist(data_dir,'dir') == 7
file_list = dir(data_dir);
fname = sprintf('%s/tables_%s_dgp%d.tex', main_dir, output_name, dgp_type_vec);
FID = fopen(fname, 'w');
fname_w_only = sprintf('%s/tables_%s_w_only_dgp%d.tex', main_dir, output_name, dgp_type_vec);
FID_w = fopen(fname_w_only, 'w');
N_files = length(file_list)-2;

fprintf(FID, '\\input{./sections/header_article.tex} \n');
fprintf(FID, '\\input{./sections/header_commands.tex} \n');
fprintf(FID, '\\begin{document} \n \n');
fprintf(FID_w, '\\input{./sections/header_article.tex} \n');
fprintf(FID_w, '\\input{./sections/header_commands.tex} \n');
fprintf(FID_w, '\\begin{document} \n \n');

sort_files = 1;

if sort_files == 1
file = 0;
for T = T_vec
    
    if T == 100 || T == 250
        innovation_type_vec = (1:7);
    elseif T == 500
        innovation_type_vec = (1:8);
    elseif T == 1000
        innovation_type_vec = (1:9);
    end
    
    for dgp_type = dgp_type_vec
        for innovation_type = innovation_type_vec
            for id_type = id_type_vec
                for f = feasible_vec

                    if f == 0
                        if T == 100 || T == 250
                            Ln_vec = unique([5, floor(T^(1/3)), floor(sqrt(T)/(log(T)/4)), floor(sqrt(T)/(log(T)/5)), floor(sqrt(T)-1), floor(.5*T/log(T))]); % Number of Lags
                        elseif T == 500 || T == 1000
                            Ln_vec = unique([5, floor(T^(1/3)), floor(sqrt(T)/(log(T)/4)), floor(sqrt(T)/(log(T)/5)), floor(sqrt(T)), floor(.5*T/log(T)), floor(T/log(T))]); % Number of Lags
                        end
                    elseif f == 1
                    end
                   
                  for Ln = Ln_vec
                        temp = sprintf('%s_dgp%d_e%d_id%d_T%d_f%d_Ln%d.tex', output_name, dgp_type, innovation_type, id_type, T, f, Ln);
                        file_in = sprintf('%s/%s/%s', main_dir, table_dir, temp);
                    if exist(file_in) == 2
                        disp(temp);
                        file = file + 1;
                        temp_file_ID = fopen(file_in, 'r+');
                        while 1
                            temp_line = fgets(temp_file_ID);
                            if ~ischar(temp_line) 
                                fprintf(FID, '\n \n');
                                break
                            end
                            fprintf(FID, '%s', temp_line);
                        end
                        fclose(temp_file_ID);
                    %fprintf(FID, '\\input{./%s/%s} \n \n', table_dir, temp);
                    end

                        temp = sprintf('%s_dgp%d_e%d_id%d_T%d_f%d_Ln%d_w_only.tex', output_name, dgp_type, innovation_type, id_type, T, f, Ln);
                        file_in = sprintf('%s/%s/%s', main_dir, table_dir, temp);
                    if exist(file_in) == 2
                        disp(temp);
                        file = file + 1;
                        temp_file_ID = fopen(file_in, 'r+');
                        while 1
                            temp_line = fgets(temp_file_ID);
                            if ~ischar(temp_line) 
                                fprintf(FID_w, '\n \n');
                                break
                            end
                            fprintf(FID_w, '%s', temp_line);
                        end
                        fclose(temp_file_ID);
                    %fprintf(FID_w, '\\input{./%s/%s} \n \n', table_dir, temp);
                    end
                  end
                end
            end
        end
    end
end
% if (file ~= N_files)
%     disp('Error: files processed != N_files');
% end

elseif sort_files == 0
  for file = 1:N_files
     if file_list(file).isdir == 0
         temp = file_list(file).name;
         if strcmp(temp(end-3:end),'.tex') == 1 && strcmp(temp,sprintf('tables_%s.tex', output_name)) == 0 && strcmp(temp,sprintf('figs_%s.tex', output_name)) == 0
             disp(temp);
             if strcmp(temp(end-9:end-4),'w_only') == 1
                 fprintf(FID_w, '\\input{./%s/%s} \n \n', table_dir, temp);
             else
                 fprintf(FID, '\\input{./%s/%s} \n \n', table_dir, temp);
             end
         end
     end
  end
end


fprintf(FID, '\\end{document}');
fprintf(FID_w, '\\end{document}');

fclose(FID);
fclose(FID_w);
end



%% Figures
data_dir = sprintf('./%s/%s', main_dir, fig_dir);
if exist(data_dir,'dir') == 7
file_list = dir(data_dir);
fname = sprintf('%s/figs_%s.tex', main_dir, output_name);
FID = fopen(fname, 'w');
N_files = length(file_list);
for file = 1:N_files
    if file_list(file).isdir == 0
        temp = file_list(file).name;
        if strcmp(temp(end-3:end),'.png') == 1
            fprintf(FID, '\centering \n \includegraphics[scale=%s]{./%s/%s} \n \n', fig_scale, fig_dir, temp);
        end
    end
end
fclose(FID);
end

end % output_name loop
