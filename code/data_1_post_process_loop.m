
clear all;
clc;

T_vec = 100; %[100 250 500 1000];

for T = T_vec
data_main_dir = sprintf('G:/Simulation_Data/max_corr_weak_id/data_ll/data_T%d',T);
data_main_dir = sprintf('./data/data_T%d',T);
%data_main_dir = sprintf('./data_T%d',T);
dir_list = dir(data_main_dir);
N_dirs = length(dir_list);

for d = 1:N_dirs
    if dir_list(d).isdir == 1 
        data_sub_dir = dir_list(d).name;
        if (strcmp(data_sub_dir,'.') == 0 && strcmp(data_sub_dir,'..') == 0)
            data_dir = sprintf('%s/%s', data_main_dir, data_sub_dir);
            file_list = dir(data_dir);
            N_files = length(file_list);
            for file = 1:N_files
                if file_list(file).isdir == 0
                    temp = file_list(file).name;
                    if temp(1:4) == 'comb'  %only use the combined data file%
                        file_name = sprintf('%s/%s',data_dir,temp);
                        disp(temp);
                        load(file_name);
                        
                        data_1_1_post_process;
                        %data_1_2_post_process_w;
                        
                    end
                end
            end
            clearvars -except dir_list N_dirs data_main_dir d;
        end
    end
end

end % T_vec