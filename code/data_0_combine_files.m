% Compile data

clear all;
clc;
T=100;
data_main_dir = sprintf('G:/Simulation_data/max_corr_weak_id/data_ll/data_T%d', T);
data_main_dir = sprintf('./data/data_T%d', T);
%data_main_dir = sprintf('./data_T%d', T);
dir_list = dir(data_main_dir);
N_dirs = length(dir_list);

for d = 1:N_dirs
    if dir_list(d).isdir == 1 
        data_sub_dir = dir_list(d).name;
        if (strcmp(data_sub_dir,'.') == 0 && strcmp(data_sub_dir,'..') == 0)
            data_dir = sprintf('%s/%s', data_main_dir, data_sub_dir);

            file_list = dir(data_dir);

            data_all = [];
            N_files = length(file_list);
            for file = 1:N_files
                if file_list(file).isdir == 0
                    temp = file_list(file).name;
                    %if str2num(temp(12:end-4)) > 10
                    if temp(1:4) ~= 'comb'
                        file_name = sprintf('%s/%s',data_dir,temp);
                        disp(temp);
                        load(file_name);
                        data_all = [data_all data];
                        clear data time
                    end
                    %end
                end
            end

            data = data_all;
            clear data_all;

            outputname=sprintf('%s/combined_%s.mat', data_dir, data_sub_dir);
            save(outputname, 'data');

        end
    end
end


