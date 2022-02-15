%Data Post Processing

if exist('data', 'var') == 1

% Check for repeated samples
J = length(data);
rng_seed = [];
for j = 1:J
    rng_seed = [rng_seed; data(j).rng_seed];
end
[C, ia, ic] = unique(rng_seed);
length(C);
if (length(C) ~= J)
    display('Error: Non-Unique Random Seeds');
end
% dr_w_2 = dr_w(:,ia);
% rf_w_2 = sum(dr_w_2,2)/ length(dr_w_2); % For H_0 true models



%%
J = length(data);
Ln_vec = data(1).Ln_vec;
len_l = length(Ln_vec);
alpha_levels = data(1).alpha_levels;
len_a = length(alpha_levels);
num_drs = size(data(1).dr_table{1},1);
num_dr_ps = size(data(1).dr_p_table{1},1);
num_drs_w_only = size(data(1).dr_table_w_only{1},1);
num_dr_ps_w_only = size(data(1).dr_p_table_w_only{1},1);
clear dr_tables dr_p_tables dr_tables_w_only dr_p_tables_w_only;
for ind_l = 1:len_l
    dr_tables{ind_l} = false(num_drs, len_a, J);
    dr_p_tables{ind_l} = false(num_dr_ps, len_a, J);
    dr_tables_w_only{ind_l} = false(num_drs_w_only, len_a, J);
    dr_p_tables_w_only{ind_l} = false(num_dr_ps_w_only, len_a, J);
    rf_tables{ind_l} = false(num_drs, len_a);
    rf_p_tables{ind_l} = false(num_dr_ps, len_a);
    rf_tables_w_only{ind_l} = false(num_drs_w_only, len_a);
    rf_p_tables_w_only{ind_l} = false(num_dr_ps_w_only, len_a);
    for j = 1:J
        % add in: check if empty or if NA first and verify that it is the
        % correct Ln
        ind_l2 = find(data(j).Ln_vec == Ln_vec(ind_l));
        if size(ind_l2,2) == 0
            dr_tables{ind_l}(:,:,j) = [];
            dr_p_tables{ind_l}(:,:,j) = [];
            dr_tables_w_only{ind_l}(:,:,j) = [];
            dr_p_tables_w_only{ind_l}(:,:,j) = [];
        elseif size(ind_l2,2) ~= 0
            dr_tables{ind_l}(:,:,j) = data(j).dr_table{ind_l2};
            dr_p_tables{ind_l}(:,:,j) = data(j).dr_p_table{ind_l2};
            dr_tables_w_only{ind_l}(:,:,j) = data(j).dr_table_w_only{ind_l2};
            dr_p_tables_w_only{ind_l}(:,:,j) = data(j).dr_p_table_w_only{ind_l2};
        end
    end
    rf_tables{ind_l} = sum(dr_tables{ind_l},3)/size(dr_tables{ind_l},3);
    rf_p_tables{ind_l} = sum(dr_p_tables{ind_l},3)/size(dr_p_tables{ind_l},3);
    rf_tables_w_only{ind_l} = sum(dr_tables_w_only{ind_l},3)/size(dr_tables_w_only{ind_l},3);
    rf_p_tables_w_only{ind_l} = sum(dr_p_tables_w_only{ind_l},3)/size(dr_p_tables_w_only{ind_l},3);
end
% dr_w = dr_w(isnan(Tn_cv(:,2))==0);



%%
outputdir = 'fig/tables';
if exist(outputdir,'dir') == 0
    mkdir(outputdir);
end;

table_name = 'RF';

clear table1 rownames colnames Title Caption;

% obj.dr_table{ind_l}  = [dr_MC_ICS{ind_l}'; dr_LBQ_ICS{ind_l}'; dr_sup_LM_ICS{ind_l}'; dr_CvM_ICS{ind_l}';
%                         dr_MC_LF{ind_l}'; dr_LBQ_LF{ind_l}'; dr_sup_LM_LF{ind_l}'; dr_CvM_LF{ind_l}';
%                         dr_MC_S{ind_l}'; dr_LBQ_S{ind_l}'; dr_sup_LM_S{ind_l}'; dr_CvM_S{ind_l}';
%                         dr_MC_NoX{ind_l}'; dr_LBQ_NoX{ind_l}'; dr_sup_LM_NoX{ind_l}'; dr_CvM_NoX{ind_l}'];
% obj.dr_p_table{ind_l}= [dr_p_MC_ICS{ind_l}'; dr_p_LBQ_ICS{ind_l}'; dr_p_sup_LM_ICS{ind_l}'; dr_p_CvM_ICS{ind_l}';
%                         dr_p_MC_LF{ind_l}'; dr_p_LBQ_LF{ind_l}'; dr_p_sup_LM_LF{ind_l}'; dr_p_CvM_LF{ind_l}';
%                         dr_p_MC_S{ind_l}'; dr_p_LBQ_S{ind_l}'; dr_p_sup_LM_S{ind_l}'; dr_p_CvM_S{ind_l}';
%                         dr_p_MC_NoX{ind_l}'; dr_p_LBQ_NoX{ind_l}'; dr_p_sup_LM_NoX{ind_l}'; dr_p_CvM_NoX{ind_l}'];
% obj.dr_table_w_only{ind_l}   = [dr_MC_W{ind_l}'; dr_LBQ_W{ind_l}'; dr_sup_LM_W{ind_l}'; dr_CvM_W{ind_l}'];
% obj.dr_p_table_w_only{ind_l} = [dr_p_MC_W{ind_l}'; dr_p_LBQ_W{ind_l}'; dr_p_sup_LM_W{ind_l}'; dr_p_CvM_W{ind_l}'];                    

    i=1;
    rownames{i} = sprintf('MC ICS'); i=i+1;
    rownames{i} = sprintf('LBQ ICS'); i=i+1;
    rownames{i} = sprintf('sup LM ICS'); i=i+1;
    rownames{i} = sprintf('CvM ICS'); i=i+1;
    rownames{i} = sprintf('MC LF'); i=i+1;
    rownames{i} = sprintf('LBQ LF'); i=i+1;
    rownames{i} = sprintf('sup LM LF'); i=i+1;
    rownames{i} = sprintf('CvM LF'); i=i+1;
    rownames{i} = sprintf('MC S'); i=i+1;
    rownames{i} = sprintf('LBQ S'); i=i+1;
    rownames{i} = sprintf('sup LM S'); i=i+1;
    rownames{i} = sprintf('CvM S'); i=i+1;
    rownames{i} = sprintf('MC NoX'); i=i+1;
    rownames{i} = sprintf('LBQ NoX'); i=i+1;
    rownames{i} = sprintf('sup LM NoX'); i=i+1;
    rownames{i} = sprintf('CvM NoX');
    i=1;
    rownames_w{i} = sprintf('MC W'); i=i+1;
    rownames_w{i} = sprintf('LBQ W'); i=i+1;
    rownames_w{i} = sprintf('sup LM W'); i=i+1;
    rownames_w{i} = sprintf('CvM W');

    
alpha_levels = data(1).alpha_levels;
for a = 1:length(alpha_levels)
    colnames{a} = sprintf('$(\\alpha=%.2f)$', alpha_levels(a));
end

%dgp_type, innovation_type, id_type, T, Ln, J
id_type = data(1).id_type; id_type_string = data(1).id_type_string;
dgp_type = data(1).dgp_type; dgp_type_string = data(1).dgp_type_string;
innovation_type = data(1).innovation_type; innovation_type_string = data(1).innovation_type_string;
T = data(1).n;
feasible_flag_string = data(1).feasible_flag_string; feasible_flag = data(1).feasible_flag;
beta_n = data(1).beta0;

for ind_l = 1:len_l
    clear table1 table2 Caption Title outputname;
    Ln = Ln_vec(ind_l);

    Caption = '';
    Title{1} = sprintf('Rejection Frequencies: %s CV based Tests', feasible_flag_string);
    Title{2} = sprintf('%s, %s, %s', dgp_type_string, innovation_type_string, id_type_string);
    Title{3} = sprintf('$T=%d$, $\\beta_n = %.3f$, $\\Le_n = %d$, $J=%d$', T, beta_n, Ln, J);

    outputname = sprintf('./%s/%s_cv_dgp%d_e%d_id%d_T%d_f%d_Ln%d', outputdir, table_name, dgp_type, innovation_type, id_type, T, feasible_flag, Ln);
    table1 = rf_tables{ind_l};
    tabletotex(table1, rownames, colnames, outputname, Title);

    Title{1} = sprintf('Rejection Frequencies: %s P-value based Tests', feasible_flag_string);
    outputname = sprintf('./%s/%s_p_dgp%d_e%d_id%d_T%d_f%d_Ln%d', outputdir, table_name, dgp_type, innovation_type, id_type, T, feasible_flag, Ln);
    table2 = rf_p_tables{ind_l};
    tabletotex(table2, rownames, colnames, outputname, Title);
end % Ln_vec loop

for ind_l = 1:len_l
    clear table1 table2 Caption Title outputname;
    Ln = Ln_vec(ind_l);

    Caption = '';
    Title{1} = sprintf('Rejection Frequencies: %s CV based Tests', feasible_flag_string);
    Title{2} = sprintf('%s, %s, %s', dgp_type_string, innovation_type_string, id_type_string);
    Title{3} = sprintf('$T=%d$, $\\beta_n = %.3f$, $\\Le_n = %d$, $J=%d$', T, beta_n, Ln, J);
    
    outputname = sprintf('./%s/%s_cv_dgp%d_e%d_id%d_T%d_f%d_Ln%d_w_only', outputdir, table_name, dgp_type, innovation_type, id_type, T, feasible_flag, Ln);
    table1 = rf_tables_w_only{ind_l};
    tabletotex(table1, rownames_w, colnames, outputname, Title);

    Title{1} = sprintf('Rejection Frequencies: %s P-value based Tests', feasible_flag_string);
    outputname = sprintf('./%s/%s_p_dgp%d_e%d_id%d_T%d_f%d_Ln%d_w_only', outputdir, table_name, dgp_type, innovation_type, id_type, T, feasible_flag, Ln);
    table2 = rf_p_tables_w_only{ind_l};
    tabletotex(table2, rownames_w, colnames, outputname, Title);
end % Ln_vec loop

end % if exist data



 