classdef class_tests < class_dgp
    properties (SetAccess=public, GetAccess=public)

        % Max Corr
        MC_test_stat
        cv_MC_ICS, cv_MC_LF, cv_MC_S, cv_MC_W, cv_MC_NoX
        p_MC_ICS, p_MC_LF, p_MC_S, p_MC_W, p_MC_NoX
        % dr_MC_ICS, dr_MC_LF, dr_MC_S, dr_MC_W, dr_MC_NoX
        
        % Ljung-Box
        LBQ_test_stat
        cv_LBQ_ICS, cv_LBQ_LF, cv_LBQ_S, cv_LBQ_W, cv_LBQ_NoX
        p_LBQ_ICS, p_LBQ_LF, p_LBQ_S, p_LBQ_W, p_LBQ_NoX
        % dr_LBQ_ICS, dr_LBQ_LF, dr_LBQ_S, dr_LBQ_W, dr_LBQ_NoX
        
        % AP sup LM
        lambda_grid
        sup_LM_test_stat
        cv_sup_LM_ICS, cv_sup_LM_LF, cv_sup_LM_S, cv_sup_LM_W, cv_sup_LM_NoX
        p_sup_LM_ICS, p_sup_LM_LF, p_sup_LM_S, p_sup_LM_W, p_sup_LM_NoX
        % dr_sup_LM_ICS, dr_sup_LM_LF, dr_sup_LM_S, dr_sup_LM_W, dr_sup_LM_NoX
        
        % CvM test
        pi_grid
        CvM_test_stat
        cv_CvM_ICS, cv_CvM_LF, cv_CvM_S, cv_CvM_W, cv_CvM_NoX
        p_CvM_ICS, p_CvM_LF, p_CvM_S, p_CvM_W, p_CvM_NoX
        % dr_CvM_ICS, dr_CvM_LF, dr_CvM_S, dr_CvM_W, dr_CvM_NoX
        
        % DWB
        Ln_vec, kn, M, n
        
        % other
        kappa_n
        b_in_vec
        pi0_in_vec
        dr_table, dr_p_table, dr_table_w_only, dr_p_table_w_only
        p_val_table
        alpha_levels
        feasible_flag, feasible_flag_string
    end
    properties (SetAccess=public, GetAccess=public, Hidden = true)
        options
    end
    methods
        function out = about(obj)
            disp('V28');
            disp('------------');
            disp('Added more lag specifications and changed parameters for feasible=1.');
            disp('------------');
            disp('------------');
            disp('V27');
            disp('------------');
            disp('Corrected bias in variance issue, and now allows for ARMA(1,1) dgp');
            disp('V24');
            disp('------------');
            disp('Allow Ln vector within the testing function to decrease computation time');
            disp('------------');
            disp('------------');
            disp('V23');
            disp('------------');
            disp('Changed to implement DWB with expansion for the other tests.');
            disp('Reorganized code');
            disp('------------');
            disp('------------');
            disp('V18');
            disp('------------');
            disp('Added in other tests');
            disp('Separated DGP class from Test class');
            disp('Currently only STAR1 DGP is allowed, need to add others');
            disp('Need to combine other tests and Max Corr Test');
            disp('Maybe Estimate bias2 instead of using asymptotic value');
            disp('------------');
            disp('------------');
            disp('Need to add other tests and simulate the asymptotic distribution of the test stat to compare to the bootstrap');
            disp('------------');
            disp('------------');
            disp('Version 16');
            disp('------------');
            disp('New to this version: ');
            disp('various updates, including adding bias terms in weak id bootstrap and adding multiple innovation types');
            disp('------------');
            disp('------------');
            disp('Version 14');
            disp('------------');
            disp('New to this version: ');
            disp('- Took beta out of d_theta to account for B^{-1}(beta) in semi-strong id');
            disp('- Changing BS computation to follow the max corr paper.');
            disp('- Hopefully this will fix the BS implementation issue.');
            disp('------------');
            disp('Version 12');
            disp('------------');
            disp('New to this version: ');
            disp('- Switched mu and pi.  Now pi is the threshold param, mu is the speed param.');
            disp('- Added in LF cv, pval, and decision_rule');
            disp('------------');
            disp('------------');
            disp('Potential Issues: implementation of the bootstrap');
            disp('------------');
        end
        function obj = class_tests(dgp_type, id_type, innovation_type, theta_in, num_params, T, seed)
            %First instantiate the DGP
            obj@class_dgp(dgp_type, id_type, innovation_type, theta_in, num_params, T, seed);
        end
        function obj = clean_up(obj)
            obj.Y = [];
            obj.Yhat = [];
            obj.ehat = [];
            obj.U = [];
            obj.pi_grid = [];
            obj.lambda_grid = [];
            obj.bias2_true = []; 
            obj.bias2_est = []; 
            obj.bias2_est_old = [];
            obj.options = [];
        end
        function obj = set_params(obj, feasible_flag)
            obj.options = optimoptions(@fmincon,'Display','off','Diagnostics','off');
            obj.n = obj.Te;  % We are testing residuals
            obj.alpha_levels = [0.01 0.05 0.10];
            obj.lambda_grid = (-.8:.005:.8); % obj.lambda_grid = (-.8:.005:.8);
            obj.pi_grid = (0:0.001:3.141); % obj.pi_grid = (0:0.001:3.141); % for CvM test, not max corr test
            obj.feasible_flag = feasible_flag;
            if obj.feasible_flag == 0
                if obj.n == 100 || obj.n == 250
                    obj.M = 1000;  %number of DWB repetitions
                    obj.Ln_vec = unique([5, floor(obj.n^(1/3)), floor(sqrt(obj.n)/(log(obj.n)/4)), floor(sqrt(obj.n)/(log(obj.n)/5)), floor(sqrt(obj.n)-1), floor(.5*obj.n/log(obj.n)), floor(obj.n/log(obj.n))]); % Number of Lags
                elseif obj.n == 500 || obj.n == 1000
                    obj.M = 500;  %number of DWB repetitions
                    obj.Ln_vec = unique([5, floor(obj.n^(1/3)), floor(sqrt(obj.n)/(log(obj.n)/4)), floor(sqrt(obj.n)/(log(obj.n)/5)), floor(sqrt(obj.n)), floor(.5*obj.n/log(obj.n)), floor(obj.n/log(obj.n))]); % Number of Lags
                end
            elseif obj.feasible_flag == 1  % Extremely computationally intensive!!!
                obj.M = 500;  %number of DWB repetitions
                if obj.n == 100
                    if (obj.innovation_type == 1 || obj.innovation_type == 2 || obj.innovation_type == 3)
                        obj.Ln_vec = unique([5, floor(sqrt(obj.n)/(log(obj.n)/5)), floor(.5*obj.n/log(obj.n))]); % Number of Lags
                    elseif (obj.innovation_type == 4 || obj.innovation_type == 5)
                        obj.Ln_vec = unique([5, floor(sqrt(obj.n)/(log(obj.n)/5)), floor(.5*obj.n/log(obj.n))]); % Number of Lags
                    elseif (obj.innovation_type == 6)
                        obj.Ln_vec = unique([floor(.5*obj.n/log(obj.n))]); % Number of Lags
                    elseif (obj.innovation_type == 7)
                        obj.Ln_vec = unique([floor(obj.n/log(obj.n))]); % Number of Lags
                    end
                elseif obj.n == 250
                    if (obj.innovation_type == 1 || obj.innovation_type == 2 || obj.innovation_type == 3)
                        obj.Ln_vec = unique([5, floor(sqrt(obj.n)/(log(obj.n)/4)), floor(sqrt(obj.n)/(log(obj.n)/5)), floor(sqrt(obj.n)-1), floor(.5*obj.n/log(obj.n))]); % Number of Lags
                    elseif (obj.innovation_type == 4 || obj.innovation_type == 5)
                        obj.Ln_vec = unique([5, floor(sqrt(obj.n)/(log(obj.n)/5)), floor(sqrt(obj.n)-1), floor(.5*obj.n/log(obj.n))]); % Number of Lags
                    elseif (obj.innovation_type == 6)
                        obj.Ln_vec = unique([floor(sqrt(obj.n)/(log(obj.n)/5)), floor(sqrt(obj.n)-1), floor(.5*obj.n/log(obj.n))]); % Number of Lags
                    elseif (obj.innovation_type == 7)
                        obj.Ln_vec = unique([floor(.5*obj.n/log(obj.n))]); % Number of Lags
                    end
                elseif obj.n == 500
                    if (obj.innovation_type == 1 || obj.innovation_type == 2 || obj.innovation_type == 3)
                        obj.Ln_vec = unique([5, floor(sqrt(obj.n)/(log(obj.n)/5)), floor(.5*obj.n/log(obj.n)), floor(obj.n/log(obj.n))]); % Number of Lags
                    elseif (obj.innovation_type == 4 || obj.innovation_type == 5)
                        obj.Ln_vec = unique([5, floor(sqrt(obj.n)/(log(obj.n)/5)), floor(.5*obj.n/log(obj.n)), floor(obj.n/log(obj.n))]); % Number of Lags
                    elseif (obj.innovation_type == 6)
                        obj.Ln_vec = unique([floor(sqrt(obj.n)/(log(obj.n)/5)), floor(.5*obj.n/log(obj.n)), floor(obj.n/log(obj.n))]); % Number of Lags
                    elseif (obj.innovation_type == 7)
                        obj.Ln_vec = unique([floor(sqrt(obj.n))]); % Number of Lags
                    elseif (obj.innovation_type == 8)
                        obj.Ln_vec = unique([floor(obj.n/log(obj.n))]); % Number of Lags
                    end
                elseif obj.n == 1000
                    if (obj.innovation_type == 1 || obj.innovation_type == 2 || obj.innovation_type == 3)
                        obj.Ln_vec = unique([5, floor(obj.n^(1/3)), floor(sqrt(obj.n)/(log(obj.n)/5)), floor(.5*obj.n/log(obj.n)), floor(obj.n/log(obj.n))]); % Number of Lags
                    elseif (obj.innovation_type == 4 || obj.innovation_type == 5)
                        obj.Ln_vec = unique([5, floor(obj.n^(1/3)), floor(sqrt(obj.n)/(log(obj.n)/5)), floor(.5*obj.n/log(obj.n)), floor(obj.n/log(obj.n))]); % Number of Lags
                    elseif (obj.innovation_type == 6)
                        obj.Ln_vec = unique([floor(sqrt(obj.n)/(log(obj.n)/5)), floor(.5*obj.n/log(obj.n)), floor(obj.n/log(obj.n))]); % Number of Lags
                    elseif (obj.innovation_type == 7)
                        obj.Ln_vec = unique([floor(sqrt(obj.n)/(log(obj.n)/5))]); % Number of Lags
                    elseif (obj.innovation_type == 8)
                        obj.Ln_vec = unique([floor(.5*obj.n/log(obj.n))]); % Number of Lags
                    elseif (obj.innovation_type == 9)
                        obj.Ln_vec = unique([floor(obj.n/log(obj.n))]); % Number of Lags
                    end
                end
            end
            obj.kn = floor(sqrt(obj.n)); % Window size for DWB
            obj.kappa_n = sqrt(log(obj.T)); % for ICS cvs
            obj.b_true = 0 * (obj.id_type == 0) + obj.beta0 * sqrt(obj.n) * (obj.id_type == 1) + obj.beta0 * (obj.id_type == 2); 
            if feasible_flag == 0
                obj.feasible_flag_string = 'Infeasible';
                obj.b_in_vec = obj.b_true;
                obj.pi0_in_vec = obj.pi0;
            elseif feasible_flag == 1  % Extremely computationally intensive!!!
                obj.feasible_flag_string = 'Feasible';
                step_size = .1;
                num_steps = 5;
                obj.b_in_vec = (obj.b_true-(num_steps*step_size):step_size:obj.b_true+(num_steps*step_size)); 
                obj.pi0_in_vec = (obj.pi0-(num_steps*step_size):step_size:obj.pi0+(num_steps*step_size));
            end
        end
        function obj = fcn_test_stats_only(obj, feasible_flag)
            % Step 0: Setup the testing parameters
            obj = obj.set_params(feasible_flag);
            len_Ln = length(obj.Ln_vec);
            
            % Step 1.1: get correlation vector:
            [rho_vec, gamma_vec] = obj.fcn_rho_Ln_no_expansion(obj.n-1);
            
            % Step 1.2: Get test statistics
            for ind_l = 1:len_Ln
                obj.MC_test_stat{ind_l} = obj.fcn_MC_Tn_hat(rho_vec(1:obj.Ln_vec(ind_l)));
                obj.LBQ_test_stat{ind_l} = obj.LBQ_fcn(rho_vec(1:obj.Ln_vec(ind_l)));
                obj.sup_LM_test_stat{ind_l} = obj.sup_LM_fcn(rho_vec);
                obj.CvM_test_stat{ind_l} = obj.CvM_fcn(gamma_vec);
            end
            clear rho_vec gamma_vec;
        end
        function obj = fcn_run_all_tests(obj, feasible_flag)
            % Step 0: Setup the testing parameters
            obj = obj.set_params(feasible_flag);
            len_Ln = length(obj.Ln_vec);
            
            % Step 1.1: get correlation vector:
            [rho_vec, gamma_vec] = obj.fcn_rho_Ln_no_expansion(obj.n-1);
            
            % Step 1.2: Get test statistics
            for ind_l = 1:len_Ln
                obj.MC_test_stat{ind_l} = obj.fcn_MC_Tn_hat(rho_vec(1:obj.Ln_vec(ind_l)));
                obj.LBQ_test_stat{ind_l} = obj.LBQ_fcn(rho_vec(1:obj.Ln_vec(ind_l)));
                obj.sup_LM_test_stat{ind_l} = obj.sup_LM_fcn(rho_vec);
                obj.CvM_test_stat{ind_l} = obj.CvM_fcn(gamma_vec);
            end
            clear rho_vec gamma_vec;
            
            % Step 2: Grab DWB Multipliers
            Z = obj.DWB_wM_fcn;
            
            % Step 3.1: Calculate correlations without expansion
            % and Step 3.2: Calculate Test Stat Distributions without expansion
            MC_DWB_NoX = zeros(obj.M,len_Ln);
            LBQ_DWB_NoX = zeros(obj.M,len_Ln);
            sup_LM_DWB_NoX = zeros(obj.M,len_Ln);
            CvM_DWB_NoX = zeros(obj.M,len_Ln);
            for m = 1:obj.M
                [rho_vec, gamma_vec] = obj.fcn_rho_Ln_no_expansion(obj.n-1, Z(:,m));
                temp_sup_LM = obj.sup_LM_fcn(rho_vec);
                temp_CvM = obj.CvM_fcn(gamma_vec);
                for ind_l = 1:len_Ln
                    MC_DWB_NoX(m,ind_l) = obj.fcn_MC_Tn_hat(rho_vec(1:obj.Ln_vec(ind_l)));
                    LBQ_DWB_NoX(m,ind_l) = obj.LBQ_fcn(rho_vec(1:obj.Ln_vec(ind_l)));
                    sup_LM_DWB_NoX(m,ind_l) = temp_sup_LM;
                    CvM_DWB_NoX(m,ind_l) = temp_CvM;
                end
            end
            clear rho_vec gamma_vec temp_sup_LM temp_CvM;
            
            % Step 3.3: Calculate no expansion cvs and p-values
            for ind_l = 1:len_Ln
                obj.cv_MC_NoX{ind_l} = obj.fcn_cv_bs(MC_DWB_NoX(:,ind_l));
                obj.cv_LBQ_NoX{ind_l} = obj.fcn_cv_bs(LBQ_DWB_NoX(:,ind_l));
                obj.cv_sup_LM_NoX{ind_l} = obj.fcn_cv_bs(sup_LM_DWB_NoX(:,ind_l));
                obj.cv_CvM_NoX{ind_l} = obj.fcn_cv_bs(CvM_DWB_NoX(:,ind_l));
                obj.p_MC_NoX{ind_l} = obj.pval_fcn(MC_DWB_NoX(:,ind_l), obj.M, obj.MC_test_stat{ind_l});
                obj.p_LBQ_NoX{ind_l} = obj.pval_fcn(LBQ_DWB_NoX(:,ind_l), obj.M, obj.LBQ_test_stat{ind_l});
                obj.p_sup_LM_NoX{ind_l} = obj.pval_fcn(sup_LM_DWB_NoX(:,ind_l), obj.M, obj.sup_LM_test_stat{ind_l});
                obj.p_CvM_NoX{ind_l} = obj.pval_fcn(CvM_DWB_NoX(:,ind_l), obj.M, obj.CvM_test_stat{ind_l});
            end
            clear MC_DWB_NoX LBQ_DWB_NoX sup_LM_DWB_NoX CvM_DWB_NoX; % Free up some space
            
            % Step 4.1: Calculate correlations under strong id
            % and Step 4.2: Use strong id correlations to calculate Strong Id Test Stat Distributions
            MC_DWB_S = zeros(obj.M,len_Ln);
            LBQ_DWB_S = zeros(obj.M,len_Ln);
            sup_LM_DWB_S = zeros(obj.M,len_Ln);
            CvM_DWB_S = zeros(obj.M,len_Ln);
            for m = 1:obj.M
                [rho_vec, gamma_vec] = obj.fcn_rho_theta_Ln(obj.n-1, Z(:,m));
                temp_sup_LM = obj.sup_LM_fcn(rho_vec);
                temp_CvM = obj.CvM_fcn(gamma_vec);
                for ind_l = 1:len_Ln
                    MC_DWB_S(m,ind_l) = obj.fcn_MC_Tn_hat(rho_vec(1:obj.Ln_vec(ind_l)));
                    LBQ_DWB_S(m,ind_l) = obj.LBQ_fcn(rho_vec(1:obj.Ln_vec(ind_l)));
                    sup_LM_DWB_S(m,ind_l) = temp_sup_LM;
                    CvM_DWB_S(m,ind_l) = temp_CvM;
                end
            end
            clear rho_vec gamma_vec temp_sup_LM temp_CvM;
            
            % Step 4.3: Calculate strong id cv's and p-values
            for ind_l = 1:len_Ln
                obj.cv_MC_S{ind_l} = obj.fcn_cv_bs(MC_DWB_S(:,ind_l));
                obj.cv_LBQ_S{ind_l} = obj.fcn_cv_bs(LBQ_DWB_S(:,ind_l));
                obj.cv_sup_LM_S{ind_l} = obj.fcn_cv_bs(sup_LM_DWB_S(:,ind_l));
                obj.cv_CvM_S{ind_l} = obj.fcn_cv_bs(CvM_DWB_S(:,ind_l));
                obj.p_MC_S{ind_l} = obj.pval_fcn(MC_DWB_S(:,ind_l), obj.M, obj.MC_test_stat{ind_l});
                obj.p_LBQ_S{ind_l} = obj.pval_fcn(LBQ_DWB_S(:,ind_l), obj.M, obj.LBQ_test_stat{ind_l});
                obj.p_sup_LM_S{ind_l} = obj.pval_fcn(sup_LM_DWB_S(:,ind_l), obj.M, obj.sup_LM_test_stat{ind_l});
                obj.p_CvM_S{ind_l} = obj.pval_fcn(CvM_DWB_S(:,ind_l), obj.M, obj.CvM_test_stat{ind_l});
            end
            clear MC_DWB_S LBQ_DWB_S sup_LM_DWB_S CvM_DWB_S; % Free up some space
            
            % Step 5.0: For feasible tests, first determine if
            % bootstrapping the weak id distribution is necessary by
            % checking the ICS statistic.  Note that for infeasible tests,
            % we still want to calculate LF cvs, so we always bootstrap the
            % distribution under weak id.
            ICS_flag = 0;
            if feasible_flag == 1
                ICS_flag = obj.fcn_ICS; % if ICS_flag == 1, then use strong id, o.w. use LF
            end
            if ICS_flag == 0 || feasible_flag == 0 % so skip bs wid if ICS_flag == 1 AND f_flag == 1
            
            % Step 5.1: Calculate correlations under weak id
            [rho_vec_mat, gamma_vec_mat, obj.bias2_est] = obj.fcn_rho_psi_Ln_loop(obj.n-1, Z);  % NOTE: For simplicity, this defines all of the correlation vectors: M x b_vec x pi0_vec

            %Calculate true bias2 for comparison
            bias2_flag = 1;
            if bias2_flag == 1
                obj.bias2_true = zeros((obj.n-1),1);
                for h = 1:(obj.n-1)
                    obj.bias2_true(h) = obj.fcn_bias2_true(h);
                end
            end
            
            % Step 5.2: Use weak id correlations to calculate Weak Id Test Stat Distributions(b,pi)
            len_b = length(obj.b_in_vec);
            len_p = length(obj.pi0_in_vec);
            MC_DWB_W = zeros(obj.M, len_b, len_p, len_Ln);
            LBQ_DWB_W = zeros(obj.M, len_b, len_p, len_Ln);
            sup_LM_DWB_W = zeros(obj.M, len_b, len_p, len_Ln);
            CvM_DWB_W = zeros(obj.M, len_b, len_p, len_Ln);
            for m = 1:obj.M
                for ind_b = 1:len_b
                    for ind_p = 1:len_p
                        rho_vec = rho_vec_mat(:, m, ind_b, ind_p);
                        gamma_vec = gamma_vec_mat(:, m, ind_b, ind_p);
                        temp_sup_LM = obj.sup_LM_fcn(rho_vec);
                        temp_CvM = obj.CvM_fcn(gamma_vec);
                        for ind_l = 1:len_Ln
                            MC_DWB_W(m, ind_b, ind_p, ind_l) = obj.fcn_MC_Tn_hat(rho_vec(1:obj.Ln_vec(ind_l)));
                            LBQ_DWB_W(m, ind_b, ind_p, ind_l) = obj.LBQ_fcn(rho_vec(1:obj.Ln_vec(ind_l)));
                            sup_LM_DWB_W(m, ind_b, ind_p, ind_l) = temp_sup_LM;
                            CvM_DWB_W(m, ind_b, ind_p, ind_l) = temp_CvM;
                        end
                    end
                end
            end
            clear rho_vec_mat gamma_vec_mat rho_vec gamma_vec temp_sup_LM temp_CvM;

            % Step 5.3: Calculate weak id cv's and p-values for each distribution(b,pi)
            len_a = length(obj.alpha_levels);
            len_b = length(obj.b_in_vec);
            len_p = length(obj.pi0_in_vec);
            cv_MC_W_mat = zeros(len_a, len_b, len_p, len_Ln);
            cv_LBQ_W_mat = zeros(len_a, len_b, len_p, len_Ln);
            cv_sup_LM_W_mat = zeros(len_a, len_b, len_p, len_Ln);
            cv_CvM_W_mat = zeros(len_a, len_b, len_p, len_Ln);
            p_MC_W_mat = zeros(1, len_b, len_p, len_Ln);
            p_LBQ_W_mat = zeros(1, len_b, len_p, len_Ln);
            p_sup_LM_W_mat = zeros(1, len_b, len_p, len_Ln);
            p_CvM_W_mat = zeros(1, len_b, len_p, len_Ln);
            for ind_b = 1:len_b
                for ind_p = 1:len_p
                    for ind_l = 1:len_Ln
                        cv_MC_W_mat(:, ind_b, ind_p, ind_l)     = obj.fcn_cv_bs(MC_DWB_W(:, ind_b, ind_p, ind_l));
                        cv_LBQ_W_mat(:, ind_b, ind_p, ind_l)    = obj.fcn_cv_bs(LBQ_DWB_W(:, ind_b, ind_p, ind_l));
                        cv_sup_LM_W_mat(:, ind_b, ind_p, ind_l) = obj.fcn_cv_bs(sup_LM_DWB_W(:, ind_b, ind_p, ind_l));
                        cv_CvM_W_mat(:, ind_b, ind_p, ind_l)    = obj.fcn_cv_bs(CvM_DWB_W(:, ind_b, ind_p, ind_l));
                        p_MC_W_mat(:, ind_b, ind_p, ind_l)     = obj.pval_fcn(MC_DWB_W(:, ind_b, ind_p, ind_l), obj.M, obj.MC_test_stat{ind_l});
                        p_LBQ_W_mat(:, ind_b, ind_p, ind_l)    = obj.pval_fcn(LBQ_DWB_W(:, ind_b, ind_p, ind_l), obj.M, obj.LBQ_test_stat{ind_l});
                        p_sup_LM_W_mat(:, ind_b, ind_p, ind_l) = obj.pval_fcn(sup_LM_DWB_W(:, ind_b, ind_p, ind_l), obj.M, obj.sup_LM_test_stat{ind_l});
                        p_CvM_W_mat(:, ind_b, ind_p, ind_l)    = obj.pval_fcn(CvM_DWB_W(:, ind_b, ind_p, ind_l), obj.M, obj.CvM_test_stat{ind_l});
                    end
                end
            end
            clear MC_DWB_W LBQ_DWB_W sup_LM_DWB_W CvM_DWB_W; % Free up some space
            
            % Step 5.4: sup over b,pi to get weak id cv's and p-values
            for ind_l = 1:len_Ln
                obj.cv_MC_W{ind_l}    = max(max(cv_MC_W_mat(:, :, :, ind_l),[],3),[],2);
                obj.cv_LBQ_W{ind_l}   = max(max(cv_LBQ_W_mat(:, :, :, ind_l),[],3),[],2);
                obj.cv_sup_LM_W{ind_l}= max(max(cv_sup_LM_W_mat(:, :, :, ind_l),[],3),[],2);
                obj.cv_CvM_W{ind_l}   = max(max(cv_CvM_W_mat(:, :, :, ind_l),[],3),[],2);
                obj.p_MC_W{ind_l}     = min(min(p_MC_W_mat(:, :, :, ind_l),[],3),[],2);
                obj.p_LBQ_W{ind_l}    = min(min(p_LBQ_W_mat(:, :, :, ind_l),[],3),[],2);
                obj.p_sup_LM_W{ind_l} = min(min(p_sup_LM_W_mat(:, :, :, ind_l),[],3),[],2);
                obj.p_CvM_W{ind_l}    = min(min(p_CvM_W_mat(:, :, :, ind_l),[],3),[],2);
            end
            clear cv_MC_W_mat cv_LBQ_W_mat cv_sup_LM_W_mat cv_CvM_W_mat p_MC_W_mat p_LBQ_W_mat p_sup_LM_W_mat p_CvM_W_mat; % Free up some space

            else %what to do if skipping weak id bs
                for ind_l = 1:len_Ln
                  obj.cv_MC_W{ind_l}    = 0;
                  obj.cv_LBQ_W{ind_l}   = 0;
                  obj.cv_sup_LM_W{ind_l}= 0;
                  obj.cv_CvM_W{ind_l}   = 0;
                  obj.p_MC_W{ind_l}     = 1;
                  obj.p_LBQ_W{ind_l}    = 1;
                  obj.p_sup_LM_W{ind_l} = 1;
                  obj.p_CvM_W{ind_l}    = 1;
                end           
            end % End the if from Step 5.0
            
            % Step 6: Calculate LF and ICS Critical Values and P-values
            for ind_l = 1:len_Ln
                obj.cv_MC_LF{ind_l}     = obj.fcn_cv_LF(obj.cv_MC_W{ind_l}, obj.cv_MC_S{ind_l});
                obj.cv_MC_ICS{ind_l}    = obj.fcn_cv_ICS(obj.cv_MC_LF{ind_l}, obj.cv_MC_S{ind_l});
                obj.cv_LBQ_LF{ind_l}    = obj.fcn_cv_LF(obj.cv_LBQ_W{ind_l}, obj.cv_LBQ_S{ind_l});
                obj.cv_LBQ_ICS{ind_l}   = obj.fcn_cv_ICS(obj.cv_LBQ_LF{ind_l}, obj.cv_LBQ_S{ind_l});
                obj.cv_sup_LM_LF{ind_l} = obj.fcn_cv_LF(obj.cv_sup_LM_W{ind_l}, obj.cv_sup_LM_S{ind_l});
                obj.cv_sup_LM_ICS{ind_l}= obj.fcn_cv_ICS(obj.cv_sup_LM_LF{ind_l}, obj.cv_sup_LM_S{ind_l});
                obj.cv_CvM_LF{ind_l}    = obj.fcn_cv_LF(obj.cv_CvM_W{ind_l}, obj.cv_CvM_S{ind_l});
                obj.cv_CvM_ICS{ind_l}   = obj.fcn_cv_ICS(obj.cv_CvM_LF{ind_l}, obj.cv_CvM_S{ind_l});

                obj.p_MC_LF{ind_l}      = obj.fcn_p_LF(obj.p_MC_W{ind_l}, obj.p_MC_S{ind_l});
                obj.p_MC_ICS{ind_l}     = obj.fcn_p_ICS(obj.p_MC_LF{ind_l}, obj.p_MC_S{ind_l});
                obj.p_LBQ_LF{ind_l}     = obj.fcn_p_LF(obj.p_LBQ_W{ind_l}, obj.p_LBQ_S{ind_l});
                obj.p_LBQ_ICS{ind_l}    = obj.fcn_p_ICS(obj.p_LBQ_LF{ind_l}, obj.p_LBQ_S{ind_l});
                obj.p_sup_LM_LF{ind_l}  = obj.fcn_p_LF(obj.p_sup_LM_W{ind_l}, obj.p_sup_LM_S{ind_l});
                obj.p_sup_LM_ICS{ind_l} = obj.fcn_p_ICS(obj.p_sup_LM_LF{ind_l}, obj.p_sup_LM_S{ind_l});
                obj.p_CvM_LF{ind_l}     = obj.fcn_p_LF(obj.p_CvM_W{ind_l}, obj.p_CvM_S{ind_l});
                obj.p_CvM_ICS{ind_l}    = obj.fcn_p_ICS(obj.p_CvM_LF{ind_l}, obj.p_CvM_S{ind_l});
            end

            % Step 7.1: Decision Rules (Using both cv's and p-values, even though theory is only for cv's right now)
            dr_MC_NoX{len_Ln} = []; dr_MC_W{len_Ln} = []; dr_MC_S{len_Ln} = []; dr_MC_LF{len_Ln} = []; dr_MC_ICS{len_Ln} = []; dr_p_MC_NoX{len_Ln} = []; dr_p_MC_W{len_Ln} = []; dr_p_MC_S{len_Ln} = []; dr_p_MC_LF{len_Ln} = []; dr_p_MC_ICS{len_Ln} = [];
            dr_LBQ_NoX{len_Ln} = []; dr_LBQ_W{len_Ln} = []; dr_LBQ_S{len_Ln} = []; dr_LBQ_LF{len_Ln} = []; dr_LBQ_ICS{len_Ln} = []; dr_p_LBQ_NoX{len_Ln} = []; dr_p_LBQ_W{len_Ln} = []; dr_p_LBQ_S{len_Ln} = []; dr_p_LBQ_LF{len_Ln} = []; dr_p_LBQ_ICS{len_Ln} = [];
            dr_sup_LM_NoX{len_Ln} = []; dr_sup_LM_W{len_Ln} = []; dr_sup_LM_S{len_Ln} = []; dr_sup_LM_LF{len_Ln} = []; dr_sup_LM_ICS{len_Ln} = []; dr_p_sup_LM_NoX{len_Ln} = []; dr_p_sup_LM_W{len_Ln} = []; dr_p_sup_LM_S{len_Ln} = []; dr_p_sup_LM_LF{len_Ln} = []; dr_p_sup_LM_ICS{len_Ln} = [];
            dr_CvM_NoX{len_Ln} = []; dr_CvM_W{len_Ln} = []; dr_CvM_S{len_Ln} = []; dr_CvM_LF{len_Ln} = []; dr_CvM_ICS{len_Ln} = []; dr_p_CvM_NoX{len_Ln} = []; dr_p_CvM_W{len_Ln} = []; dr_p_CvM_S{len_Ln} = []; dr_p_CvM_LF{len_Ln} = []; dr_p_CvM_ICS{len_Ln} = [];
            for ind_l = 1:len_Ln
                dr_MC_NoX{ind_l} = obj.fcn_dec_rule(obj.MC_test_stat{ind_l}, obj.cv_MC_NoX{ind_l});
                dr_MC_W{ind_l}   = obj.fcn_dec_rule(obj.MC_test_stat{ind_l}, obj.cv_MC_W{ind_l});
                dr_MC_S{ind_l}   = obj.fcn_dec_rule(obj.MC_test_stat{ind_l}, obj.cv_MC_S{ind_l});
                dr_MC_LF{ind_l}  = obj.fcn_dec_rule(obj.MC_test_stat{ind_l}, obj.cv_MC_LF{ind_l});
                dr_MC_ICS{ind_l} = obj.fcn_dec_rule(obj.MC_test_stat{ind_l}, obj.cv_MC_ICS{ind_l});
                dr_p_MC_NoX{ind_l}   = obj.fcn_dec_rule_p(obj.p_MC_NoX{ind_l});
                dr_p_MC_W{ind_l}     = obj.fcn_dec_rule_p(obj.p_MC_W{ind_l});
                dr_p_MC_S{ind_l}     = obj.fcn_dec_rule_p(obj.p_MC_S{ind_l});
                dr_p_MC_LF{ind_l}    = obj.fcn_dec_rule_p(obj.p_MC_LF{ind_l});
                dr_p_MC_ICS{ind_l}   = obj.fcn_dec_rule_p(obj.p_MC_ICS{ind_l});
            
                dr_LBQ_NoX{ind_l} = obj.fcn_dec_rule(obj.LBQ_test_stat{ind_l}, obj.cv_LBQ_NoX{ind_l});
                dr_LBQ_W{ind_l}   = obj.fcn_dec_rule(obj.LBQ_test_stat{ind_l}, obj.cv_LBQ_W{ind_l});
                dr_LBQ_S{ind_l}   = obj.fcn_dec_rule(obj.LBQ_test_stat{ind_l}, obj.cv_LBQ_S{ind_l});
                dr_LBQ_LF{ind_l}  = obj.fcn_dec_rule(obj.LBQ_test_stat{ind_l}, obj.cv_LBQ_LF{ind_l});
                dr_LBQ_ICS{ind_l} = obj.fcn_dec_rule(obj.LBQ_test_stat{ind_l}, obj.cv_LBQ_ICS{ind_l});
                dr_p_LBQ_NoX{ind_l}   = obj.fcn_dec_rule_p(obj.p_LBQ_NoX{ind_l});
                dr_p_LBQ_W{ind_l}     = obj.fcn_dec_rule_p(obj.p_LBQ_W{ind_l});
                dr_p_LBQ_S{ind_l}     = obj.fcn_dec_rule_p(obj.p_LBQ_S{ind_l});
                dr_p_LBQ_LF{ind_l}    = obj.fcn_dec_rule_p(obj.p_LBQ_LF{ind_l});
                dr_p_LBQ_ICS{ind_l}   = obj.fcn_dec_rule_p(obj.p_LBQ_ICS{ind_l});
            
                dr_sup_LM_NoX{ind_l} = obj.fcn_dec_rule(obj.sup_LM_test_stat{ind_l}, obj.cv_sup_LM_NoX{ind_l});
                dr_sup_LM_W{ind_l}   = obj.fcn_dec_rule(obj.sup_LM_test_stat{ind_l}, obj.cv_sup_LM_W{ind_l});
                dr_sup_LM_S{ind_l}   = obj.fcn_dec_rule(obj.sup_LM_test_stat{ind_l}, obj.cv_sup_LM_S{ind_l});
                dr_sup_LM_LF{ind_l}  = obj.fcn_dec_rule(obj.sup_LM_test_stat{ind_l}, obj.cv_sup_LM_LF{ind_l});
                dr_sup_LM_ICS{ind_l} = obj.fcn_dec_rule(obj.sup_LM_test_stat{ind_l}, obj.cv_sup_LM_ICS{ind_l});
                dr_p_sup_LM_NoX{ind_l}   = obj.fcn_dec_rule_p(obj.p_sup_LM_NoX{ind_l});
                dr_p_sup_LM_W{ind_l}     = obj.fcn_dec_rule_p(obj.p_sup_LM_W{ind_l});
                dr_p_sup_LM_S{ind_l}     = obj.fcn_dec_rule_p(obj.p_sup_LM_S{ind_l});
                dr_p_sup_LM_LF{ind_l}    = obj.fcn_dec_rule_p(obj.p_sup_LM_LF{ind_l});
                dr_p_sup_LM_ICS{ind_l}   = obj.fcn_dec_rule_p(obj.p_sup_LM_ICS{ind_l});
            
                dr_CvM_NoX{ind_l} = obj.fcn_dec_rule(obj.CvM_test_stat{ind_l}, obj.cv_CvM_NoX{ind_l});
                dr_CvM_W{ind_l}   = obj.fcn_dec_rule(obj.CvM_test_stat{ind_l}, obj.cv_CvM_W{ind_l});
                dr_CvM_S{ind_l}   = obj.fcn_dec_rule(obj.CvM_test_stat{ind_l}, obj.cv_CvM_S{ind_l});
                dr_CvM_LF{ind_l}  = obj.fcn_dec_rule(obj.CvM_test_stat{ind_l}, obj.cv_CvM_LF{ind_l});
                dr_CvM_ICS{ind_l} = obj.fcn_dec_rule(obj.CvM_test_stat{ind_l}, obj.cv_CvM_ICS{ind_l});
                dr_p_CvM_NoX{ind_l}   = obj.fcn_dec_rule_p(obj.p_CvM_NoX{ind_l});
                dr_p_CvM_W{ind_l}     = obj.fcn_dec_rule_p(obj.p_CvM_W{ind_l});
                dr_p_CvM_S{ind_l}     = obj.fcn_dec_rule_p(obj.p_CvM_S{ind_l});
                dr_p_CvM_LF{ind_l}    = obj.fcn_dec_rule_p(obj.p_CvM_LF{ind_l});
                dr_p_CvM_ICS{ind_l}   = obj.fcn_dec_rule_p(obj.p_CvM_ICS{ind_l});
            end
            
            % Step 7.2: Consolidate decision rules for table
            % Note: this table is num drs x num alphas = 16 x 3
            for ind_l = 1:len_Ln
                obj.dr_table{ind_l}  = [dr_MC_ICS{ind_l}'; dr_LBQ_ICS{ind_l}'; dr_sup_LM_ICS{ind_l}'; dr_CvM_ICS{ind_l}';
                                        dr_MC_LF{ind_l}'; dr_LBQ_LF{ind_l}'; dr_sup_LM_LF{ind_l}'; dr_CvM_LF{ind_l}';
                                        dr_MC_S{ind_l}'; dr_LBQ_S{ind_l}'; dr_sup_LM_S{ind_l}'; dr_CvM_S{ind_l}';
                                        dr_MC_NoX{ind_l}'; dr_LBQ_NoX{ind_l}'; dr_sup_LM_NoX{ind_l}'; dr_CvM_NoX{ind_l}'];
                obj.dr_p_table{ind_l}= [dr_p_MC_ICS{ind_l}'; dr_p_LBQ_ICS{ind_l}'; dr_p_sup_LM_ICS{ind_l}'; dr_p_CvM_ICS{ind_l}';
                                        dr_p_MC_LF{ind_l}'; dr_p_LBQ_LF{ind_l}'; dr_p_sup_LM_LF{ind_l}'; dr_p_CvM_LF{ind_l}';
                                        dr_p_MC_S{ind_l}'; dr_p_LBQ_S{ind_l}'; dr_p_sup_LM_S{ind_l}'; dr_p_CvM_S{ind_l}';
                                        dr_p_MC_NoX{ind_l}'; dr_p_LBQ_NoX{ind_l}'; dr_p_sup_LM_NoX{ind_l}'; dr_p_CvM_NoX{ind_l}'];
                obj.dr_table_w_only{ind_l}   = [dr_MC_W{ind_l}'; dr_LBQ_W{ind_l}'; dr_sup_LM_W{ind_l}'; dr_CvM_W{ind_l}'];
                obj.dr_p_table_w_only{ind_l} = [dr_p_MC_W{ind_l}'; dr_p_LBQ_W{ind_l}'; dr_p_sup_LM_W{ind_l}'; dr_p_CvM_W{ind_l}'];                    
            end
        end
        
        
        %%%% Max Corr Test
        %%%% Test Stat
        function MC_test_stat = fcn_MC_Tn_hat(obj, rho_vec)
            MC_test_stat = sqrt(obj.n) .* max(abs(rho_vec));
        end

        
        % Ljung-Box Test
        function LBQ = LBQ_fcn(obj, rho_vec)
            len_r = length(rho_vec);
            for h = 1:len_r
                rho_vec(h) = ((obj.n-2)/(obj.n-h)) * (obj.n * (rho_vec(h)^2) - 1);
            end
            LBQ = ((2 * len_r)^(-1/2)) * sum(rho_vec);
        end
        
        % AP Sup-LM Test
        function sup_LM = sup_LM_fcn(obj, rho_vec)
            len_r = length(rho_vec);
            len_l = length(obj.lambda_grid);
            index = (1:len_r)';
            index = repmat(index, 1, len_l);
            lambdas = repmat(obj.lambda_grid, len_r, 1);
            lambdas = lambdas .^(index-1);
            rhos = repmat(rho_vec, 1, len_l);
            lambda_rho = lambdas .* rhos;
            LM = obj.n * (1-obj.lambda_grid.^2) .* sum(lambda_rho,1).^2;
            sup_LM = max(LM);
        end
        
        % CvM Test
        function CvM = CvM_fcn(obj, gamma_vec)
            len = length(obj.pi_grid);
            S = zeros(len,1);
            for l = 1:len
                p = obj.pi_grid(l);
                psi_vec = zeros(length(gamma_vec),1);
                for h = 1:(obj.n-1)  % should h start at 1-obj.n or 1???
                    if h ~= 0
                        psi_vec(h) = (1/(h*pi))*sin(h*p);
                    elseif h == 0
                        psi_vec(h) = p / (2 * pi);
                    end
                end
                S(l) = sqrt(obj.n) * gamma_vec' * psi_vec;
            end
            CvM = sum(S.^2) / floor(len/pi);
        end
        
        %%%%%
        % Correlation and Covariance Calculations
        % Weak Id Expansion
        function [rho_psi_Ln, gamma_vec, bias2] = fcn_rho_psi_Ln_loop(obj, Ln, Z)  % NOTE: For simplicity, this defines all of the correlation vectors: M x b_vec x pi0_vec
            rho_psi_Ln = zeros(Ln, obj.M, length(obj.b_in_vec), length(obj.pi0_in_vec));
            gamma_vec = zeros(Ln, obj.M, length(obj.b_in_vec), length(obj.pi0_in_vec));
            for m = 1:obj.M
                Z1 = Z(:,m);
                for ind_b = 1:length(obj.b_in_vec)
                    b_in = obj.b_in_vec(ind_b);
                    for ind_p = 1:length(obj.pi0_in_vec)
                        pi0_in = obj.pi0_in_vec(ind_p);
                        [rho_psi_Ln(:, m, ind_b, ind_p), gamma_vec(:, m, ind_b, ind_p), bias2] = fcn_rho_psi_Ln(obj, Z1, Ln, b_in, pi0_in);
                    end
                end
            end
        end
        function [rho_psi_Ln, gamma_vec, bias2_vec] = fcn_rho_psi_Ln(obj, Z1, Ln, b_in, pi0_in)
            pi_star_bs = obj.bs_pi_star_fcn(Z1, b_in, pi0_in);
            theta_in = [0, obj.zeta_hat, pi_star_bs];
            d_psi = obj.d_psi_fcn(pi_star_bs); % d_psi x T
            e_hat_0n = obj.fcn_e_hat(theta_in);
            [m_psi_t, H_inv, K_n] = obj.bs_GHK_fcn(theta_in, pi0_in);
            bias1 = (H_inv * K_n * b_in)';
            rho_psi_Ln = zeros(Ln,1);
            gamma_vec = zeros(Ln,1);
            bias2_vec = zeros(Ln,1);
            for h = 1:Ln
                bias2 = obj.fcn_bias2(b_in, pi0_in, h); %obj.bias2_true(h); %obj.fcn_bias2(b_in, h); %calculated in dgp
                bias2_vec(h) = bias2;
                D_n_bs = (1/obj.n) * ( (d_psi(:,1:end-h) * e_hat_0n(1+h:end)) + (d_psi(:,1+h:end) * e_hat_0n(1:end-h)) );
                X_n_bs = e_hat_0n(1:end-h) .* e_hat_0n(1+h:end); 
                [rho_psi_Ln(h), gamma_vec(h)] = obj.rho_psi_n_bs_fcn(h, m_psi_t, H_inv, D_n_bs, X_n_bs, Z1, bias1, bias2);
            end
        end
        function [rho_psi_n_bs, numerator] = rho_psi_n_bs_fcn(obj, h, m_psi_t, H_inv, D_n_bs, X_n_bs, Z1, bias1, bias2)
            % E_t = e_t e_{t-h} - (H_inv * m_psi)' * D_n
            Z2 = Z1(1+h:end);
            G_n_t = m_psi_t(:,1+h:end) - repmat(((1/obj.n) .* sum(m_psi_t,2)),1,(obj.n-h));
            E_t = X_n_bs - (H_inv * G_n_t)' * D_n_bs - repmat( (1 / sqrt(obj.n)) * bias2, (obj.n-h), 1);
            numerator = (1/obj.n) * Z2' * ( E_t - ((1/obj.n) * sum(E_t)) ) + (1 / sqrt(obj.n)) * bias1 * D_n_bs + (1 / sqrt(obj.n)) * bias2;
            rho_psi_n_bs = numerator / obj.sigma2_hat;
        end
        function Xi_bs = bs_Xi_fcn(obj, Z1, pi_in, b_in, pi0_in)
            theta_in = [obj.beta_hat, obj.zeta_hat, pi_in];
            [m_psi_t, H_inv, K_n] = obj.bs_GHK_fcn(theta_in, pi0_in);
            k_psi = sum(obj.num_params(1:2));
            G_n = zeros(k_psi, obj.n);
            for t = 1:obj.n
            G_n(:,t) = m_psi_t(:,t) .* repmat(Z1(t),k_psi,1);
            end
            tempG = G_n - repmat(mean(G_n,2),1,obj.n);
            G_n_temp = (1/sqrt(obj.n)) .* sum(tempG,2);
            Xi_bs = -(1/2) .* (G_n_temp + K_n * b_in)' * (H_inv) * (G_n_temp + K_n * b_in);
        end
        function pi_star_bs = bs_pi_star_fcn(obj, Z1, b_in, pi0_in)
            LB_pi = obj.LB(sum(obj.num_params(1:2))+1:sum(obj.num_params(1:3)));
            UB_pi = obj.UB(sum(obj.num_params(1:2))+1:sum(obj.num_params(1:3)));
            step_size = .1;
            Xi_temp = inf;
            % First get a starting point
            for pi_temp = (LB_pi(1):step_size:UB_pi(1))
                % pi_temp = i; %[i, j];
                Xi_temp2 = obj.bs_Xi_fcn(Z1, pi_temp, b_in, pi0_in);
                if (Xi_temp2 < Xi_temp)
                    Xi_temp = Xi_temp2;
                    pi_init = pi_temp;
                end
            end
            [pi_star_bs, fval_temp, exit_flag_temp] = fmincon(@(pi_in) obj.bs_Xi_fcn(Z1, pi_in, b_in, pi0_in), pi_init, [], [], [], [], LB_pi, UB_pi, [], obj.options);
            if exit_flag_temp < 0
                pi_star_bs = NaN;
            end
        end
        %%%% Strong Id Expansion
        function [rho_theta_Ln, gamma_vec] = fcn_rho_theta_Ln(obj, Ln, Z1)
            d_theta = obj.d_theta_fcn();
            [m_theta_t, J_theta_inv] = obj.bs_GJ_theta(d_theta);
            rho_theta_Ln = zeros(Ln,1);
            gamma_vec = zeros(Ln,1);
            for h = 1:Ln
                [rho_theta_Ln(h), gamma_vec(h)] = obj.rho_theta_n_bs_fcn(h, m_theta_t, J_theta_inv, Z1, d_theta);
            end
        end
        function [rho_theta_n_bs, numerator] = rho_theta_n_bs_fcn(obj, h, m_theta_t, J_inv, Z1, d_theta_in)
            % X_n - J^{-1} * G * D_n
            Z2 = Z1(1+h:end);
            X_n_bs = obj.ehat(1:end-h) .* obj.ehat(1+h:end);
            D_theta_n_bs = (1/obj.n) * ( (d_theta_in(:,1:end-h) * obj.ehat(1+h:end)) + (d_theta_in(:,1+h:end) * obj.ehat(1:end-h)) );
            %%%%
            m_theta_t = m_theta_t(:,1+h:end);
            E_t = X_n_bs - (J_inv * m_theta_t)' * D_theta_n_bs;
            numerator = (1/obj.n) * Z2' * ( E_t - ((1/obj.n) * sum(E_t)) );
            rho_theta_n_bs = numerator / obj.sigma2_hat;
        end
        %%%%% No Expansion
        function [rho_n_no_expansion, gamma_vec] = fcn_rho_Ln_no_expansion(obj, Ln, Z1)
            rho_n_no_expansion = zeros(Ln,1);
            gamma_vec = zeros(Ln,1);
            if nargin <= 2
                for h = 1:Ln
                    [rho_n_no_expansion(h), gamma_vec(h)] = obj.corr_fcn(h);
                end
            elseif nargin == 3
                for h = 1:Ln
                    [rho_n_no_expansion(h), gamma_vec(h)] = obj.corr_fcn(h, Z1);
                end
            end
        end
        function [rho, gamma_h] = corr_fcn(obj, h, w)
            if nargin <= 2
                gamma_h = obj.covar_fcn(h);    
            elseif nargin == 3
                gamma_h = obj.covar_fcn(h, w);
            end
            gamma_0 = obj.covar_fcn(0);
            rho = gamma_h / gamma_0;
        end
        function gamma = covar_fcn(obj, h, w)
            h = abs(h);
            e = obj.ehat; % - mean(obj.ehat);
            E = e(1:(end - h)) .* e((1 + h):end);
            if nargin <= 2
                gamma = (1 / obj.n) * sum(E);
            elseif nargin == 3
                gamma = (1 / obj.n) * sum( w((1 + h):end) .* (E - (1 / obj.n) * sum(E)) );
            end
        end
        
        
        %%%% Functions for Critical Values, P-values, and decision rules
        function cv_LF = fcn_cv_LF(obj, cv_w, cv_s)
            cv_LF = max(cv_w, cv_s);
        end
        function cv_ICS = fcn_cv_ICS(obj, cv_LF, cv_s)
            d_theta = obj.d_theta_fcn();
            [m_theta_t, J_theta_n_inv] = obj.bs_GJ_theta(d_theta);
            V_n = (m_theta_t * m_theta_t')  / obj.n;
            Sigma = J_theta_n_inv * V_n * J_theta_n_inv;
            Sigma_beta = Sigma(1:obj.num_params(1), 1:obj.num_params(1));
            Sigma_beta_inv = Sigma_beta \ eye(obj.num_params(1));
            A_n = sqrt(obj.n * obj.beta_hat' * Sigma_beta_inv * obj.beta_hat / obj.num_params(1));
            if (A_n <= obj.kappa_n)
                cv_ICS = cv_LF;
            elseif (A_n > obj.kappa_n)
                cv_ICS = cv_s;
            end
        end
        function ICS_flag = fcn_ICS(obj) % if ICS_flag == 1, then use strong id, o.w. use LF
            d_theta = obj.d_theta_fcn();
            [m_theta_t, J_theta_n_inv] = obj.bs_GJ_theta(d_theta);
            V_n = (m_theta_t * m_theta_t')  / obj.n;
            Sigma = J_theta_n_inv * V_n * J_theta_n_inv;
            Sigma_beta = Sigma(1:obj.num_params(1), 1:obj.num_params(1));
            Sigma_beta_inv = Sigma_beta \ eye(obj.num_params(1));
            A_n = sqrt(obj.n * obj.beta_hat' * Sigma_beta_inv * obj.beta_hat / obj.num_params(1));
            if (A_n <= obj.kappa_n)
                ICS_flag = 0;
            elseif (A_n > obj.kappa_n)
                ICS_flag = 1;
            end
        end
        function cv_bs_out = fcn_cv_bs(obj, distr_in)
            temp = distr_in(isnan(distr_in)==0);
            temp = sort(temp);
            cv_bs_out = zeros(length(obj.alpha_levels),1);
            for a = 1:length(obj.alpha_levels)
                alpha_level = obj.alpha_levels(a);
                A = floor(length(temp)*(1-alpha_level));
                cv_bs_out(a) = temp(A);
            end
        end
        function p_LF = fcn_p_LF(obj, p_w, p_s)
            p_LF = min(p_w, p_s);
        end
        function p_ICS = fcn_p_ICS(obj, p_LF, p_s)
            d_theta = obj.d_theta_fcn();
            [m_theta_t, J_theta_n_inv] = obj.bs_GJ_theta(d_theta);
            V_n = (m_theta_t * m_theta_t')  / obj.n;
            Sigma = J_theta_n_inv * V_n * J_theta_n_inv;
            Sigma_beta = Sigma(1:obj.num_params(1), 1:obj.num_params(1));
            Sigma_beta_inv = Sigma_beta \ eye(obj.num_params(1));
            A_n = sqrt(obj.n * obj.beta_hat' * Sigma_beta_inv * obj.beta_hat / obj.num_params(1));
            if (A_n <= obj.kappa_n)
                p_ICS = p_LF;
            elseif (A_n > obj.kappa_n)
                p_ICS = p_s;
            end
        end
        function p = pval_fcn(obj, T, M, statistic)
            p = (1/M) * sum((T>statistic));
        end
        function dr = fcn_dec_rule(obj, T, cv)
            dr = (T > cv);
        end
        function dr = fcn_dec_rule_p(obj, p)
            dr = false(length(obj.alpha_levels),1); % false preallocates a boolean array with 0s
            for a = 1:length(obj.alpha_levels)
                alpha_level = obj.alpha_levels(a);
                dr(a) = (p < alpha_level);
            end
        end

        
        % Dependent Wild Bootstrap
        function w = DWB_w_fcn(obj)
            num_windows = ceil((obj.n) / obj.kn);
            w1 = randn(1,num_windows);
            w1 = repmat(w1,obj.kn,1);
            w1 = w1(:);
            w = w1((end-obj.n+1):end);
        end
        function wM = DWB_wM_fcn(obj) 
            wM = zeros(obj.n, obj.M);
            for m = 1:obj.M
                wM(:,m) = obj.DWB_w_fcn;
            end
        end
        function test(obj, a, b) %Just to check argument numbers
            display(nargin);
        end

    end
end
