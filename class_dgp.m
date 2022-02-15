classdef class_dgp < class_innovations
    properties (SetAccess=public, GetAccess=public)
        Y
        Yhat
        ehat
        theta0
        theta_hat
        sigma2_hat
        sigma_hat
        dgp_type
        dgp_type_string
        id_type
        id_type_string
    end
    properties (SetAccess=public, GetAccess=public, Hidden = true)
        b_true
        Te
        beta0
        pi0
        mu0
        zeta0
        beta_hat
        pi_hat
        mu_hat
        zeta_hat
        chi_hat
	    LB
        UB
        %lag_num
        num_params
        num_lags
        theta_init
        fval
        exit_flag
        bias2_true, bias2_est, bias2_est_old
    end
    methods
        function out = about(obj)
            display('This class instantiates the dgp');
        end
        function obj = class_dgp(dgp_type, id_type, innovation_type, theta_in, num_params, T, seed)
            %First instantiate the innovations
            obj@class_innovations(innovation_type, T, seed);
            %Next get the thetas in order
            obj.theta0 = theta_in;
            obj.num_params = num_params;
            obj.num_lags = max(obj.num_params);
            obj.Te = obj.T - obj.num_lags;
            % tell us what things are
            obj.dgp_type = dgp_type;
            obj.id_type = id_type;
            if id_type == 0
                obj.id_type_string = 'No Id';
            elseif id_type == 1
                obj.id_type_string = 'Weak Id';
            elseif id_type == 2
                obj.id_type_string = 'Strong Id';
            end
            % Now generate the data
            switch obj.dgp_type
                case 1
                    obj.dgp_type_string = 'STAR1';
                    obj = obj.dgp_Y_STAR1;
                case 2
                    obj.dgp_type_string = 'STAR2';
                    obj = obj.dgp_Y_STAR2;
                case 3
                    obj.dgp_type_string = 'ARMA';
                    obj = obj.dgp_Y_ARMA11;
                otherwise
                    error('Please provide appropriate DGP');
            end
            obj.Yhat = zeros(obj.T,1);
            obj.ehat = zeros(obj.Te,1);
        end
        function obj = estimation(obj, LB, UB)
            options = optimoptions(@fmincon,'Display','off', 'Diagnostics','off');
            obj.LB = LB';
            obj.UB = UB';
            switch obj.dgp_type
                case 1
                    temp = sort(obj.Y);
                    obj.LB(3) = temp(floor(.15*length(temp)));
                    obj.UB(3) = temp(floor(.85*length(temp)));
                    obj.theta_init = (obj.UB - obj.LB) .* rand(1,length(obj.LB))' + obj.LB;
                    [obj.theta_hat, obj.fval, obj.exit_flag] = fmincon(@(theta) obj.loss_fcn_STAR1(theta), obj.theta_init, [1, 1, 0], 1, [], [], obj.LB, obj.UB, [], options);
                    obj.theta_hat = obj.theta_hat';
                case 2
                    temp = sort(obj.Y);
                    obj.LB(3) = temp(floor(.15*length(temp)));
                    obj.UB(3) = temp(floor(.85*length(temp)));
                    obj.theta_init = (obj.UB - obj.LB) .* rand(1,length(obj.LB))' + obj.LB;
                    [obj.theta_hat, obj.fval, obj.exit_flag] = fmincon(@(theta) obj.loss_fcn_STAR1(theta), obj.theta_init, [1, 1, 0], 1, [], [], obj.LB, obj.UB, [], options);
                    obj.theta_hat = obj.theta_hat';
                case 3
                    obj.theta_init = (obj.UB - obj.LB) .* rand(1,length(obj.LB))' + obj.LB;
                    [obj.theta_hat, obj.fval, obj.exit_flag] = fmincon(@(theta) obj.loss_fcn_ARMA11(theta), obj.theta_init, [1, 0, 1], 1, [], [], obj.LB, obj.UB, [], options);
                    obj.theta_hat = obj.theta_hat';
            end
            %%%%
            obj.Yhat = obj.fcn_Yhat(obj.theta_hat);
            obj.ehat = obj.fcn_e_hat(obj.theta_hat);
            obj.sigma2_hat = (1/obj.Te) * (obj.ehat' * obj.ehat);
            obj.sigma_hat = sqrt(obj.sigma2_hat); %1; 
            %%%%
            obj.beta_hat = obj.theta_hat(1:obj.num_params(1));
            obj.zeta_hat = obj.theta_hat((obj.num_params(1)+1):sum(obj.num_params(1:2)));
            obj.pi_hat = obj.theta_hat((sum(obj.num_params(1:2))+1):sum(obj.num_params(1:3)));
            %obj.mu_hat = obj.theta0((sum(obj.num_params(1:3))+1):sum(obj.num_params(1:4)));
        end
        function Yhat = fcn_Yhat(obj, theta_in)
            switch obj.dgp_type
                case 1
                    Yhat = obj.fcn_Yhat_STAR1(theta_in);
                case 2
                    Yhat = obj.fcn_Yhat_STAR1(theta_in);
                case 3
                    Yhat = obj.fcn_Yhat_ARMA11(theta_in);
            end
        end
        function ehat = fcn_e_hat(obj, theta_in) %Note that the input was changed from Yhat to theta_in.  Change the test class code accordingly.
            switch obj.dgp_type
                case 1
                    ehat = obj.fcn_e_hat_STAR1(theta_in);
                case 2
                    ehat = obj.fcn_e_hat_STAR1(theta_in);
                case 3
                    ehat = obj.fcn_e_hat_ARMA11(theta_in);
            end
        end
        %%%% Derivatives
        function d_psi = d_psi_fcn(obj, pi_in)
            switch obj.dgp_type
                case 1
                    d_psi = obj.d_psi_fcn_STAR1(pi_in);
                case 2
                    d_psi = obj.d_psi_fcn_STAR1(pi_in);
                case 3
                    d_psi = obj.d_psi_fcn_ARMA11(pi_in);
            end
        end
        function d_theta = d_theta_fcn(obj)
            switch obj.dgp_type
                case 1
                    d_theta = obj.d_theta_fcn_STAR1;
                case 2
                    d_theta = obj.d_theta_fcn_STAR1;
                case 3
                    d_theta = obj.d_theta_fcn_ARMA11;
            end
        end
        function [m_psi_t, H_inv, K_n] = bs_GHK_fcn(obj, theta_in, pi0_in)
            switch obj.dgp_type
                case 1
                    [m_psi_t, H_inv, K_n] = obj.bs_GHK_fcn_STAR1(theta_in, pi0_in);
                case 2
                    [m_psi_t, H_inv, K_n] = obj.bs_GHK_fcn_STAR1(theta_in, pi0_in);
                case 3
                    [m_psi_t, H_inv, K_n] = obj.bs_GHK_fcn_ARMA11(theta_in, pi0_in);
            end
        end
        function [m_theta_t, J_theta_n_inv] = bs_GJ_theta(obj, d_theta)
            switch obj.dgp_type
                case 1
                    [m_theta_t, J_theta_n_inv] = obj.bs_GJ_theta_STAR1(d_theta);
                case 2
                    [m_theta_t, J_theta_n_inv] = obj.bs_GJ_theta_STAR1(d_theta);
                case 3
                    [m_theta_t, J_theta_n_inv] = obj.bs_GJ_theta_ARMA11(d_theta);
            end
        end
        function bias2 = fcn_bias2(obj, b_in, pi_in, h)
            switch obj.dgp_type
                case 1
                    bias2 = obj.fcn_bias2_STAR1(b_in, pi_in, h);
                case 2
                    bias2 = obj.fcn_bias2_STAR1(b_in, pi_in, h);
                case 3
                    bias2 = obj.fcn_bias2_ARMA11(b_in, pi_in, h);
            end
        end
        function bias2_true = fcn_bias2_true(obj, h)
            switch obj.dgp_type
                case 1
                    bias2_true = obj.fcn_bias2_true_STAR1(h);
                case 2
                    bias2_true = obj.fcn_bias2_true_STAR1(h);
                case 3
                    bias2_true = obj.fcn_bias2_true_ARMA11(h);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% STAR1
        function Q = loss_fcn_STAR1(obj, theta_in)
            theta_in = theta_in';
            %Yhat_temp = obj.fcn_Yhat_STAR1(theta_in);
            ehat_temp = obj.fcn_e_hat(theta_in);
            Q = (1/(obj.Te)) * sum((ehat_temp.^2));
            %dQ = 0; %later
        end
        function ehat = fcn_e_hat_STAR1(obj, theta_in)
            temp_Yhat = obj.fcn_Yhat_STAR1(theta_in);
            ehat = obj.Y - temp_Yhat;
            ehat = ehat(1+obj.num_lags:obj.T);
        end
        function Yhat = fcn_Yhat_STAR1(obj, theta_in)
            beta = theta_in(1:obj.num_params(1));
            zeta = theta_in((obj.num_params(1)+1):sum(obj.num_params(1:2)));
            pi = theta_in((sum(obj.num_params(1:2))+1):sum(obj.num_params(1:3)));
            mu = obj.theta0((sum(obj.num_params(1:3))+1):sum(obj.num_params(1:4)));
            %%%%
            Yhat = zeros(obj.T,1);
            Yhat(1:obj.num_lags) = obj.Y(1:obj.num_lags);
            for t = (obj.num_lags+1):obj.T
                XB   = beta' * obj.Y(t-obj.num_params(1):t-1,:);
                Xpi  = obj.Y(t-obj.num_params(3):t-1,:);
                Gt   = 1 ./ (1 + exp(- mu' * (Xpi-pi)));
                Zxi  = zeta' * obj.Y(t-obj.num_params(2):t-1,:);
                Yhat(t) = XB .* (Gt) + Zxi;
            end;
        end
        function obj = dgp_Y_STAR1(obj)
            % yt = Gt(pi, Xpi) * Xt' * B + Zt' * xi + et
            obj.beta0 = obj.theta0(1:obj.num_params(1));
            obj.zeta0 = obj.theta0((obj.num_params(1)+1):sum(obj.num_params(1:2)));
            obj.pi0 = obj.theta0((sum(obj.num_params(1:2))+1):sum(obj.num_params(1:3)));
            obj.mu0 = obj.theta0((sum(obj.num_params(1:3))+1):sum(obj.num_params(1:4)));
            %%%
            Y0 = zeros(obj.init_T, 1);
            for t = (obj.num_lags+1):obj.init_T
                XB  = obj.beta0' * Y0(t-obj.num_params(1):t-1);
                Xpi = Y0(t-obj.num_params(3):t-1);
                Gt  = 1 ./ (1 + exp(- obj.mu0' * (Xpi-obj.pi0)));
                Zzeta = obj.zeta0' * Y0(t-obj.num_params(2):t-1);
                Y0(t) = XB .* (Gt) + Zzeta + obj.U(t);
            end
            obj.Y = Y0(end-obj.T+1:end); % remove burn-in values
        end
        
        %%%% Now the derivatives and other objects that will be used in the
        %%%% max corr test
        
        %Bias2 term
        function bias2_true = fcn_bias2_true_STAR1(obj, h)
            theta_in = obj.theta0(1:sum(obj.num_params(1:3)));
            theta_in(1) = 0;
            e_hat_0n = fcn_e_hat_STAR1(obj, theta_in);
            e = obj.U(end-obj.Te+1:end);
            bias2_true = (e_hat_0n(1:end-h) .* e_hat_0n(1+h:end)) - (e(1:end-h) .* e(1+h:end));
            bias2_true = (1/sqrt(obj.Te)) * sum(bias2_true);
        end
%         function bias2 = fcn_bias2_old(obj, b_in, h)
%             bias2 = b_in * obj.zeta_hat^(h) * obj.sigma2_hat;
%         end
        function bias2 = fcn_bias2_STAR1(obj, b_in, pi0_in, h)
            g_t = zeros(obj.Te, 1);
            for t = (obj.num_lags+1):obj.Te
                z_t    = obj.Y(t-obj.num_params(3):t-1,:);
                g_t(t) = 1 ./ (1 + exp(- obj.mu0' * (z_t - pi0_in)));
            end
            temp = sum(obj.ehat(1:end-h) .* obj.ehat(1:end-h) .* g_t(1+h:end)) / obj.Te;
            bias2 = b_in * obj.zeta_hat^(h-1) * temp;
        end
        %Weak ID
        function d_psi = d_psi_fcn_STAR1(obj, pi_in)
            %d_psi,t = -[x_t' g(z_t, \pi), x_t']'
            pi = pi_in;
            mu = obj.mu0;
            %%%%
            d_psi = zeros(sum(obj.num_params(1:2)), obj.T);
            for t = (obj.num_lags+1):obj.T
                x_t   = obj.Y(t-obj.num_params(1):t-1,:);
                z_t   = obj.Y(t-obj.num_params(3):t-1,:);
                g_t   = 1 ./ (1 + exp(- mu' * (z_t - pi)));
                d_psi(:,t) = -[x_t' * g_t, x_t']';
            end
            d_psi = d_psi(:,obj.num_lags+1:end);
        end
        function [m_psi_t, H_inv, K_n] = bs_GHK_fcn_STAR1(obj, theta_in, pi0_in)
            k_psi = sum(obj.num_params(1:2));
            pi_in = theta_in((k_psi+1):end);
            S_beta = [eye(obj.num_params(1)); zeros(obj.num_params(2),obj.num_params(1))];
            d_beta_at_pi0 = S_beta' * obj.d_psi_fcn(pi0_in); % d_beta x T
            d_psi = obj.d_psi_fcn(pi_in); % d_psi x T
            H_n = (1/obj.Te) * (d_psi * d_psi'); % d_psi x d_psi
            H_inv = H_n \ eye(size(H_n,1));
            %H_inv = pinv(H_n);
            K_n = -(1/obj.Te) * d_psi * d_beta_at_pi0'; % d_psi x d_beta
            %y_hat = obj.fcn_Yhat(theta_in);
            e_hat = obj.fcn_e_hat(theta_in);
            m_psi_t = zeros(k_psi, obj.Te);
            for t = 1:obj.Te
                m_psi_t(:,t) = d_psi(:,t) .* repmat(e_hat(t),k_psi,1);
            end
        end

        % Strong ID
        function d_theta = d_theta_fcn_STAR1(obj)
            %d_theta,t = -[x_t' g(z_t, \pi), x_t', (\beta' x_t g_{\pi}(z_t, \pi))']'
            beta = obj.beta_hat;
            pi = obj.pi_hat;
            mu = obj.mu0;
            %%%%
            d_theta = zeros(sum(obj.num_params(1:3)), obj.T);
            for t = (obj.num_lags+1):obj.T
                x_t   = obj.Y(t-obj.num_params(1):t-1,:);
                z_t   = obj.Y(t-obj.num_params(3):t-1,:);
                g_t   = 1 ./ (1 + exp(- mu' * (z_t - pi)));
                g_pit = -(1 + exp(- mu' * (z_t - pi)))^(-2) * exp(- mu' * (z_t - pi)) * (mu);
                if isnan(g_pit) == 1
                    g_pit = 0;
                end
                % (1 + exp(- mu' * (z_t - pi)))^{-1}
                % -(1 + exp(- mu' * (z_t - pi)))^{-2} * exp(- mu' * (z_t - pi)) * (--mu)
                % -(1 + exp(- mu' * (z_t - pi)))^{-2} * exp(- mu' * (z_t - pi)) * (mu)
                d_theta(:,t) = -[x_t' * g_t, x_t', (x_t * g_pit)']'; %,(beta' * x_t * g_pit)'
            end
            d_theta = d_theta(:,obj.num_lags+1:end);
        end
        function [m_theta_t, J_theta_n_inv] = bs_GJ_theta_STAR1(obj, d_theta)
            %d_theta = obj.d_theta_fcn();
    %         B = [ones(obj.num_params(1),1); ones(obj.num_params(2),1); obj.beta_hat];
    %         B = diag(B);
    %         B_inv = pinv(B);
    %         J_theta_n = B_inv * (1/obj.Te) * (d_theta * d_theta') * B_inv;
    %         G_theta_n = B_inv * sigma_in * (1/sqrt(obj.Te)) * d_theta * Z';
            J_theta_n = (1/obj.Te) * (d_theta * d_theta');
            % J_inv = pinv(J_in);
            J_theta_n_inv = J_theta_n \ eye(size(J_theta_n,1));
            p = sum(obj.num_params(1:3));
            m_theta_t = zeros(p, obj.Te);
            for t = 1:obj.Te
                m_theta_t(:,t) = d_theta(:,t) .* repmat(obj.ehat(t),p,1);
            end
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% STAR2
        function obj = dgp_Y_STAR2(obj)
            % yt = Gt(pi, Xpi) * Xt' * B + Zt' * xi + et
            obj.beta0 = obj.theta0(1:obj.num_params(1));
            obj.zeta0 = obj.theta0((obj.num_params(1)+1):sum(obj.num_params(1:2)));
            obj.pi0 = obj.theta0((sum(obj.num_params(1:2))+1):sum(obj.num_params(1:3)));
            obj.mu0 = obj.theta0((sum(obj.num_params(1:3))+1):sum(obj.num_params(1:4)));
            zeta_temp = [obj.zeta0 (.15 / sqrt(obj.Te))];
            %%%
            Y0 = zeros(obj.init_T, 1);
            for t = (max(obj.num_lags,length(zeta_temp))+1):obj.init_T
                XB  = obj.beta0' * Y0(t-obj.num_params(1):t-1);
                Xpi = Y0(t-obj.num_params(3):t-1);
                Gt  = 1 ./ (1 + exp(- obj.mu0' * (Xpi-obj.pi0)));
                Zzeta = zeta_temp * Y0((t-length(zeta_temp)):t-1);
                Y0(t) = XB .* (Gt) + Zzeta + obj.U(t);
            end
            obj.Y = Y0(end-obj.T+1:end); % remove burn-in values
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% ARMA11
        function Q = loss_fcn_ARMA11(obj, theta_in)
            beta = theta_in(1:obj.num_params(1));
            zeta = theta_in((obj.num_params(1)+1):sum(obj.num_params(1:2)));
            pi = theta_in((sum(obj.num_params(1:2))+1):sum(obj.num_params(1:3)));
            Q_out = zeros(obj.T,1);
            lag_num = 1;
            for t = (lag_num+1):obj.T
                clear temp ind;
                ind = -(-t+lag_num:lag_num:-1);
                pi_vec = zeros(size(ind,2),1);
                for j = 1:size(ind,2)
                    pi_vec(j) = pi^(j-1);
                end;
                temp = pi_vec' * obj.Y(ind);
                Q_out(t) = (1/2) * log(zeta) + (1/(2*zeta)) * (obj.Y(t) - beta * temp)^2;
            end;
            Q = (1/(obj.T-lag_num)) * sum(Q_out);
            %%%%
            %[rho_beta, rho_zeta, rho_pi] = obj.fcn_DQ_ARMA11(theta_in);
            %dQ = (1/(obj.T-lag_num)) * [sum(rho_beta), sum(rho_zeta), sum(rho_pi)]';
        end
        function ehat = fcn_e_hat_ARMA11(obj, theta_in) % Fix this!!! need zeta, also see if anything else if missing...
            beta = theta_in(1:obj.num_params(1));
            zeta = theta_in((obj.num_params(1)+1):sum(obj.num_params(1:2)));
            pi = theta_in((sum(obj.num_params(1:2))+1):sum(obj.num_params(1:3)));
            ehat = zeros(obj.T,1);
            lag_num = 1;
            for t = (lag_num+1):obj.T
                clear temp ind;
                ind = -(-t+lag_num:lag_num:-1);
                pi_vec = zeros(size(ind,2),1);
                for j = 1:size(ind,2)
                    pi_vec(j) = pi^(j-1);
                end
                temp = pi_vec' * obj.Y(ind);
                ehat(t) = (zeta^(-1/2)) .* (obj.Y(t) - beta * temp);
            end    
            ehat = ehat(1+obj.num_lags:obj.T);
        end
        function Yhat_out = fcn_Yhat_ARMA11(obj, theta_in)
            beta = theta_in(1:obj.num_params(1));
            zeta = theta_in((obj.num_params(1)+1):sum(obj.num_params(1:2)));
            pi = theta_in((sum(obj.num_params(1:2))+1):sum(obj.num_params(1:3)));
            temp_ehat = obj.fcn_e_hat_ARMA11(theta_in);
            temp_ehat = [0; temp_ehat];
            Yhat_out = zeros(obj.T, 1);
            Yhat_out(1) = obj.Y(1);
            for t = 2:obj.T
                Yhat_out(t) = (pi + beta) .* Yhat_out(t-1) + (zeta^(1/2)) .* temp_ehat(t) - pi .* (zeta^(1/2)) .* temp_ehat(t-1);
            end
            Yhat_out = Yhat_out(1+obj.num_lags:end);
        end
        function obj = dgp_Y_ARMA11(obj)  %Note: Need to make sure that obj.U has var = 1;
            obj.beta0 = obj.theta0(1:obj.num_params(1));
            obj.zeta0 = obj.theta0((obj.num_params(1)+1):sum(obj.num_params(1:2)));
            obj.pi0 = obj.theta0((sum(obj.num_params(1:2))+1):sum(obj.num_params(1:3)));
            Y0 = zeros(obj.init_T, 1);
            e = (obj.zeta0^(1/2)) .* obj.U;
            for t = 2:obj.init_T
                Y0(t) = (obj.pi0 + obj.beta0) .* Y0(t-1) + e(t) - obj.pi0 .* e(t-1);
            end
            obj.Y = Y0(end-obj.T+1:end);
        end
        %%%% Now the derivatives and other objects that will be used in the
        %%%% max corr test
        %%%% See AC 2012 page 50 for these quantities
        function [d_b, d_z, d_p] = fcn_D_eps_ARMA11(obj, theta_in)
            beta = theta_in(1:obj.num_params(1));
            zeta = theta_in((obj.num_params(1)+1):sum(obj.num_params(1:2)));
            pi = theta_in((sum(obj.num_params(1:2))+1):sum(obj.num_params(1:3)));
            %B = [ones(obj.num_params(1),1); ones(obj.num_params(2),1); obj.beta_hat];
            %B = diag(B);
            %B_inv = B \ eye(size(B,1));
            d_b = zeros(obj.T,1);
            d_z = zeros(obj.T,1);
            d_p = zeros(obj.T,1);
            lag_num = 1;
            for t = (lag_num+1):obj.T
                clear temp ind;
                ind = -(-t+lag_num:lag_num:-1);
                pi_vec = zeros(size(ind,2),1);
                pi_vec_d = zeros(size(ind,2),1);
                for j = 1:size(ind,2)
                    pi_vec(j) = pi^(j-1);
                    pi_vec_d(j) = (j-1) * pi^(j-2);
                end;
                temp = pi_vec' * obj.Y(ind);
                temp2 = pi_vec_d' * obj.Y(ind);
                d_b(t) = -zeta^(-1/2) * temp;
                d_z(t) = -(1/2) * zeta^(-3/2) * (obj.Y(t) - beta * temp);
                d_p(t) = -zeta^(-1/2) * temp2; % -zeta^(-1/2) * beta * temp2;
            end;
            d_b = d_b((lag_num+1):end);
            d_z = d_z((lag_num+1):end);
            d_p = d_p((lag_num+1):end);
        end
        function [rho_beta, rho_zeta, rho_pi] = fcn_DQ_ARMA11(obj, theta_in)
            beta = theta_in(1:obj.num_params(1));
            zeta = theta_in((obj.num_params(1)+1):sum(obj.num_params(1:2)));
            pi = theta_in((sum(obj.num_params(1:2))+1):sum(obj.num_params(1:3)));
            rho_beta = zeros(obj.T,1);
            rho_zeta = zeros(obj.T,1);
            rho_pi = zeros(obj.T,1);
            lag_num = 1;
            for t = (lag_num+1):obj.T
                clear temp ind;
                ind = -(-t+lag_num:lag_num:-1);
                pi_vec = zeros(size(ind,2),1);
                pi_vec_d = zeros(size(ind,2),1);
                for j = 1:size(ind,2)
                    pi_vec(j) = pi^(j-1);
                    pi_vec_d(j) = (j-1) * pi^(j-2);
                end;
                temp = pi_vec' * obj.Y(ind);
                temp2 = pi_vec_d' * obj.Y(ind);
                rho_beta(t) = -zeta^(-1) * (obj.Y(t) - beta * temp) * temp;
                rho_zeta(t) = -(1/2) * zeta^(-2) * ((obj.Y(t) - beta * temp)^2 - zeta);
                rho_pi(t)   = -zeta^(-1) * (obj.Y(t) - beta * temp) * temp2; % -zeta^(-1) * (obj.Y(t) - beta * temp) * beta * temp2; 
            end;
            rho_beta = rho_beta((lag_num+1):end);
            rho_zeta = rho_zeta((lag_num+1):end);
            rho_pi = rho_pi((lag_num+1):end);
        end
        function [rho_bb, rho_bz, rho_zz, rho_bp, rho_zp, rho_pp] = fcn_D2Q_ARMA11(obj, theta_in)
            beta = theta_in(1:obj.num_params(1));
            zeta = theta_in((obj.num_params(1)+1):sum(obj.num_params(1:2)));
            pi = theta_in((sum(obj.num_params(1:2))+1):sum(obj.num_params(1:3)));
            rho_bb = zeros(obj.T,1);
            rho_bz = zeros(obj.T,1);
            rho_zz = zeros(obj.T,1);
            rho_bp = zeros(obj.T,1);
            rho_zp = zeros(obj.T,1);
            rho_pp = zeros(obj.T,1);
            lag_num = 1;
            for t = (lag_num+1):obj.T
                clear temp ind;
                ind = -(-t+lag_num:lag_num:-1);
                pi_vec = zeros(size(ind,2),1);
                pi_vec_d = zeros(size(ind,2),1);
                pi_vec_d2 = zeros(size(ind,2),1);
                for j = 1:size(ind,2)
                    pi_vec(j) = pi^(j-1);
                    pi_vec_d(j) = (j-1) * pi^(j-2);
                    pi_vec_d2(j) = (j-1) * (j-2) * pi^(j-3);
                end;
                temp = pi_vec' * obj.Y(ind);
                temp2 = pi_vec_d' * obj.Y(ind);
                temp3 = pi_vec_d2' * obj.Y(ind);
                rho_bb(t) = zeta^(-1) * temp * temp;
                rho_bz(t) = zeta^(-2) * (obj.Y(t) * temp);
                rho_zz(t) = -(1/2) * zeta^(-2) + zeta^(-3) * obj.Y(t)^2;
                rho_bp(t) = zeta^(-1) * temp2 * (2 * beta * temp - obj.Y(t));
                rho_zp(t) = zeta^(-2) * (obj.Y(t) - beta * temp) * beta * temp2;
                rho_pp(t) = zeta^(-1) * (temp2)^2; % zeta^(-1) * ((beta * temp2)^2 - ((obj.Y(t) - beta * temp) * beta * temp3));
            end;
            rho_bb = rho_bb((lag_num+1):end);
            rho_bz = rho_bz((lag_num+1):end);
            rho_zz = rho_zz((lag_num+1):end);
            rho_bp = rho_bp((lag_num+1):end);
            rho_zp = rho_zp((lag_num+1):end);
            rho_pp = rho_pp((lag_num+1):end);
        end
        %Weak ID
        function bias2_test = fcn_bias2_true_ARMA11(obj, h)
            theta_in = obj.theta0(1:sum(obj.num_params(1:3)));
            theta_in(1) = 0;
            e_hat_0n = fcn_e_hat_ARMA11(obj, theta_in);
            e = obj.U(end-obj.Te+1:end);
            bias2_test = (e_hat_0n(1:end-h) .* e_hat_0n(1+h:end)) - (e(1:end-h) .* e(1+h:end));
            bias2_test = (1/sqrt(obj.Te)) * sum(bias2_test);
        end
        function bias2 = fcn_bias2_ARMA11(obj, b_in, pi0_in, h)
            bias2 = b_in * pi0_in^(h-1) * obj.zeta_hat;
        end
        function d_psi = d_psi_fcn_ARMA11(obj, pi_in)
            theta_in = [obj.beta_hat; obj.zeta_hat; pi_in];
            [d_b, d_z, ~] = obj.fcn_D_eps_ARMA11(theta_in);
            %d_psi = zeros(sum(obj.num_params(1:2)),length(d_b));
            d_psi = [d_b'; d_z'];
        end
        function [m_psi_t, H_inv, K_n] = bs_GHK_fcn_ARMA11(obj, theta_in, pi0_in)
            %beta = theta_in(1:obj.num_params(1));
            %zeta = theta_in((obj.num_params(1)+1):sum(obj.num_params(1:2)));
            pi = theta_in((sum(obj.num_params(1:2))+1):sum(obj.num_params(1:3)));
            [rho_beta, rho_zeta, ~] = obj.fcn_DQ_ARMA11(theta_in);
            [rho_bb, rho_bz, rho_zz, ~, ~, ~] = obj.fcn_D2Q_ARMA11(theta_in);
            m_psi_t = [rho_beta, rho_zeta]';
            r_bb = mean(rho_bb); r_bz = mean(rho_bz); r_zz = mean(rho_zz);
            H_n = [r_bb r_bz; r_bz r_zz];
            H_inv = H_n \ eye(size(H_n,1));
            K_n = [-(1 - (pi0_in * pi))^(-1) ; 0]; % Add in estimator later
        end
        % Strong ID
        function d_theta = d_theta_fcn_ARMA11(obj)
            [d_b, d_z, d_p] = obj.fcn_D_eps_ARMA11(obj.theta_hat);
            d_theta = [d_b'; d_z'; d_p'];
        end
        function [m_theta_t, J_theta_n_inv] = bs_GJ_theta_ARMA11(obj, d_theta)
            theta_in = obj.theta_hat;
            [rho_beta, rho_zeta, rho_pi] = obj.fcn_DQ_ARMA11(theta_in);
            [rho_bb, rho_bz, rho_zz, rho_bp, rho_zp, rho_pp] = obj.fcn_D2Q_ARMA11(theta_in);
            m_theta_t = [rho_beta, rho_zeta, rho_pi]';
            J_theta_n = zeros(sum(obj.num_params(1:3)),sum(obj.num_params(1:3)),length(rho_bb));
            J_theta_n(1,1,:) = rho_bb; J_theta_n(1,2,:) = rho_bz; J_theta_n(1,3,:) = rho_bp;
            J_theta_n(2,1,:) = rho_bz; J_theta_n(2,2,:) = rho_zz; J_theta_n(2,3,:) = rho_zp;
            J_theta_n(3,1,:) = rho_bp; J_theta_n(3,2,:) = rho_zp; J_theta_n(3,3,:) = rho_pp;
            J_theta_n = mean(J_theta_n,3);
            J_theta_n_inv = J_theta_n \ eye(sum(obj.num_params(1:3)));
        end

    end
 end
