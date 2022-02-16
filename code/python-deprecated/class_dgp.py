import numpy as np
from scipy.optimize import minimize
from class_innovations import class_innovations

class class_dgp(class_innovations):
    """This class instantiates the dgp"""
    
    def __init__(self, dgp_type, id_type, innovation_type, theta_in, num_params, T, seed_input = 'shuffle'):
        # First instantiate the innovations
        super(class_dgp,self).__init__(innovation_type, T, seed_input);
        # Next get the thetas in order
        self.theta0 = theta_in;
        self.num_params = num_params;
        self.num_lags = max(num_params);
        self.Te = self.T - self.num_lags;  # For residuals, we work with Te rather than T
        # tell us what things are
        self.dgp_type = dgp_type;
        self.id_type = id_type;
        if id_type == 0:
            self.id_type_string = 'No Id';
        elif id_type == 1:
            self.id_type_string = 'Weak Id';
        elif id_type == 2:
            self.id_type_string = 'Strong Id';
        # Now generate the data
        if self.dgp_type == 1:
            self.dgp_type_string = 'STAR1';
            self.dgp_Y_STAR1();
        elif self.dgp_type == 2:
            self.dgp_type_string = 'STAR2';
            self.dgp_Y_STAR2();
        elif self.dgp_type == 3:
            self.dgp_type_string = 'ARMA';
            self.dgp_Y_ARMA11();
        else:
            print('Please provide appropriate DGP');

    ## Pointer functions
    def fcn_Yhat(self, theta_in):
        if self.dgp_type == 1:
            Yhat = self.fcn_Yhat_STAR1(theta_in);
        elif self.dgp_type == 2:
            Yhat = self.fcn_Yhat_STAR1(theta_in);
        elif self.dgp_type == 3:
            Yhat = self.fcn_Yhat_ARMA11(theta_in);
        return Yhat;

    def fcn_e_hat(self, theta_in): #Note that the input was changed from Yhat to theta_in.  Change the test class code accordingly.
        if self.dgp_type == 1:
            ehat = self.fcn_e_hat_STAR1(theta_in);
        elif self.dgp_type == 2:
            ehat = self.fcn_e_hat_STAR1(theta_in);
        elif self.dgp_type == 3:
            ehat = self.fcn_e_hat_ARMA11(theta_in);
        return ehat;
    
    def d_psi_fcn(self, pi_in):
        if self.dgp_type == 1:
            d_psi = self.d_psi_fcn_STAR1(pi_in);
        elif self.dgp_type == 2:
            d_psi = self.d_psi_fcn_STAR1(pi_in);
        elif self.dgp_type == 3:
            d_psi = self.d_psi_fcn_ARMA11(pi_in);
        return d_psi;

    def d_theta_fcn(self):
        if self.dgp_type == 1:
            d_theta = self.d_theta_fcn_STAR1;
        elif self.dgp_type == 2:
            d_theta = self.d_theta_fcn_STAR1;
        elif self.dgp_type == 3:
            d_theta = self.d_theta_fcn_ARMA11;
        return d_theta;

    def bs_GHK_fcn(self, theta_in, pi0_in):
        if self.dgp_type == 1:
            mHK = self.bs_GHK_fcn_STAR1(theta_in, pi0_in);
        elif self.dgp_type == 2:
            mHK = self.bs_GHK_fcn_STAR1(theta_in, pi0_in);
        elif self.dgp_type == 3:
            mHK = self.bs_GHK_fcn_ARMA11(theta_in, pi0_in);
        return mHK; # mHK = { 'm_psi_t':m_psi_t, 'H_inv':H_inv, 'K_n':K_n };

    def bs_GJ_theta(self, d_theta):
        if self.dgp_type == 1:
            mJ = self.bs_GJ_theta_STAR1(d_theta);
        elif self.dgp_type == 2:
            mJ = self.bs_GJ_theta_STAR1(d_theta);
        elif self.dgp_type == 3:
            mJ = self.bs_GJ_theta_ARMA11(d_theta);
        return mJ; # mJ = { 'm_theta_t':m_theta_t, 'J_theta_n_inv':J_theta_n_inv };

    def fcn_bias2(self, b_in, pi_in, h):
        if self.dgp_type == 1:
            bias2 = self.fcn_bias2_STAR1(b_in, pi_in, h);
        elif self.dgp_type == 2:
            bias2 = self.fcn_bias2_STAR1(b_in, pi_in, h);
        elif self.dgp_type == 3:
            bias2 = self.fcn_bias2_ARMA11(b_in, pi_in, h);
        return bias2;

    def fcn_bias2_true(self, h):
        if self.dgp_type == 1:
            bias2_true = self.fcn_bias2_true_STAR1(h);
        elif self.dgp_type == 2:
            bias2_true = self.fcn_bias2_true_STAR1(h);
        elif self.dgp_type == 3:
            bias2_true = self.fcn_bias2_true_ARMA11(h);
        return bias2_true;

    def estimation(self, LB, UB, optns = {'disp': True}):
        self.LB = LB;
        self.UB = UB;
        self.theta_init = (self.UB - self.LB) * np.random.randn(len(self.LB)) + self.LB;
        if self.dgp_type in [1, 2]:
            temp = np.sort(self.Y);
            ind = int(np.floor(.15*len(temp)));
            self.LB[2] = temp[ind];
            ind = int(np.floor(.85*len(temp)));
            self.UB[2] = temp[ind];
            bnds = tuple(zip(LB,UB));
            cons = ({'type': 'ineq', 'fun' : lambda x: np.array(1 - x[0] - x[1]), 'jac' : lambda x: np.array([-1.0, -1.0, 0.0])}, {'type': 'ineq', 'fun' : lambda x: np.array(-1 + x[0] + x[1]), 'jac' : lambda x: np.array([1.0, 1.0, 0.0])});
            soln = minimize(self.loss_fcn_STAR1, self.theta_init, jac = self.dloss_fcn_STAR1, method='SLSQP', bounds=bnds, constraints=cons, options = optns);
        elif self.dgp_type == 3:
            bnds = tuple(zip(LB,UB));
            cons = ({'type': 'ineq', 'fun' : lambda x: np.array(1 - x[0] - x[2]), 'jac' : lambda x: np.array([-1.0, 0.0, -1.0])}, {'type': 'ineq', 'fun' : lambda x: np.array(-1 + x[0] + x[2]), 'jac' : lambda x: np.array([1.0, 0.0, 1.0])});
            soln = minimize(self.loss_fcn_ARMA11, self.theta_init, jac = self.dloss_fcn_ARMA11, method='SLSQP', bounds=bnds, constraints=cons, options = optns);
        ####
        self.theta_hat = soln.x;
        self.Yhat = self.fcn_Yhat(self.theta_hat);
        self.ehat = self.fcn_e_hat(self.theta_hat);
        self.sigma2_hat = (1/self.Te) * sum(self.ehat * self.ehat);
        self.sigma_hat = np.sqrt(self.sigma2_hat); 
        ####
        self.beta_hat = self.theta_hat[0:self.num_params[0]];
        self.zeta_hat = self.theta_hat[(self.num_params[0]):sum(self.num_params[0:2])];
        self.pi_hat = self.theta_hat[(sum(self.num_params[0:2])):sum(self.num_params[0:3])];
        
    ## STAR(1)
    def dgp_Y_STAR1(self):
        # yt = Gt(pi, Xpi) * Xt' * B + Zt' * xi + et
        self.beta0 = self.theta0[0:self.num_params[0]];
        self.zeta0 = self.theta0[(self.num_params[0]):sum(self.num_params[0:2])];
        self.pi0 = self.theta0[(sum(self.num_params[0:2])):sum(self.num_params[0:3])];
        self.mu0 = self.theta0[(sum(self.num_params[0:3])):sum(self.num_params[0:4])];
        Y0 = np.zeros(self.init_T);
        for t in range((self.num_lags),self.init_T):
            XB  = np.dot(self.beta0, Y0[t-self.num_params[0]:t]);
            Xpi = Y0[t-self.num_params[2]:t];
            Gt  = 1 / (1 + np.exp(- np.dot(self.mu0 , (Xpi-self.pi0) )));
            Zzeta = np.dot(self.zeta0 , Y0[t-self.num_params[1]:t]);
            Y0[t] = (XB * Gt) + Zzeta + self.U[t];
        self.Y = Y0[len(Y0)-self.T:len(Y0)]; # remove burn-in values

    def loss_fcn_STAR1(self, theta_in):
        ehat_temp = self.fcn_e_hat(theta_in);
        Q = (1/(self.Te)) * sum((ehat_temp**2));
        return Q;
 
    def dloss_fcn_STAR1(self, theta_in): #later
        temp_ehat = self.fcn_e_hat_STAR1(theta_in);
        temp_d_theta = self.d_theta_fcn_STAR1(theta_in);
        p = sum(self.num_params[0:3]);
        m_theta_t = np.zeros((p, self.Te));
        for t in range(self.Te):
            m_theta_t[:,t] = temp_ehat[t] * temp_d_theta[:,t];
        dQ = np.sum(m_theta_t,1) / self.Te;
        #####
        #dQ = np.array((0,0,0));
        return dQ;

    def fcn_e_hat_STAR1(self, theta_in):
        temp_Yhat = self.fcn_Yhat_STAR1(theta_in);
        ehat = self.Y - temp_Yhat;
        ehat = ehat[self.num_lags:self.T];
        return ehat;
    
    def fcn_Yhat_STAR1(self, theta_in):
        beta = theta_in[0:self.num_params[0]];
        zeta = theta_in[(self.num_params[0]):sum(self.num_params[0:2])];
        pi = theta_in[(sum(self.num_params[0:2])):sum(self.num_params[0:3])];
        mu = self.theta0[(sum(self.num_params[0:3])):sum(self.num_params[0:4])];
        ####
        Yhat = np.zeros(self.T);
        Yhat[0:self.num_lags] = self.Y[0:self.num_lags];
        for t in range((self.num_lags),self.T):
            XB  = np.dot(beta, self.Y[t-self.num_params[0]:t]);
            Xpi = self.Y[t-self.num_params[2]:t];
            Gt  = 1 / (1 + np.exp(- np.dot(mu , (Xpi - pi) )));
            Zzeta = np.dot(zeta , self.Y[t-self.num_params[1]:t]);
            Yhat[t] = (XB * Gt) + Zzeta;
        return Yhat;

    #Bias2 term
    def fcn_bias2_true_STAR1(self, h):
        theta_in = self.theta0[0:sum(self.num_params[0:3])];
        theta_in[0] = 0;
        e_hat_0n = self.fcn_e_hat_STAR1(theta_in);
        e = self.U[len(self.U)-self.Te:len(self.U)];
        bias2_true = (e_hat_0n[0:len(e_hat_0n)-h] * e_hat_0n[1+h:len(e_hat_0n)]) - (e[1:len(e)-h] * e[1+h:len(e)]);
        bias2_true = (1/np.sqrt(self.Te)) * sum(bias2_true);
        return bias2_true;
        
    def fcn_bias2_STAR1(self, b_in, pi0_in, h):
        g_t = np.zeros(self.Te);
        for t in range(self.num_lags, self.Te):
            z_t = self.Y[t-self.num_params[2]:t];
            g_t[t] = 1 / (1 + np.exp(- self.mu0 * (z_t - pi0_in)));
        temp = sum(self.ehat[1:len(self.ehat)-h] * self.ehat[1:len(self.ehat)-h] * g_t[1+h:len(g_t)]) / self.Te;
        bias2 = b_in * (self.zeta_hat**(h-1)) * temp;
        return bias2;

    # Weak ID
    def d_psi_fcn_STAR1(self, pi_in):
        #d_psi,t = -[x_t' g(z_t, \pi), x_t']'
        pi = pi_in;
        mu = self.mu0;
        ####
        d_psi = np.zeros((sum(self.num_params[0:2]), self.T));
        for t in range(self.num_lags, self.T):
            x_t   = self.Y[t-self.num_params[0]-1:t-1];
            z_t   = self.Y[t-self.num_params[2]-1:t-1];
            g_t   = 1 / (1 + np.exp(- mu * (z_t - pi)));
            d_psi[:,t] = - np.array((x_t * g_t, x_t));
        d_psi = d_psi[:,self.num_lags:self.T];
        return d_psi;

    def bs_GHK_fcn_STAR1(self, theta_in, pi0_in):
        k_psi = sum(self.num_params[0:2]);
        pi_in = theta_in[(k_psi+1):len(theta_in)];
        S_beta = np.vstack((np.eye(self.num_params[0]), np.zeros((self.num_params[1],self.num_params[0]))));
        d_beta_at_pi0 = np.dot(np.transpose(S_beta), self.d_psi_fcn(pi0_in)); # d_beta x T
        d_psi = self.d_psi_fcn(pi_in); # d_psi x T
        H_n = (1/self.Te) * np.dot(d_psi, np.transpose(d_psi)); # d_psi x d_psi              
        H_inv = np.linalg.inv(H_n);
        K_n = -(1/self.Te) * np.dot(d_psi, np.transpose(d_beta_at_pi0)); # d_psi x d_beta  
        e_hat = self.fcn_e_hat(theta_in);
        m_psi_t = np.zeros((k_psi, self.Te));
        for t in range(self.Te):
            m_psi_t[:,t] = d_psi[:,t] * np.tile(e_hat[t],(k_psi,1));
        mHK = { 'm_psi_t':m_psi_t, 'H_inv':H_inv, 'K_n':K_n };
        return mHK;

    # Strong ID
    def d_theta_fcn_STAR1(self, theta_in = None):
        #d_theta,t = -[x_t' g(z_t, \pi), x_t', (\beta' x_t g_{\pi}(z_t, \pi))']'
        if theta_in is None:
        #beta = self.beta_hat;
            pi = self.pi_hat;
            mu = self.mu0;
        else:
            beta = theta_in[0:self.num_params[0]];
            #zeta = theta_in[(self.num_params[0]):sum(self.num_params[0:2])];
            pi = theta_in[(sum(self.num_params[0:2])):sum(self.num_params[0:3])];
            mu = self.theta0[(sum(self.num_params[0:3])):sum(self.num_params[0:4])];
        ####
        d_theta = np.zeros((sum(self.num_params[0:3]), self.T));
        for t in range(self.num_lags,self.T):
            x_t   = self.Y[t-self.num_params[0]:t];
            z_t   = self.Y[t-self.num_params[2]:t];
            g_t   = 1 / (1 + np.exp(- mu * (z_t - pi)));
            g_pit = -(1 + np.exp(- mu * (z_t - pi)))**(-2) * np.exp(- mu * (z_t - pi)) * (mu);
            if np.isnan(g_pit) == 1:
                g_pit = 0;
            # (1 + exp(- mu' * (z_t - pi)))^{-1}
            # -(1 + exp(- mu' * (z_t - pi)))^{-2} * exp(- mu' * (z_t - pi)) * (--mu)
            # -(1 + exp(- mu' * (z_t - pi)))^{-2} * exp(- mu' * (z_t - pi)) * (mu)
            if theta_in is None:
                d_theta[:,t] = - np.hstack( (np.dot(x_t, g_t), x_t, np.dot(x_t, g_pit) ) );  #,(beta' * x_t * g_pit)'
            else:
                d_theta[:,t] = - np.hstack( ( np.dot(x_t, g_t), x_t, np.dot(beta,np.dot(x_t, g_pit)) ) );  #,(beta' * x_t * g_pit)'
        d_theta = d_theta[:,self.num_lags:self.T];
        return d_theta;    
        
    def bs_GJ_theta_STAR1(self, d_theta):
        ehat_temp = self.ehat;
        J_theta_n = (1/self.Te) * np.dot(d_theta, np.transpose(d_theta));
        J_theta_n_inv = np.linalg.inv(J_theta_n);
        p = sum(self.num_params[0:3]);
        m_theta_t = np.zeros((p, self.Te));
        for t in range(self.Te):
            m_theta_t[:,t] = ehat_temp[t] * d_theta[:,t];
        mJ = { 'm_theta_t':m_theta_t, 'J_theta_n_inv':J_theta_n_inv};
        return mJ;

    #### STAR2
    def dgp_Y_STAR2(self):
        # yt = Gt(pi, Xpi) * Xt' * B + Zt' * xi + et
        self.beta0 = self.theta0[0:self.num_params[0]];
        self.zeta0 = self.theta0[(self.num_params[0]):sum(self.num_params[0:2])];
        self.pi0   = self.theta0[(sum(self.num_params[0:2])):sum(self.num_params[0:3])];
        self.mu0   = self.theta0[(sum(self.num_params[0:3])):sum(self.num_params[0:4])];
        zeta_temp  = np.append(self.zeta0, .15 / np.sqrt(self.Te));
        ###
        Y0 = np.zeros(self.init_T);
        for t in range( max(self.num_lags,len(zeta_temp)), self.init_T):
            XB  = np.dot(self.beta0, Y0[t-self.num_params[0]:t]);
            Xpi = Y0[t-self.num_params[2]:t];
            Gt  = 1 / (1 + np.exp(- np.dot(self.mu0 , (Xpi-self.pi0) )));
            Zzeta = np.dot(self.zeta0 , Y0[t-self.num_params[1]:t]);
            Y0[t] = (XB * Gt) + Zzeta + self.U[t];
        self.Y = Y0[len(Y0)-self.T:len(Y0)]; # remove burn-in values
                
    #### ARMA11
    def loss_fcn_ARMA11(self, theta_in):
        beta = theta_in[0:self.num_params[0]];
        zeta = theta_in[(self.num_params[0]):sum(self.num_params[0:2])];
        pi   = theta_in[(sum(self.num_params[0:2])):sum(self.num_params[0:3])];
        Q_out = np.zeros(self.T);
        lag_num = 1;
        for t in range((lag_num),self.T):
            ind = range(t-lag_num, 1, lag_num);
            pi_vec = np.zeros(len(ind));
            for j in range(len(ind)):
                pi_vec[j] = pi**(j-1);
            temp = np.dot(pi_vec, self.Y[ind]);
            Q_out[t] = (1/2) * np.log(zeta) + (1/(2*zeta)) * (self.Y[t] - beta * temp)**2;
            del(temp, ind);
        Q = (1/(self.T-lag_num)) * sum(Q_out);
        return Q;

    def dloss_fcn_ARMA11(self, theta_in):
        rho = self.fcn_DQ_ARMA11(theta_in, 1);
        rho_b = np.mean(rho['rho_beta']);
        rho_z = np.mean(rho['rho_zeta']);
        rho_p = np.mean(rho['rho_pi']);
        dQ = np.vstack((rho_b, rho_z, rho_p));
        return dQ;
        
    def fcn_e_hat_ARMA11(self, theta_in):
        beta = theta_in[0:self.num_params[0]];
        zeta = theta_in[(self.num_params[0]):sum(self.num_params[0:2])];
        pi   = theta_in[(sum(self.num_params[0:2])):sum(self.num_params[0:3])];
        ehat = np.zeros(self.T);
        lag_num = 1;
        for t in range((lag_num),self.T):
            ind = range(t-lag_num, 1, lag_num);
            pi_vec = np.zeros(len(ind));
            for j in range(len(ind)):
                pi_vec[j] = pi**(j-1);
            temp = np.dot(pi_vec, self.Y[ind]);
            ehat[t] = (zeta**(-1/2)) * (self.Y[t] - np.dot(beta, temp));
            del(temp, ind);
        ehat = ehat[self.num_lags:self.T];
        return ehat;
        
    def fcn_Yhat_ARMA11(self, theta_in):
        beta = theta_in[0:self.num_params[0]];
        zeta = theta_in[(self.num_params[0]):sum(self.num_params[0:2])];
        pi   = theta_in[(sum(self.num_params[0:2])):sum(self.num_params[0:3])];
        temp_ehat_0 = self.fcn_e_hat_ARMA11(theta_in);
        temp_ehat = np.zeros(self.T);
        temp_ehat[1:self.T] = temp_ehat_0;
        Yhat_out = np.zeros(self.T);
        Yhat_out[1] = self.Y[1];
        for t in range(1,self.T):
            Yhat_out[t] = (pi + beta) * Yhat_out[t-1] + (zeta**(1/2)) * temp_ehat[t] - pi * (zeta**(1/2)) * temp_ehat[t-1];
        Yhat_out = Yhat_out[self.num_lags:self.T];
        return Yhat_out;

    def dgp_Y_ARMA11(self):   #Note: Need to make sure that obj.U has var = 1;
        self.beta0 = self.theta0[0:self.num_params[0]];
        self.zeta0 = self.theta0[(self.num_params[0]):sum(self.num_params[0:2])];
        self.pi0   = self.theta0[(sum(self.num_params[0:2])):sum(self.num_params[0:3])];
        Y0 = np.zeros(self.init_T);
        e = (self.zeta0**(1/2)) * self.U;
        for t in range(1,self.init_T):
            Y0[t] = (self.pi0 + self.beta0) * Y0[t-1] + e[t] - self.pi0 * e[t-1];
        self.Y = Y0[(self.init_T-self.T):self.init_T];
 
        #### Now the derivatives and other objects that will be used in the max corr test
        #### See AC 2012 page 50 for these quantities
    def fcn_D_eps_ARMA11(self, theta_in):
        beta = theta_in[0:self.num_params[0]];
        zeta = theta_in[(self.num_params[0]):sum(self.num_params[0:2])];
        pi = theta_in[(sum(self.num_params[0:2])):sum(self.num_params[0:3])];
        #B = [ones(obj.num_params(1),1); ones(obj.num_params(2),1); obj.beta_hat];
        #B = diag(B);
        #B_inv = B \ eye(size(B,1));
        d_b = np.zeros(self.T);
        d_z = np.zeros(self.T);
        d_p = np.zeros(self.T);
        lag_num = 1;
        for t in range(lag_num, self.T):
            ind = range(t-lag_num, 1, lag_num);
            pi_vec = np.zeros(len(ind));
            pi_vec_d = np.zeros(len(ind));
            for j in range(len(ind)):
                pi_vec[j] = pi**(j-1);
                pi_vec_d[j] = (j-1) * pi**(j-2);
            temp = np.dot(pi_vec, self.Y[ind]);
            temp2 = np.dot(pi_vec_d, self.Y[ind]);
            d_b[t] = -zeta**(-1/2) * temp;
            d_z[t] = -(1/2) * zeta**(-3/2) * (self.Y[t] - beta * temp);
            d_p[t] = -zeta**(-1/2) * temp2; # -zeta**(-1/2) * beta * temp2;
            del(temp, temp2, ind);
        d_b = d_b[lag_num:self.T];
        d_z = d_z[lag_num:self.T];
        d_p = d_p[lag_num:self.T];
        derivatives = { 'd_b':d_b, 'd_z':d_z, 'd_p':d_p};
        return derivatives;

    def fcn_DQ_ARMA11(self, theta_in, loss_fcn_flag = None):
        beta = theta_in[0:self.num_params[0]];
        zeta = theta_in[(self.num_params[0]):sum(self.num_params[0:2])];
        pi = theta_in[(sum(self.num_params[0:2])):sum(self.num_params[0:3])];
        rho_beta = np.zeros(self.T);
        rho_zeta = np.zeros(self.T);
        rho_pi = np.zeros(self.T);
        lag_num = 1;
        for t in range(lag_num, self.T):
            ind = range(t-lag_num, 1, lag_num);
            pi_vec = np.zeros(len(ind));
            pi_vec_d = np.zeros(len(ind));
            for j in range(len(ind)):
                pi_vec[j] = pi**(j-1);
                pi_vec_d[j] = (j-1) * pi**(j-2);
            temp = np.dot(pi_vec, self.Y[ind]);
            temp2 = np.dot(pi_vec_d, self.Y[ind]);
            rho_beta[t] = -zeta**(-1) * (self.Y[t] - beta * temp) * temp;
            rho_zeta[t] = -(1/2) * zeta**(-2) * ((self.Y[t] - beta * temp)**2 - zeta);
            if loss_fcn_flag is None:
                rho_pi[t] = -zeta**(-1) * (self.Y[t] - beta * temp) * temp2; 
            else:
                rho_pi[t] = -zeta**(-1) * (self.Y[t] - beta * temp) * beta * temp2; 
            del(temp, temp2, ind);
        rho_beta = rho_beta[lag_num:self.T];
        rho_zeta = rho_zeta[lag_num:self.T];
        rho_pi = rho_pi[lag_num:self.T];
        rho = { 'rho_beta':rho_beta, 'rho_zeta':rho_zeta, 'rho_pi':rho_pi};
        return rho;

    def fcn_D2Q_ARMA11(self, theta_in):
        beta = theta_in[0:self.num_params[0]];
        zeta = theta_in[(self.num_params[0]):sum(self.num_params[0:2])];
        pi = theta_in[(sum(self.num_params[0:2])):sum(self.num_params[0:3])];
        rho_bb = np.zeros(self.T);
        rho_bz = np.zeros(self.T);
        rho_zz = np.zeros(self.T);
        rho_bp = np.zeros(self.T);
        rho_zp = np.zeros(self.T);
        rho_pp = np.zeros(self.T);
        lag_num = 1;
        for t in range(lag_num, self.T):
            ind = range(t-lag_num, 1, lag_num);
            pi_vec = np.zeros(len(ind));
            pi_vec_d = np.zeros(len(ind));
            pi_vec_d2 = np.zeros(len(ind));
            for j in range(len(ind)):
                pi_vec[j] = pi**(j-1);
                pi_vec_d[j] = (j-1) * pi**(j-2);
                pi_vec_d2[j] = (j-1) * (j-2) * pi**(j-3);
            temp = np.dot(pi_vec, self.Y[ind]);
            temp2 = np.dot(pi_vec_d, self.Y[ind]);
            temp3 = np.dot(pi_vec_d2, self.Y[ind]);
            rho_bb[t] = zeta**(-1) * temp * temp;
            rho_bz[t] = zeta**(-2) * (self.Y[t] * temp);
            rho_zz[t] = -(1/2) * zeta**(-2) + zeta**(-3) * self.Y[t]**2;
            rho_bp[t] = zeta**(-1) * temp2 * (2 * beta * temp - self.Y[t]);
            rho_zp[t] = zeta**(-2) * (self.Y[t] - beta * temp) * beta * temp2;
            rho_pp[t] = zeta**(-1) * (temp2)**2; # zeta**(-1) * ((beta * temp2)**2 - ((self.Y[t] - beta * temp) * beta * temp3));
            del(temp, temp2, temp3, ind);
        rho_bb = rho_bb[lag_num:self.T];
        rho_bz = rho_bz[lag_num:self.T];
        rho_zz = rho_zz[lag_num:self.T];
        rho_bp = rho_bp[lag_num:self.T];
        rho_zp = rho_zp[lag_num:self.T];
        rho_pp = rho_pp[lag_num:self.T];
        rho = { 'rho_bb':rho_bb, 'rho_bz':rho_bz, 'rho_zz':rho_zz, 'rho_bp':rho_bp, 'rho_zp':rho_zp, 'rho_pp':rho_pp};
        return rho;

    #Weak ID
    def fcn_bias2_true_ARMA11(self, h):
        theta_in = self.theta0[0:sum(self.num_params[0:3])];
        theta_in[0] = 0;
        e_hat_0n = self.fcn_e_hat_ARMA11(theta_in);
        e = self.U[len(self.U)-self.Te:len(self.U)];
        bias2_test = (e_hat_0n[0:(self.Te-h)] * e_hat_0n[h:self.Te]) - (e[0:(self.Te-h)] * e[h:self.Te]);
        bias2_test = (1/(self.Te**(1/2))) * sum(bias2_test);
        return bias2_test;

    def fcn_bias2_ARMA11(self, b_in, pi0_in, h):
        bias2 = b_in * (pi0_in**(h-1)) * self.zeta_hat;
        return bias2;

    def d_psi_fcn_ARMA11(self, pi_in):
        theta_in = np.vstack((self.beta_hat, self.zeta_hat, pi_in));
        deriv = self.fcn_D_eps_ARMA11(theta_in);
        d_b = deriv['d_b'];
        d_z = deriv['d_z'];
        d_psi = np.vstack((d_b, d_z));
        return d_psi;

    def bs_GHK_fcn_ARMA11(self, theta_in, pi0_in):
        pi = theta_in[(sum(self.num_params[0:2])):sum(self.num_params[0:3])];
        rho = self.fcn_DQ_ARMA11(theta_in);
        rho_beta = rho['rho_beta'];
        rho_zeta = rho['rho_zeta'];
        d_rho = self.fcn_D2Q_ARMA11(theta_in);
        rho_bb = d_rho['rho_bb'];
        rho_bz = d_rho['rho_bz'];
        rho_zz = d_rho['rho_zz'];
        m_psi_t = np.vstack((rho_beta, rho_zeta));
        r_bb = np.mean(rho_bb); 
        r_bz = np.mean(rho_bz); 
        r_zz = np.mean(rho_zz);
        H_n = np.vstack(( np.hstack((r_bb, r_bz)), np.hstack((r_bz, r_zz)) ));
        H_inv = np.linalg.inv(H_n);
        K_n = np.vstack(( (-(1 - (pi0_in * pi))**(-1)), 0 )); # Add in estimator later
        mHK = { 'm_psi_t':m_psi_t, 'H_inv':H_inv, 'K_n':K_n };
        return mHK;

    # Strong ID
    def d_theta_fcn_ARMA11(self):
        deriv = self.fcn_D_eps_ARMA11(self.theta_hat);
        d_b = deriv['d_b'];
        d_z = deriv['d_z'];
        d_p = deriv['d_p'];
        d_theta = np.vstack((d_b, d_z, d_p));
        return d_theta;

    def bs_GJ_theta_ARMA11(self, d_theta):
        theta_in = self.theta_hat;
        rho = self.fcn_DQ_ARMA11(theta_in);
        rho_beta = rho['rho_beta'];
        rho_zeta = rho['rho_zeta'];
        rho_pi = rho['rho_pi'];
        d_rho = self.fcn_D2Q_ARMA11(theta_in);
        rho_bb = np.mean(d_rho['rho_bb']);
        rho_bz = np.mean(d_rho['rho_bz']);
        rho_zz = np.mean(d_rho['rho_zz']);
        rho_bp = np.mean(d_rho['rho_bp']);
        rho_zp = np.mean(d_rho['rho_zp']);
        rho_pp = np.mean(d_rho['rho_pp']);
        m_theta_t = np.vstack(( rho_beta, rho_zeta, rho_pi ));
        J_theta_n = np.vstack(( np.hstack(( rho_bb, rho_bz, rho_bp )), np.hstack(( rho_bz, rho_zz, rho_zp )), np.hstack(( rho_bp, rho_zp, rho_pp )) ));
        J_theta_n_inv = np.linalg.inv(J_theta_n);
        mJ = { 'm_theta_t':m_theta_t, 'J_theta_n_inv':J_theta_n_inv };
        return mJ;
        