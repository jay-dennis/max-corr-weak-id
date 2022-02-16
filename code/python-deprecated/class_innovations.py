import numpy as np

class class_innovations:
    """This class instantiates the innovations"""

    def __init__(self, innovation_type, T, seed_input = 'shuffle'):
        T = T + 1;
        init_T = T + T;
        self.innovation_type = innovation_type;
        self.T = T;
        self.seed = seed_input;
        if type(seed_input).__name__ == 'int':
            np.random.seed(self.seed); #, version=2)
        elif type(seed_input).__name__ != 'int':
            np.random.seed(); #(version=2)
        self.rng_seed = np.random.get_state(); # for reproducibility
        # Now generate the data
        self.init_T = init_T;
        eps = np.random.randn(self.init_T)
        if innovation_type == 1:
            self.U = self.dgp_iid(eps);
            self.innovation_type_string = 'iid';
        elif innovation_type == 2:
            self.U = self.dgp_GARCH(eps);
            self.innovation_type_string = 'GARCH(1,1)';
        elif innovation_type == 3:
            self.U = self.dgp_bilinear(eps);
            self.innovation_type_string = 'Bilinear';
        elif innovation_type == 4:
            p=2;
            self.U = self.dgp_ARp(eps,p);
            self.innovation_type_string = 'AR(%d)' % (p);
        elif innovation_type == 5:
            q=1;
            self.U = self.dgp_MAq(eps,q);
            self.innovation_type_string = 'MA(%d)' % (q);
        elif innovation_type == 6:
            q=10;
            self.U = self.dgp_MAq(eps,q);
            self.innovation_type_string = 'MA(%d)' % (q);
        elif innovation_type == 7:
            q=21;
            self.U = self.dgp_MAq(eps,q);
            self.innovation_type_string = 'MA(%d)' % (q);
        elif innovation_type > 7:
            print('Error.  Enter a valid innovation type.');

    def dgp_iid(self, eps):
        return eps;
        
    def dgp_GARCH(self, eps, beta = [.3, .6]):
        y = np.zeros(self.init_T);
        sigma2 = np.ones(self.init_T);
        y[0] = (sigma2[0]**(.5)) * eps[0];
        for t in range(1, self.init_T):
            sigma2[t] = 1 + beta[0] * y[t-1]**2 + beta[1] * sigma2[t-1];
            y[t] = (sigma2[t]**.5) * eps[t];
        return y;
            
    def dgp_bilinear(self, eps, beta = .5):
        y = np.zeros(self.init_T);
        y[1] = eps[1];
        for t in range(1, self.init_T):
            y[t] = beta * y[t-1] * eps[t-1] + eps[t];
        return y;
        
    def dgp_ARp(self, eps, p = 1, beta = .5):
        y = np.zeros(self.init_T);
        for t in range(p, self.init_T):
            y[t] = beta * y[t-p] + eps[t];
        return y;
    
    def dgp_MAq(self, eps, q = 1, beta = .5):
        y = np.zeros(self.init_T);
        for t in range(q, self.init_T):
            y[t] = beta * eps[t-q] + eps[t];
        return y;
