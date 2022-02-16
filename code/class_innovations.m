classdef class_innovations
    properties (SetAccess=public, GetAccess=public)
	 T
     init_T
	 seed
     U
     innovation_type_string
     innovation_type
     rng_seed
    end
    methods
       function obj = class_innovations(innovation_type, T, seed_input)
            obj.innovation_type = innovation_type;
            obj.T = T+1;
            obj.seed = seed_input;
            rng(obj.seed); %Note: matlab doc says not to seed when
            % running in parallel, bc default behavior is to make sims
            % independent across workers.  Seeding 'shuffle' will result in
            % repeated seeds.  When using array jobs on LongLeaf, control
            % the seed by passing '$SLURM_ARRAY_TASK_ID' as the seed_input.
            temp = rng;
            obj.rng_seed = temp.Seed;
            clear temp;
            % Now generate the data
            obj.init_T = obj.T + obj.T;
            eps = randn(obj.init_T,1);
            switch innovation_type
                case 1
                    obj.U = obj.dgp_iid(eps);
                    obj.innovation_type_string = 'iid';
                case 2
                    obj.U = obj.dgp_GARCH(eps);
                    obj.innovation_type_string = 'GARCH(1,1)';
                case 3
                    obj.U = obj.dgp_bilinear(eps);
                    obj.innovation_type_string = 'Bilinear';
                case 4
                    p=2;
                    obj.U = obj.dgp_ARp(eps,p);
                    obj.innovation_type_string = sprintf('AR(%d)',p);
                case 5
                    q=1;
                    obj.U = obj.dgp_MAq(eps,q);
                    obj.innovation_type_string = sprintf('MA(%d)',q);
                case 6
                    q=10;
                    obj.U = obj.dgp_MAq(eps,q);
                    obj.innovation_type_string = sprintf('MA(%d)',q);
                case 7
                    q=21;
                    obj.U = obj.dgp_MAq(eps,q);
                    obj.innovation_type_string = sprintf('MA(%d)',q);
                case 8
                    q=50;
                    obj.U = obj.dgp_MAq(eps,q);
                    obj.innovation_type_string = sprintf('MA(%d)',q);
                case 9
                    q=100;
                    obj.U = obj.dgp_MAq(eps,q);
                    obj.innovation_type_string = sprintf('MA(%d)',q);
            end
            %obj.U = obj.U((obj.init_T - obj.T + 1):end); % remove burn-in values
        end
        function y = dgp_iid(obj, eps)
            y = eps;
        end
        function y = dgp_GARCH(obj, eps, beta)
            beta = [.3, .6];
            sigma2 = ones(obj.init_T,1);
            y = zeros(obj.init_T,1);
            y(1) = sqrt(sigma2(1)) * eps(1);
            for t = 2:obj.init_T
                sigma2(t) = 1 + beta(1) * y(t-1)^2 + beta(2) * sigma2(t-1);
                y(t) = sqrt(sigma2(t)) * eps(t);
            end
        end
        function y = dgp_bilinear(obj, eps, beta)
            beta = .5;
            y = zeros(obj.init_T,1);
            y(1) = eps(1);
            for t = 2:obj.init_T
                y(t) = beta * y(t-1) * eps(t-1) + eps(t);
            end
        end
        function y = dgp_ARp(obj, eps, p, beta)
            beta = .5;
            y = zeros(obj.init_T,1);
            for t = (p+1):obj.init_T
                y(t) = beta * y(t-p) + eps(t);
            end
        end
        function y = dgp_MAq(obj, eps, q, beta)
            beta = .5;
            y = zeros(obj.init_T,1);
            for t = (q+1):obj.init_T
                y(t) = beta * eps(t-q) + eps(t);
            end
        end
    end
end
