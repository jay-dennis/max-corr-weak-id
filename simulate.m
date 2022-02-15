function simulate(dgp_type, id_type, innovation_type, feasible_flag, J, T, seed)

% simulation
%clear all;
%clc;

%% Setup

temp_b = .3;  % value of b kept outside of all loops for re-initialization

warning('off','all');

% dgp_type_vec = (1:3);
% T_vec = [100 250 500 1000];
% id_type_vec = [0 1 2];
% innovation_type_vec = (1:9);
if nargin == 0
    dgp_type = 1;
    id_type = 0;
    innovation_type = 6;
    feasible_flag = 0;
    J = 2;
    T = 100;
    seed = 'shuffle';
elseif nargin < 7
    seed = 'shuffle';    
end

rng(seed);
temp = rng;
sim_number = temp.Seed;
clear temp;


b=temp_b;
if id_type == 1  % Weak id
    beta_in = (b/sqrt(T));
elseif id_type == 2 % Strong id
    beta_in = b;
elseif id_type == 0 % No id
    beta_in = 0; b = 0;
end

if dgp_type == 1
    zeta_in = .5;
    pi_in = 0;
    mu_in = 10;
    num_params = [1 1 1 1]; % For vector beta, zeta, etc.
    LB = [-.95 -.95 -2]; % For Parameter Estimation (mu is excluded)
    UB = -LB;
    theta0 = [beta_in zeta_in pi_in mu_in];     % For STAR1 [beta, zeta, pi, mu]
elseif dgp_type == 2
    zeta_in = .5;
    pi_in = 0;
    mu_in = 10;
    num_params = [1 1 1 1]; % For vector beta, zeta, etc.
    LB = [-.95 -.95 -2]; % For Parameter Estimation (mu is excluded)
    UB = -LB;
    theta0 = [beta_in zeta_in pi_in mu_in];      % For STAR2 bc zeta_2 = .85 - zeta_1
elseif dgp_type == 3
    zeta_in = 1;
    pi_in = .5;
    num_params = [1 1 1]; % For vector beta, zeta, etc.
    LB = [-.95 0 -.95]; % For Parameter Estimation
    UB = [.95  2  .95];
    theta0 = [beta_in zeta_in pi_in];      % For STAR2 bc zeta_2 = .85 - zeta_1
end

fprintf('%d   DGP: %d   Innov: %d   ID: %d   T: %d  J: %d  Feasible?=%d \n', sim_number, dgp_type, innovation_type, id_type, T, J, feasible_flag);

%% Procedures
flag_run_tests = 1;
flag_clean_up = 1;

% DGP
clear data time;
%data = class_tests.empty(J,0);
tic;
data = repmat(class_tests(dgp_type, id_type, innovation_type, theta0, num_params, T, seed),1,J);
for j=1:J
    %data = [data class_tests(dgp_type, id_type, innovation_type, theta0, num_params, T, seed)];
    data(j) = class_tests(dgp_type, id_type, innovation_type, theta0, num_params, T, seed);
end
time.dgp = toc/60;

% Estimation
tic;
for j=1:J
    data(j) = data(j).estimation(LB, UB);
end
time.est = toc/60;

if flag_run_tests == 0
tic;
for j=1:J
    data(j) = data(j).fcn_test_stats_only(feasible_flag);
end
time.test_stat = toc/60;
disp(time.test_stat);
end

% Bootstrap for Test Stat    
if flag_run_tests == 1
tic;
for j=1:J
    data(j) = data(j).fcn_run_all_tests(feasible_flag);
end
time.bs = toc/60;
disp(time.bs);
end

if flag_clean_up == 1
for j=1:J
    data(j) = data(j).clean_up();
end
end

%% Save
%output_main_dir = sprintf('G:/Simulation_data/max_test_many_zeros/9July/data_T%d', T);
output_main_dir = sprintf('./data_T%d', T);

outputdir = sprintf('%s/output_dgp%d_innov%d_id%d_T%d_f%d_beta%d_zeta%d', output_main_dir, dgp_type, innovation_type, id_type, T, feasible_flag, floor(1000*beta_in), floor(100*zeta_in));
if exist(outputdir,'dir') == 0
    mkdir(outputdir);
end;
outputname=sprintf('%s/data_dgp%d_innov%d_id%d_T%d_f%d_J%d_%d', outputdir, dgp_type, innovation_type, id_type, T, feasible_flag, J, sim_number);

%display(outputname);

while exist(outputname,'file') == 2
    disp('Error - File Name Exists; iterating forward');
    sim_number = sim_number + 1;
    outputname=sprintf('%s/data_dgp%d_innov%d_id%d_T%d_f%d_J%d_%d', outputdir, dgp_type, innovation_type, id_type, T, feasible_flag, J, sim_number);
end
save(outputname, 'data', 'time');

%    end % innovation_type loop
%end % id_type loop


end



