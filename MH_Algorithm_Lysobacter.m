%Script made by Adline Vouillamoz and Isaline Guex
clear;
close all;

%Save or Not
save_data = 1; %1 if save, 0 otherwise
Name_file = 'Lysobacter_Fitting_v1';

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','48:69', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','118:139', 'Format','auto');
Parameters_Senka_Lag_time = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','142:163', 'Format','auto');
Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 6, 'Range','1:22', 'Format','auto');
Data_Evol_test = Data_Evol;
Time_step = [0 1 3 7 10 21]*24; %Measured time step in hours
tspan = [0, max(Time_step)]; %Time interval in hours
S = height(Data_Evol);

yield_vect = table2array(Parameters_set(1:S,5))'; %Change it according to the desired mu_max
name = string(table2array(Parameters_set(1:S,1)));%string(table2array(table(2:22,3)));

mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8)); %Normal Parameters_Senka_mu_max(:,2:3). Log-normal Parameters_Senka_mu_max(:,7:8)
mu_max_vect = mu_max_dist(:,1);
Res_Percentage = table2array(Parameters_Senka_mu_max(:,4)); %Percentage of resources that can be consumed
Parameters_Senka_Lag_time = table2array(Parameters_Senka_Lag_time(:,7:8)); %table2array(Parameters_Senka_Lag_time(:,2:3));
mean_param_LN = mu_max_dist(:,1);
var_param_LN = mu_max_dist(:,2);

Data_Evol_temp = table2array(Data_Evol(:, 2:end));
Data_Evol_temp = Data_Evol_temp(:,~ismember(mod(1:length(Data_Evol_temp(1,:)), 4), 2));%Data_Evol_temp(:,mod(1:length(Data_Evol_temp(1,:)),4) ~= 0);%80% of the data to train. Hold-out method.
nb_obs = length(Data_Evol_temp(1,:));
nb_time_step = length(Time_step);
nb_rep = nb_obs/nb_time_step;
Measured_Abund = zeros(length(mu_max_vect), nb_time_step, nb_rep); %Number species, number times, number replicates.
for i = 1:nb_rep
    Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
end
mean_y_0 = mean(Measured_Abund(:,1,:), 3);
std_y_0 = std(Measured_Abund(:,1,:), 1, 3);
StackPlot_Meas = Measured_Abund./sum(Measured_Abund);

%Fitted parameters
mu_max_vect_temp = max(lognrnd(mu_max_vect, 0*mu_max_dist(:,2)), 0.001);
yield = max(normrnd(yield_vect, 0.0), 1e-05);
ind_Lyso = 6;
%Load parameters
Name_file_load = 'Death_CF_Model_V2'; %Name of the load data for initial setting
kappa_mat = load(strcat('Data/', Name_file_load, '_Kappa_mat.mat')); kappa_mat = kappa_mat.kappa_mat;
CrossFeed_Mat_Lyso = load(strcat('Data/', Name_file_load, '_CrossFeed_Mat.mat')); CrossFeed_Mat_Lyso = CrossFeed_Mat_Lyso.CrossFeed_Mat_Temp;
Resource_Matrix_Lyso = load(strcat('Data/', Name_file_load, '_Resource_Matrix.mat')); Resource_Matrix_Lyso = Resource_Matrix_Lyso.Resource_Matrix;
Death_Mat = load(strcat('Data/', Name_file_load, '_Pred_Mat.mat')); Death_Mat = Death_Mat.Death_Mat_Temp;%Pred_Mat_Temp
Threshold_CF = load(strcat('Data/', Name_file_load, '_Threshold.mat')); Threshold_CF = Threshold_CF.Threshold_CF;
Threshold_death = load(strcat('Data/', Name_file_load, '_Threshold_Pred.mat')); Threshold_death = Threshold_death.Threshold_death;%Threshold_pred.Threshold_pred;
Lag_time_Cons_Lyso = load(strcat('Data/', Name_file_load, '_Lag_time_Cons.mat')); Lag_time_Cons_Lyso = Lag_time_Cons_Lyso.Lag_time_Cons;
R = load(strcat('Data/', Name_file_load, '_R_mat.mat')); R = R.R;
nb_Res = length(R);
kappa_mat_Lyso = Increased_Mat(kappa_mat, ind_Lyso, S, 4, 0);
kappa_mat_Lyso(ind_Lyso,:) = [2.5e+05 mu_max_vect_temp(6) (mu_max_vect_temp(ind_Lyso)/yield(ind_Lyso) - mu_max_vect_temp(ind_Lyso)) 2.5e+05];
Death_Mat_Lyso = Increased_Mat(Death_Mat, ind_Lyso, S, S, 1); 
Death_Mat_Lyso(ind_Lyso, ind_Lyso) = 0.0214;
CrossFeed_Mat_Lyso = Increased_Mat(CrossFeed_Mat_Lyso, ind_Lyso, S, S, 1);
CrossFeed_Mat_Lyso(ind_Lyso, :) = mean(CrossFeed_Mat_Lyso([17 19], :));
CrossFeed_Mat_Lyso(:, ind_Lyso) = mean(CrossFeed_Mat_Lyso(:,[17 19]), 2);
CrossFeed_Mat_Lyso(ind_Lyso, ind_Lyso) = 0;
Lag_time_Cons_Lyso = Increased_Mat(Lag_time_Cons_Lyso, ind_Lyso, S, nb_Res, 0);
Lag_time_Cons_Lyso(ind_Lyso,:) = lognrnd(Parameters_Senka_Lag_time(ind_Lyso,1), Parameters_Senka_Lag_time(ind_Lyso,2));
x = rand(1, nb_Res);
rand_indice = x > 0.3;
Lag_time_Rand = normrnd(50, 10, 1, nb_Res);
Lag_time_Cons_Lyso(ind_Lyso, rand_indice) = max(Lag_time_Cons_Lyso(ind_Lyso, rand_indice), Lag_time_Rand(rand_indice));
Var_Resource_Matrix = lognrnd(zeros(S,nb_Res), 0*repmat(mu_max_dist(:,2), 1, nb_Res), S, nb_Res);
x = rand(1, nb_Res);
rand_inidice = x > 0.24;
Resource_Matrix_Lyso = Increased_Mat(Resource_Matrix_Lyso, ind_Lyso, S, nb_Res, 0); 
Resource_Matrix_Lyso = Resource_Matrix_Lyso.*Var_Resource_Matrix;
Resource_Matrix_Lyso(ind_Lyso, rand_inidice) = lognrnd(zeros(1, sum(rand_inidice)), mu_max_dist(ind_Lyso,2).*ones(1, sum(rand_inidice)));
Mat_kappa_3_Lyso = kappa_mat_Lyso(:,3).*CrossFeed_Mat_Lyso./kappa_mat_Lyso(:,2);
Lag_time_Pred_Lyso = 0*max(normrnd(150, 0, S, S), 0);
Pred_Mat_Lyso = zeros(S,S);
Prey_num = [3 7 11 13 15 17 20];
Pred_Mat_Lyso(ind_Lyso, Prey_num) = kappa_mat_Lyso(ind_Lyso,2)*ones(1, length(Prey_num));

%Initialization of the model parameters fixed for all replicates 
t_0 = 0; %Time 0
CrossFeed_Mat_Lyso_Temp = CrossFeed_Mat_Lyso; %Temporal predation matrix
Death_Mat_Lyso_Temp = Death_Mat_Lyso;
yield_Pred = 0.2; 

%Setting for Matlab ODE solver
%nb_tot_Species = S*3 + nb_Res;
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.

%Initialization of the colors
colors = distinguishable_colors(S + 6);
colors(7:21,:) = colors(6:20,:);
colors(6,:) = colors(27,:);

num_fig = 1;
%Initial concentrations using a normal distribution
mat_y_0 = mean_y_0;
mat_y_0 = [mat_y_0 zeros(S,1) zeros(S,1)];

%Fit one unique threshold value
Threshold_Pred = max(normrnd(1.0e-04, 0*0.1e-05, S, 1),0);
mat_y_0 = reshape(mat_y_0', 1, []);
mat_y_0 = [mat_y_0 R];
kappa_mat_Lyso_init = kappa_mat_Lyso;    
Mat_kappa_3_Lyso_init = Mat_kappa_3_Lyso;
    
beta = 0.5*1e2; 
num_steps = 500;
nb_loops = 4; %Number of different loops (CF, predation. resource, nb of resources)
     
n = 1; r = 1; l = 1; p = 1; u = 1; %Total energy
max_tot_iter = 15;
liste_nrj_CF = zeros(1, num_steps*max_tot_iter);
Struct_CrossFeed_Mat = cell(num_steps*max_tot_iter, 1); %Creation d'une structure pour sauver les matrices générées
Struct_CrossFeed_Mat{1} = CrossFeed_Mat_Lyso_Temp; CrossFeed_Mat_init = CrossFeed_Mat_Lyso_Temp;
CrossFeed_Mat_Lyso_Temp_cand = CrossFeed_Mat_Lyso_Temp; %Initilization candidat
Threshold_CF_cand = Threshold_CF; Threshold_CF_init = Threshold_CF;
acceptance_ratio_CF = zeros(1, num_steps*max_tot_iter);
liste_nrj_pred = zeros(1, num_steps*max_tot_iter);
Struct_Death_Mat = cell(num_steps*max_tot_iter, 1); %Creation d'une structure pour sauver les matrices générées
Struct_Death_Mat{1} = Death_Mat_Lyso_Temp; Death_Mat_Lyso_init = Death_Mat_Lyso_Temp;
Death_Mat_Lyso_Temp_cand = Death_Mat_Lyso_Temp; %Initilization candidat
Threshold_death_cand = Threshold_death; Threshold_death_init = Threshold_death_cand;
Threshold_Pred_cand = Threshold_Pred; Threshold_Pred_init = Threshold_Pred;
acceptance_ratio_pred = zeros(1, num_steps*max_tot_iter);
liste_nrj_res = zeros(1, num_steps*max_tot_iter);
Struct_Res_Mat = cell(num_steps*max_tot_iter, 1); %Creation d'une structure pour sauver les matrices générées
Struct_Res_Mat{1} = Resource_Matrix_Lyso; Resource_Matrix_init = Resource_Matrix_Lyso;
Resource_Matrix_Lyso_cand = Resource_Matrix_Lyso; %Initilization candidat
Struct_Pred_Mat = cell(num_steps*max_tot_iter, 1);
Struct_Pred_Mat{1} = Pred_Mat_Lyso; Pred_Mat_Lyso_init = Pred_Mat_Lyso;
Pred_Mat_Lyso_cand = Pred_Mat_Lyso;
acceptance_ratio_res = zeros(1, num_steps*max_tot_iter);
liste_nrj_tot = zeros(1, nb_loops*num_steps*max_tot_iter);
acceptance_ratio_tot = zeros(1, nb_loops*num_steps*max_tot_iter);
Lag_time_Cons_Lyso_cand = Lag_time_Cons_Lyso; Lag_time_Cons_Lyso_init = Lag_time_Cons_Lyso;
Lag_time_Pred_Lyso_cand = Lag_time_Pred_Lyso; Lag_time_Pred_Lyso_init = Lag_time_Pred_Lyso;
mat_y_0_new = mat_y_0; 
liste_nrj_nb_res = zeros(1, num_steps*max_tot_iter);
acceptance_ratio_nb_res = zeros(1, num_steps*max_tot_iter);
Evol_nb_Res = zeros(1, num_steps*max_tot_iter);

var_theta = 0.1;
acceptance_ratio_temp = 10000;
liste_nrj_tot_temp = 10000;
N = nb_loops*num_steps*max_tot_iter;
max_val = [10./mu_max_vect_temp 10*ones(S,1) 10./mu_max_vect_temp];%[10*ones(S,1) 10*ones(S,1) 10./mu_max_vect_temp];
Hill_CF_coeff = 10; Hill_Pred_coeff = 10;
Hill_CF_coeff_cand = Hill_CF_coeff; Hill_Pred_coeff_cand = Hill_Pred_coeff;
delta_dist = 1;
mu_Likelihood = reshape(repmat(zeros(S,1), 1, nb_Res), 1, []);
var_Likelihood = reshape(repmat(5*mu_max_dist(:,2), 1, nb_Res), 1, []);%reshape(repmat(max_val(:,3), 1, nb_Res), 1, []);%
theta_res = var_theta*mean(mu_max_dist(:,2));
delta_dist_vect = [];
LT_CF_Death = zeros(S, S); %No lag time for CF and Death. Depent only on the threshold.
Cov_Lt_eye = eye(S*S);
while (liste_nrj_tot_temp > 3.5 || (acceptance_ratio_temp < 0.21 || acceptance_ratio_temp > 0.29) || p < 5000) &&  p < N && num_steps > 0%5.6, 0.22:0.27

    %%%%%%%%%% Cross-feeding loop %%%%%%%%

    Mat_kappa_3_Lyso = kappa_mat_Lyso(:,3).*CrossFeed_Mat_Lyso_Temp_cand./mu_max_vect_temp;%CrossFeed_Mat_Temp_cand./yield' - CrossFeed_Mat_Temp_cand;
    theta_T = 0; %Same threshold as the fitted one
    Cov_Mat = 0.01*eye(S*S);
    for k = 1:num_steps
        sol = ode45(@(t, y) fun_CF_Death_Pred(t, y, kappa_mat_Lyso, CrossFeed_Mat_Lyso_Temp_cand, Mat_kappa_3_Lyso, Resource_Matrix_Lyso, Threshold_CF_cand, Threshold_death, Threshold_Pred, Death_Mat_Lyso_Temp, Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons_Lyso, Lag_time_Pred_Lyso, nb_Res, Hill_CF_coeff_cand, Hill_Pred_coeff), tspan,  mat_y_0, opts_1); %Multiple resource groups
        z_temp = deval(sol, Time_step);
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
        X = X + P;
        nu = 2*k^(-3/4);
        StackPlotTot = X./sum(X);

        CrossFeed_Mat_Lyso_Temp = CrossFeed_Mat_Lyso_Temp./(mu_max_vect_temp);
        CrossFeed_Mat_Lyso_Temp_cand = CrossFeed_Mat_Lyso_Temp_cand./(mu_max_vect_temp);

        [CrossFeed_Mat_Lyso_Temp_cand, CrossFeed_Mat_Lyso_Temp, ~, Threshold_CF, Struct_CrossFeed_Mat, liste_nrj_tot, liste_nrj_CF, acceptance_ratio_tot, acceptance_ratio_CF, p, n, ~, theta_T, ~, nu, Cov_Mat] = ...
            fun_MH_Candidate_Rob(CrossFeed_Mat_Lyso_Temp_cand, CrossFeed_Mat_Lyso_Temp,  Threshold_CF_cand, Threshold_CF, Struct_CrossFeed_Mat, p, n, liste_nrj_tot, liste_nrj_CF,...
            acceptance_ratio_tot, acceptance_ratio_CF, S, S, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, 0, theta_T, 0, nu, Cov_Mat, Cov_Lt_eye, max_val(:,1), nb_rep, LT_CF_Death, LT_CF_Death);
        
        temp = CrossFeed_Mat_Lyso_Temp;
        CrossFeed_Mat_Lyso_Temp_cand(1:1+size(CrossFeed_Mat_Lyso_Temp_cand,1):end) = 0;
        temp(ind_Lyso,:) = CrossFeed_Mat_Lyso_Temp_cand(ind_Lyso,:);
        temp(:,ind_Lyso) = CrossFeed_Mat_Lyso_Temp_cand(:,ind_Lyso);
        CrossFeed_Mat_Lyso_Temp_cand = temp;
        Threshold_CF_cand = Threshold_CF;

        CrossFeed_Mat_Lyso_Temp = CrossFeed_Mat_Lyso_Temp.*(mu_max_vect_temp);
        CrossFeed_Mat_Lyso_Temp_cand = CrossFeed_Mat_Lyso_Temp_cand.*(mu_max_vect_temp);

        Mat_kappa_3_Lyso = kappa_mat_Lyso(:,3).*CrossFeed_Mat_Lyso_Temp_cand./mu_max_vect_temp;%CrossFeed_Mat_Temp_cand./yield' - CrossFeed_Mat_Temp_cand;

        if mod(p, 500) == 0
            disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
        end          
    end
    Cov_Mat_CF = Cov_Mat;
        
    %%%%%%%%%% Death loop %%%%%%%%

    Mat_kappa_3_Lyso = kappa_mat_Lyso(:,3).*CrossFeed_Mat_Lyso_Temp./mu_max_vect_temp;
    theta_T_death = 0;%Same threshold as the fitted one
    Cov_Mat = 0.01*eye(S);
    for k = 1:num_steps
        sol = ode45(@(t, y) fun_CF_Death_Pred(t, y, kappa_mat_Lyso, CrossFeed_Mat_Lyso_Temp, Mat_kappa_3_Lyso, Resource_Matrix_Lyso, Threshold_CF, Threshold_death_cand, Threshold_Pred, Death_Mat_Lyso_Temp_cand, Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons_Lyso, Lag_time_Pred_Lyso, nb_Res, Hill_CF_coeff_cand, Hill_Pred_coeff), tspan,  mat_y_0, opts_1); %Multiple resource groups
        z_temp = deval(sol, Time_step);
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
        X = X + P;
        nu = 2*k^(-3/4); 
        StackPlotTot = X./sum(X);
        
        [Death_Mat_Lyso_Temp_cand, Death_Mat_Lyso_Temp, ~, Threshold_death, Struct_Death_Mat, liste_nrj_tot, liste_nrj_pred, acceptance_ratio_tot, acceptance_ratio_pred, p, r, ~, theta_T_death, ~, nu, Cov_Mat] = ...
            fun_MH_Candidate_Rob_Death(Death_Mat_Lyso_Temp_cand, Death_Mat_Lyso_Temp, Threshold_death_cand, Threshold_death, Struct_Death_Mat, p, r, liste_nrj_tot, liste_nrj_pred,...
            acceptance_ratio_tot, acceptance_ratio_pred, S, S, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, 0, theta_T_death, 0, nu, Cov_Mat, Cov_Lt_eye, max_val(:,2), nb_rep, LT_CF_Death, LT_CF_Death);
 
        temp = Death_Mat_Lyso_Temp;
        temp(ind_Lyso, ind_Lyso) = Death_Mat_Lyso_Temp_cand(ind_Lyso, ind_Lyso);
        Death_Mat_Lyso_Temp_cand = temp;
        Threshold_death_cand = Threshold_death;
        
        ind = eye(S,S);
        Death_Mat_Lyso_Temp_cand(~ind) = 0; %Changement pour diagonale zeros

        if mod(p, 500) == 0
            disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
        end
    end
    Cov_Mat_Death = Cov_Mat;

    %%%%%%%%%% Predation matrix loop %%%%%%%%

    Cov_Mat = 0.01*eye(S*S);
    theta_Pred = 0.01*var_theta*mean(Threshold_CF);
    Cov_Mat_LT = 0.01*eye(S*S);
    for k = 1:num_steps
        sol = ode45(@(t, y) fun_CF_Death_Pred(t, y, kappa_mat_Lyso, CrossFeed_Mat_Lyso_Temp, Mat_kappa_3_Lyso, Resource_Matrix_Lyso, Threshold_CF, Threshold_death, Threshold_Pred_cand, Death_Mat_Lyso_Temp, Pred_Mat_Lyso_cand, yield_Pred, S, Lag_time_Cons_Lyso, Lag_time_Pred_Lyso_cand, nb_Res, Hill_CF_coeff, Hill_Pred_coeff), tspan,  mat_y_0, opts_1); %Multiple resource groups
        z_temp = deval(sol, Time_step);
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
        X = X + P;
        nu = 2*k^(-3/4); 
        StackPlotTot = X./sum(X);
        
        [Pred_Mat_Lyso_cand, Pred_Mat_Lyso, Threshold_Pred_cand, Threshold_Pred, Struct_Res_Mat, liste_nrj_tot, ~, acceptance_ratio_tot, ~, p, u, ~, theta_Pred, ~, nu, Cov_Mat] = ...
            fun_MH_Candidate_Rob(Pred_Mat_Lyso_cand, Pred_Mat_Lyso, Threshold_Pred_cand, Threshold_Pred, Struct_Pred_Mat, p, u, liste_nrj_tot, 0,...
            acceptance_ratio_tot, 0, S, S, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, 0, theta_Pred, 0, nu, Cov_Mat, Cov_Mat_LT, max_val(:,2), nb_rep, LT_CF_Death, LT_CF_Death);
        
        temp = Pred_Mat_Lyso;
        temp(ind_Lyso,:) = Pred_Mat_Lyso_cand(ind_Lyso,:);
        temp(ind_Lyso, ind_Lyso) = 0;
        Pred_Mat_Lyso_cand = temp;
        
        if mod(p, 500) == 0
            disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
        end
    end
    Cov_Mat_Pred = Cov_Mat;

    %%%%%%%%%% Resource matrix loop %%%%%%%%

    Cov_Mat = 0.01*eye(S*nb_Res);
    Cov_Mat_LT = 0.01*eye(S*nb_Res);
    for k = 1:num_steps
        sol = ode45(@(t, y) fun_CF_Death_Pred(t, y, kappa_mat_Lyso, CrossFeed_Mat_Lyso_Temp, Mat_kappa_3_Lyso, Resource_Matrix_Lyso_cand, Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Lyso_Temp, Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons_Lyso_cand, Lag_time_Pred_Lyso, nb_Res, Hill_CF_coeff_cand, Hill_Pred_coeff), tspan,  mat_y_0, opts_1); %Multiple resource groups
        z_temp = deval(sol, Time_step);
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
        X = X + P;
        nu = 2*k^(-3/4); 
        StackPlotTot = X./sum(X);

        [Resource_Matrix_Lyso_cand, Resource_Matrix_Lyso, T_cand, T, Struct_Res_Mat, liste_nrj_tot, liste_nrj_res, acceptance_ratio_tot, acceptance_ratio_res, p, l, ~, ~, ~, nu, Cov_Mat, Cov_Mat_LT, Lag_time_Cons_cand, Lag_time_Cons] = ...
            fun_MH_Candidate_Rob(Resource_Matrix_Lyso_cand, Resource_Matrix_Lyso, zeros(nb_Res,1), zeros(nb_Res,1), Struct_Res_Mat, p, l, liste_nrj_tot, liste_nrj_res,...
            acceptance_ratio_tot, acceptance_ratio_res, S, nb_Res, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, 0, 0, 0, nu, Cov_Mat, Cov_Mat_LT, max_val(:,3), nb_rep, Lag_time_Cons_Lyso_cand, Lag_time_Cons_Lyso);
        
        temp = Resource_Matrix_Lyso;
        temp(ind_Lyso,:) = Resource_Matrix_Lyso_cand(ind_Lyso,:);
        Resource_Matrix_Lyso_cand = temp;
        
        if mod(p, 500) == 0
            disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
        end
    end
    Cov_Mat_Res = Cov_Mat;
    Resource_Matrix_Lyso_cand = Resource_Matrix_Lyso; Lag_time_Cons_Lyso_cand = Lag_time_Cons_Lyso;
    mat_y_0_new = mat_y_0;
    acceptance_ratio_temp = sum(acceptance_ratio_tot)/p;
    liste_nrj_tot_temp = liste_nrj_tot(p - 1);
    delta_dist = 1;
    mu_Likelihood = reshape(repmat(zeros(S,1), 1, nb_Res), 1, []);
    var_Likelihood = reshape(repmat(mu_max_dist(:,2), 1, nb_Res), 1, []);
    Temp = 1/beta; Temp = Temp*(1-u/(3*N)); beta = 1/Temp;
end     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Data to test on half of the remaining data
Time_step = [0 1 3 7 10 21]*24; %Measured time step in hours
tspan = [0, max(Time_step)]; %Time interval in hours
Data_Evol_temp = table2array(Data_Evol_test(:, 2:end)); %table2array(Data_Evol(:, 2:end));
Data_Evol_temp = Data_Evol_temp(:,mod(1:length(Data_Evol_temp(1,:)), 4) == 1);%20% of the data to test
nb_time_step = length(Time_step);
nb_obs = length(Data_Evol_temp(1,:));
nb_rep = nb_obs/nb_time_step;
Measured_Abund = zeros(length(mu_max_vect), nb_time_step, nb_rep); %Number species, number times, number replicates.
for i = 1:nb_rep
    Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
end
StackPlot_Meas = Measured_Abund./sum(Measured_Abund);
gamma_Stack = std(StackPlot_Meas, 1, 3);
Threshold_Surviving = 1e-10;
nb_Surv_Obs = sum(mean(Measured_Abund,3) > Threshold_Surviving);
mean_y_0 = mean(Measured_Abund(:,1,:), 3);%table2array(Data_Evol(1:20, 2));
mat_y_0 = [mean_y_0 zeros(S,1) zeros(S,1)];
mat_y_0 = reshape(mat_y_0', 1, []);
mat_y_0 = [mat_y_0 R];

sol = ode45(@(t, y) fun_CF_Death_Pred(t, y, kappa_mat_Lyso, CrossFeed_Mat_Lyso_Temp, Mat_kappa_3_Lyso, Resource_Matrix_Lyso, Threshold_CF, Threshold_death, Threshold_Pred, Death_Mat_Lyso_Temp, Pred_Mat_Lyso, yield_Pred, S, Lag_time_Cons_Lyso, Lag_time_Pred_Lyso, nb_Res, Hill_CF_coeff, Hill_Pred_coeff), tspan,  mat_y_0, opts_1);
z_temp = deval(sol, Time_step);
X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
X = X + P;

CrossFeed_Mat_Temp_rates = CrossFeed_Mat_Lyso_Temp.*mu_max_vect_temp; Resource_Matrix_rates = Resource_Matrix_Lyso.*mu_max_vect_temp;
        
%%% Figure generation %%%

acceptance_ratio_fin = sum(acceptance_ratio_tot)/(nb_loops*num_steps*max_tot_iter);
ratio_fin_abund = sum(X(:,end))/sum(Measured_Abund(:,end));

% figure(num_fig) %figure 1 : Simulated absolute abundances
for j = 1:S
    plot(Time_step, X(j,:), '-o', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
    hold on
end
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
axis([0 600 0 3e-04])

num_fig = num_fig + 1;
figure(num_fig) %figure 2 : Observed absolute abundances
for j = 1:S
    plot(Time_step, mean(Measured_Abund(j,:,:), 3), '-*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
    hold on
end
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
axis([0 600 0 3e-04])

num_fig = num_fig + 1;
figure(num_fig)
mean_col = mean(Measured_Abund, 3);
std_col = std(Measured_Abund, 1, 3);
for i = 1:S
 errorbar(Time_step, mean_col(i,:), std_col(i,:), '-*', 'Color', colors(j,:));
 hold on 
end

num_fig = num_fig + 1; 
z_temp = z_temp(1:(end-nb_Res), end);
z_temp = reshape(z_temp', 3, S);
z_temp = z_temp';
z = sum(z_temp(:,1:2), 2);%Total biomass of all species after 1 week
StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
[Shann_sim, Simp_sim] = Shannon_Simpson_Indices(S,StackPlot);
[Shann_obs, Simp_obs] = Shannon_Simpson_Indices(S, mean(StackPlot_Meas,3));

StackPlotTot = X./sum(X);

z_fin_sim = z(1:S,end); %Absolute stationary abundances
z_fin_obs = mean(Measured_Abund(1:S, end, :), 3); %Absolute stationary abundances
nb_Surv_Sim = sum(StackPlotTot > Threshold_Surviving);
sum(z_fin_sim)/sum(z_fin_obs)
acceptance_ratio_temp
liste_nrj_tot_temp

figure(num_fig); %figure 3 abondances relatives simulations
bar(StackPlotTot', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked all replicates')
num_fig = num_fig + 1;

figure(num_fig); %figure 4 Relative observed abundances
bar(mean(StackPlot_Meas,3)', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked observed')
num_fig = num_fig + 1;

[diff_vect, I, ~] = Sort_diff(StackPlotTot(:,end),mean(StackPlot_Meas(:,end,:),3));
name_order = name(I);
figure(num_fig);  %figure 5 Difference between observed and simulated resluts at the end of the experiment
stem(diff_vect);
xtickangle(90)
set(gca,'xtick',1:21,'xticklabel',name_order)
ylabel('Sim - Obs')

num_fig = num_fig + 1;

figure(num_fig); %figure 7 
for j = 1:S
    scatter(mean(StackPlot_Meas(j, 2, :), 3), StackPlotTot(j,2), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors(j,:), 'MarkerFaceColor', colors(S+1-j,:));%col(S+1-j,:))      
    hold on
end
axis([0 0.5 0 0.5]);
axis square
reflin = refline(1,0);
reflin.Color = 'r';
xlabel('Experiment'); 
ylabel('Simulation');
legend(name);
title('Scatter Tot');

num_fig = num_fig + 1;

figure(num_fig); %figure 7
for j = 1:S
    scatter(z_fin_obs(j), z_fin_sim(j), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors(j,:), 'MarkerFaceColor', colors(S+1-j,:));%col(S+1-j,:))      
    hold on
end
axis([0 3*10^(-4) 0 3*10^(-4)]);
reflin = refline(1,0);
axis square
reflin.Color = 'r';
xlabel('Experiment'); 
ylabel('Simulation');
legend(name);
title('Scatter absolute abundance');
num_fig = num_fig + 1;

figure(num_fig); %figure 8
plot(1:length(Time_step), Simp_obs, 'b--o')
hold on
plot(1:length(Time_step), Simp_sim, 'r--o')
num_fig = num_fig + 1;

figure(num_fig); %figure 9
plot(1:length(Time_step), nb_Surv_Obs, 'b--o')
hold on
plot(1:length(Time_step), nb_Surv_Sim, 'r--o')
num_fig = num_fig + 1;

num_fig = num_fig+1; %figure 14
figure(num_fig);
plot(1:length(liste_nrj_tot), liste_nrj_tot)

if save_data == 1
    FolderName = strcat(cd, '/Figures/');
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName   = num2str(get(FigHandle, 'Number'));
        FigName = strcat('Fig', FigName);
        FigName = strcat(FolderName, FigName, Name_file);
        set(0, 'CurrentFigure', FigHandle);
        saveas(FigHandle, FigName, 'pdf');
    end
    save(strcat(cd, '/Data/', Name_file,'_CrossFeed_Mat.mat'), 'CrossFeed_Mat_Lyso_Temp');
    save(strcat(cd, '/Data/', Name_file,'_Threshold.mat'), 'Threshold_CF');
    save(strcat(cd, '/Data/', Name_file,'_Threshold_Pred.mat'), 'Threshold_death');
    save(strcat(cd, '/Data/', Name_file,'_Threshold_Pred_Lyso.mat'), 'Threshold_Pred');
    save(strcat(cd, '/Data/', Name_file,'_Pred_Mat.mat'), 'Death_Mat_Lyso_Temp');
    save(strcat(cd, '/Data/', Name_file,'_Pred_Mat_Lyso.mat'), 'Pred_Mat_Lyso');
    save(strcat(cd, '/Data/', Name_file,'_Resource_Matrix.mat'), 'Resource_Matrix_Lyso');
    save(strcat(cd, '/Data/', Name_file,'_Lag_time_Pred.mat'), 'Lag_time_Pred_Lyso');
    save(strcat(cd, '/Data/', Name_file,'_Lag_time_Cons.mat'), 'Lag_time_Cons_Lyso');
    save(strcat(cd, '/Data/', Name_file,'_Mat_kappa_3.mat'), 'Mat_kappa_3_Lyso');
    save(strcat(cd, '/Data/', Name_file,'_Kappa_mat.mat'), 'kappa_mat_Lyso');
    save(strcat(cd, '/Data/', Name_file,'_R_mat.mat'), 'R');
    save(strcat(cd, '/Data/', Name_file,'_CrossFeed_Mat_init.mat'), 'CrossFeed_Mat_init');
    save(strcat(cd, '/Data/', Name_file,'_Threshold_init.mat'), 'Threshold_CF_init');
    save(strcat(cd, '/Data/', Name_file,'_Threshold_Pred_init.mat'), 'Threshold_death_init');
    save(strcat(cd, '/Data/', Name_file,'_Threshold_Pred_Lyso_init.mat'), 'Threshold_Pred_init');
    save(strcat(cd, '/Data/', Name_file,'_Pred_Mat_init.mat'), 'Death_Mat_Lyso_init');
    save(strcat(cd, '/Data/', Name_file,'_Pred_Mat_Lyso_init.mat'), 'Pred_Mat_Lyso_init');
    save(strcat(cd, '/Data/', Name_file,'_Resource_Matrix_init.mat'), 'Resource_Matrix_init');
    save(strcat(cd, '/Data/', Name_file,'_Lag_time_Pred_init.mat'), 'Lag_time_Pred_Lyso_init');
    save(strcat(cd, '/Data/', Name_file,'_Lag_time_Cons_init.mat'), 'Lag_time_Cons_Lyso_init');
    save(strcat(cd, '/Data/', Name_file,'_Mat_kappa_3_init.mat'), 'Mat_kappa_3_Lyso_init');
    save(strcat(cd, '/Data/', Name_file,'_Kappa_mat_init.mat'), 'kappa_mat_Lyso_init');
end