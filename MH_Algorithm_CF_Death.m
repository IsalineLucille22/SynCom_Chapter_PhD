%Script made by Adline Vouillamoz and Isaline Guex
clear;
close all;

%Save or Not
save_data = 1; %1 if save, 0 otherwise
Name_file = 'Death_CF_Model_V5';

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','25:46', 'Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','71:91', 'Format','auto');
Parameters_Senka_Lag_time = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','94:114', 'Format','auto');
Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 4, 'Range','1:21', 'Format','auto');%readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 4, 'Range','1:21', 'Format','auto');%
% Data_Evol_test = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 4, 'Range','1:21', 'Format','auto');%readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 4, 'Range','1:21', 'Format','auto');
Data_Evol_test = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 3, 'Range','1:21', 'Format','auto'); %Test with the new SynCom 20
Time_step = [0 1 3 7 10 15 21 22]*24; %[0 1 3 7 10 15 21 22]*24; %[0 1 3 7 10 15 21 22]*24; %Measured time step in hours
tspan = [0, max(Time_step)]; %Time interval in hoursx
S = height(Data_Evol);

yield_vect = table2array(Parameters_set(1:S,5)); %Change it according to the desired mu_max
name = string(table2array(Parameters_set(1:S,1)));%string(table2array(table(2:22,3)));

mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8)); %Normal Parameters_Senka_mu_max(:,2:3). Log-normal Parameters_Senka_mu_max(:,7:8)
mu_max_vect = mu_max_dist(:,1);
Res_Percentage = table2array(Parameters_Senka_mu_max(:,12)); %table2array(Parameters_Senka_mu_max(:,4)); %Percentage of resources that can be consumed
Parameters_Senka_Lag_time = table2array(Parameters_Senka_Lag_time(:,7:8)); %table2array(Parameters_Senka_Lag_time(:,2:3));
mean_param_LN = mu_max_dist(:,1);
var_param_LN = mu_max_dist(:,2);

Data_Evol_temp = table2array(Data_Evol(:, 2:end));
Data_Evol_temp = Data_Evol_temp(:,~ismember(mod(1:length(Data_Evol_temp(1,:)), 8),[1, 2, 3]));%Data_Evol_temp(:,mod(1:length(Data_Evol_temp(1,:)),4) ~= 0);%80% of the data to train. Hold-out method.
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

liste_nrj_tot_temp = 10000;
while liste_nrj_tot_temp > 3.9
%Initialization of the model parameters fixed for all replicates 
t_0 = 0; %Time 0
nb_Res = 1; %Initial number of resources (group of resources)
nb_Res_max = 20; %Upper bound for the number of resources
CrossFeed_Mat = lognrnd((log(0.8) + mu_max_dist(:,1)).*ones(S,S), 0*mu_max_dist(:,2).*ones(S,S));%.val val% of zeros. ConsumerxProducer
% In Philip's experiment Lysobacter is not present
Death_Mat = lognrnd((log(0.8) + mean(mu_max_dist(:,1))).*ones(S,S), 0*mu_max_dist(:,2).*ones(S,S));%PredatorsxPreys %Rates matrix does not depend on mu_max.
Lag_time_Pred = zeros(S, S);
Lag_time_Pred_init = Lag_time_Pred;

CrossFeed_Mat_Temp = zeros(S,S); %Temporal predation matrix
rand_indice = rand(S,S) > 0.8; %Percentage of zeros. Larger is the value smaller is the number of interaction
CrossFeed_Mat_Temp(rand_indice) = CrossFeed_Mat(rand_indice);%Put some element to 0
Death_Mat_Temp = zeros(S,S);
Death_Mat_Temp(16,16) = Death_Mat(16,16);
Death_Mat_Temp(1:1+size(Death_Mat_Temp,1):end) = Death_Mat(1:1+size(Death_Mat,1):end);
CrossFeed_Mat_Temp(1:1+size(CrossFeed_Mat_Temp,1):end) = 0;
yield_Pred = 0.2; 

% Create a resource matrix
%Load parameters
Name_file_load = 'Death_CF_Model_V6';
R_0 = 0.5;
x = rand(S,nb_Res);
Res_Percentage_old = Res_Percentage;
Res_Percentage = min(max(Res_Percentage, 1), 1);%min(max(Res_Percentage, 0.24), 0.24);
rand_indice = x > Res_Percentage;
Resource_Matrix = ones(S,nb_Res);
Resource_Matrix(rand_indice) = 0;
Resource_Matrix = Resource_Matrix.*lognrnd(zeros(S, nb_Res), 0*mu_max_dist(:,2).*ones(S, nb_Res));%Resource_Matrix.*max(normrnd(ones(S, nb_Res), 0.05.*ones(S, nb_Res)), 0.01);%Resource_Matrix.*max(normrnd(ones(S, nb_Res), mu_max_dist(:,2)./mu_max_dist(:,1).*ones(S, nb_Res)), 0.01);%Centered on mu_max because multiplied mby mu_max so depend on the measured mu_max.
R = max(R_0*15*10^(-3)/nb_Res*min(normrnd(1, 0.3, 1, nb_Res), 1.1),0);%Concentration of carbon resources in g/ml;
Lag_time_Cons = lognrnd(Parameters_Senka_Lag_time(:,1).*ones(S,nb_Res), Parameters_Senka_Lag_time(:,2).*ones(S,nb_Res)); %max(normrnd(Parameters_Senka_Lag_time(:,1).*ones(S,nb_Res), Parameters_Senka_Lag_time(:,2).*ones(S,nb_Res)), 0); %Lag time for resource consumption
x = rand(S,nb_Res);
rand_indice = x > 0.3;% Res_Percentage_Nb_Cons ;
Lag_time_Rand = normrnd(50, 10, S, nb_Res);
Lag_time_Cons(rand_indice) = max(Lag_time_Cons(rand_indice), Lag_time_Rand(rand_indice));
Res_Percentage = Res_Percentage_old;

%Setting for Matlab ODE solver
%nb_tot_Species = S*3 + nb_Res;
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.

%Initialization of the colors

colors = distinguishable_colors(S);

num_fig = 1;
%Initial concentrations using a normal distribution
mat_y_0 = mean_y_0;
mat_y_0 = [mat_y_0 zeros(S,1) zeros(S,1)];

%Fit one unique threshold value
Threshold_CF = max(normrnd(1.0e-03, 0*0.1e-05, S, 1),0); 
Threshold_death = max(normrnd(1.0e-04, 0*0.1e-05, S, 1),0);%max(normrnd(1.0e-03, 0*0.1e-05, S, 1),0);
mat_y_0 = reshape(mat_y_0', 1, []);
mat_y_0 = [mat_y_0 R];
mu_max_vect_temp = exp(mu_max_vect + mu_max_dist(:,2).^2./2);%lognrnd(mu_max_vect, 0*mu_max_dist(:,2));%
yield = max(normrnd(yield_vect, 0.0), 0);%max(normrnd(yield_vect, 0.0), 1e-05);
kappa_mat = [2.5e+05*ones(S,1) mu_max_vect_temp (mu_max_vect_temp./yield - mu_max_vect_temp) 2.5e+05*ones(S,1)];
kappa_mat_init = kappa_mat; 
Mat_kappa_3_tot = repmat(kappa_mat(:,3), 1, S);%repmat((kappa_temp./yield' - kappa_temp), 1, S);%Kappa_3 rates defined by species not byproduct.
Mat_kappa_3 = Mat_kappa_3_tot;
Mat_kappa_3(CrossFeed_Mat_Temp == 0) = 0;Mat_kappa_3_init = Mat_kappa_3;
    
beta = 0.5*1e2; 
num_steps = 500;
nb_loops = 4; %Number of different loops (CF, predation. resource, nb of resources)
     
n = 1; r = 1; l = 1; p = 1; u = 1; %Total energy
max_tot_iter = 18;
liste_nrj_CF = zeros(1, num_steps*max_tot_iter);
Struct_CrossFeed_Mat = cell(num_steps*max_tot_iter, 1); %Creation d'une structure pour sauver les matrices générées
Struct_CrossFeed_Mat{1} = CrossFeed_Mat_Temp; CrossFeed_Mat_init = CrossFeed_Mat_Temp;
CrossFeed_Mat_Temp_cand = CrossFeed_Mat_Temp; %Initilization candidat
Threshold_CF_cand = Threshold_CF; Threshold_CF_init = Threshold_CF;
acceptance_ratio_CF = zeros(1, num_steps*max_tot_iter);
liste_nrj_pred = zeros(1, num_steps*max_tot_iter);
Struct_Death_Mat = cell(num_steps*max_tot_iter, 1); %Creation d'une structure pour sauver les matrices générées
Struct_Death_Mat{1} = Death_Mat_Temp; Death_Mat_init = Death_Mat_Temp;
Death_Mat_Temp_cand = Death_Mat_Temp; %Initilization candidat
Threshold_death_cand = Threshold_death; Threshold_death_init = Threshold_death_cand;
acceptance_ratio_pred = zeros(1, num_steps*max_tot_iter);
liste_nrj_res = zeros(1, num_steps*max_tot_iter);
Struct_Res_Mat = cell(num_steps*max_tot_iter, 1); %Creation d'une structure pour sauver les matrices générées
Struct_Res_Mat{1} = Resource_Matrix; Resource_Matrix_init = Resource_Matrix;
Resource_Matrix_cand = Resource_Matrix; %Initilization candidat
acceptance_ratio_res = zeros(1, num_steps*max_tot_iter);
liste_nrj_tot = zeros(1, nb_loops*num_steps*max_tot_iter);
acceptance_ratio_tot = zeros(1, nb_loops*num_steps*max_tot_iter);
nb_Res_new = nb_Res;
Lag_time_Cons_cand = Lag_time_Cons; Lag_time_Cons_init = Lag_time_Cons;
R_new = R; R_init = R; mat_y_0_new = mat_y_0; 
liste_nrj_nb_res = zeros(1, num_steps*max_tot_iter);
acceptance_ratio_nb_res = zeros(1, num_steps*max_tot_iter);
Evol_nb_Res = zeros(1, num_steps*max_tot_iter);

var_theta = 0.1;
acceptance_ratio_temp = 10000;
N = nb_loops*num_steps*max_tot_iter;
max_val = [10./mu_max_vect_temp 10*ones(S,1) 10./mu_max_vect_temp];
Hill_CF_coeff = 10; Hill_Pred_coeff = 10;
Hill_CF_coeff_cand = Hill_CF_coeff; Hill_Pred_coeff_cand = Hill_Pred_coeff;
delta_dist = 1;
mu_Likelihood = reshape(repmat(zeros(S,1), 1, nb_Res), 1, []);
var_Likelihood = reshape(repmat(5*mu_max_dist(:,2), 1, nb_Res), 1, []);
theta_res = var_theta*mean(mu_max_dist(:,2));
delta_dist_vect = [];
LT_CF_Death = zeros(S, S);
Cov_Lt_eye = eye(S*S);
while (liste_nrj_tot_temp > 3.25 || (acceptance_ratio_temp < 0.21 || acceptance_ratio_temp > 0.29) || p < 5000) &&  p < N && num_steps > 0%5.6, 0.22:0.27

    %%%%%%%%%% Resource matrix loop %%%%%%%%

    Cov_Mat_Res = 0.01*eye(S*nb_Res);
    Cov_Mat_LT = 0.01*eye(S*nb_Res);
    for k = 1:num_steps
        
        sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat, CrossFeed_Mat_Temp, Mat_kappa_3, Resource_Matrix_cand, Threshold_CF, Threshold_death, Death_Mat_Temp, yield_Pred, S, Lag_time_Cons_cand, Lag_time_Pred, nb_Res, Hill_CF_coeff, Hill_Pred_coeff), tspan,  mat_y_0, opts_1); %Multiple resource groups
        z_temp = deval(sol, Time_step);
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
        X = X + P;
        nu = k^(-3/4); 
        StackPlotTot = X./sum(X);
        
        [Resource_Matrix_cand, Resource_Matrix, T_cand, T, Struct_Res_Mat, liste_nrj_tot, liste_nrj_res, acceptance_ratio_tot, acceptance_ratio_res, p, l, ~, ~, ~, nu, Cov_Mat_Res, Cov_Mat_LT, Lag_time_Cons_cand, Lag_time_Cons] = ...
            fun_MH_Candidate_Rob(Resource_Matrix_cand, Resource_Matrix, zeros(nb_Res,1), zeros(nb_Res,1), Struct_Res_Mat, p, l, liste_nrj_tot, liste_nrj_res,...
            acceptance_ratio_tot, acceptance_ratio_res, S, nb_Res, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, 0, 0, 0, nu, Cov_Mat_Res, Cov_Mat_LT, max_val(:,3), nb_rep, Lag_time_Cons_cand, Lag_time_Cons);  
        
        if mod(p, 500) == 0
            disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
        end
    end

    %%%%%%%%%% Number of resources loop %%%%%%%%
    Lag_time_Cons_cand = Lag_time_Cons;
    Resource_Matrix_cand = Resource_Matrix;
    for k = 1:num_steps
        
        sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat, CrossFeed_Mat_Temp, Mat_kappa_3, Resource_Matrix_cand, Threshold_CF, Threshold_death, Death_Mat_Temp, yield_Pred, S, Lag_time_Cons_cand, Lag_time_Pred, nb_Res_new, Hill_CF_coeff, Hill_Pred_coeff), tspan,  mat_y_0_new, opts_1); %Multiple resource groups
        z_temp = deval(sol, Time_step);
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
        X = X + P;
        nu = k^(-3/4); 
        StackPlotTot = X./sum(X);
                  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Energy determination %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        somme_abs_tot_nb_res = sum(sum(sum((X - Measured_Abund).^2)))/nb_rep;
        somme_abs_tot_nb_res = somme_abs_tot_nb_res/sum(mean(Measured_Abund(:,end,:), 3).^2) + 2/nb_rep*sum(sum(sum((StackPlotTot - StackPlot_Meas).^2)))/sum(mean(StackPlot_Meas(:,end,:), 3).^2);
        energie_nb_res = somme_abs_tot_nb_res;

        rand_val = unifrnd(0,1);
        liste_nrj_nb_res(u) = energie_nb_res;
        liste_nrj_tot(p) = energie_nb_res;
    
        if p > 1
            delta_nrj_nb_res = liste_nrj_tot(p) - liste_nrj_tot(p-1);
        else
            delta_nrj_nb_res = 0;
        end

        ratio = min(exp(-beta*delta_nrj_nb_res), 1);
        if rand_val < ratio 
            Resource_Matrix = Resource_Matrix_cand;
            Lag_time_Cons = Lag_time_Cons_cand;
            acceptance_ratio_nb_res(u) = 1;
            acceptance_ratio_tot(p) = 1;
            Evol_nb_Res(u) = nb_Res_new;
            nb_Res = nb_Res_new;
            R = R_new; mat_y_0 = mat_y_0_new;
        else
            liste_nrj_tot(p) = liste_nrj_tot(p - 1);
            Evol_nb_Res(u) = nb_Res;
        end
        
        u = u + 1;
        p = p + 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Matrix modification %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %Modification matrice
        Num_add = [-1,0,1]; 
        Num_add = Num_add(randperm(numel(Num_add),1)); %0 if we don't want to change the number of resources
        nb_Res_new = nb_Res + Num_add; 
        Conc_new = sum(R);
        if Num_add > 0 && nb_Res_new <= nb_Res_max
            add_col = ones(S, Num_add);
            x = rand(S,Num_add);
            rand_indice = x > Res_Percentage;
            add_col(rand_indice) = 0;
            add_col = add_col.*lognrnd(zeros(S, Num_add), mu_max_dist(:,2).*ones(S, Num_add));%Centered on mu_max because multiplied mby mu_max so depend on the measured mu_max.
            Resource_Matrix_cand = [Resource_Matrix add_col];
            Lag_time_Cons_cand = [Lag_time_Cons lognrnd(Parameters_Senka_Lag_time(:,1).*ones(S, Num_add), Parameters_Senka_Lag_time(:,2).*ones(S, Num_add))]; %Lag time for resource consumption
            x = rand(S, Num_add);
            rand_indice = x > 0.3;%Res_Percentage_Nb_Cons;
            Lag_time_Rand = normrnd(50, 10, S, Num_add);
            Lag_time_Cons_cand(rand_indice, end) = max(Lag_time_Cons_cand(rand_indice, end), Lag_time_Rand(rand_indice));
            R_new = max(R_0*15*10^(-3)/nb_Res*min(normrnd(1, 0.3, 1, Num_add), 1.1),0);%Concentration of the new resource not normalized.
            R_new = [R R_new];
            R_new = R_new./(sum(R_new)/Conc_new); %Renormalization
            mat_y_0_new = [mat_y_0(1:(S*3)) R_new]; %Update initial concentrations
        elseif nb_Res_new > 0 && nb_Res_new <= nb_Res_max
            Resource_Matrix_cand = Resource_Matrix(:, 1:nb_Res_new);
            Lag_time_Cons_cand = Lag_time_Cons(:, 1:nb_Res_new);
            R_new = R(1:nb_Res_new); 
            mat_y_0_new = [mat_y_0(1:(S*3)) R_new]; %Update initial concentrations
        elseif nb_Res_new > 0
            Num_add = [-1,0];
            Num_add = Num_add(randperm(numel(Num_add),1));
            nb_Res_new = nb_Res + Num_add; 
            Resource_Matrix_cand = Resource_Matrix(:, 1:nb_Res_new);%Resource_Matrix(:, ind_temp);%
            Lag_time_Cons_cand = Lag_time_Cons(:, 1:nb_Res_new);%Lag_time_Cons(:, ind_temp);%
            R_new = R(1:nb_Res_new);
            mat_y_0_new = [mat_y_0(1:(S*3)) R_new]; %Update initial concentrations
        else
            nb_Res_new = nb_Res;
            Resource_Matrix_cand = Resource_Matrix;
            Lag_time_Cons_cand = Lag_time_Cons;
            R_new = R;
            mat_y_0_new = [mat_y_0(1:(S*3)) R_new]; %Update initial concentrations
        end
     
        if mod(p, 500) == 0
            disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
        end
    end
       
    
    %%%%%%%%%% Cross-feeding loop %%%%%%%%

    theta_T = 0.01*var_theta*mean(Threshold_CF);
    Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat_Temp_cand./mu_max_vect_temp;%CrossFeed_Mat_Temp_cand./yield' - CrossFeed_Mat_Temp_cand;
    Cov_Mat = 0.01*eye(S*S);
    for k = 1:num_steps
        sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat, CrossFeed_Mat_Temp_cand, Mat_kappa_3, Resource_Matrix, Threshold_CF_cand, Threshold_death, Death_Mat_Temp, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, Hill_CF_coeff_cand, Hill_Pred_coeff), tspan,  mat_y_0, opts_1); %Multiple resource groups
        z_temp = deval(sol, Time_step);
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
        X = X + P;
        nu = k^(-3/4);
        StackPlotTot = X./sum(X);

        CrossFeed_Mat_Temp = CrossFeed_Mat_Temp./(mu_max_vect_temp);
        CrossFeed_Mat_Temp_cand = CrossFeed_Mat_Temp_cand./(mu_max_vect_temp);

        [CrossFeed_Mat_Temp_cand, CrossFeed_Mat_Temp, Threshold_CF_cand, Threshold_CF, Struct_CrossFeed_Mat, liste_nrj_tot, liste_nrj_CF, acceptance_ratio_tot, acceptance_ratio_CF, p, n, ~, theta_T, ~, nu, Cov_Mat] = ...
            fun_MH_Candidate_Rob(CrossFeed_Mat_Temp_cand, CrossFeed_Mat_Temp, Threshold_CF_cand, Threshold_CF, Struct_CrossFeed_Mat, p, n, liste_nrj_tot, liste_nrj_CF,...
            acceptance_ratio_tot, acceptance_ratio_CF, S, S, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, 0, theta_T, 0, nu, Cov_Mat, Cov_Lt_eye, max_val(:,1), nb_rep, LT_CF_Death, LT_CF_Death);

        CrossFeed_Mat_Temp_cand(1:1+size(CrossFeed_Mat_Temp_cand,1):end) = 0;

        CrossFeed_Mat_Temp = CrossFeed_Mat_Temp.*(mu_max_vect_temp);
        CrossFeed_Mat_Temp_cand = CrossFeed_Mat_Temp_cand.*(mu_max_vect_temp);

        Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat_Temp_cand./mu_max_vect_temp;

        if mod(p, 500) == 0
            disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
        end          
    end
    Cov_Mat_CF = Cov_Mat;
        
    %%%%%%%%%% Death loop %%%%%%%%
    
    theta = var_theta*mean(mean(Death_Mat_Temp_cand)); 
    theta_T_death = 0.01*var_theta*mean(Threshold_death);
    Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat_Temp./mu_max_vect_temp;
    Cov_Mat = 0.01*eye(S);
    for k = 1:num_steps
        
        sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat, CrossFeed_Mat_Temp, Mat_kappa_3, Resource_Matrix, Threshold_CF, Threshold_death_cand, Death_Mat_Temp_cand, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, Hill_CF_coeff, Hill_Pred_coeff_cand), tspan,  mat_y_0, opts_1); %Multiple resource groups
        z_temp = deval(sol, Time_step);
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
        X = X + P;
        nu = k^(-3/4); 
        StackPlotTot = X./sum(X);
        
        [Death_Mat_Temp_cand, Death_Mat_Temp, Threshold_death_cand, Threshold_death, Struct_Death_Mat, liste_nrj_tot, liste_nrj_pred, acceptance_ratio_tot, acceptance_ratio_pred, p, r, theta, theta_T_death, ~, nu, Cov_Mat] = ...
            fun_MH_Candidate_Rob_Death(Death_Mat_Temp_cand, Death_Mat_Temp, Threshold_death_cand, Threshold_death, Struct_Death_Mat, p, r, liste_nrj_tot, liste_nrj_pred,...
            acceptance_ratio_tot, acceptance_ratio_pred, S, S, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, theta, theta_T_death, 0, nu, Cov_Mat, Cov_Lt_eye, max_val(:,2), nb_rep, LT_CF_Death, LT_CF_Death);
 
        ind = eye(S,S);
        Death_Mat_Temp_cand(~ind) = 0; %Changement pour diagonale zeros

        if mod(p, 500) == 0
            disp(strcat("Iteration ", num2str(p), "/", num2str(N)))
        end
    end
    Cov_Mat_Death = Cov_Mat;
    Resource_Matrix_cand = Resource_Matrix; nb_Res_new = nb_Res; Lag_time_Cons_cand = Lag_time_Cons; R_new = R; mat_y_0_new = mat_y_0; %Re-initialization if the last candidate is not accepted. Otherwise issues with the resource loop.
    acceptance_ratio_temp = sum(acceptance_ratio_tot)/p;
    liste_nrj_tot_temp = liste_nrj_tot(p - 1);
    delta_dist = 1;
    mu_Likelihood = reshape(repmat(zeros(S,1), 1, nb_Res), 1, []);
    var_Likelihood = reshape(repmat(mu_max_dist(:,2), 1, nb_Res), 1, []);
    Temp = 1/beta; Temp = Temp*(1-u/(3*N)); beta = 1/Temp;
    if mod(p, 12001) == 0 && liste_nrj_tot_temp > 5
        plot(1:length(liste_nrj_tot), liste_nrj_tot)
        break      
    end
end 
acceptance_ratio_temp
liste_nrj_tot_temp
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Data to test on half of the remaining data
Time_step = [0 1 3 7 10 21]*24; %Measured time step in hours
tspan = [0, max(Time_step)]; %Time interval in hours
Data_Evol_temp = table2array(Data_Evol_test(:, 2:end)); %table2array(Data_Evol(:, 2:end));
Data_Evol_temp = Data_Evol_temp(:,ismember(mod(1:length(Data_Evol_temp(1,:)), 4), [1, 2, 3]));%20% of the data to test
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

sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat, CrossFeed_Mat_Temp, Mat_kappa_3, Resource_Matrix, Threshold_CF, Threshold_death, Death_Mat_Temp, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, Hill_CF_coeff, Hill_Pred_coeff), tspan,  mat_y_0, opts_1); %Multiple resource groups
z_temp = deval(sol, Time_step);
X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
X = X + P;

CrossFeed_Mat_Temp_rates = CrossFeed_Mat_Temp.*mu_max_vect_temp; Resource_Matrix_rates = Resource_Matrix.*mu_max_vect_temp;
        
%%% Figure generation %%%

acceptance_ratio_fin = sum(acceptance_ratio_tot)/(nb_loops*num_steps*max_tot_iter);
ratio_fin_abund = sum(X(:,end))/sum(Measured_Abund(:,end));

% figure(num_fig) %figure 1 : Simulated absolute abundances
for j = 1:S
    plot(Time_step, X(j,:), '-o', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
    hold on
end
%legend(name, 'Orientation', 'vertical', 'Location', 'southeast')

num_fig = num_fig + 1;
figure(num_fig) %figure 2 : Observed absolute abundances
for j = 1:S
    plot(Time_step, mean(Measured_Abund(j,:,:), 3), '-*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
    hold on
end
%legend(name, 'Orientation', 'vertical', 'Location', 'southeast')

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
z_temp = reshape(z_temp',3, S);
z_temp = z_temp';
z = sum(z_temp(:,1:2), 2);%Total biomass of all species after 1 week
StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
[Shann_sim, Simp_sim] = Shannon_Simpson_Indices(S,StackPlot);
[Shann_obs, Simp_obs] = Shannon_Simpson_Indices(S, mean(StackPlot_Meas,3));
rand_prop = max(normrnd(1, 0.3, S, 1), 0);
mat_y_0 = (z_temp(:,1:2)/10).*rand_prop; %Take the 10% of the system
mat_y_0 = [mat_y_0 zeros(S,1)];

StackPlotTot = X./sum(X);

z_fin_sim = z(1:S,end); %Absolute stationary abundances
z_fin_obs = mean(Measured_Abund(1:S, end, :), 3); %Absolute stationary abundances
nb_Surv_Sim = sum(StackPlotTot > Threshold_Surviving);
sum(z_fin_sim)/sum(z_fin_obs)
acceptance_ratio_temp
liste_nrj_tot_temp

figure(num_fig); 
bar(StackPlotTot', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked all replicates')
num_fig = num_fig + 1;

figure(num_fig); 
bar(mean(StackPlot_Meas,3)', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked observed')
num_fig = num_fig + 1;

[diff_vect, I, ~] = Sort_diff(StackPlotTot(:,end),mean(StackPlot_Meas(:,end,:),3));
name_order = name(I);
figure(num_fig);
stem(diff_vect);
xtickangle(90)
set(gca,'xtick',1:21,'xticklabel',name_order)
ylabel('Sim - Obs')

num_fig = num_fig + 1;

figure(num_fig);
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

figure(num_fig);
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

figure(num_fig);
plot(1:length(Time_step), Simp_obs, 'b--o')
hold on
plot(1:length(Time_step), Simp_sim, 'r--o')
num_fig = num_fig + 1;

figure(num_fig);
plot(1:length(Time_step), nb_Surv_Obs, 'b--o')
hold on
plot(1:length(Time_step), nb_Surv_Sim, 'r--o')
num_fig = num_fig + 1;

num_fig = num_fig+1;
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
    save(strcat(cd, '/Data/', Name_file,'_CrossFeed_Mat.mat'), 'CrossFeed_Mat_Temp');
    save(strcat(cd, '/Data/', Name_file,'_Threshold.mat'), 'Threshold_CF');
    save(strcat(cd, '/Data/', Name_file,'_Threshold_Pred.mat'), 'Threshold_death');
    save(strcat(cd, '/Data/', Name_file,'_Pred_Mat.mat'), 'Death_Mat_Temp');
    save(strcat(cd, '/Data/', Name_file,'_Resource_Matrix.mat'), 'Resource_Matrix');
    save(strcat(cd, '/Data/', Name_file,'_Lag_time_Pred.mat'), 'Lag_time_Pred');
    save(strcat(cd, '/Data/', Name_file,'_Lag_time_Cons.mat'), 'Lag_time_Cons');
    save(strcat(cd, '/Data/', Name_file,'_Mat_kappa_3.mat'), 'Mat_kappa_3');
    save(strcat(cd, '/Data/', Name_file,'_Kappa_mat.mat'), 'kappa_mat');
    save(strcat(cd, '/Data/', Name_file,'_R_mat.mat'), 'R');
    save(strcat(cd, '/Data/', Name_file,'_CrossFeed_Mat_init.mat'), 'CrossFeed_Mat_init');
    save(strcat(cd, '/Data/', Name_file,'_Threshold_init.mat'), 'Threshold_CF_init');
    save(strcat(cd, '/Data/', Name_file,'_Threshold_Pred_init.mat'), 'Threshold_death_init');
    save(strcat(cd, '/Data/', Name_file,'_Pred_Mat_init.mat'), 'Death_Mat_init');
    save(strcat(cd, '/Data/', Name_file,'_Resource_Matrix_init.mat'), 'Resource_Matrix_init');
    save(strcat(cd, '/Data/', Name_file,'_Lag_time_Pred_init.mat'), 'Lag_time_Pred_init');
    save(strcat(cd, '/Data/', Name_file,'_Lag_time_Cons_init.mat'), 'Lag_time_Cons_init');
    save(strcat(cd, '/Data/', Name_file,'_Mat_kappa_3_init.mat'), 'Mat_kappa_3_init');
    save(strcat(cd, '/Data/', Name_file,'_Kappa_mat_init.mat'), 'kappa_mat_init');
    save(strcat(cd, '/Data/', Name_file,'_R_mat_init.mat'), 'R_init');
    %Structures saved
    save(strcat(cd, '/Data/', Name_file,'_Struct_Res_Mat.mat'), 'Struct_Res_Mat');
    save(strcat(cd, '/Data/', Name_file,'_Struct_Death_Mat.mat'), 'Struct_Death_Mat');
    save(strcat(cd, '/Data/', Name_file,'_Struct_CrossFeed_Mat.mat'), 'Struct_CrossFeed_Mat');
    %Covariance estimators saved
    save(strcat(cd, '/Data/', Name_file,'_Cov_Mat_CF.mat'), 'Cov_Mat_CF');
    save(strcat(cd, '/Data/', Name_file,'_Cov_Mat_Death.mat'), 'Cov_Mat_Death');
    save(strcat(cd, '/Data/', Name_file,'_Cov_Mat_Res.mat'), 'Cov_Mat_Res');
    save(strcat(cd, '/Data/', Name_file,'_Cov_Mat_LT.mat'), 'Cov_Mat_LT');
end