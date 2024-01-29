%Simulations subsytems with fitted interspecific interactions
clear
close all


%Save or Not
save_data = 0; %1 if save, 0 otherwise
Name_file = 'Death_CF_Model_V15';

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','25:46','Format','auto');
% Data_Evol = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 4, 'Range','24:44', 'Format','auto');
Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 3, 'Range','1:21', 'Format','auto'); %Test with the new SynCom 20
% Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 4, 'Range','1:21', 'Format','auto'); %Mixed data

%Load parameters
kappa_mat = load(strcat('Data/', Name_file, '_Kappa_mat.mat')); kappa_mat = kappa_mat.kappa_mat;
CrossFeed_Mat = load(strcat('Data/', Name_file, '_CrossFeed_Mat.mat')); CrossFeed_Mat = CrossFeed_Mat.CrossFeed_Mat_Temp;
Resource_Matrix = load(strcat('Data/', Name_file, '_Resource_Matrix.mat')); Resource_Matrix = Resource_Matrix.Resource_Matrix;
Pred_Mat = load(strcat('Data/', Name_file, '_Pred_Mat.mat')); Pred_Mat = Pred_Mat.Death_Mat_Temp;
Threshold = load(strcat('Data/', Name_file, '_Threshold.mat')); Threshold = Threshold.Threshold_CF;%Threshold.Threshold; For Philip data
Threshold_Pred = load(strcat('Data/', Name_file, '_Threshold_Pred.mat')); Threshold_Pred = Threshold_Pred.Threshold_death;
Lag_time_Cons = load(strcat('Data/', Name_file, '_Lag_time_Cons.mat')); Lag_time_Cons = Lag_time_Cons.Lag_time_Cons;
Lag_time_Pred = load(strcat('Data/', Name_file, '_Lag_time_Pred.mat')); Lag_time_Pred = Lag_time_Pred.Lag_time_Pred;
Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
R = load(strcat('Data/', Name_file, '_R_mat.mat')); R = R.R;
name = string(table2array(Parameters_set(1:20,1)));
Time_step = 0:1:500;%[0 1 3 7 10 21]*24;%Time step in hours

Data_Evol_temp = table2array(Data_Evol(:, 2:end));
std_y_0 = 0;
Data_Evol_temp = Data_Evol_temp(:,mod(1:length(Data_Evol_temp(1,:)), 4) == 2);%20% of the data to train. Hold-out method.
mean_y_0 = Data_Evol_temp(:,1);

% Measured_Abund = table2array(Data_Evol(1:20, 2:7));
%Number of surviving species after 8 weeks
Threshold_Surviving = 1e-10;

%Initialization of the model parameters fixed for all replicates 
t_0 = 0; %Time 0
nb_Res = length(R); %Number of resources (group of resources)
n_exp = 1; %No transfer in Philip experiment 
tspan = [0, max(Time_step)]; %[0, 22*24]; %Time interval in hours
S_sim = 1:20; %Species present into the subsystem
S = length(S_sim); %Number of species

kappa_mat = kappa_mat(S_sim,:);
CrossFeed_Mat = CrossFeed_Mat(S_sim, S_sim);
Resource_Matrix = Resource_Matrix(S_sim,:);
Pred_Mat = Pred_Mat(S_sim,S_sim);
Threshold = Threshold(S_sim);
Threshold_Pred = Threshold_Pred(S_sim);
Lag_time_Cons = Lag_time_Cons(S_sim,:);
Lag_time_Pred = 0;%Lag_time_Pred(S_sim,S_sim);
name = name(S_sim);
yield_Pred = 0;% 20% of yield for predation
mean_y_0 = mean_y_0(S_sim);

%Setting for Matlab ODE solver
nb_tot_Species = S*3 + nb_Res;
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.

%Initialization of the colors

colors = distinguishable_colors(60);

nb_replicates = 1;

num_fig = 1;
for i = 1:nb_replicates
    %Initial concentrations using a normal distribution
    mat_y_0 = mean_y_0;
    
    mat_y_0 = [mat_y_0 zeros(S,1) zeros(S,1)];

    for k = 1:n_exp
        y_0 = sum(mat_y_0(:,1:2),2);
        mat_y_0 = reshape(mat_y_0', 1, []);
        mat_y_0 = [mat_y_0 R];
   
%         sol = ode45(@(t, y) fun_Monod_Mult_Res(t, y, kappa_mat, CrossFeed_Mat, Resource_Matrix, Threshold, Threshold_Pred, Pred_Mat, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res), tspan,  mat_y_0, opts_1); %Multiple resource groups
%         sol = ode45(@(t, y) fun_Monod_Inhibit_Mult_Res(t, y, kappa_mat, CrossFeed_Mat, Resource_Matrix, Threshold, Threshold_Pred, Pred_Mat, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1);
%         sol = ode45(@(t, y) fun_Death(t, y, kappa_mat, CrossFeed_Mat, Resource_Matrix, Threshold, Threshold_Pred, Pred_Mat, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1); %Multiple resource groups
        sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat, CrossFeed_Mat, Mat_kappa_3, Resource_Matrix, Threshold, Threshold_Pred, Pred_Mat, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1); 
        z_temp = deval(sol, Time_step);
        sum(z_temp(:,1))
        sum(z_temp(:,4))
        sum(z_temp(:,end))
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
        R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass
        figure(num_fig)
        for j = 1:S
            plot(Time_step, X(j,:), '-*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
            hold on
        end
        legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
        num_fig = num_fig + 1;
        
        figure(num_fig)
        plot(Time_step, sum(R_temp), '-*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
        num_fig = num_fig + 1;

        figure(num_fig)
        for j = 1:S
            plot(Time_step, W(j,:), '-*', 'Color', colors(j,:));%, Time_step, R_temp, 'o');
            hold on
        end
        num_fig = num_fig + 1;


        z_temp = z_temp(1:(end-nb_Res), end);
        z_temp = reshape(z_temp',3, S);
        z_temp = z_temp';
        z = sum(z_temp(:,1:2), 2);%Total biomass of all species after 1 week
        StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
        [Shann_sim, Simp_sim] = Shannon_Simpson_Indices(S,StackPlot);
    end
    StackPlotTot = X./sum(X);
end

StackPlotTot = StackPlotTot./nb_replicates;
z_fin_sim = z(1:S,end); %Absolute stationary abundances
nb_Surv_Sim = sum(StackPlotTot > Threshold_Surviving);


figure(num_fig);
bar(StackPlotTot', 'stacked');
axis([0 11.5 0 1])
legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
title('Stacked all replicates')
num_fig = num_fig + 1;


figure(num_fig);
plot(1:length(Time_step), Simp_sim, 'r--o')
num_fig = num_fig + 1;

figure(num_fig);
plot(1:length(Time_step), nb_Surv_Sim, 'r--o')
num_fig = num_fig + 1;