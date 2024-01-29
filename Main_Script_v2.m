%Simulations with fitted interspecific interactions
clear
close all


%Save or Not
save_data = 0; %1 if save, 0 otherwise
Name_file = 'Death_CF_Model_V6';%'Rand_Parameters';%
Name_file_Saved = 'Rand_Parameters';

%Loadind data
Parameters_set = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 8,'Range','25:46','Format','auto');
Parameters_Senka_mu_max = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 9, 'Range','71:91', 'Format','auto');
% Data_Evol = readtable(strcat('Data/','MergedData.xlsx'), 'Sheet', 4, 'Range','24:44', 'Format','auto');
Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 3, 'Range','1:21', 'Format','auto'); %Test with the new SynCom 20
% Data_Evol = readtable(strcat('Data/','S20_S21_abs_abund_cfu_Senka.xlsx'), 'Sheet', 4, 'Range','1:21', 'Format','auto'); %Mixed data
mu_max_dist = table2array(Parameters_Senka_mu_max(:,7:8));
S = height(Data_Evol);

%Load parameters
kappa_mat = load(strcat('Data/', Name_file, '_Kappa_mat.mat')); kappa_mat = kappa_mat.kappa_mat;
% kappa_mat(:,2) = lognrnd(log(0.32).*ones(S,1), 0.5*ones(S,1));
CrossFeed_Mat = load(strcat('Data/', Name_file, '_CrossFeed_Mat.mat')); CrossFeed_Mat = CrossFeed_Mat.CrossFeed_Mat_Temp;%0.8*repmat(kappa_mat(:,2), 1, 20);%CrossFeed_Mat.CrossFeed_Mat;%lognrnd(log(0.15).*ones(S,S), 1.7.*ones(S,S));%
Resource_Matrix = load(strcat('Data/', Name_file, '_Resource_Matrix.mat')); Resource_Matrix = Resource_Matrix.Resource_Matrix;%lognrnd(log(0.4).*ones(S,S), ones(S,S));%
Pred_Mat = load(strcat('Data/', Name_file, '_Pred_Mat.mat')); Pred_Mat = Pred_Mat.Death_Mat_Temp;%Pred_Mat.Pred_Mat;%diag(lognrnd(log(0.1).*ones(1,S), 0.1.*ones(1,S)));%
Threshold = load(strcat('Data/', Name_file, '_Threshold.mat')); Threshold = Threshold.Threshold_CF;%Threshold.Threshold_CF;%Threshold.Threshold; For Philip data
Threshold_Pred = load(strcat('Data/', Name_file, '_Threshold_Pred.mat')); Threshold_Pred = Threshold_Pred.Threshold_death; %Threshold_Pred.Threshold_Pred;%
Lag_time_Cons = load(strcat('Data/', Name_file, '_Lag_time_Cons.mat')); Lag_time_Cons = Lag_time_Cons.Lag_time_Cons;
Lag_time_Pred = load(strcat('Data/', Name_file, '_Lag_time_Pred.mat')); Lag_time_Pred = Lag_time_Pred.Lag_time_Pred;
R = load(strcat('Data/', Name_file, '_R_mat.mat')); R = R.R;
nb_Res = length(R); %Number of resources (group of resources)
Mat_kappa_3 = kappa_mat(:,3).*CrossFeed_Mat./kappa_mat(:,2);
% CrossFeed_Mat = CrossFeed_Mat./kappa_mat(:,2);
name = string(table2array(Parameters_set(1:20,1)));
Time_step = [0 1 3 7 10 21]*24; %[0 1 3 7 15 22]*24;%Measured time step in hours
% Time_step = [0 1 3 7 10 15 21 22]*24; %Measured time step in hours
Var_Resource_Matrix = lognrnd(zeros(S,nb_Res), 0*repmat(mu_max_dist(:,2), 1, nb_Res), S, nb_Res);
Resource_Matrix = Resource_Matrix.*Var_Resource_Matrix;

% x = rand(S,nb_Res);
% rand_indice = x > 0.3;% Res_Percentage_Nb_Cons ;
% Lag_time_Rand = normrnd(50, 10, S, nb_Res);%lognrnd(log(50), 1, S, nb_Res);
% Lag_time_Cons(rand_indice) = max(Lag_time_Cons(rand_indice), Lag_time_Rand(rand_indice));

% CrossFeed_Mat = 0.8*repmat(kappa_mat(:,2), 1, S);
% CrossFeed_Mat(1:1+size(CrossFeed_Mat,1):end) = 0;
% CrossFeed_Mat_Temp = zeros(S,S); %Temporal predation matrix
% rand_indice = rand(S,S) > 0.8; %Percentage of zeros. Larger is the value smaller is the number of interaction
% CrossFeed_Mat(rand_indice) = 0;%Put some element to 0


Data_Evol_temp = table2array(Data_Evol(:, 2:end));
std_y_0 = 0;
Data_Evol_temp = Data_Evol_temp(:, ismember(mod(1:length(Data_Evol_temp(1,:)),4), [1, 2, 3]));%20% of the data to train. Hold-out method.
nb_obs = length(Data_Evol_temp(1,:));
nb_time_step = length(Time_step);
nb_rep = nb_obs/nb_time_step;
Measured_Abund = zeros(S, nb_time_step, nb_rep); %Number species, number times, number replicates.
for i = 1:nb_rep
    Measured_Abund(:,:,i) = Data_Evol_temp(:, mod(1:nb_obs, nb_rep) == (i - 1));
end
StackPlot_Meas = Measured_Abund./sum(Measured_Abund);
rand = unifrnd(0,1);
mean_y_0 = Measured_Abund(:, 1, 2) +  (1 - rand)*Measured_Abund(:, 1, 1);
mean_y_0 = mean_y_0 + normrnd(0, 0.05*mean(mean_y_0));
Measured_Abund = mean(Measured_Abund,3);


% Measured_Abund = table2array(Data_Evol(1:20, 2:7));
%Number of surviving species after 8 weeks
 Threshold_Surviving = 1e-10;
nb_Surv_Obs = sum(Measured_Abund > Threshold_Surviving);

%Initialization of the model parameters fixed for all replicates 
t_0 = 0; %Time 0
n_exp = 1; %No transfer in Philip experiment 
tspan = [0, max(Time_step)*24]; %Time interval in hours

% R = max(0.5*15*10^(-3)/nb_Res*min(normrnd(1, 0.1, 1, nb_Res), 1.1),0);%[0.00149 0.001449 0.0016 0.00127 0.00172];%Concentration of carbon resources in g/ml;
yield_Pred = 0;% 20% of yield for predation

%Setting for Matlab ODE solver
nb_tot_Species = S*3 + nb_Res;
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-9);%,'NonNegative',1:nb_tot_Species); %To smooth the curves obtained using ode45.

%Initialization of the colors

%Initialization of the colors
colors_ind = [4, 1, 10, 15, 12, 3, 8, 6, 20, 18, 19, 2, 5, 21, 13, 9, 14, 7, 16, 17];
colors_init = {'#B35806', '#E08214', '#D53E4F', '#B2182B', '#B6604B', '#C51B7D', ...
               '#DE77AE', '#F1B6DA', '#FDAE61', '#FEE090', '#A6D96A', '#5AAE61', ...
               '#01665E', '#35978F', '#1B7837', '#C2A5CF', '#9970AB', '#762A83', ...
               '#80CDC1', '#C7EAE5', '#2166AC', '#4393C3'};%distinguishable_colors(S);
colors = {};
for i = 1:S
    colors{i} = colors_init{colors_ind(i)};
end

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
        sol = ode45(@(t, y) fun_CF_Death(t, y, kappa_mat, CrossFeed_Mat, Mat_kappa_3, Resource_Matrix, Threshold, Threshold_Pred, Pred_Mat, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res, 10, 10), tspan,  mat_y_0, opts_1); %Multiple resource groups
        z_temp = deval(sol, Time_step);
        sum(z_temp(:,1))
        sum(z_temp(:,4))
        sum(z_temp(:,end))
        X = z_temp(mod(1:S*3,3) == 1,:); %Species'biomass
        P = z_temp(mod(1:S*3,3) == 2,:); %Complexes'biomass
        W = z_temp(mod(1:S*3,3) == 0,:); %Byproducts'biomass
        R_temp = z_temp(S*3 + 1: end,:); %Byproducts'biomass
        X = X + P;
        figure(num_fig)
        for j = 1:S
            plot(Time_step, X(j,:), '-o', 'Color', colors{j});%, Time_step, R_temp, 'o');
            hold on
        end
        %legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
        axis([0 600 0 4e-04])
        num_fig = num_fig + 1;
        figure(num_fig)
        for j = 1:S
            plot(Time_step, Measured_Abund(j,:), '-*', 'Color', colors{j});%, Time_step, R_temp, 'o');
            hold on
        end
        %legend(name, 'Orientation', 'vertical', 'Location', 'southeast')
        axis([0 600 0 4e-04])

        num_fig = num_fig + 1;
        z_temp = z_temp(1:(end-nb_Res), end);
        z_temp = reshape(z_temp',3, S);
        z_temp = z_temp';
        z = sum(z_temp(:,1:2), 2);%Total biomass of all species after 1 week
        StackPlot = X./sum(X); %Proportion of each species into the system at the end of the cycle
        [Shann_sim, Simp_sim] = Shannon_Simpson_Indices(S,StackPlot);
        [Shann_obs, Simp_obs] = Shannon_Simpson_Indices(S, mean(StackPlot_Meas,3));
    end
    StackPlotTot = X./sum(X);
end

StackPlotTot = StackPlotTot./nb_replicates;
z_fin_sim = z(1:S,end); %Absolute stationary abundances
z_fin_obs = Measured_Abund(1:S, end); %Absolute stationary abundances
nb_Surv_Sim = sum(StackPlotTot > Threshold_Surviving);
sum(z_fin_sim)/sum(z_fin_obs)

h_vect = zeros(1,S);
p_vect = zeros(1,S);
for i = 1:S
    [h_vect(i),p_vect(i)] = ttest2(StackPlot_Meas(i,:)', StackPlotTot(i,:)');
end
name_diff = name(logical(h_vect));

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

[diff_vect, I, fact] = Sort_diff(StackPlotTot(:,end),StackPlot_Meas(:,end));
name_order = name(I);
h = figure(num_fig);
stem(diff_vect);
xtickangle(90)
set(gca,'xtick',1:21,'xticklabel',name_order)
ylabel('Sim - Obs')
num_fig = num_fig + 1;

figure(num_fig);
for j = 1:S
    scatter(StackPlot_Meas(j, n_exp+1), StackPlotTot(j,n_exp+1), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors{j}, 'MarkerFaceColor', colors{S+1-j});%col(S+1-j,:))      
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
    scatter(z_fin_obs(j), z_fin_sim(j), 100, 'd', 'LineWidth', 5, 'MarkerEdgeColor', colors{j}, 'MarkerFaceColor', colors{S+1-j});%col(S+1-j,:))      
    hold on
end
axis([0 3*10^(-4) 0 4*10^(-4)]);
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


if save_data == 1
    save(strcat(cd, '/Data/', Name_file_Saved,'_CrossFeed_Mat.mat'), 'CrossFeed_Mat');
    save(strcat(cd, '/Data/', Name_file_Saved,'_Threshold.mat'), 'Threshold');
    save(strcat(cd, '/Data/', Name_file_Saved,'_Threshold_Pred.mat'), 'Threshold_Pred');
    save(strcat(cd, '/Data/', Name_file_Saved,'_Pred_Mat.mat'), 'Pred_Mat');
    save(strcat(cd, '/Data/', Name_file_Saved,'_Resource_Matrix.mat'), 'Resource_Matrix');
    save(strcat(cd, '/Data/', Name_file_Saved,'_Lag_time_Pred.mat'), 'Lag_time_Pred');
    save(strcat(cd, '/Data/', Name_file_Saved,'_Lag_time_Cons.mat'), 'Lag_time_Cons');
    save(strcat(cd, '/Data/', Name_file_Saved,'_Mat_kappa_3.mat'), 'Mat_kappa_3');
    save(strcat(cd, '/Data/', Name_file_Saved,'_Kappa_mat.mat'), 'kappa_mat');
    save(strcat(cd, '/Data/', Name_file_Saved,'_Pred_Mat_Lyso.mat'), 'Pred_Mat');
    save(strcat(cd, '/Data/', Name_file_Saved,'_R_mat.mat'), 'R');
end