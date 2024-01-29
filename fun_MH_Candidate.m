function [Loop_Mat_Temp_cand, Loop_Mat_Temp, Threshold_Loop_cand, Threshold_Loop, Struct_Loop_Mat, liste_nrj_tot, liste_nrj_loop, acceptance_ratio_tot, acceptance_ratio_loop, p, n_loop, theta, theta_T_Loop, theta_LT_Loop, nu, LT_Loop_cand, LT_Loop] = fun_MH_Candidate(Loop_Mat_Temp_cand, Loop_Mat_Temp, Threshold_Loop_cand, Threshold_Loop, Struct_Loop_Mat, p, n_loop, liste_nrj_tot, liste_nrj_loop, acceptance_ratio_tot, acceptance_ratio_loop, S_consumer, S_consumed, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, theta, theta_T_Loop, theta_LT_Loop, nu, max_val, nb_rep, LT_Loop_cand, LT_Loop)
%Function Metropolis Hasting. Compute the energy of the new candidate and
%compare it to the present state.      
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Energy determination %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    somme_abs_tot = sum(sum(sum((X - Measured_Abund).^2)))/nb_rep; %sum(sum(sum(abs((X - Measured_Abund)))))/4;%
    somme_abs_tot = somme_abs_tot/sum(mean(Measured_Abund(:,end,:), 3).^2) + 2/nb_rep*sum(sum(sum((StackPlotTot - StackPlot_Meas).^2)))/sum(mean(StackPlot_Meas(:,end,:), 3).^2);%somme_abs_tot/sum(mean(Measured_Abund(:,end,:), 3)) + 3/4*sum(sum(sum((StackPlotTot - StackPlot_Meas))));%somme_abs_tot/sum(mean(Measured_Abund(:,end,:), 3).^2) + 3*sum(sum(sum((StackPlotTot - StackPlot_Meas).^2)))/sum(mean(StackPlot_Meas(:,end,:), 3).^2);%2*3*sum(sum(sum((StackPlotTot - StackPlot_Meas).^2)));
    energie = somme_abs_tot;

    rand_val = unifrnd(0,1);
    liste_nrj_loop(n_loop) = energie; 
    liste_nrj_tot(p) = energie;

    if p > 1
        delta_nrj = liste_nrj_tot(p) - liste_nrj_tot(p-1); 
    else
        delta_nrj = 0;
    end

    ratio = min(exp(-beta*delta_nrj), 1);%min(exp(delta_nrj),1);%
    if rand_val < ratio 
        Loop_Mat_Temp = Loop_Mat_Temp_cand;
        Threshold_Loop = Threshold_Loop_cand;
        LT_Loop = LT_Loop_cand;
        acceptance_ratio_tot(p) = 1;
        acceptance_ratio_loop(n_loop) = 1;
    else
%         liste_nrj_loop(n_loop) = liste_nrj_loop(n_loop - 1);
        liste_nrj_tot(p) = liste_nrj_tot(p - 1);
    end
    
    Struct_Loop_Mat{n_loop + 1} = Loop_Mat_Temp;
    n_loop = n_loop + 1;
    p = p + 1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Matrix modification %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Loop_Mat_Temp_cand = Loop_Mat_Temp; %Création du mouveau candidat à partir de la matrice sauvée.
    rand_inidice = rand(S_consumer, S_consumed) > 0.95;
    Loop_Mat_Temp_cand(rand_inidice) = max(Loop_Mat_Temp(rand_inidice) + normrnd(0, theta, sum(sum(rand_inidice)), 1),0); %min(max(Loop_Mat_Temp(rand_inidice) + normrnd(0, theta, sum(sum(rand_inidice)), 1),0), max_val); %
    High_val = find(Loop_Mat_Temp_cand > max_val);
    Loop_Mat_Temp_cand(High_val) = Loop_Mat_Temp(High_val); %Upper boundary
    logtheta = log(theta) + nu*(ratio - 0.234);
    theta = exp(logtheta);

    % Modification thresholds
%     rand_inidice = (sum(rand_inidice) > 0)';
%     Threshold_Loop_cand = Threshold_Loop;
%     zeros_ind = unique(randi(S_consumed, 15, 1));
%     Threshold_Loop_cand(rand_inidice) = max(Threshold_Loop(rand_inidice) + normrnd(0, theta_T_Loop, sum(rand_inidice), 1),0);
    Threshold_Loop_cand(1) = max(Threshold_Loop(1) + normrnd(0, theta_T_Loop),0);
%     Threshold_Loop_cand(zeros_ind) = Threshold_Loop(zeros_ind);
    logtheta = log(theta_T_Loop) + nu*(ratio - 0.234);
    theta_T_Loop = exp(logtheta);

    % Modification lag times %In loop either threshold fitted or lag times
    % but not both (perhaps for model with predation)% No upper bound for
    % the lag time value
    LT_Loop_cand = LT_Loop;
    rand_inidice = rand(S_consumer, S_consumed) > 0.98;
    LT_Loop_cand(rand_inidice) = max(LT_Loop(rand_inidice) + normrnd(0, theta_LT_Loop, sum(sum(rand_inidice)), 1),0);
%     LT_Loop_cand(rand_inidice) = max(LT_Loop(rand_inidice) + normrnd(zeros(sum(sum(rand_inidice)), 1), theta_LT_Loop(rand_inidice), sum(sum(rand_inidice)), 1),0);
    logtheta = log(theta_LT_Loop) + nu*(ratio - 0.234);%
%     logtheta = log(theta_LT_Loop(rand_inidice)) + nu*(ratio - 0.234);
    theta_LT_Loop = exp(logtheta);%theta_LT_Loop(rand_inidice) = exp(logtheta);
end