function [Loop_Mat_Temp_cand, Loop_Mat_Temp, Threshold_Loop_cand, Threshold_Loop, Hill_coeff_cand, Hill_coeff, Struct_Loop_Mat, liste_nrj_tot, liste_nrj_loop, acceptance_ratio_tot, acceptance_ratio_loop, p, n_loop, theta, theta_T_Loop, theta_Hill_Loop, nu] = fun_MH_Candidate_Hill(Loop_Mat_Temp_cand, Loop_Mat_Temp, Threshold_Loop_cand, Threshold_Loop, Hill_coeff_cand, Hill_coeff, Struct_Loop_Mat, p, n_loop, liste_nrj_tot, liste_nrj_loop, acceptance_ratio_tot, acceptance_ratio_loop, S_consumer, S_consumed, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, theta, theta_T_Loop, theta_Hill_Loop, nu, max_val, nb_rep)
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
        Hill_coeff = Hill_coeff_cand;
        acceptance_ratio_tot(p) = 1;
        acceptance_ratio_loop(n_loop) = 1;
    else
        liste_nrj_loop(n_loop) = liste_nrj_loop(n_loop - 1);
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
    Threshold_Loop_cand = Threshold_Loop;
if mod(p,2) == 0
    Loop_Mat_Temp_cand(rand_inidice) = max(Loop_Mat_Temp(rand_inidice) + normrnd(0, theta, sum(sum(rand_inidice)), 1),0); %min(max(Loop_Mat_Temp(rand_inidice) + normrnd(0, theta, sum(sum(rand_inidice)), 1),0), max_val); %
    High_val = find(Loop_Mat_Temp_cand > max_val);
    Loop_Mat_Temp_cand(High_val) = Loop_Mat_Temp(High_val); %Upper boundary
    logtheta = log(theta) + nu*(ratio - 0.234);
    theta = exp(logtheta);
else
    % Modification thresholds
    rand_inidice = (sum(rand_inidice) > 0)';
    zeros_ind = unique(randi(S_consumed, 15, 1));
    Threshold_Loop_cand(rand_inidice) = max(Threshold_Loop(rand_inidice) + normrnd(0, theta_T_Loop, sum(rand_inidice), 1),0);
    Threshold_Loop_cand(zeros_ind) = Threshold_Loop(zeros_ind);
    logtheta = log(theta_T_Loop) + nu*(ratio - 0.234);
    theta_T_Loop = exp(logtheta);
end

    % Modification Hill factor %On the index?
    Num_add = 1 - sign(rand(1))*poissrnd(1);%[-1,0,1]; 
    Hill_coeff_cand = min(max(Hill_coeff + Num_add, 1), 20); 
    logtheta = log(theta_Hill_Loop) + nu*(ratio - 0.234);
    theta_Hill_Loop = exp(logtheta);
end