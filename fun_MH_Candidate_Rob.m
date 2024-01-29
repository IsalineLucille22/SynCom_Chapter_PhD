function [Loop_Mat_Temp_cand, Loop_Mat_Temp, Threshold_Loop_cand, Threshold_Loop, Struct_Loop_Mat, liste_nrj_tot, liste_nrj_loop, acceptance_ratio_tot, acceptance_ratio_loop, p, n_loop, theta, theta_T_Loop, theta_LT_Loop, nu, Cov_Mat, Cov_Mat_LT, LT_Loop_cand, LT_Loop] = fun_MH_Candidate_Rob(Loop_Mat_Temp_cand, Loop_Mat_Temp, Threshold_Loop_cand, Threshold_Loop, Struct_Loop_Mat, p, n_loop, liste_nrj_tot, liste_nrj_loop, acceptance_ratio_tot, acceptance_ratio_loop, S_consumer, S_consumed, X, Measured_Abund, StackPlotTot, StackPlot_Meas, beta, theta, theta_T_Loop, theta_LT_Loop, nu, Cov_Mat, Cov_Mat_LT, max_val, nb_rep, LT_Loop_cand, LT_Loop)
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
        liste_nrj_tot(p) = liste_nrj_tot(p - 1);
    end
    
    Struct_Loop_Mat{n_loop + 1} = Loop_Mat_Temp;
    n_loop = n_loop + 1;
    p = p + 1;

    % Modification matrix coefficients
    nu_2 = min(1, 2*nu);
    N = S_consumer*S_consumed;
    Loop_Mat_Temp = reshape(Loop_Mat_Temp,[],1);
    W_n = normrnd(zeros(N,1), ones(N,1));
    prod = Cov_Mat*W_n;
    temp = Loop_Mat_Temp + prod;
    Loop_Mat_Temp_cand = max(temp, 0); %Cov_Mat*W_n only change a fraction of indices
    W_n(temp<0) = prod(temp<0) - temp(temp<0);
    Loop_Mat_Temp = reshape(Loop_Mat_Temp, S_consumer, S_consumed);
    Loop_Mat_Temp_cand = reshape(Loop_Mat_Temp_cand, S_consumer, S_consumed);
    W_n = reshape(W_n, S_consumer, S_consumed);
    High_val = find(Loop_Mat_Temp_cand > max_val);
    Loop_Mat_Temp_cand(High_val) = Loop_Mat_Temp(High_val); %Upper boundary
    W_n(High_val) = 0;
    W_n = reshape(W_n,[],1);
    Cov_Mat = Cov_Mat*(eye(N) + nu_2*(ratio - 0.234)*(W_n*W_n')/norm(W_n)^2)*Cov_Mat';
    Cov_Mat = chol(Cov_Mat,'lower');

    % Modification thresholds
    Threshold_Loop_cand(1) = max(Threshold_Loop(1) + normrnd(0, theta_T_Loop),0);
    logtheta = log(theta_T_Loop) + nu*(ratio - 0.234);
    theta_T_Loop = exp(logtheta);
    
    % Modification lag times
    LT_Loop = reshape(LT_Loop, [], 1);
    W_n = normrnd(zeros(N,1), ones(N,1));
    prod = Cov_Mat_LT*W_n;
    temp = LT_Loop + prod;
    LT_Loop_cand = max(temp, 0); 
    W_n(temp<0) = prod(temp<0) - temp(temp<0);
    LT_Loop = reshape(LT_Loop, S_consumer, S_consumed);
    LT_Loop_cand = reshape(LT_Loop_cand, S_consumer, S_consumed);
    Cov_Mat_LT = Cov_Mat_LT*(eye(N) + nu_2*(ratio - 0.234)*(W_n*W_n')/norm(W_n)^2)*Cov_Mat_LT';
    Cov_Mat_LT = chol(Cov_Mat_LT,'lower');
end