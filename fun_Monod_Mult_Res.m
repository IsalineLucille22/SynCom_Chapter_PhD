function dz_vect = fun_Monod_Mult_Res(t, z, kappa, CrossFeed_Mat, Resource_Matrix, Threshold, Threshold_pred, Pred_Mat, yield_Pred, S, Lag_time_Cons, Lag_time_Pred, nb_Res)
z(z<0) = 0;
r = z((end - nb_Res + 1):end);
z = z(1:(end - nb_Res));
z = reshape(z',3, S); 
z = z';%[S_i, P_i, W_i]
k = 10;
T = (1./(1 + (Threshold./z(:,3)).^k));
T(isnan(T)) = 0;
k = 10;
T_pred = (1./(1 + (Threshold_pred./z(:,1)).^k));
T_pred(isnan(T_pred)) = 0;
K_S_res = (kappa(:,2) + kappa(:,3))./kappa(:,1).*Resource_Matrix;%(kappa(:,2).*Resource_Matrix + kappa(:,3))./kappa(:,1);%
K_S_res(K_S_res==0) = 1; %To avoid nan, if 0 then multiplied by the zeros of the matrix, so we can put a 1.
K_S_Waste = (kappa(:,2) + kappa(:,3))./kappa(:,1).*CrossFeed_Mat;%(kappa(:,2).*CrossFeed_Mat + kappa(:,3))./kappa(:,1);%
K_S_Waste(K_S_Waste==0) = 1; %To avoid nan, if 0 then multiplied by the zeros of the matrix, so we can put a 1.
K_S_Pred = (Pred_Mat + kappa(:,3))./kappa(:,1);%(kappa(:,2) + kappa(:,3))./kappa(:,1).*Pred_Mat;%Normal %Comment
% K_S_Pred = (kappa(:,2) + kappa(:,3))./kappa(:,1).*Pred_Mat; %Uncomment
K_S_Pred(K_S_Pred==0) = 1; %To avoid nan, if 0 then multiplied by the zeros of the matrix, so we can put a 1.
% dx_vect =  z(:,1).*sum(kappa(:,2).*((t > Lag_time_Cons).*Resource_Matrix).*(r'./(r' + K_S_res)),2)...
%     + z(:,1).*sum(kappa(:,2).*CrossFeed_Mat.*((T.*z(:,3))'./(z(:,3)' + K_S_Waste)),2)...
%     + yield_Pred.*(t > Lag_time_Pred).*z(:,1).*sum(Pred_Mat.*z(:,1)'./(z(:,1)' + K_S_Pred), 2) - sum((t > Lag_time_Pred).*z(:,1).*Pred_Mat.*z(:,1)'./(z(:,1)' + K_S_Pred))';%Species
dx_vect =  z(:,1).*sum(kappa(:,2).*((t > Lag_time_Cons).*Resource_Matrix).*(r'./(r' + K_S_res)),2)...
    + z(:,1).*sum(kappa(:,2).*CrossFeed_Mat.*((T.*z(:,3))'./(z(:,3)' + K_S_Waste)),2)...
    + yield_Pred.*z(:,1).*sum(((t > Lag_time_Pred).*Pred_Mat).*(T_pred.*z(:,1))'./(z(:,1)' + K_S_Pred), 2)...
    - sum(z(:,1).*((t > Lag_time_Pred).*Pred_Mat).*(T_pred.*z(:,1))'./(z(:,1)' + K_S_Pred))';
dy_vect = zeros(S,1); %No complex
dw_vect = z(:,1).*sum(kappa(:,3).*((t > Lag_time_Cons).*Resource_Matrix).*(r'./(r' + K_S_res)),2)...
    + z(:,1).*sum(kappa(:,3).*CrossFeed_Mat.*((T.*z(:,3))'./(z(:,3)' + K_S_Waste)),2)...
    - sum(z(:,1).*(kappa(:,2) + kappa(:,3)).*CrossFeed_Mat.*((T.*z(:,3))'./(z(:,3)' + K_S_Waste)))';%Waste
dr = - sum(z(:,1).*(kappa(:,2) + kappa(:,3)).*((t > Lag_time_Cons).*Resource_Matrix).*(r'./(r' + K_S_res)));
dz_vect = [dx_vect dy_vect dw_vect];
dz_vect = reshape(dz_vect', 1, []);
dz_vect = [dz_vect dr]';
end
