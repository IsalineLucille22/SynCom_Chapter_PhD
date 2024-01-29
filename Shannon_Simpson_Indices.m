function [Shann, Simp] = Shannon_Simpson_Indices(S,props_vect)
%Compute the Shannon index. S is the total number of species, props_vect
%contains the proportion of each species in the sample.
Shann = props_vect.*log(props_vect);
Shann = - sum(Shann(~isnan(Shann)))/log(S); %Shannon equitability index. Takes values between 0 (no richness) to 1 (high richness)
Simp = sum(props_vect.^2); %simpson index. The probability to randomly choose the same species twice, 0 (high richness) 1 (no richness)
end