function [Sample_rate, LN_fit, nb_neg_val] = Rates_Distribution(ind_species, Mat_Struct, num_fig)
%Generate the histrogram and log-normal fit for the rates of the
%interaction matrices of Mat_Struct.
N = length(Mat_Struct);
mat_0 = Mat_Struct{1};
m = length(mat_0(1,:));
Sample_rate = zeros(1, N*m);
for i = 1:N
    mat_temp = Mat_Struct{i};
    Sample_rate((i - 1)*m + 1: i*m) = mat_temp(ind_species,:);
end
non_neg_val = Sample_rate(Sample_rate>0);
nb_neg_val = N*m - length(non_neg_val);
figure(num_fig)
LN_fit = lognfit(non_neg_val);
% histfit(non_neg_val, [], 'lognormal')
% axis([0 20 0 10]);
histogram(Sample_rate, 50)
% [f, x] = hist(Sample_rate);
% bar(x, f / trapz(x, f));
% dim_int = 1.5*max(Sample_rate);
% y = lognpdf(0:dim_int/100:dim_int, LN_fit(1), LN_fit(2));
% plot(0:dim_int/100:dim_int, y)
end