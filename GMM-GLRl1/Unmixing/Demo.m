clear all;close all;clc
%% -----------------Load data---------------------
dataset = '001';
switch dataset
    case '001'
        load Samson_lwd.mat;
        R = 3;
    case '002'
        load JasperRidge_lwd.mat;
        R = 4;
    case '003'
        load SalinasA_jqw.mat;
        R = 4;
end
beta1 = 0.01;
beta2 = 0.1;
max_num_comp = 4;
max_iter = 200;
A_gt_ori = A_gt;
[Nr, Nc, B_ori] = size(I);
D_ori = 0 * eye(B_ori);
%% -----------------Unmixing-------------------
[A, A_gt, w_jk, mu_jk, sigma_jk] = Unmixing(I, R, A_gt_ori, D_ori, beta1, beta2, max_num_comp, max_iter);
%% ---------------Show results-----------------
[A, ~] = alignresults(R_gt', R_gt', A', R);
A = A';
[error_avg,error_A] = calc_abundance_error(A_gt,A);
show_abundances(A, Nr, Nc, '' ,1 , R);