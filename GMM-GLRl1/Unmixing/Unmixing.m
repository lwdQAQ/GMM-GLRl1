function [A, A_gt, w_jk, mu_jk, sigma_jk] = Unmixing(I, R, A_gt, D_ori, beta1, beta2, max_num_comp, max_iter)
start_t = tic;
%% -----------------PCA--------------------
I_ori = I;
[Nr, Nc, B_ori] = size(I_ori);
N = Nr * Nc;
Y_ori = reshape(I_ori, [N, B_ori]);
[Y, mapping] = pca(Y_ori, 5);
[~, B] = size(Y);
D = mapping.M'*D_ori*mapping.M;
%% -----------------Get GMM parameters---------------------
A_gt_ori = A_gt;
A_gt = reshape(A_gt, [N, R]);
endmembers = cell(1, R);
for j = 1:R            
    X = Y(A_gt(:,j)>0.99, :);            
    endmembers{j} = X;
end
A = [];X = [];
for j = 1:R
    N1 = size(endmembers{j},1);
    X1 = endmembers{j};
    A1 = zeros(N1,R);
    A1(:,j) = 1;
    A = cat(1,A,A1);
    X = cat(1,X,X1);
end
sizes = [N, R]; shrink_size = 0;
[K, w_jk, mu_jk, sigma_jk] = Estimate_num_comp(X, A, sizes, shrink_size, max_num_comp);
%% -----------------Abundance Initialization---------------------
A = calc_A_from_mus(Y, mu_jk);
eta = 0.05;
[L,H] = calc_Laplacians_GL1_1(I, A, eta, R);
WA = calc_sparsity_weight(A, Nr, Nc, R);
%% -----------------EM Algorithm---------------------
K_all = K2K_all(K);
w_k = w_jk2w_k(w_jk,K_all);
%%% learning rate initialization
delta_t0 = 1e-12;
delta_t_A = delta_t0;
%%% grad initialization

eval_totals = [];

fix_A = 0;

params = [];
params.mu_jk = mu_jk;
params.sigma_jk = sigma_jk;
params.K_all = K_all;
params.w_k = w_k;
params.D = D;
for iter = 1:max_iter
    %%% E Step
    tstart = tic;
    N_nk = calc_log_gaussians(A, D, mu_jk, sigma_jk, K_all, Y);
    gamma_nk = calc_gamma_E_step(N_nk, w_k);
    params.gamma_nk = gamma_nk;
    params.L = L;
    params.H = H;
    params.WA = WA;
    params.beta1 = beta1;
    params.beta2 = beta2;
    %%% M step
    % update A, L, KL, H
    if ~fix_A
        der_A = calc_der_A(A, Y, params);
        delta_t_A = calc_time_step_adaptive(@eval_obj_fun_A, @update_A, ...
            A, params, der_A, delta_t_A, delta_t0, Y);
        A = update_A(A, der_A, delta_t_A);
        [L,H] = calc_Laplacians_GL1_1(I, A, eta, R);
        WA = calc_sparsity_weight(A, Nr, Nc, R);
        WA = ones(size(A));
    end
    
    telapsed = toc(tstart);
    disp(['EM iteration ', num2str(iter), ' lasts ', num2str(telapsed)]);
end
w_jk = w_k2w_jk(w_k, K_all);
[mu_jk,sigma_jk] = restore_from_projection(mu_jk,sigma_jk,[],mapping.mean,mapping.M);

elapsed_t = toc(start_t);
disp(['Total algorithm execution time is ',num2str(elapsed_t)]);

function eval_total = calc_obj_fun(A, Y, params)
if isstruct(params)
    arg_set = fieldnames(params);
    for i = 1:length(arg_set)
        eval([arg_set{i},'=params.',arg_set{i},';']);
    end
end
[N,~] = size(Y);
K = length(w_k);
N_nk = calc_log_gaussians(A, D, mu_jk, sigma_jk, K_all, Y);
% Need to avoid too large negative logarithm
N_nk_max = max(N_nk,[],2);
temp = log(sum(repmat(w_k,N,1) .* exp(N_nk - repmat(N_nk_max,1,K)), 2));
value = N_nk_max + temp;

eval_N = -sum(value);
eval_A = (beta1/2) * trace(A'*L*A) + (beta2/2) * norm_p(WA.*A,0.5);

eval_total = eval_N + eval_A;

function N_nk = calc_log_gaussians(A, D, mu_jk, sigma_jk, K_all, Y)
[N,B] = size(Y);
K1 = size(K_all,1);

[mu_nk,sigma_nk] = calc_mu_sigma_nk(A, D, mu_jk, sigma_jk, K_all);

N_nk = zeros(N,K1);
for k = 1:K1
    N_nk(:,k) = mvnpdf(Y, mu_nk(:,:,k), sigma_nk(:,:,:,k));
end
N_nk = log(N_nk + 0.0001);

function [mu_nk,sigma_nk] = calc_mu_sigma_nk(A, D, mu_jk, sigma_jk, K_all)
[N,M] = size(A);
K1 = size(K_all,1);
B = size(mu_jk{1},2);

[mu_all,sigma_all] = calc_mu_sigma_all(mu_jk, sigma_jk, K_all);

sigma2_I = D;

mu_nk = zeros(N,B,K1);
sigma_nk = zeros(B*B,N,K1);
for k = 1:K1
    mu_nk(:,:,k) = A * mu_all(:,:,k)';
    sigma_nk(:,:,k) = sigma_all(:,:,k) * (A.^2)' + repmat(sigma2_I(:),1,N);
end

sigma_nk = reshape(sigma_nk,[B,B,N,K1]);

function [mu_all,sigma_all] = calc_mu_sigma_all(mu_jk, sigma_jk, K_all)
mu_all = calc_mu_all(mu_jk, K_all);
sigma_all = calc_sigma_all(sigma_jk, K_all);

function w_k = w_jk2w_k(w_jk, K_all)
w_k = ones(1, size(K_all,1));
for i = 1:size(K_all,1)
    k = K_all(i,:);
    for j = 1:length(w_jk)
        w_k(i) = w_k(i)*w_jk{j}(k(j));
    end
end

function w_jk = w_k2w_jk(w_k, K_all)
M = size(K_all,2);
w_jk = cell(1,M);
K = K_all2K(K_all);

for j = 1:M
    w_jk{j} = zeros(1,K(j));
end

for j = 1:M
    for l = 1:K(j)
        w = w_k(K_all(:,j)==l);
        w_jk{j}(l) = sum(w);
    end
end

function der_A = calc_der_A(A, Y, params)
if isstruct(params)
    arg_set = fieldnames(params);
    for i = 1:length(arg_set)
        eval([arg_set{i},'=params.',arg_set{i},';']);
    end
end
[N,K1] = size(gamma_nk);
M = length(mu_jk);
B = size(mu_jk{1},2);

[lambda_nk,psi_nk] = calc_lambda_psi_nk(A, D, mu_jk, sigma_jk, K_all, gamma_nk, Y);

[mu_all,sigma_all] = calc_mu_sigma_all(mu_jk, sigma_jk, K_all);
psi_k = reshape(psi_nk, [B*B,N,K1]);

der_A = zeros(N,M);
temp1 = zeros(N,M);
for k = 1:K1
    der_A = der_A - lambda_nk(:,:,k) * mu_all(:,:,k);
    temp1 = temp1 - psi_k(:,:,k)' * sigma_all(:,:,k);
end
der_A = der_A + 2 * A .* temp1;

% der_A = der_A + beta1 * KL * A;
der_A = der_A + beta1 * L * A + beta2/4 * WA .* (A.^(-0.5));

function [lambda_nk,psi_nk] = calc_lambda_psi_nk(A, D, mu_jk, ...
    sigma_jk, K_all, gamma_nk, Y)
[mu_nk,sigma_nk,prec] = calc_mu_sigma_prec(A, D, mu_jk, sigma_jk, K_all);

[N,B,K1] = size(mu_nk);

lambda_nk = zeros(N,B,K1);
psi_nk = zeros(B,B,N,K1);

for k = 1:K1
    Y1 = Y - mu_nk(:,:,k);
    sigma_y_n_mu_nk = multiprod(prec(:,:,:,k), Y1', [1 2], [1]);
    lambda_nk(:,:,k) = (repmat(gamma_nk(:,k)', B, 1) .* sigma_y_n_mu_nk)';
    
    tmp1 = reshape(sigma_y_n_mu_nk, [B 1 N]);
    tmp2 = reshape(sigma_y_n_mu_nk, [1 B N]);
    tmp3 = multiprod(tmp1, tmp2, [1 2], [1 2]) - prec(:,:,:,k);
    psi_nk(:,:,:,k) = 0.5 * multiprod(reshape(gamma_nk(:,k), [1 1 N]), tmp3);
end

function [mu_nk,sigma_nk,prec] = calc_mu_sigma_prec(A, D, mu_jk, sigma_jk, K_all)
[mu_nk,sigma_nk] = calc_mu_sigma_nk(A, D, mu_jk, sigma_jk, K_all);

% prec1 = multinv(sigma_nk);
% To avoid out of memory. The running time is similar to the above single
% statement.
prec = zeros(size(sigma_nk));
for k = 1:size(prec,4)
% parfor k = 1:size(prec,4)
    prec(:,:,:,k) = multinv(sigma_nk(:,:,:,k));
end

function A_new = update_A(A, der_A, delta_t)
A_new = A - delta_t * der_A;
A_new = max(1e-4,A_new);
A_new = project_to_simplex(A_new);

function val = eval_obj_fun_A(A, Y, params)
if isstruct(params)
    arg_set = fieldnames(params);
    for i = 1:length(arg_set)
        eval([arg_set{i},'=params.',arg_set{i},';']);
    end
end
N_nk = calc_log_gaussians(A, D, mu_jk, sigma_jk, K_all, Y);
val = -sum(sum(0.5 * gamma_nk .* N_nk)) + beta1/2 * trace(A'*L*A) + beta2/2 * norm_p(WA.*A, 0.5);