function [L,H] = calc_Laplacians_GL1_1(I, A, eta, R) 
[W,~] = image2graphGL1(I,A,eta,1e-9);  
D = diag(sum(W,2)); 
L = D - W; 
L = sparse(L); 
% KL = L + beta2/beta1 * speye(N);  

W = ones(R,R); 
D = diag(sum(W,2)); 
H = D - W;
end