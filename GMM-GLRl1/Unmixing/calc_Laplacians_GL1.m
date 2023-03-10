function [L,H] = calc_Laplacians_GL1(I, A, eta, R) 
[Nr, Nc, B] = size(I);
N = Nr*Nc;
Y = reshape(I, [N, B]);
epsilon = 0.001;
W = zeros(N);
for i = 1:N
    for j = 1:N
        Nci = ceil(i/Nr); Nri = i - (Nci - 1)*Nr;
        Ncj = ceil(i/Nr); Nrj = j - (Ncj - 1)*Nr;
        if abs(Nci - Ncj)==1||abs(Nri - Nrj)==1
            W(i,j) = exp(-sum((Y(i,:)-Y(j,:)).^2)/(2*B*eta));
            W(i,j) = W(i,j)*(sum(abs(A(i,:)-A(j,:)))+sqrt(R)*epsilon)/...
                (sum((A(i,:)-A(j,:)).^2)+epsilon^2);
        end
    end
end

D = diag(sum(W, 2));
L = D - W;
% L = sparse(L);

W = ones(R, R);
D = diag(sum(W, 2));
H = D - W;
% [rows,cols,B] = size(I); N = rows * cols;  
% [W,Neighbors] = image2graphGL1(I,A,eta,1e-9);  
% D = diag(sum(W,2)); 
% L = D - W; 
% L = sparse(L); 
% % KL = L + beta2/beta1 * speye(N);  
% 
% W = ones(M,M); 
% D = diag(sum(W,2)); 
% H = D - W;
end