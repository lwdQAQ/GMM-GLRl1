function [s_est, A_est] = alignresults(M, A_est1, s_est1, p)

CRD = corrcoef([M A_est1]);
D = abs(CRD(p+1:2*p,1:p));
% permute results
perm_mtx = zeros(p,p);
aux=zeros(p,1);
for i=1:p
    [ld cd]=find(max(D(:))==D);
    ld=ld(1);cd=cd(1); % in the case os more than one maximum
    perm_mtx(ld,cd)=1;
    D(:,cd)=aux; D(ld,:)=aux';
end
A_est = A_est1 * perm_mtx;
s_est = s_est1' * perm_mtx;
s_est = s_est';

