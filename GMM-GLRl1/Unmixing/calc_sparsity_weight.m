function WA = calc_sparsity_weight(A,rows,cols,M)
WA = zeros(size(A));
A1 = reshape(A,rows,cols,M);
B = [0 0.25 0;0.25 0 0.25;0 0.25 0];
for m = 1:M
    temp = conv2(A1(:,:,m), B, 'same');
    WA(:,m) = temp(:);
end
WA = (1./(1e-5+abs(WA)));
end

