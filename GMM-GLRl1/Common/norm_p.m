function result = norm_p(X,p)
% p > 0
if ndims(X) == 1
    result = sum(abs(X).^p);
else
    result = sum(sum(abs(X).^p));
end
end