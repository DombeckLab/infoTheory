function Ix = selfInformation(x,X,Y)
%SELFINFORMATION The estimated self-information at x for the expSpline
%   defined by X and Y
lambda = expSpline(x,X,Y);

Ix = zeros(size(lambda));

Ix(lambda~=0) = lambda(lambda~=0).*log(lambda(lambda~=0))/log(2);

end
