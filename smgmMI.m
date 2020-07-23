function [I,lambda_] = smgmMI(X,Y,varargin)
%SMGMMI The Skaggs, McNaughton Gothard mutual information estimate
%   [I,lambda_] = SMGMMIRATE(X,Y,...) Returns the Skaggs, McNaughton, 
%   Gothard and Markus information rate of Y across X. X and Y are N long 
%   vectors, where N is the number of samples. I is the information rate in 
%   bits per sample. lambda_ is the average Y per sample.
%
%   Optional Parameters
%   x_min (min(X)) - The lower bounds for the binning of X
%   x_max (max(X)) - The upper bounds for the binning of X
%   n_bin (60) - The number of bins across X
%
% See Skaggs, W.E., McNaughton, B.L., and Gothard, K.M. (1993). An 
%   Information-Theoretic Approach to Deciphering the Hippocampal Code. In 
%   Advances in Neural Information Processing Systems 5, S.J. Hanson, J.D. 
%   Cowan, and C.L. Giles, eds. (Morgan-Kaufmann), pp. 1030–1037.
%
% Written by Jason Climer, PhD (jrclimer@northwestern.edu)

% Parse input
if ~exist('Y','var')
   if iscell(X)
       Y = X{2};
       X = X{1};
   else
    Y = X(:,2);
    X = X(:,1);
   end
end

if numel(X)~=numel(Y)
   throw(MException('infoTheory:sizeMismatch','X and Y must be the same size')); 
end

ip = inputParser;
ip.addParameter('x_min',min(X));
ip.addParameter('x_max',max(X));
ip.addParameter('n_bin',60);
ip.parse(varargin{:});
for j=fields(ip.Results)'
   eval(sprintf('%s=ip.Results.%s;',j{1},j{1})); 
end

% Make maps
[oc,~,bin] = histcounts(X,linspace(x_min,x_max,n_bin+1));
lambdai = sum((bin==1:n_bin).*Y)./oc;
oc = oc/sum(oc);% Occupancy
lambda_ = nansum(lambdai.*oc);% Mean value
I = nansum(lambdai.*oc.*(log(lambdai)-log(lambda_)))/log(2);% Information

end

