function [lambda,d] = expSpline(x,X,Y,varargin)
%EXPSPLINE Evaluates at x the exponentiated spline, normalized by the mean
%
% INPUTS
%   x - The values to evaluate the spline, bounded by [0 1]
%   X - The sorted control points for the spline, bounded by [0 1]
%   Y - The locations for the control points in log space
%
% RETURNS
%   lambda - The values of the exponentiated spline at x
%   d - The normalization factor
global Z D
MAXSTORE = 10;

if ~nargin
    Z = {};
    D = {};
    return
end

ip = inputParser;
ip.addParameter('trackLength',0);
ip.addParameter('trackStart',0);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

if isequal(trackStart,'auto')
    trackStart = min(x);
end

if isequal(trackLength,'auto')
   trackLength = range(x);
end

x = x-trackStart;
if trackLength~=0
   x = x/trackLength; 
end

if any(X(:)<=0|X(:)>=1)
    throw(MException('infoTheory:badControlPoints','Spline X positions out of range'));
end

if numel(Y)-numel(X)~=2
     throw(MException('infoTheory:badControlPoints','There must be two more Y positions than X positions'));
end

Y = Y-max(Y)+3;
[X,i] = sort(X);
Y = Y(:)';X = X(:)';
Y = Y([1 i(:)'+1 numel(X)+2:end]);

% Note:
f=@(x)exp(interp1([0;X(:);1],Y(:),x,'spline'));
if numel(Z)<numel(Y)
   D = cat(2,D,cell(1,numel(Y)-numel(Z)));
   Z = cat(2,Z,cell(1,numel(Y)-numel(Z)));   
end

if isempty(Z{numel(Y)})
    k=[];
else   
    try
        k=find(all(Z{numel(Y)}==[X(:)' Y(:)'],2));
    catch err
        save(['/projects/p30593/' datestr(now,'YYmmDDHHMMSS') '_err.mat']);
        rethrow err;
    end
end

if isempty(k)
%     keyboard
    d=quadgk(f,0,1);    
    
    if size(Z{numel(Y)},1)>=2*MAXSTORE
       k = sort(randsample(1:size(Z{numel(Y)},1),MAXSTORE,false));
       Z{numel(Y)}=Z{numel(Y)}(k,:);
       D{numel(Y)}=D{numel(Y)}(k);
    end
    
    D{numel(Y)} = [D{numel(Y)};d];
    Z{numel(Y)} = [Z{numel(Y)};X(:)' Y(:)'];
   
%     [Z{numel(Y)},i] = sortrows(Z{numel(Y)});
%     D{numel(Y)} = D{numel(Y)}(i);    
elseif numel(k)>1
    error('expSpline:badPreStore','Could not recall partial solutions');    
else
%     keyboard
    d = D{numel(Y)}(k);
end

try
    lambda = reshape(f(x(:))/d,size(x));
catch err
   keyboard; 
end


end

