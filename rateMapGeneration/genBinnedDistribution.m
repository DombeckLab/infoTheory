function [z,Ihat] = genBinnedDistribution(I0,lambda_,varargin)
%GENBINNEDDISTRIBUTION Generate a binned probablility map that matches a spike information rate
%
%   REQUIRED INPUT
%       I0 - The target bit-rate (bits/second)
%       lambda_ - The mean firing rate assuming uniform occupancy
%
%   OPTIONAL INPUT
%       Z0 ([]) - The initial guess for [X Y]. If empty is random.
%       npoints (5) - The number of control points for the continuous initial guess.
%       n (60) - The number of spatial bins
%       m (10) - The maximim number of spikes per frame
%       Fs (1e3) - The sampling rate
%
%   RETURNS
%       z - The log probability distribution
%       Ihat - The actual theorhetical value for the map defined by X and Y
%
% See also: expSpline
ip = inputParser;
ip.addParameter('Z0',[]);
ip.addParameter('npoints',5);
ip.addParameter('n',60);
ip.addParameter('m',10);
ip.addParameter('Fs',1e3);
ip.parse(varargin{:});
for j=fields(ip.Results)'
   eval([j{1} '=ip.Results.' j{1} ';']); 
end


I = I0/lambda_;
[X,Y,~] = genExpSpline(I,'npoints',5);
lambda = expSpline(linspace(0,1,n),X,Y);
lambda = lambda/mean(lambda)*lambda_;

p = poisspdf(repmat((0:m)',[1 n]),repmat(lambda/Fs,[m+1 1]));
z = log(p);
p = @(z)exp(z)./sum(exp(z),1)/n;
rshp = @(z)reshape(z,m+1,n);
lambdahat = @(p)sum(sum(p.*(0:m)'))*Fs
Ihat = @(p)sum(sum(p.*log2(p./(sum(p,1).*sum(p,2)))))*Fs
fun = @(z)sum([I-Ihat(p(rshp(z))) lambda_-lambdahat(p(rshp(z)))].^2);


zhat = fminunc(fun,z(:),optimoptions(@fminunc,'display','none','MaxFunctionEvaluations',inf));
Ihat = Ihat(p(rshp(zhat)));
z = rshp(zhat);
z = z-log(sum(exp(z),1))-log(n);


end

