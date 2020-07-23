function [X,Y,Ihat] = genExpSpline(I,varargin)
%GENEXPSPLINE Generate a map that matches a SMGM spike information
%
%   REQUIRED INPUT
%       I - The target bit-rate (bits/spike)
%
%   OPTIONAL INPUT
%       Z0 ([]) - The initial guess for [X Y]. If empty is random.
%       npoints (5) - The number of control points.
%
%   RETURNS
%       X - The control points across the space
%       Y - The control points in rate
%       Ihat - The actual theorhetical value for the map defined by X and Y
%
% See also: expSpline

ip = inputParser;
ip.addParameter('Z0',[]);
ip.addParameter('npoints',5);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.(j{1});']);
end

if ~ismember('Z0',ip.UsingDefaults)    
    npoints = numel(Z0)/2+1;
end

if I<=0
   throw(MException('infoTheory:badInfo','Information target must be >=0')); 
end

% The information in the map defined by Z, where Z = [X Y]
inffun = @(Z)quadgk(@(x)selfInformation(x,Z(1:numel(Z)/2-1),Z(numel(Z)/2:end)),0,1);

% The squared residual between the target information and the information
%   in the map defined by Z
fun = @(Z)(inffun(Z)-I).^2;

if ~ismember('Z0',ip.UsingDefaults)    
    e = fun(Z0);
    npoints = npoints-1;
else
    e = inf;
end


while e>1e-10% Cannot converge with this few points
    npoints=npoints+1;% Increment the number of points
    
    if ~ismember('Z0',ip.UsingDefaults)
        Z0 = ip.Results.Z0;
    else
        e = inf;
        Z0 = NaN(1,2*npoints-2);
        while any([...
                isnan([e;Z0(:)])...
                ;isinf([e;Z0(:)])...
                ])
            Z0 = [sort(rand(1,npoints-2)) normrnd(0,1,1,npoints)];% Randomly select a starting point
            e = fun(Z0);% Ensure the starting point is valis
        end
    end
    
    [Z,e] = fmincon(...
        fun ...% Minimize the squared residual
        ,Z0...% Starting point
        ,[[eye(npoints-3) zeros(npoints-3,1)]-[zeros(npoints-3,1) eye(npoints-3)],zeros(npoints-3,npoints)],zeros(npoints-3,1)...The x points must be in order
        ,[],[]...
        ,[zeros(1,npoints-2) -inf(1,npoints)]...X points are bounded on [0 1], y points are unbounded
        ,[ones(1,npoints-2) inf(1,npoints)]...
        ,[],optimoptions(@fmincon,'OptimalityTolerance',0,'display','none','PlotFcn',[])...@fig1plotfun)...
        );
    
    if any(~ismember({'Z0','npoints'},ip.UsingDefaults))&&e>1e-10
        warning('Could not fit with given inputs');
        [X,Y,Ihat] = genExpSpline(I);
        e = (I-Ihat)^2;
    end
end

% Format output
Ihat = inffun(Z);
X = Z(1:npoints-2);Y=Z(npoints-1:end);

end
