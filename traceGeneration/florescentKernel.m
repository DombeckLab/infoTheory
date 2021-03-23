function [ a,b,riseTauHat,fallTauHat] = fluorescentKernel( riseTau,fallTau,asab )
%FLUORESCENTKERNEL Solve for the double-exponential based on rise and fall times
%
% INPUT
%   riseTau - The time to the peak of the kernel (sec)
%   fallTau - The half fall time of the kernel (sec)
%
% OPTIONAL INPUT
%   asab (false) - If set to true, the syntax is
%   flurorescentKenrel(a,b,true) 
%
% RETURNS
%   a - The a parameter for doubleExp
%   b - The b parameter for doubleExp
%   riseTauHat - The numerically solved resulting rise time
%   fallTauHat - The numerically solved resulting fall time
%
% See also: doubleExp
if ~exist('riseTau','var')
    riseTau = 45e-3;
    fallTau = 142e-3;
    
    % Precalculated values, in hex for full double precision
    a = hex2num('40166d8cd99cf4b0');
    b = hex2num('404c9d940a3e7103');
    riseTauHat = hex2num('3fa70a3d70a3d70b');
    fallTauHat = hex2num('3fc22d0e56041893');
    
    return
end
if ~exist('asab','var'), asab = false; end

if asab
    a = riseTau;
    b = fallTau;
    riseTauHat = (log(a)-log(b))/(a-b);    
    fallTauHat = exp(fzero(@(t)doubleExp(a,b,riseTauHat+exp(t))-0.5,-10,optimset('display','none')));
 else
    b = @(a)-lambertw(-1,-a*exp(-a*riseTau)*riseTau)/riseTau;% Theorhetical solution for b given a
    l = @(a)doubleExp(a,b(a),riseTau+fallTau)-0.5;% The difference between the theorhetical and target maps
    
    [a,~,exitflag]=fzero(l,1/riseTau/2,optimset('display','none'));% Find a theorhetical a
    if exitflag~=1% Could not solve        
        try
            [a,~,exitflag]=fzero(l,2/riseTau,optimset('display','none'));% Try a higher starting point
        catch err
            exitflag = -10;
        end
        
        if exitflag~=1% Could not solve            
            a=fminsearch(@(a)sum((doubleExp(exp(a(:,1)),exp(a(:,2)),[riseTau riseTau+fallTau])-[1 0.5]).^2,2)...
                ,log([1/2 2]/riseTau),optimset('display','none'));% Minimize the squared residuals
            b=exp(a(2));
            a=exp(a(1));
            if ~isnan(a)&&~isnan(b)&&~isinf(a)&&~isinf(b)&&a~=b
               exitflag=1; 
            end
        end
    end
    
    if exitflag==1&&(~isnumeric(b)||b~=a)% We found it
        if ~isnumeric(b), b = b(a);end
        riseTauHat = (log(a)-log(b))/(a-b);   
        fallTauHat = fzero(@(t)doubleExp(a,b,riseTauHat+t)-0.5,4*fallTau,optimset('display','none'));
        if fallTauHat<0
            fallTauHat = exp(fzero(@(t)doubleExp(a,b,riseTauHat+exp(t))-0.5,4+log(fallTauHat),optimset('display','none')));
        end
    else
        a = NaN;
        b = NaN;
        riseTauHat = NaN;
        fallTauHat = NaN;
    end
end

if isnan(a)||isnan(b)||isnan(riseTauHat)||isnan(fallTauHat)
   a = nan;b=nan;riseTauHat=nan;fallTauHat=nan; 
end
