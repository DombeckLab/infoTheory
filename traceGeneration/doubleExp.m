function [ z ] = doubleExp(a,b,t)
%DOUBLEEXP Evaluates the double exponential at t
%   The double exponential is (e^(-a t)-e^(-b t))*c, where c is chosen so
%   the peak of the function is 1.
if ~exist('b','var')% only given t
    t = a;
    % Precalculated values, in hex for full double precision with riseTau =
    % 0.0450 s and fall tau = 0.142 s, as reported for gCamp6f
    a = hex2num('40166d8cd99cf4b0');
    b = hex2num('404c9d940a3e7103');
    
end

t = t(:)';
a = a(:);
b = b(:);

z=(exp(-a*t)-exp(-b*t))./((a./b).^(a./(b-a))-(a./b).^(b./(b-a)));
z(t<=0) = 0;
z(:,a==b)=0;

end

