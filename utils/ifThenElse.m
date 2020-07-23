function out = ifThenElse(cond,do,other)
% IFTHENELSE Conditional programatic evaluation
%
% INPUTS
%    cond - The condition, true of false
%    do - A function handle of what to run if true
% 
% OPTIONAL INPUTS
%    other (@()false) - A function handle to run if false
%
% Jason Climer (jrclimer@northwestern.edu)
if ~exist('other','var'), other = @()false; end

if cond
    out = do();
else
    out = other();
end
end