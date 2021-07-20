function [ x,behaviorIs ] = loadBehaviorT( time,varargin )
%LOADBEHAIORT Loads behavior files
%
% INPUTS
%   time - The duration of the session in minutes
%
% OPTIONAL PARAMETERS
%   behaviorIs ([]) - The indices of the behaviors chosen. If left empty,
%       uses
%   behaviorPath ([]) - The location of the behavior file. If empty, is
%       chosen via a user interface
%
% RETURNS
% x - The position of the animal, sampled at 1 kH
% behaviorIs - The indices of the behavior, used to recreate the behavior
% global behavior;
% global nBehavior;

ip = inputParser;
ip.addParameter('behaviorIs',[]);
ip.addParameter('behaviorPath',[]);

ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

getBehaviorIs = isempty(behaviorIs);

% if isempty(behavior)
if isempty(behaviorPath)
    [behavior{1},behavior{2}]=uigetfile('*.mat','Pick behavior file');
    behaviorPath = [behavior{2} behavior{1}];
end
behavior = matfile(behaviorPath);

if ~isequal(behaviorPath(end-3:end),'.mat')
    throw(MException('infoTheory:badPath','behaviorPath must target a .mat file'));
end
nBehavior = prod(size(behavior,'x'));
% end

if ~exist('time','var')
    x = [];
    behaviorIs = [];
    return;
end

x = [];

i = 0;
while numel(x)/1e3/60<time
    i = i+1;
    if getBehaviorIs
        behaviorIs = [behaviorIs;randi(nBehavior)];
    end
    
    z = behavior.x(1,behaviorIs(i));x=[x;z{1}];
end
x = x(1:floor(time*60*1e3));

end

