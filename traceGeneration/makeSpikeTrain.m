function [spk,spkTs] = makeSpikeTrain(x,X,Y,varargin)
%MAKESPIKETRAIN Generates a spike train as an inhomogenous Poisson process
%
% REQUIRED INPUTS
%   x - The position of the animal, sampled at Fs
%   X - The x position of nodes in the rate map
%   Y - The y position of nodes in the rate map
%
% OPTIONAL PARAMETERS
%   trackLength (300 cm) - The length of the track
%   goodEpochs ('longRuns') - Can be a n x 2 matrix of start and stop
%       times. If set to the string 'longRuns' finds the long running
%       epochs.
%   Fs (1e3 Hz) - The sample frequency
%   FsVid (30 Hz) - The sample frequency of the output
%   meanRate (1 Hz) - The mean firing rate during the good epochs
%
% RETURNS
%   spk - The spike counts sampled at Fs
%   spkTs - The times of the APs
ip = inputParser();
ip.addParameter('trackLength',300);
ip.addParameter('goodEpochs','longRuns');
ip.addParameter('Fs',1e3);
ip.addParameter('FsVid',30);
ip.addParameter('meanRate',1);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

[x,xVid,speed,ts,tsVid,~] = formatBehavior( x , 'Fs', Fs, 'FsVid', FsVid, 'trackLength',trackLength);
speed(speed<-10) = 0;
speedVid = speed;
speed = interp1(tsVid,speed,ts,'pchip');

lambda = expSpline(x/trackLength,X,Y);

if isequal(goodEpochs,'longRuns')
    goodEpochs = longRunningEpochs(xVid,speedVid);
end

if isempty(goodEpochs)
    lambda = lambda./mean(lambda);% Make means 1 across whole session
else% Make means 1 during supplied epochs
    lambda = lambda./mean(lambda(any(1:size(lambda,1)>=goodEpochs(:,1)&1:size(lambda,1)<=goodEpochs(:,2),1)',:));
end

spk = poissrnd(lambda*meanRate/Fs);% Generate spikes

% Find spike times
spkTs = arrayfun(@(i)find(spk>=i),1:max(spk),'UniformOutput',false);
spkTs = cat(1,spkTs{:});
spkTs = sort(spkTs+rand(numel(spkTs),1))/Fs;
