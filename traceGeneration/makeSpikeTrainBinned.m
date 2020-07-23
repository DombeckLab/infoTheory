function [spk,spkTs] = makeSpikeTrainBinned(x,z,varargin)
%MAKESPIKETRAINBINNED Generates a spike train from a 2D distribution
%
% REQUIRED INPUTS
%   x - The position of the animal
%   z - The log probability map
%
% OPTIONAL PARAMETERS
%   trackLength (300 cm) - The length of the track
%   Fs (1e3 Hz) - The sample frequency
%   FsVid (30 Hz) - The sample frequency of the output
%
% RETURNS
%   spk - The spike histogram, sampled at Fs
%   spkTs - The spike times
ip = inputParser();
ip.addParameter('trackLength',300);
ip.addParameter('Fs',1e3);
ip.addParameter('FsVid',30);
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

n = size(z,2);
m = size(z,1);

[x,xVid,speed,ts,tsVid,~] = formatBehavior( x , 'Fs', Fs, 'FsVid', FsVid, 'trackLength',trackLength);
speed(speed<-10) = 0;
speedVid = speed;

[oc,xi] = histc(x,linspace(0,300,n+1));
xi(xi==n+1)=n;

P = cumsum(exp(z+log(n)));
spk = rand(size(xi));
spk = sum(P(:,xi)<=spk');
spkTs = arrayfun(@(i)find(spk>=i),1:max(spk),'UniformOutput',false);
spkTs = cat(2,spkTs{:})/Fs;


