function [epochs] = longRunningEpochs(xVid,velVid,varargin)
%LONGRUNNINGEPOCHS Return epochs for long running periods
% Defines long-running epochs as in Sheffield 2017.
%
%   REQUIRED INPUT
%       xVid - A t long vector of positions (cm)
%       velVid - A t long vector of speeds (cm/sec)
%
%   RETURNS
%       epochs - A Nx2 matrix of the onset and offset times for each
%       long-running epoch, where N is the number of running epochs
%
%   OPTIONAL PARAMETERS
%       speedThresh (7) - The minimum speed (cm/sec) for inclusion
%       runLength (40) - The minimum continuous run length (cm)
%
%   See	Sheffield, M. E. J., Adoff, M. D. & Dombeck, D. A. Increased 
%       Prevalence of Calcium Transients across the Dendritic Arbor during 
%       Place Field Formation. Neuron 96, 490–504.e5 (2017)

%Parse inputs
ip = inputParser;
ip.addParameter('speedThresh',7);
ip.addParameter('runLength',40);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% Get epochs where the animal was running fast enough
epochs = threshEpochs(velVid>speedThresh);

% Keep the ones where the animal ran far enough
epochs = epochs(diff(xVid(epochs),[],2)>runLength,:);
end

