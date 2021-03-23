function [cel,params] = randomCell(varargin)
%RANDOMCELL Generate physiology from modeled neural responses to behavior
%
% RETURNS
%   cel - A struct contianing
%       spkTs: The times of spikes
%       spkVid: The spike train binned by the fluorescence timestamps
%       Fc: The fluorescence trace
%   params - A struct containing all the optional parameters for the cell
%
% OPTIONAL PARAMETERS
%   - BEHAVIOR PARAMETERS -
%   x ([]) - The behavior, sampled at Fs. If empty, randomly generates
%       using the parameters below.
%   time (rand()*60) - the duration of the session in minutes
%   behaviorIs ([]) - The indices of the behaviors chosen. If left empty,
%       uses
%   behaviorPath ([]) - The location of the behavior file. If empty, is
%       chosen via a user interface
%
%   - MAP PARAMETERS -
%   spatialInformationRateTarget (rand()*6) - Spatial information rate
%       (bits/AP)
%   X ([]) - X control points for spatial map. If empty, determined
%       from spatialInformationRateTarget
%   Y ([]) - Y control points for spatial map. If empty, determined
%       from spatialInformationRateTarget
%
%   - FIRING RATE PARAMETERS -
%   meanRate (exprnd(1)) - The expected firing rate (Hz)
%   goodEpochs ([]) - If supplied, used to normalize the expected mean
%       rate. Can be a [n x 2] matrix of start (left) and stop (right)
%       frames, or the string 'longRuns' which uses the long running
%       epochs. If empty uses entire session.
%
%   - ΔF/F parameters -
%   riseTau (45e-3) - The rise time for the kernel (sec)
%   fallTau (142e-3) - The 1/2 fall time for the kernel
%       (sec)
%   transientHeight (0.19) - The average height of each AP-evoked
%       transient (ΔF/F)
%   transientStd (0.028) - The standard deviation of the transient height
%       (ΔF/F)
%   noiseAmount (0.03) - The amplitude of the florescence noise
%   (ΔF/F)
%   params.Fs (1e3) - The sample frequency (Hz)
%   params.FsVid (30) - The video sample frequency (Hz)
%
% Jason Climer, jrclimer@northwestern.edu
params.rng = rng;% Save rng status
warning('off','MATLAB:quadgk:NonFiniteValue')
% Parse inputs
ip = inputParser;
ip.CaseSensitive = 1;
ip.addParameter('x',[]);
ip.addParameter('time',rand*60);
ip.addParameter('behaviorIs',[]);
ip.addParameter('behaviorPath',[]);

ip.addParameter('spatialInformationRateTarget',rand()*6);
ip.addParameter('X',[]);
ip.addParameter('Y',[]);

ip.addParameter('meanRate',unifrnd(0.1,30));
ip.addParameter('goodEpochs','longRuns');

ip.addParameter('riseTau',45e-3);
ip.addParameter('fallTau',142e-3);
ip.addParameter('transientHeight', 0.19);%19+-2.8 single ap, Chen et al, 2014
ip.addParameter('transientStd',0.028);
ip.addParameter('noiseAmount', 0.03);

ip.addParameter('Fs',1e3);
ip.addParameter('FsVid',30);
ip.addParameter('trackLength',300);

ip.addParameter('verbose',false);

ip.parse(varargin{:});
for j=fields(ip.Results)'
    params.(j{1})=ip.Results.(j{1});
end

verbose = ip.Results.verbose;

if verbose, fprintf('%s validating parameters...\n',datestr(now)); end

if isempty(params.X)||isempty(params.Y)% spatial map not made
    if verbose, fprintf('%s Making spatial map...\n',datestr(now)); end
    % Make spatial map
    [params.X,params.Y,params.spatialIhat] = genExpSpline(params.spatialInformationRateTarget);    
elseif ~ismember('spatialIhat',fields(params))||isempty(params.spatialIhat)
    if verbose, fprintf('%s Using user defined spatial map...\n',datestr(now)); end
    
    if ~ismember('spatialIhat',fields(params))||isempty(params.spatialIhat)
    % Store spatialIhat if not already set
        params.spatialIhat = quadgk(@(x)selfInformation(x,params.X,params.Y),0,1);
    end
end

if verbose, fprintf('%s Formatting behavior...\n',datestr(now)); end
% Format the behavior & build smoothed speed
if isempty(params.x)
    [x,params.behaviorIs] = loadBehaviorT(params.time,'behaviorIs',params.behaviorIs,'behaviorPath',params.behaviorPath);
else
   x = params.x; 
end

[x,xVid,speed,ts,tsVid,~] = formatBehavior( x , 'Fs', params.Fs, 'FsVid', params.FsVid, 'trackLength',params.trackLength);
speed(speed<-10) = 0;
speed = interp1(tsVid,speed,ts,'pchip');

if verbose, fprintf('%s Building CIFs...\n',datestr(now)); end

if verbose, fprintf('%s Building spike train...\n',datestr(now)); end


[spk,cel.spkTs] = makeSpikeTrain(x,params.X,params.Y...
    ,'trackLength',params.trackLength...
    ,'goodEpochs',params.goodEpochs...
    ,'Fs',params.Fs...
    ,'FsVid',params.FsVid...
    ,'meanRate',params.meanRate...
    );

% Bin spike times by the video
cel.spkVid = histc(cel.spkTs,tsVid);

% Solve for kernel
[params.a,params.b,params.riseTauHat,params.fallTauHat] = florescentKernel(params.riseTau,params.fallTau);

% Make fluorescent trace
cel.Fc = spk2F(spk...
    ,'a',params.a...
    ,'b',params.b...
    ,'Fs',params.Fs...
    ,'FsVid',params.FsVid...
    ,'transientHeight',params.transientHeight...
    ,'transientStd',params.transientStd...
    ,'noiseAmount',params.noiseAmount...
    );

if verbose, fprintf('%s Cell created!\n',datestr(now)); end

end

