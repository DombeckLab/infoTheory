function [ Fc ] = spk2F(spk,varargin)
%SPK2F Generates a simulated fluorescence trace
%
%   INPUTS
%       spk - A time trace of bined spike counts
%
%   PARAMETERS
%       a - The a parameter for doubleExp. Defaults to gCamp6f kernel.
%       b - The b parameter for doubleExp. Defaults to gCamp6f kernel.
%       riseTau - The rise time for the kernel (sec). Defaults to gCamp6f
%           kernel.
%       fallTau - The fall time for the kernel (sec). Defaults to gCamp6f
%           kernel.
%       Fs (1000) - The sampling frequency for the spike time trace
%       FsVid (30) - The sampling frequency for the video
%       transientHeight (0.19) - The average height of transients
%       transientStd (0.028) - The variance of the height of the transients
%       noiseAmount (0.03) - The size of Gaussian shot noise
%
%   RETURNS
%       Fc - The simulated fluorescence trace

ip = inputParser;
ip.addParameter('a',hex2num('40166d8cd99cf4b0'));
ip.addParameter('b',hex2num('404c9d940a3e7103'));
ip.addParameter('riseTau',[]);
ip.addParameter('fallTau',[]);

ip.addParameter('Fs',1000);
ip.addParameter('FsVid',30);
ip.addParameter('transientHeight',0.19);
ip.addParameter('transientStd',0.028);
ip.addParameter('noiseAmount',0.03);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

% Find the minimum size for a precise kernel
if ~isempty(riseTau)&&~isempty(fallTau)
   [a,b,riseTauHat,fallTauHat] = florescentKernel(riseTau,fallTau); 
else
   [~,~,riseTauHat,fallTauHat] = florescentKernel(a,b,true); 
end

kernelt = fzero(@(t)doubleExp(a,b,exp(t))-1e-6,min(log(riseTauHat+fallTauHat)+log(100),log(300)));
kernelt = ((1:round(exp(kernelt)*Fs))-1)/Fs;% Make time array for kernel
kernel = doubleExp(a,b,kernelt);% Make kernel
kernel = kernel/max(kernel);% Ensure kernel peak is 1
kernel = [zeros(1,numel(kernel)-1) kernel];% Make kernel causal

ts = (0:numel(spk)-1)/Fs;
tsVid = 0:(1/FsVid):max(ts);

Fc = conv(max(...
    normrnd(transientHeight*spk,sqrt(transientStd^2*spk))...Scale heights as sums of normal random variable
    ,0),kernel,'same');% convolve with single spike kernel
Fc = interp1(ts,Fc,tsVid,'nearest');% Subsample (intentionally with aliasing)
Fc = Fc+normrnd(0,noiseAmount,size(Fc));% Add noise



end

