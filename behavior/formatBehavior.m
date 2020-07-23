function [x,xVid,vel,ts,tsVid,rewardTimes] = formatBehavior( x,varargin )
%FORMATBEHAVIOR Summary of this function goes here
%   Detailed explanation goes here
ip = inputParser;
ip.addParameter('Fs',1e3);
ip.addParameter('FsVid',30);
ip.addParameter('trackLength',300);
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

x = (x-min(x))/range(x)*trackLength;

ts = (0:numel(x)-1)'/Fs;
tsVid = (0:1/FsVid:max(ts))';tsVid=tsVid(1:end-1);
[A,B] = butter(6,FsVid/Fs,'low');

for i=1:2
rewardTimes = find(diff(x)<-20/300*trackLength);
temp=[0;rewardTimes;numel(ts)];
xVid=zeros(size(tsVid));
for j=1:numel(temp)-1
    if temp(j+1)-temp(j)-temp(j)>18
    xVid(round(ts(temp(j)+1)*FsVid+1):min(round(ts(temp(j+1))*FsVid+1),numel(tsVid)))=...
        interp1(ts(temp(j)+1:temp(j+1))...
        ,filtfilt(A,B,x(temp(j)+1:temp(j+1),1))...
        ,tsVid(round(ts(temp(j)+1)*FsVid+1):min(round(ts(temp(j+1))*FsVid+1),numel(tsVid))),'pchip');
    else
        kernel = ts(temp(j)+1:temp(j+1));
        kernel = kernel-mean(kernel);
        kernel = normpdf(kernel,0,1/FsVid);
        kernel = kernel/sum(kernel);
        kernel = conv([repmat(x(temp(j)+1),[100 1]);x(temp(j)+1:temp(j+1),1);repmat(x(temp(j+1)),[100 1])],kernel,'same');
        kernel = kernel(101:end-100);
        xVid(round(ts(temp(j)+1)*FsVid+1):min(round(ts(temp(j+1))*FsVid+1),numel(tsVid)))=...
             interp1(ts(temp(j)+1:temp(j+1))...
             ,kernel....
             ,tsVid(round(ts(temp(j)+1)*FsVid+1):min(round(ts(temp(j+1))*FsVid+1),numel(tsVid))),'pchip');
    end
end

if i==1
   x=min(max((x-min(xVid))/range(xVid)*trackLength,0),trackLength);
end
end

vel = [0;diff(xVid)]*FsVid;
vel(vel<-100)=0;

end

