function out = threshEpochs(in)
%%THRESHEPOCHS Make epochs from boolean vector
ons = find(diff(in)==1);
if in(1), ons = [1;ons]; end
offs = find(diff(in)==-1);
if in(end), offs = [offs;numel(in)];end
out = [ons offs];
end