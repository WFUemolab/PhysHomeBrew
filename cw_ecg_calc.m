function [stat, raw] = cw_ecg_calc(dat, conds, sr, PHYS)
% function stat = cw_ecg_calc(dat, conds, sr, trig)
%
% This function calculates the initial IBI deceleration and subsequent IBI 
% acceleration from some pre-stimulus baseline.
%
% dat is a (m x n) matrix with m trials and n timepts per trial
% conds is a (m x 1) matrix with the condition identifiers for each trial
% sr is the sampling rate (e.g. 20 hz)
% trig is a struct:
%   trig.pre = pre-stimulus baseline time (in seconds)
%   trig.trig = trig time (how long the trigger occupies)
%   trig.post = post-stimulus time
%
% output is a (1 x k*2) stat matrix with k numbers of conditions and 2
% stats per condition (IBI deceleration and acceleration, respectively)


respWindow = PHYS.EMG.responseWindow * sr; %the window(s) within which we are going to average EMG
baseWindow = PHYS.EMG.baseWindow * sr; %the window within which we are going to get the proximal baseline
trialstart = PHYS.EMG.trialstart * sr +1;
datcrit = zeros(size(dat,1), length(respWindow(1):respWindow(2)));
decper = 3; %the window(s) within which we are going to id IBI deceleration
accper = 4; %the window(s) within which we are going to id IBI acceleration 
sampleraw = 2; %sample at which to save the raw data (for graphing purposes)

allstat = zeros(size(dat,1), 2);
allraw = zeros(size(dat,1), length(respWindow(1):respWindow(2)));
preraw = zeros(size(dat,1), length(respWindow(1):respWindow(2))/sr*sampleraw);

% First, get the stats of interest

for i = 1:size(dat,1)
    datcrit = dat(i,trialstart+respWindow(1):trialstart+respWindow(2)) - mean(dat(i,trialstart+baseWindow(1):trialstart+baseWindow(2)));
    t = min(datcrit(1:decper*sr)); %find IBI deceleration
    tpt = find(datcrit(1:decper*sr) == t);
    if ~isempty(tpt)
        t2 = max(datcrit(tpt:accper*sr)); %find IBI acceleration
    else
        t2 = NaN;
    end
    allstat(i,1) = t;
    allstat(i,2) = t2;
end

%Now Determine the averages for each condition

[NConds,ConVec] = anslabfinddiffel(conds);
ConVec = sort(ConVec);
prestat = zeros(NConds,2);

for j = 1:NConds
    tmp = allstat(conds ==ConVec(j),:);
    prestat(j,:) = nanmean(tmp);
end

% Create a stat matrix with one row - the first NCond columns are IBI dec
% the second NCond columns are IBI acc
stat(1:NConds) = prestat(:,1)';
stat(NConds+1:NConds+NConds) = prestat(:,2)';

% Next, we do something similar for the raw data

for i = 1:size(dat,1)
    allraw(i,:) = dat(i,trialstart+respWindow(1):trialstart+respWindow(2)) - mean(dat(i,trialstart+baseWindow(1):trialstart+baseWindow(2)));
    preraw(i,:) = resample(allraw(i,:), sampleraw, sr);
end

raw = zeros(1,NConds*size(preraw,2));

for j = 1:NConds
    tmp = preraw(conds == ConVec(j),:);
    raw(1,(j-1)*size(preraw,2)+1:j*size(preraw,2)) = nanmean(tmp);
end

end
