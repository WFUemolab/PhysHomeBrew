function [stat, allstat, binDat] = cw_scl_calc(dat, conds, sr, PHYS, sclbase)
% function stat = cw_scl_calc(dat, conds, sr, trig)
%
% This function calculates the post-stimulus change in SCR.
%
% dat is a (m x n) matrix with m trials and n timepts per trial
% conds is a (m x 1) matrix with the condition identifiers for each trial
% sr is the sampling rate (e.g. 20 hz)
% trig is a struct:
%   trig.pre = pre-stimulus baseline time (in seconds)
%   trig.trig = trig time (how long the trigger occupies)
%   trig.post = post-stimulus time
%
% output is a (1 x k) stat matrix with k numbers of conditions 


respWindow = PHYS.SCL.responseWindow * sr; %the window(s) within which we are going to average EMG
baseWindow = PHYS.SCL.baseWindow * sr; %the window within which we are going to get the proximal baseline
trialstart = PHYS.SCL.trialstart * sr;
respDat = zeros(size(dat,1), length(respWindow(1):respWindow(2)));
baseDat = zeros(size(dat,1), length(baseWindow(1):baseWindow(2)));
datcrit = respDat;
if isfield(PHYS.SCL, 'binsize') 
   binsize = PHYS.SCL.binsize;
else
   binsize = 1;
end


allstat = zeros(size(dat,1),1);
sclstd = nanstd(sclbase);
sclmean = nanmean(sclbase);
bin = 1*sr;

binDat = zeros(size(dat,1), size([baseDat respDat],2)/(sr*binsize));

for i = 1:size(dat,1)
     respDat = dat(i,trialstart+respWindow(1)+1:trialstart+respWindow(2));
     baseDat = dat(i,trialstart+baseWindow(1)+1:trialstart+baseWindow(2));
    if strcmp(PHYS.SCL.scrtype, 'proximal')
       datcrit = respDat - mean(baseDat);
        t = max(datcrit); %find SCR
        if t < 0
            t = 0;
        end
        if isfield(PHYS.SCL, 'transform')
            switch PHYS.SCL.transform
                case 'log'
                    t = log(t+1);
                case 'sqrt'
                    t = sqrt(t);
            end
        end
                
    % the above removes proximal baseline - the bottom doesn't
    elseif strcmp(PHYS.SCL.scrtype, 'z')
        datcrit = (respDat  - repmat(sclmean, 1, length(respWindow(1)+1:respWindow(2))))./repmat(sclstd, 1, length(respWindow(1)+1:respWindow(2)));
        t = nanmean(datcrit);
    %the above standardizes from an initial baseline period
    elseif strcmp(PHYS.SCL.scrtype, 'minmax')
        datcrit = dat(i,trialstart+1:trialstart+respWindow(2));
        peak = max(datcrit(respWindow(1)+1:end));
        latency = find(datcrit == peak);
        trough = min(datcrit(1:latency));
        t = max(peak-trough, 0);
        
        if isfield(PHYS.SCL, 'transform')
            switch PHYS.SCL.transform
                case 'log'
                    t = log(t+1);
                case 'sqrt'
                    t = sqrt(t);
            end
        end
    end
        allstat(i) = t;
        binDat(i,:) = binavg([baseDat respDat], binsize*sr);
end

%Now Determine the averages for each condition

[NConds,ConVec] = anslabfinddiffel(conds);
ConVec = sort(ConVec);
prestat = zeros(NConds,1);
alltrials = [];

for j = 1:NConds
    tmp = allstat(conds ==ConVec(j));
    prestat(j) = nanmean(tmp);
    alltrials = [alltrials tmp'];
end

% Create a stat matrix with one row - the first NCond columns are IBI dec
% the second NCond columns are IBI acc
stat(1:NConds) = prestat';

end