function [stat, allstat, binDat, emgmean, emgstd] = cw_emg_calc(dat, conds, sr, PHYS, emgbase)
% function stat = cw_emg_calc(dat, conds, sr, trig)
%
% This function calculates the post-stimulus change in EMG activity.
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


respWindow = PHYS.EMG.responseWindow * sr; %the window(s) within which we are going to average EMG
baseWindow = PHYS.EMG.baseWindow * sr; %the window within which we are going to get the proximal baseline
trialstart = PHYS.EMG.trialstart * sr;
binsize = PHYS.EMG.binsize;

respDat = zeros(size(dat,1), length(respWindow(1):respWindow(2)));
baseDat = zeros(size(dat,1), length(baseWindow(1):baseWindow(2)));
datcrit = respDat;
allstat = zeros(size(dat,1),1);
bin = 1*sr; % bin for moving average analyses; default is 1s
mavg = zeros(length(dat) - bin,1); %initialize moving average variable

% get baseline STD info for standardizing data
emgstd = nanstd(emgbase);
emgmean = nanmean(emgbase);

%get extended trial data (binned)
binDat = zeros(size(dat,1), size([baseDat respDat],2)/(sr*binsize));

for i = 1:size(dat,1)
     respDat = dat(i,trialstart+respWindow(1)+1:trialstart+respWindow(2));
     baseDat = dat(i,trialstart+baseWindow(1)+1:trialstart+baseWindow(2));
     if strcmp(PHYS.EMG.basetype, 'proximal')
        datcrit = respDat - mean(baseDat);
        if strcmp(PHYS.EMG.meanmax, 'max')
            t = max(datcrit); %find 
        else 
            t = nanmean(datcrit);
        end
        %if t < 0
         %   t = 0;
        %end
    % the above removes proximal baseline - the bottom doesn't
    elseif strcmp(PHYS.EMG.basetype, 'z')
        datcrit = (respDat  - repmat(emgmean, 1, size(respDat,2))./repmat(emgstd, 1, size(respDat,2)));
         if strcmp(PHYS.EMG.meanmax, 'max')
            t = max(datcrit); %find 
        else 
            t = nanmean(datcrit);
        end
    %the above standardizes from an initial baseline period
     elseif strcmp(PHYS.EMG.basetype, 'nobase')
         datcrit = dat(i,trialstart+respWindow(1)+1:trialstart+respWindow(2));
          if strcmp(PHYS.EMG.meanmax, 'max')
            t = max(datcrit); %find 
        else 
            t = nanmean(datcrit);
        end
    end
    allstat(i) = t;
    binDat(i,:) = binavg([baseDat respDat], binsize*sr);
end

%Now Determine the averages for each condition

%[NConds,ConVec] = anslabfinddiffel(conds);
ConVec = unique(conds); %workaround when anslab is unavailable
NConds = size(ConVec,2);

ConVec = sort(ConVec);
prestat = zeros(NConds,1);

for j = 1:NConds
    tmp = allstat(conds ==ConVec(j));
    prestat(j) = nanmean(tmp);
end

% Create a stat matrix with one row - the first NCond columns are IBI dec
% the second NCond columns are IBI acc
stat(1:NConds) = prestat';

end