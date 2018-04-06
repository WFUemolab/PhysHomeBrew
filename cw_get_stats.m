function cw_get_stats(batch, PHYS, varargin)
% function cw_get_stats(batch, PHYS, 'ecg', 'scl', 'emg', 'startle')
%
% Need to run this function from the top directory (e.g. 'LAB/')
%
% PHYS is a structure of information and parameters created by createPHYS
%
% This function gathers information from the .mat files that ANSLAB creates
% and puts in the various channel folders. This package assumes that you
% have already run the 'CW edited version' of the 'Extended trial extraction' from ANSLAB's menu.
% 
% batch = the name of a batch file with the subjects of interest (raw .acq files)
% OR a string of filenames to search in the current directory (e.g.
% 'LAB*.acq'
%
% IMPORTANT: You can use this function to collect startle stats, however,
% you must remember to set the subject list to the appropriate name (e.g.
% 'LAB02*.acq')
%
% Then need to enter in the various channels of interest. The output text 
% files will be put into the 'stats' folder with the name of each channel
% appended to each subject's file name:
% e.g. 'LAB01000_ecg.txt'
%
% These files can then be uploaded into a stats package

%suffix = '';
suffix = PHYS.suffix;
numchan = length(varargin);
chans = varargin;
base = PHYS.base; %baseline time for EMG (s)

[pathname,name,ext] = fileparts(batch);

if ~isdir(pathname) && strcmp(ext, '.m')
    batch = ['batch/' batch];
end

expdir = pwd;

% Get subject info from the batch file or from the search string
if exist(batch)
    files = textread(batch, '%s', 'delimiter', '\n', 'whitespace','');
else
    filea = ls(['raw/' batch]);
    if ~isempty(filea)
    for i = 1:size(filea,1)
        files{i,1} = [expdir, filesep, 'raw', filesep, filea(i,:)];
    end
    else
    disp('NO FILES LOADED! Try going to the top directory');
    return
    end
end
numsubs = length(files);
for i = 1:numsubs
    [FileDir{i}, SubName{i}] = fileparts(files{i});
end


% For each subject:

for i = 1:numsubs
    
    fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(['\nGetting data from ' SubName{i} '\n']);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    
    % Load subject's Events file that was extracted from their trials
    if ~strcmp(chans{1}, 'startle') && ~strcmp(chans{1}, 'rating')
        ALL = load(['raw/' SubName{i} '_Events.mat']);
    end
    
    T = anslabreadtiming([expdir '/raw/' SubName{i} '.m']);
    if isfield(PHYS, 'baseTrig')
        baseLoc = find(T(:,1) == PHYS.baseTrig);
        baseWindow = [T(baseLoc,2), T(baseLoc,3)];
    end
    % for each channel (e.g. 'ecg') of interest
    
    for j = 1:numchan
        switch chans{j}
            case 'ecg'
                % find where the IBI channel is in the DATA file
                ibiloc = [];
                for k = 1:size(ALL.HDR.ChLabelCell,1)
                    if strcmp(deblank(ALL.HDR.ChLabelCell(k,:)), 'ibi')
                    ibiloc = k;
                    end
                end
                if isempty(ibiloc)
                    disp('NO ECG CHANNEL IN DATA FILE!');
                    return
                end
                
                % get the data from the IBI channel in the DATA file
                ecgdat = zeros(length(ALL.Conditions), size(ALL.DATA{1},2));
                for l = 1:size(ecgdat,1)
                    ecgdat(l,:) = ALL.DATA{l}(ibiloc,:);
                end
                
                % get the stats for this channel
                [ecgstat(i,:), ecgraw(i,:)] = cw_ecg_calc(ecgdat, ALL.Conditions, ALL.HDR.SR{ibiloc}, PHYS);
                
            case 'scl'
                BASE = load(['eda/' PHYS.SCL.channame '/' SubName{i}, '00.mat']);
                
                % find where the SCL channel is in the DATA file
                
                 sclloc = [];
                for k = 1:size(ALL.HDR.ChLabelCell,1)
                    if strcmp(deblank(ALL.HDR.ChLabelCell(k,:)), 'SCL') % should be 'SCL'
                    sclloc = k;
                    end
                end
                if isempty(sclloc)
                    disp('NO SCL CHANNEL IN DATA FILE!');
                    %return
                end
                
                 if strcmp(base, 'Trial')
                        eval(['sclbase = BASE.sc0(baseWindow(1)*BASE.ep:baseWindow(2)*BASE.ep);']);
                        eval(['scltask = resample(BASE.sc0,1,BASE.ep);']);
                    else
                        eval(['sclbase = BASE.sc0(1:base*BASE.ep);']);
                        eval(['scltask = resample(BASE.sc0(base*BASE.ep+1:end),1,BASE.ep);']);
                 end
                                  
                % get the data from the SCL channel in the DATA file
                scldat = zeros(length(ALL.Conditions), size(ALL.DATA{2},2));
               for l = 1:size(scldat,1)
                        if size(scldat(l,:),2) == size(ALL.DATA{l}(sclloc,:),2)
                        scldat(l,:) = ALL.DATA{l}(sclloc,:);
                        else
                            scldat(l,2:end) = ALL.DATA{l}(sclloc,:);
                        end
               end
                  if isfield(PHYS.SCL, 'binsize') 
                        binsize = PHYS.SCL.binsize;
                    else
                        binsize = 1;
                    end
               if i == 1
                        %for getting extended trial by trial data (in bins)
                        windowSize = ((PHYS.SCL.responseWindow(2) - PHYS.SCL.responseWindow(1)) + (PHYS.SCL.baseWindow(2) - PHYS.SCL.baseWindow(1)))/binsize;
                        SCLBinTrials = zeros(length(ALL.Conditions), windowSize + 3);
                        SCLBinTrials(:,1) = repmat(str2num(SubName{i}(5:8)), length(ALL.Conditions),1);
                        SCLBinTrials(:,2) = (1:length(ALL.Conditions))';
                        SCLBinTrials(:,3) = ALL.Conditions;
                        
                        SCLAllTrials = zeros(length(ALL.Conditions), 4);
                        SCLAllTrials(:,1) = repmat(str2num(SubName{i}(5:8)), length(ALL.Conditions),1);
                        SCLAllTrials(:,2) = (1:length(ALL.Conditions))';
                        SCLAllTrials(:,3) = ALL.Conditions;
               else
                        SCLBinTrials(end+1:end+length(ALL.Conditions),1:3) = [repmat(str2num(SubName{i}(5:8)), length(ALL.Conditions),1), (1:length(ALL.Conditions))', ALL.Conditions'];
                        SCLAllTrials(end+1:end+length(ALL.Conditions),1:3) = [repmat(str2num(SubName{i}(5:8)), length(ALL.Conditions),1), (1:length(ALL.Conditions))', ALL.Conditions'];
               end
                % get the stats for this channel
                [sclstatTMP, SCLAllTrials(end-length(ALL.Conditions)+1:end, 4), SCLBinTrials(end-length(ALL.Conditions)+1:end, 4:windowSize+3)] = cw_scl_calc(scldat, ALL.Conditions, ALL.HDR.SR{sclloc}, PHYS, sclbase);
                %[sclstat(i,:), sclalltrials(i,:)] = cw_scl_calc(scldat, ALL.Conditions, ALL.HDR.SR{sclloc}, PHYS, sclbase);
                SCLAllTrials(end-length(ALL.Conditions)+1:end, 5) = T(:,4);
                SCLBinTrials(end-length(ALL.Conditions)+1:end, windowSize+4) = T(:,4);
                
                if i == 1
                    sclstat(i,:) = sclstatTMP;
                    else
                        sclstat(i,:) = padarray(sclstatTMP, [0 size(sclstat(1,:),2) - size(sclstatTMP,2)], NaN, 'post');
                    end
                sclbrms(i,1) = rmssd(resample(sclbase,1,BASE.ep));
                sclrms(i,1) = rmssd(scltask);
            case 'emg'
                % Get Raw data to calculate baseline
                BASE = load(['emg/' PHYS.EMG.channame '\' SubName{i}, '00.mat']);
                numemg = length(PHYS.EMG.names); % number of emg channels of interest (e.g. zyg, orb, cor)
                emgnames = PHYS.EMG.names;
                % Need to get all EMG channels for the trials
                for m = 1:numemg
                    emgloc = [];
                    for k = 1:size(ALL.HDR.ChLabelCell,1)
                        if strcmp(deblank(ALL.HDR.ChLabelCell(k,:)), ['EMG ' int2str(m)]) % should be ['EMG ' int2str(m)]
                            emgloc = k;
                        end
                    end
                    if isempty(emgloc)
                        disp('NO EMG CHANNEL IN DATA FILE!');
                        return
                    end
                    
                    % get the data from the EMG channel in the DATA file
                    emgdat = zeros(length(ALL.Conditions), size(ALL.DATA{2},2));
                    if isfield(PHYS.EMG, 'binsize') 
                        binsize = PHYS.EMG.binsize;
                    else
                        binsize = 1;
                    end
                    
                    if (i == 1 && m == 1)
                        %for getting average trial by trial data
                        emgAvgTrials = zeros(length(ALL.Conditions), numemg+3);
                        emgAvgTrials(:,1) = repmat(str2num(SubName{i}(5:8)), length(ALL.Conditions),1);
                        emgAvgTrials(:,2) = (1:length(ALL.Conditions))';
                        emgAvgTrials(:,3) = ALL.Conditions;
                        
                      
                        
                    elseif m ==1
                        emgAvgTrials(end+1:end+length(ALL.Conditions),1:3) = [repmat(str2num(SubName{i}(5:8)), length(ALL.Conditions),1), (1:length(ALL.Conditions))', ALL.Conditions'];
                        
                    end
                    
                    if i == 1
                        %for getting extended trial by trial data (in bins)
                        windowSize = ((PHYS.EMG.responseWindow(2) - PHYS.EMG.responseWindow(1)) + (PHYS.EMG.baseWindow(2) - PHYS.EMG.baseWindow(1)))/binsize;
                        emgBinTrials{m} = zeros(length(ALL.Conditions), windowSize + 3);
                        emgBinTrials{m}(:,1) = repmat(str2num(SubName{i}(5:8)), length(ALL.Conditions),1);
                        emgBinTrials{m}(:,2) = (1:length(ALL.Conditions))';
                        emgBinTrials{m}(:,3) = ALL.Conditions;
                    else
                        emgBinTrials{m}(end+1:end+length(ALL.Conditions),1:3) = [repmat(str2num(SubName{i}(5:8)), length(ALL.Conditions),1), (1:length(ALL.Conditions))', ALL.Conditions'];
                    end
                    
                    for l = 1:size(emgdat,1)
                        if size(emgdat(l,:),2) == size(ALL.DATA{l}(emgloc,:),2)
                        emgdat(l,:) = ALL.DATA{l}(emgloc,:);
                        else
                            emgdat(l,2:end) = ALL.DATA{l}(emgloc,:);
                        end
                    end
                    % get base data for this emg channel
                    if strcmp(base, 'trial')
                        eval(['emgbase = BASE.N' int2str(m) '_0;']);
                        eval(['emgtask = resample(BASE.N' int2str(m) '_0,1,BASE.ep);']);
                    else
                        eval(['emgbase = BASE.N' int2str(m) '_0(1:base*BASE.ep);']);
                        eval(['emgtask = resample(BASE.N' int2str(m) '_0(base*BASE.ep+1:end),1,BASE.ep);']);
                    end
                    % get the stats for this channel
                    [emgstatTMP, emgAvgTrials(end-length(ALL.Conditions)+1:end, m+3), emgBinTrials{m}(end-length(ALL.Conditions)+1:end, 4:windowSize+3), emgbaseStats(i,m*2-1), emgbaseStats(i,m*2)] = cw_emg_calc(emgdat, ALL.Conditions, ALL.HDR.SR{emgloc}, PHYS, emgbase, binsize);
                    if i == 1
                    emgstat{m}(i,:) = emgstatTMP;
                    else
                        emgstat{m}(i,:) = padarray(emgstatTMP, [0 size(emgstat{m}(1,:),2) - size(emgstatTMP,2)], NaN, 'post');
                    end
                    emgbrms{m}(i,:) = rmssd(resample(emgbase, 1, BASE.ep));
                    emgrms{m}(i,:) = rmssd(emgtask);
                end
            case 'startle'
                
                %first, need to get the data and the condition numbers
                %[trialnum,time,amp,latency,onslate,bamp,condx] = textread(['startle/' PHYS.startle.channame '/' SubName{i} '00.txt'], '%n%n%n%n%n%n%n', 'headerlines', 1);
                amp2 = load(['startle/' PHYS.startle.channame '/' SubName{i} '00.mat']);
                amp = amp2.sa;
                T = anslabreadtiming(['raw/' SubName{i} '.m']);
                conds = T(:,1);
                
                % Need to get two variables, one for amplitude (no 0's) and
                % one with magnitude (with 0')
                mag = amp;
                amp(amp == -9999 | amp == 0) = NaN;
                mag(mag == -9999) = NaN;
                
                
                % This will standardize the startle data for each subject
                % and set the mean = 50 and STd = 10
                amp = tscale(amp);
                mag = tscale(mag);
                
                % This will winsorize the data (set data > 3sds = 3sds)
                amp = winsor(amp,3);
                mag = winsor(mag,3);
                
                %try
                 %   startletot(i,:) = amp';
                %catch
                 %   startletot(i,:) = [amp', repmat(-9999, 1, 96-size(amp,1))];
                %end
                
                %Set up condition information
                [NConds,ConVec] = anslabfinddiffel(conds);
                ConVec = sort(ConVec);
                prestat = zeros(NConds,2);
                
                %Mean across trials within each condition
                for k = 1:NConds
                    tmpamp = amp(conds ==ConVec(k));
                    tmpmag = mag(conds ==ConVec(k));
                    prestat(k,1) = nanmean(tmpamp);
                    prestat(k,2) = nanmean(tmpmag);
                    prestat(k,3) = nansum(tmpamp) / nanmean(tmpamp);
                end
                
                % Create a stat matrix with two sets of cols - the first
                % NCond columns are Amp the second NCond columns are Mag
                startstat(i,1:NConds) = prestat(:,1)';
                startstat(i,NConds+1:NConds+NConds) = prestat(:,2)';
                startnums(i,1:NConds) = prestat(:,3)';
            case 'rating'
                
                    [tmp(:,1), tmp(:,2)] = textread(['rating/' SubName{i} '.txt'], '%d%f', 'headerlines', 1);
                    rate(1,i*2-1:i*2) = repmat(str2num(SubName{i}(4:8)),1,2);
                    rate(2:size(tmp,1)+1, i*2-1:i*2) = tmp;
                    clear tmp;
               
        end
    end
    
end

cd([expdir '\stats']);

% Now write out the data to .mat and .txt formats
for j = 1:numchan
    switch chans{j}
        case 'ecg'
            save('ecgstats', 'ecgstat', 'ecgraw');
            csvwrite(['ecgstats' suffix '.txt'], ecgstat);
            csvwrite(['ecgraw' suffix '.txt'], ecgraw);
        case 'scl'
            save('sclstats', 'sclstat');
           % save('sclalltrials', 'sclalltrials');
            csvwrite(['sclstats' suffix '.txt'], sclstat);
            %csvwrite(['sclalltrials' suffix '.txt'], sclalltrials);
            save('sclrms', 'sclrms');
            csvwrite(['sclrms' suffix '.txt'], sclrms);
            csvwrite(['sclbrms' suffix '.txt'], sclbrms);
            save('sclraw', 'scldat');
            csvwrite(['sclraw' suffix '.txt'], scldat);
            save('SCLBinTrials', 'SCLBinTrials');
            csvwrite(['SCLBinTrials' suffix '.txt'], SCLBinTrials);
            save('SCLAllTrials', 'SCLAllTrials');
            csvwrite(['SCLAllTrials' suffix '.txt'], SCLAllTrials);
        case 'emg'
            for m = 1:numemg
                emg = emgstat{m};
                emgr = emgrms{m};
                emgb = emgbrms{m};
                emgBin = emgBinTrials{m};
                save(['emgstat_' emgnames{m}], 'emg');
                save(['emgrms_' emgnames{m}], 'emgr');
                save(['emgBinTrials_' emgnames{m}], 'emgBin');
                csvwrite(['emgstat_' emgnames{m} suffix '.txt'], emg);
                csvwrite(['emgrms_' emgnames{m} suffix '.txt'], emgr);
                csvwrite(['emgbrms_' emgnames{m} suffix '.txt'], emgb);
                csvwrite(['emgBinTrials_' emgnames{m} suffix '.txt'], emgBin);
            end
            save('emgAvgTrials', 'emgAvgTrials');
            save('emgbaseStats', 'emgbaseStats');
            csvwrite(['emgAvgTrials' suffix '.txt'], emgAvgTrials);
            csvwrite(['emgbaseStats' suffix '.txt'], emgbaseStats);
        case 'startle'
            save('startle_stats', 'startstat');
            save('startle_nums', 'startnums');
            %save('startle_rawamps', 'startletot');
            csvwrite(['startle_stats' suffix '.txt'], startstat);
            csvwrite(['startle_nums' suffix '.txt'], startnums);
        case 'rating'
            save('ratings', 'rate');
            csvwrite(['rateDataAll' suffix '.txt'], rate);
    end
end
   