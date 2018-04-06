function cw_getbasestats(batch, PHYS, varargin)
% function getstats(subindex, 'ecg', 'scl') is used to get stats from mat files that are
% sprinkled throughout the directories. Used mainly for getting baseline
% stats. But can also get ecg recovery 'ecgrec'. If so, need to specify in
% the script the baseline length and critical window for calculating
% recovery.

% Make sure to call in the top directory i.e. 'LAB'
% 
% batch should either be a batch file name or a string of filenames to 
% search for (e.g. 'CS201*.acq') 
%
% channel corresponds to the channel of interest, e.g.:
% 'ecg' - get mean HR for an interval
% 'ecghrv' - get mean HR and spectral info for an interval
% 'scl' - get mean SCL, initial SCL, and slope over an interval
% 'emg' - get mean EMG and SD of EMG over an interval
% 'bprec' - get BP recovery info
% 'bpicg' - get bp info with icg - like TPR and MAP


numchan = length(varargin);
expdir = pwd;
base = PHYS.base; %length of baseline in seconds


% baserec = ; %length of baseline for calculating ecg recovery
% recwin = 60; %window of ecg data in which to calculate recovery
% %numsegs = 4; %number of segments throughout your session
% newbins = 60; % in seconds, bins to be resampled to for easy analysis

chans = varargin;

expdir = pwd;

% Get subject info from the batch file or from the search string
if exist(['batch/' batch])
    files = textread(['batch/' batch], '%s', 'delimiter', '\n', 'whitespace','');
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
numsubs = size(files,1);
for i = 1:numsubs
    [FileDir{i}, SubName{i}] = fileparts(files{i});
end

for i = 1:numsubs
    
    fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(['\nGetting data from ' SubName{i} '\n']);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    
     %T = anslabreadtiming([expdir '/raw/' SubName{i} '.m']);
     T = anslabreadtiming([expdir '/raw/' SubName{i} '.spectral.m']);
     numsegs = size(T,1); 
     
    for j = 1:numchan
        switch chans{j}
            case 'ecg'
                load(['ecg/' PHYS.ECG.channame '/' SubName{i}, '00.mat']);
                binsize = PHYS.ECG.binsize;
                stats.ecg(i,1) = nanmean(ibi0(1:base*4)); %mean IBI for baseline at 4hz
                stats.ecg(i,2) = 60000/stats.ecg(i,1); %mean HR for baseline
                Tsamp = round(T * ep);
               Tnew = Tsamp / (binsize*ep);
                for k = 1:size(T,1)
                     if k == 1
                        range = [1,ceil(Tnew(1,4))];
                     else
                    range(k,1) = range(k-1,2)+1;
                    range(k,2) = range(k,1) + Tnew(k,4) -1;
                     end
                     try
                     %stats.sclall(i,range(k,1):range(k,2)) = resample(sc0(Tsamp(k,2)+1:Tsamp(k,3)), 1, newbins*ep,0);
                      stats.ecgall(i,range(k,1):range(k,2)) = binavg(ibi0(Tsamp(k,2)+1:Tsamp(k,3)), binsize*ep);
                     catch
                      if k == size(T,1)
                          try
                          %stats.sclall(i,Tnew(k,2)+1:Tnew(k,3)) = resample(sc0(Tsamp(k,2)+1:end), 1, newbins*ep,0);
                          stats.ecgall(i,Tnew(k,2)+1:Tnew(k,3)) = binavg(ibi0(Tsamp(k,2)+1:end), binsize*ep);
                          catch
                              disp(['Cannot get data from ' SubName{i}]);
                          end
                      end
                  end
                end
                stats.ecgall(i,size(T,1)+1) = str2num(SubName{i}(5:8));
            case 'ecghrv'
                load(['spectral/ibi0/' SubName{i}, '00.mat']);
                numsegs = size(meanHR,2);
%                 if size(meanIBI,2) > numsegs
%                     meanIBI = meanIBI(size(meanIBI,2)-numsegs+1:size(meanIBI,2));
%                     meanHR = meanHR(size(meanHR,2)-numsegs+1:size(meanHR,2));
%                     VLFBpower = VLFBpower(size(LFBpower,2)-numsegs+1:size(VLFBpower,2));
%                     LFBpower = LFBpower(size(LFBpower,2)-numsegs+1:size(LFBpower,2));
%                     HFBpower = HFBpower(size(HFBpower,2)-numsegs+1:size(HFBpower,2));
%                     RMSsd = RMSsd(size(RMSsd,2)-numsegs+1:size(RMSsd,2));
%                 end
                stats.ecghrv(i,(1*numsegs)-(numsegs-1):1*numsegs) = meanIBI;
                stats.ecghrv(i,(2*numsegs)-(numsegs-1):2*numsegs) = meanHR;
                stats.ecghrv(i,(3*numsegs)-(numsegs-1):3*numsegs) = VLFBpower;
                stats.ecghrv(i,(4*numsegs)-(numsegs-1):4*numsegs) = LFBpower;
                stats.ecghrv(i,(5*numsegs)-(numsegs-1):5*numsegs) = HFBpower;
                stats.ecghrv(i,(6*numsegs)-(numsegs-1):6*numsegs) = log(HFBpower);
                stats.ecghrv(i,(7*numsegs)-(numsegs-1):7*numsegs) = RMSsd;
                stats.ecghrv(i,(8*numsegs)-(numsegs-1)) = str2num(SubName{i}(5:8));
                stats.labels = {'MeanIBI' 'MeanHR' 'LF power' 'VLF power' 'HF power' 'RSA' 'RMSSD' 'SubID'};
                
                
            case 'bpicg'
                havebp = 1;
                haveicg = 1;
                try
                    load(['bp/' PHYS.BP.channame '/' SubName{i}, '00.mat']);
                catch
                    havebp = 0;
                end
                try
                    load(['icg/' PHYS.ICG.channame '/' SubName{i}, '00.mat']);
                catch
                    haveicg = 0;
                end
%                 if size(fiDIA,2) > numsegs
%                     fiDIA = fiDIA(size(fiDIA,2)-numsegs+1:size(fiDIA,2));
%                     fiSYS = fiSYS(size(fiSYS,2)-numsegs+1:size(fiSYS,2));
%                     end
%                 stats.bp(i,(1*numsegs)-(numsegs-1):1*numsegs) = fiDIA;
%                 stats.bp(i,(2*numsegs)-(numsegs-1):2*numsegs) = fiSYS;
%                 stats.bp(i,(3*numsegs)-(numsegs-1)) = str2num(SubName{i}(5:8));
                   if havebp == 1
                    Tnew = round(T * ep);
                    timediff = diff(time);
                    badtimes = find(timediff > 200);
                    badspikes = sort([find((fiSYS < (nanmean(fiSYS)-3*nanstd(fiSYS)))) find(fiSYS > (nanmean(fiSYS)+3*nanstd(fiSYS)))]);
                    newSYS = fiSYS;
                    newDIA = fiDIA;
                    newSYS(badtimes) = NaN; newSYS(badspikes) = NaN;
                    newDIA(badtimes) = NaN; newDIA(badspikes) = NaN;
                    newSYS0 = epoch(time, newSYS, 25);
                    newDIA0 = epoch(time, newDIA, 25);
                    MAP_0 = newDIA0 + 1/3*(newSYS0 - newDIA0);
                    SYS = zeros(1,numsegs);
                    DIA = zeros(1,numsegs);
                    MAP = zeros(1,numsegs);
                    BAD = zeros(1,numsegs);
                                        
                    for k = 1:size(Tnew,1)
                        if k == size(Tnew,1)
                            SYS(k) = nanmean(newSYS0(Tnew(k,2)+1:end));
                            DIA(k) = nanmean(newDIA0(Tnew(k,2)+1:end));
                            MAP(k) = nanmean(MAP_0(Tnew(k,2)+1:end));
                        else
                            SYS(k) = nanmean(newSYS0(Tnew(k,2)+1:Tnew(k,3)));
                            DIA(k) = nanmean(newDIA0(Tnew(k,2)+1:Tnew(k,3)));
                            MAP(k) = nanmean(MAP_0(Tnew(k,2)+1:Tnew(k,3)));
                        end
                    end
                    
                    stats.bp(i,(1*numsegs)-(numsegs-1):1*numsegs) = SYS; 
                    stats.bp(i,(2*numsegs)-(numsegs-1):2*numsegs) = DIA;
                    stats.bp(i,(3*numsegs)-(numsegs-1):3*numsegs) = MAP;
                   end
                   
                   if haveicg == 1
                        if size(PEP_en,2) > numsegs
                    PEP_en = PEP_en(size(PEP_en,2)-numsegs+1:size(PEP_en,2));
                    STV_en = STV_en(size(STV_en,2)-numsegs+1:size(STV_en,2));
                    CO_en = CO_en(size(CO_en,2)-numsegs+1:size(CO_en,2));
                    HI_en = HI_en(size(HI_en,2)-numsegs+1:size(HI_en,2));
                    LVET_en = LVET_en(size(LVET_en,2)-numsegs+1:size(LVET_en,2));
                        end
                    stats.icg(i,(1*numsegs)-(numsegs-1):1*numsegs) = PEP_en;
                    stats.icg(i,(2*numsegs)-(numsegs-1):2*numsegs) = STV_en;
                    stats.icg(i,(3*numsegs)-(numsegs-1):3*numsegs) = CO_en;
                    stats.icg(i,(4*numsegs)-(numsegs-1):4*numsegs) = HI_en;
                    stats.icg(i,(5*numsegs)-(numsegs-1):5*numsegs) = LVET_en;
                    stats.icg(i,(6*numsegs)-(numsegs-1)) = str2num(SubName{i}(5:8));
                   end
                   
                  if havebp == 1 && haveicg == 1
                       
                    CO_0 = epoch(rtbbb, CO_bbb, 250);
                    TPR_0 = MAP_0(1:min(size(MAP_0,1), size(CO_0,1))) ./ CO_0(1:min(size(MAP_0,1), size(CO_0,1))) *80;
                    
                    for k = 1:size(Tnew,1)
                        if k == size(Tnew,1)
                            
                            TPR(k) = nanmean(TPR_0(Tnew(k,2)+1:end));
                            CO(k) = nanmean(CO_0(Tnew(k,2)+1:end));
                            BAD(k) = numel(find(isnan(TPR_0(Tnew(k,2)+1:end)))) ./ size(TPR_0(Tnew(k,2)+1:end),1);
                        else
                            
                            TPR(k) = nanmean(TPR_0(Tnew(k,2)+1:Tnew(k,3)));
                            CO(k) = nanmean(CO_0(Tnew(k,2)+1:Tnew(k,3)));
                            BAD(k) = numel(find(isnan(TPR_0(Tnew(k,2)+1:Tnew(k,3))))) ./ size(TPR_0(Tnew(k,2)+1:Tnew(k,3)),1);
                        end
                    end
                                     
                  
                  stats.bp(i,(4*numsegs)-(numsegs-1):4*numsegs) = TPR;
                  stats.bp(i,(5*numsegs)-(numsegs-1):5*numsegs) = CO;
                  stats.bp(i,(6*numsegs)-(numsegs-1):6*numsegs) = BAD;
                   end
                   
                  stats.bp(i,(7*numsegs)-(numsegs-1)) = str2num(SubName{i}(5:8));


           case 'icg'
                T = anslabreadtiming([expdir '/raw/' SubName{i} '.icg.m']);
                numsegs = size(T,1); 
                eval(['cd ' expdir '\icg\dzdt']);
                eval(['load ' SubName{i} '00.mat']);
                if size(PEP_en,2) > numsegs
                    PEP_en = PEP_en(size(PEP_en,2)-numsegs+1:size(PEP_en,2));
                    STV_en = STV_en(size(STV_en,2)-numsegs+1:size(STV_en,2));
                    CO_en = CO_en(size(CO_en,2)-numsegs+1:size(CO_en,2));
                    HI_en = HI_en(size(HI_en,2)-numsegs+1:size(HI_en,2));
                    LVET_en = LVET_en(size(LVET_en,2)-numsegs+1:size(LVET_en,2));
                end
                stats.icg(i,(1*numsegs)-(numsegs-1):1*numsegs) = PEP_en;
                stats.icg(i,(2*numsegs)-(numsegs-1):2*numsegs) = STV_en;
                stats.icg(i,(3*numsegs)-(numsegs-1):3*numsegs) = CO_en;
                stats.icg(i,(4*numsegs)-(numsegs-1):4*numsegs) = HI_en;
                stats.icg(i,(5*numsegs)-(numsegs-1):5*numsegs) = LVET_en;
                stats.icg(i,(6*numsegs)-(numsegs-1)) = str2num(SubName{i}(5:8));
            case 'ecgrec'
                ecg = load(['ecg/' PHYS.ECG.channame '/' SubName{i}, '00.mat']);
                T = anslabreadtiming([expdir '/raw/' SubName{i} '.spectral.m']);
                Tsamp = round(T * ecg.ep);
                ecgbase.ibi0 = ecg.ibi0(Tsamp(1,2)+1:Tsamp(1,3));
                ecgrec.ibi0 = ecg.ibi0(Tsamp(4,2)+1:Tsamp(4,3));
                stats.ecgrec(i,1) = nanmean(ecgbase.ibi0);
                stats.ecgrec(i,2) = nanmean(ecgrec.ibi0(1:5*ecg.ep));
                stats.ecgrec(i,3) = nanmean(ecgrec.ibi0(60*ecg.ep:65*ecg.ep));
                stats.ecgrec(i,4) = nanmean(ecgrec.ibi0(120*ecg.ep:125*ecg.ep));
                [stats.ecgrec(i,5), stats.ecgrec(i,6)] = reccalc(ecgbase.ibi0, ecgrec.ibi0, ecg.ep, recwin, 1);
                if length(ecgrec.ibi0) < recwin*ecg.ep
                    stats.ecgrecraw(i,1:length(resample(ecgrec.ibi0, 1, 4))) = resample(ecgrec.ibi0,1,4)';
                else
                    stats.ecgrecraw(i,:) = resample(ecgrec.ibi0(1:recwin*ecg.ep), 1, 4)';
                               
                end
               % catch
                   % disp(['No file found for ' SubName{i}])
                %end
                
            case 'scl'
                load(['eda/' PHYS.SCL.channame '/' SubName{i}, '00.mat']);
                binsize = PHYS.SCL.binsize;
%                 sclbase = sc0(1:base*25);
%                 if size(sclbase,1) < size(sclbase,2)
%                     sclbase = sclbase';
%                 end
%                 Xlin(:,1) = ones(base*25,1);
%                 Xlin(:,2) = (1:base*25)';
%                 Blin = regress(sclbase, Xlin);
%                 stats.scl(i,1) = Blin(1); % intercept
%                 stats.scl(i,2) = Blin(2)*25*60; %linear slope (transformed to us per minute)
%                 Xlog(:,1) = ones(base*25,1);
%                 Xlog(:,2) = log(1:base*25)';
%                 Blog = regress(sclbase, Xlog);
%                 stats.scl(i,3) = Blog(1); % intercept
%                 stats.scl(i,4) = Blog(2)*25*60; % log slope ;
%                 stats.scl(i,5) = nanmean(sclbase);
%                 stats.scl(i,6) = nanstd(sclbase);
                Tsamp = round(T * ep);
                Tnew = Tsamp / (binsize*ep);
                for k = 1:size(T,1)
                     if k == 1
                        range = [1,ceil(Tnew(1,4))];
                     else
                    range(k,1) = range(k-1,2)+1;
                    range(k,2) = range(k,1) + Tnew(k,4) -1;
                     end
                     try
                     %stats.sclall(i,range(k,1):range(k,2)) = resample(sc0(Tsamp(k,2)+1:Tsamp(k,3)), 1, newbins*ep,0);
                      stats.sclall(i,range(k,1):range(k,2)) = binavg(sc0(Tsamp(k,2)+1:Tsamp(k,3)), binsize*ep);                      
                     catch
                      if k == size(T,1)
                          try
                          %stats.sclall(i,Tnew(k,2)+1:Tnew(k,3)) = resample(sc0(Tsamp(k,2)+1:end), 1, newbins*ep,0);
                          stats.sclall(i,Tnew(k,2)+1:Tnew(k,3)) = binavg(sc0(Tsamp(k,2)+1:end), binsize*ep);
                          catch
                              disp(['Cannot get data from ' SubName{i}]);
                          end
                      end
                     end                            
                end
                
            case 'emg'
                numemg = PHYS.EMG.num; %number of emg channels
                load(['emg/' PHYS.EMG.channame '/' SubName{i}, '00.mat']);
                binsize = PHYS.EMG.binsize;
                for m = 1:numemg
                    %eval(['dat = detrend(N' int2str(m) '_0);']);
                    eval(['dat = N' int2str(m) '_0;']);
                    stats.emgbase(i,m) = mean(dat(1:ep*base));
                    stats.emgbase(i,m+numemg) = std(dat(1:ep*base));
                     Tsamp = round(T * ep);
                    %Tnew = Tsamp / (binsize*ep);
                    numseg = size(T,1);
                    
                    if strcmp(PHYS.EMG.analysis, 'frequency')
                    thr = stats.emgbase(i,m) + 1*stats.emgbase(i,m+numemg);
                        for k = 1:numseg
                        tempdat = binavg(dat(Tsamp(k,2)+1:Tsamp(k,3)), binsize*ep);
                        stats.emg(i,k+(m-1)*numsegs) = numel(find(tempdat > thr));
                       
                        end
                    else
                        for k = 1:numseg
                        tempdat = nanmean(dat(Tsamp(k,2)+1:Tsamp(k,3)));
                        stats.emg(i,k+(m-1)*numsegs) = tempdat;
                        end
                     stats.subs(i,:) = SubName(:,i);
                    end
                end           
            case 'bprec'
                eval(['cd ' expdir '\bp\BP'' - NIBP100D - Non-invasive Blood Pres''']);
                try
                bpbase = load([SubName{i} '02.mat']);
                bprec = load([SubName{i} '04.mat']);
                stats.bprec(i,1) = nanmean(bpbase.fiSYS0);
                stats.bprec(i,2) = nanmean(bprec.fiSYS0(1:5*bprec.ep));
                stats.bprec(i,3) = nanmean(bprec.fiSYS0(60*bprec.ep:65*bprec.ep));
                stats.bprec(i,4) = nanmean(bprec.fiSYS0(120*bprec.ep:125*bprec.ep));
                [stats.bprec(i,5), stats.bprec(i,6)] = reccalc(bpbase.fiSYS0, bprec.fiSYS0, bprec.ep, recwin);
                if length(bprec.fiSYS0) < recwin*bprec.ep
                    stats.bprecraw(i,1:length(resample(bprec.fiSYS0, 1, 5*bprec.ep,0))) = resample(bprec.fiSYS0,1,5*bprec.ep,0)';
                else
                    stats.bprecraw(i,:) = resample(bprec.fiSYS0(1:recwin*bprec.ep), 1, 5*bprec.ep,0)';
                end
                catch
                    disp(['No file found for ' SubName{i}])
                end
                
            case 'resp'
                load(['resp/' PHYS.RESP.channame '/' SubName{i}, '00.mat']);
                binsize = PHYS.RESP.binsize;
%                 sclbase = sc0(1:base*25);
%                 if size(sclbase,1) < size(sclbase,2)
%                     sclbase = sclbase';
%                 end
%                 Xlin(:,1) = ones(base*25,1);
%                 Xlin(:,2) = (1:base*25)';
%                 Blin = regress(sclbase, Xlin);
%                 stats.scl(i,1) = Blin(1); % intercept
%                 stats.scl(i,2) = Blin(2)*25*60; %linear slope (transformed to us per minute)
%                 Xlog(:,1) = ones(base*25,1);
%                 Xlog(:,2) = log(1:base*25)';
%                 Blog = regress(sclbase, Xlog);
%                 stats.scl(i,3) = Blog(1); % intercept
%                 stats.scl(i,4) = Blog(2)*25*60; % log slope ;
%                 stats.scl(i,5) = nanmean(sclbase);
%                 stats.scl(i,6) = nanstd(sclbase);
                Tsamp = round(T * ep);
                Tnew = Tsamp / (binsize*ep);
                for k = 1:size(T,1)
                     if k == 1
                        range = [1,ceil(Tnew(1,4))];
                     else
                    range(k,1) = range(k-1,2)+1;
                    range(k,2) = range(k,1) + Tnew(k,4) -1;
                     end
                     try
                     %stats.sclall(i,range(k,1):range(k,2)) = resample(sc0(Tsamp(k,2)+1:Tsamp(k,3)), 1, newbins*ep,0);
                      stats.respall(i,range(k,1):range(k,2)) = binavg(RR0(Tsamp(k,2)+1:Tsamp(k,3)), binsize*ep);                      
                     catch
                      if k == size(T,1)
                          try
                          %stats.sclall(i,Tnew(k,2)+1:Tnew(k,3)) = resample(sc0(Tsamp(k,2)+1:end), 1, newbins*ep,0);
                          stats.respall(i,Tnew(k,2)+1:Tnew(k,3)) = binavg(RR0(Tsamp(k,2)+1:end), binsize*ep);
                          catch
                              disp(['Cannot get data from ' SubName{i}]);
                          end
                      end
                     end                            
                end
            case 'bp2'
                eval(['cd ' expdir '\bp\BP'' - NIBP100D - Non-invasive Blood Pres'''])
                eval(['load ' SubName{i} '00.mat']);
                T = anslabreadtiming([expdir '/raw/' SubName{i} '.spectral.m']);
                stats.bp(i,1) = nanmean(fiSYS0(1:base*4)); %mean IBI for baseline at 4hz
                stats.bp(i,2) = 60000/stats.bp(i,1); %mean HR for baseline
                Tsamp = round(T * ep);
                Tnew = Tsamp / (newbins*ep);
                for k = 1:size(T,1)
                     if k == 1
                        range = [1,Tnew(1,3)];
                     else
                    range(k,1) = range(k-1,2)+1;
                    range(k,2) = range(k,1) + Tnew(k,4) -1;
                     end
                     try
                     stats.bpall(i,range(k,1):range(k,2)) = resample(fiSYS0(Tsamp(k,2)+1:Tsamp(k,3)), 1, 20,0);
                   catch
                      if k == size(T,1)
                          try
                          stats.bpall(i,Tnew(k,2)+1:Tnew(k,3)) = resample(fiSYS0(Tsamp(k,2)+1:end), 1, 20,0);
                          catch
                              disp(['Cannot get data from ' SubName{i}]);
                          end
                      end
                    end
                end
        end
    end
end

eval(['cd ' expdir '\stats']);
save('basestats', 'stats');

for k = 1:numchan
    switch chans{k}
        case 'ecg'
            csvwrite(['ecgstats_base' PHYS.suffix '.txt'], stats.ecg);
            csvwrite(['ecgstats_raw' PHYS.suffix '.txt'], stats.ecgall);
        case 'scl'
%             csvwrite(['sclstats_base' PHYS.suffix '.txt'], stats.scl);
            csvwrite(['sclstats_raw' PHYS.suffix '.txt'], stats.sclall);
        case 'resp'
%             csvwrite(['sclstats_base' PHYS.suffix '.txt'], stats.scl);
            csvwrite(['respstats_raw' PHYS.suffix '.txt'], stats.respall);
        case 'emg'
            csvwrite(['emgstats_base' PHYS.suffix '.txt'] , stats.emg);
            csvwrite(['emgstats_subnames' PHYS.suffix '.txt'], stats.subs);
        case 'ecghrv'
            csvwrite(['ecghrvstats_base' PHYS.suffix '.txt'], stats.ecghrv);
            strcell2file(['ecghrvstats_labels' PHYS.suffix '.csv'], stats.labels);
        case 'bpicg'
             havebp = 1;
             haveicg = 1;
             try 
                csvwrite(['bpstats' PHYS.suffix '.txt'], stats.bp);
             catch
                havebp = 0;
             end
             try
                csvwrite(['icgstats' PHYS.suffix '.txt'], stats.icg);
             catch
                haveicg = 0;
             end
        case 'icg'
            csvwrite(['icgstats_base' PHYS.suffix '.txt'], stats.icg);
        case 'ecgrec'      
            csvwrite('ecgrecstats.txt', stats.ecgrec);
            csvwrite('ecgrecraw.txt', stats.ecgrecraw);
        case 'bp2'
            csvwrite('bpstats.txt', stats.bp);
            csvwrite('bpraw.txt', stats.bpall);
            
    end
end









