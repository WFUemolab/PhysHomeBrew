function [stimtime, stimcat] = CardiacCycleStim(batch, varargin)
% function [stimtime, stimcat] = CardiacCycleStim('Batch.m.m', [ECG
% directory name])

if isempty(dir('anslabdef.m'))
    cd ..
end


expdir = pwd;

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

if ~isempty(varargin)
    ecgdir = varargin{1};
else
    ecgdir = 'ECG';
end

for i = 1:numsubs
    T = anslabreadtiming([expdir '/raw/' SubName{i} '.m']);
    TOrig = T(:,2) .* 1000;

    load(['ecg/' ecgdir '/' SubName{i}, '00.mat']);
    rtOrig = rt .* 2.5;


    for j = 1:size(TOrig,1)

        stimtmp = max(find(rtOrig <= TOrig(j)));

        stimtime(j,i) = TOrig(j) - rtOrig(stimtmp);

        if stimtime(j,i) < 200
                stimcat(j,i) = 0;
        elseif stimtime(j,i) <= 400
                stimcat(j,i) = 1;
        elseif stimtime(j,i) < 450
                stimcat(j,i) = 0;
        elseif stimtime(j,i) <= 800
                stimcat(j,i) = 2;
        else
                stimcat(j,i) = 0;
        end

    end
end

