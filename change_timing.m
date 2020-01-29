function change_timing(batch, inc, varargin)
% function change_timing('STR01201.m', 60, 'rebin')
% function change_timing('Batch.m', 60, 'rebin')
% function change_timing('STR01*.acq', 60, 'rebin')
    
% This function changes the timing file to reflect
% the initial timing + some increment that you specify
% also includes other functions that are study-specific
% to add your own, just add a function at the end of the script and list it as a 'case' below
%
% Example with inc = .5:
% Old: T = [...
%            1  204  210  6
% New: T = [...
%            1  204  204.5  .5
%
% batch = either batch file or a string of subj names to search 
%(e.g. 'STR02*.acq')
% inc = new increment



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


for i = 1:numsubs
    T = anslabreadtiming([expdir '/raw/' SubName{i} '.m']);
    switch varargin{1}
        case 'aer'
            T = aer(T);
        case 'addtrials'
        %T = AddTrials(T);
        T = lasttwotrials(T);       
        case 'change_conds'
        T = change_conds(T, inc, varargin{i+1});
        case 'rebin'
        T = rebin(T, inc);
        case 'newinc'
        T(:,3) = T(:,2) + inc; T(:,4) = inc;
        case 'fpfrebin'
        T = fpfrebin(T);
        case 'bartdur'
           anslabwritetiming([expdir '/raw/' SubName{i} '.old.m'], T);
            Tnew = T;
            for j = 1:size(T,1)
                 if j ~= size(T,1)
                     Tnew(j,3) = T(j+1,2);
                     Tnew(j,4) = T(j+1,2) - T(j,2);
                 else
                     Tnew(j,3) = T(j,2)+4;
                     Tnew(j,4) = 4;
                 end
                switch T(j,1)
                    case 4
                        Tnew(j,1) = 0;
                    case 1
                        newend = Tnew(j,3) - Tnew(j,2) - 4;
                        
                        Tnew(end+1,:) = [2, Tnew(j,2) + 4, Tnew(j,3), newend];
                        
                    case 2
                        Tnew(j,1) = 4;
                        Tnew(end+1,:) = [3, Tnew(j,2) - 3, Tnew(j,2), Tnew(j-1,4)-4];
                end
            end
            T = sortrows(Tnew, 2);
    end   
    anslabwritetiming([expdir '/raw/' SubName{i} '.m'], T);
    anslabwritetiming([expdir '/raw/' SubName{i} '.spectral.m'], T);
    anslabwritetiming([expdir '/raw/' SubName{i} '.icg.m'], T);
end

% function Tnew = AddTrials(T)
% iter = 0;
% Tnew = T(1,:); 
% for i = 2:size(T,1)
%     if T(i-1,1) < 11 
%         Tnew(i+iter,:) = T(i,:); 
%     else
%         
%         Tnew(i+iter,:) = [T(i-1,1), T(i-1,2)+16.1, T(i-1,3), 4.05];
%         Tnew(i+iter+1,:) = T(i,:);
%         iter = iter + 1; 
%        
%             Tnew(i+iter,:) = [T(i-1,1), T(i-1,2)+16.1, T(i-1,3), 4.05]; 
%         Tnew(i+iter+1,:) = [T(i-1,1), T(i-1,2)+32.2, T(i-1,3)+16.1, 4.05]; 
%         Tnew(i+iter+2,:) = T(i,:);
%         iter = iter + 2; 
%         end
%     end
% end


function T = lasttwotrials(T)
T(11,:) = [11, 0, 60, 60];
T(12,:) = [12, 60, 120, 60];
T(13,:) = [13, 120, 180, 60];
% T(67,:) = [4, T(66,2)+16.1, T(66,3), 16.1];
% T(68,:) = [5, T(66,3), T(66,3)+16.1, 16.1];
end

function T = change_conds(T, inc, numsegs)
    segs = inc:inc:numsegs*inc;
    
    for j = 1:numsegs
        condrange{j} = find(T(:,2) < segs(j));
    end
    
    for k = 1:numsegs
        if k == 1
            T(condrange{k},1) = 1;
        else
            T(~ismember(condrange{k}, condrange{k-1}),1) = k;
        end
    end
end

function Tnew = aer(T)
    c = 1;
    Tnew(1,1) = c;
    Tnew(1,2) = T(1,2);
    Tnew(1,3) = T(1,2) + 60;
    Tnew(1,4) = 60;
    for k = 1:2
        c = c + 1;
        Tnew(c,1) = c;
        Tnew(c,2) = Tnew(c-1,2) + 60;
        Tnew(c,3) = Tnew(c,2) + 60;
        Tnew(c,4) = 60;
    end
    c = c+1;
    Tnew(c,1) = c;
    Tnew(c,2) = T(2,2);
    Tnew(c,3) = T(2,2) + 60;
    Tnew(c,4) = 60;
    for k = 1:9
        c = c + 1;
        Tnew(c,1) = c;
        Tnew(c,2) = Tnew(c-1,2) + 60;
        Tnew(c,3) = Tnew(c,2) + 60;
        Tnew(c,4) = 60;
    end
    c = c+1;
    Tnew(c,1) = c;
    Tnew(c,2) = T(3,2);
    Tnew(c,3) = T(3,2) + 60;
    Tnew(c,4) = 60;
     for k = 1:2
        c = c + 1;
        Tnew(c,1) = c;
        Tnew(c,2) = Tnew(c-1,2) + 60;
        Tnew(c,3) = Tnew(c,2) + 60;
        Tnew(c,4) = 60;
     end
end
    

function Tnew = rebin(T, binsize)
        for i = 1:size(T,1)
            newrows = T(i,4)/binsize;
            for j = 1:newrows
                if i == 1 && j == 1
                    Tnew(1,1) = T(1,1);
                    Tnew(1,2) = T(1,2);
                    Tnew(1,3) = T(1,2) + binsize;
                    Tnew(1,4) = binsize;
                else
                    currrow = i*newrows-(newrows-j);
                    Tnew(currrow,1) = currrow;
                    Tnew(currrow,2) = Tnew(currrow-1,3);
                    Tnew(currrow,3) = Tnew(currrow,2) + binsize;
                    Tnew(currrow,4) = binsize;
                end
            end
        end
 end 
        
        function Tnew = fpfrebin(T)
            
            k = 0;
            for i = 1:size(T,1)
                if T(i,1) > 20
                    j = i+k;
                    Tnew(j,1) = T(i,1)*10;
                    Tnew(j,2) = T(i,2);
                    Tnew(j,3) = T(i,2) + 60;
                    Tnew(j,4) = 60;
                    
                    Tnew(j+1, 1) = Tnew(j,1) + 1;
                    Tnew(j+1, 2) = Tnew(j,2) + 60;
                    Tnew(j+1, 3) = Tnew(j+1, 2) + 60;
                    Tnew(j+1, 4) = 60;
                    
                    Tnew(j+2, 1) = Tnew(j+1,1) + 1;
                    Tnew(j+2, 2) = Tnew(j+1,2) + 60;
                    Tnew(j+2, 3) = Tnew(j+2, 2) + 60;
                    Tnew(j+2, 4) = 60;
                    
                    k = k + 2;
                else
                    j=i+k;
                    
                    Tnew(j,1) = T(i,1);
                    Tnew(j,2) = T(i,2);
                    Tnew(j,3) = T(i,2) + 60;
                    Tnew(j,4) = 60;
                    
                    
                end
            end
      end
      
function Tnew = addur(T)
                
for i = 1:(size(T,1) - 1)
    Tnew = T;
    Tnew(i,3) = T(i+1,2);
    Tnew(i,4) = T(i+1,2) - T(i,2);
end
end            
     
   
end    
   
    
            
    

    
