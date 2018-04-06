% edits to transacqhisteventext.m necessary to run PhysHomeBrew files
% ANSLAB ver 2.51
% its purpose is to save all the trial information in a _Events file to be
% read in later


lines 161-163 should be:

 %CellFilePath = [InFilePath,'.Event.Cell.ans']; CW Edits
 [FileDir,FilePath] = fileparts(InFilePath); % CW Edits
 CellFilePath = [FileDir, filesep, FilePath, '_Events']; % CW Edits
 
 
 lines 2441-2445 should be:
 
 %==================================================================
%% segmentationtype 2: save trial data file (.Event.Cell.ans)
%================================================================== 
                %save(CellFilePath,'DATA','HDR');
                save(CellFilePath,'DATA'); %CW edits