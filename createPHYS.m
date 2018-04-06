function PHYS = createPHYS(varargin)
%function PHYS = createPHYS('EMG', 'SCL', 'ECG')
%
% creates a PHYS struct that can be inputted into cw_get_stats and
% cw_getbasestats
%
% You can input (as an option) only those signals that you care about, EMG
% or SCL or ECG (including HRV). To process other signals, like BP or RESP,
% please contact Christian Waugh, waughce@wfu.edu
%
% To start:
% 
% 1.	Enter a suffix for this dataset (make sure it’s different from previous sets, or will overwrite them)
% a.	>> Suffix: _pilot
% 2.	Then enter the number of seconds for your baseline period (e.g., 60) if not part of your timing file, or ‘Trial’ if part of your timing file (this is the only time you need to enter quotes for a string)
% a.	>> Baseline: ‘Trial’
% 3.    If you chose 'Trial' you will then be asked to input the number in
% your timing file that corresponds to your baseline condition.
%       >> 0
%
% For ECG:
% 
% 1.	
% 4.	Enter next the precise channel name used for the ECG directory (can find by looking in ecg\CHANNELNAME)
% a.	Channel name: ECG
% 5.	Enter the type of analysis desired: either Trial for trial by trial analyses (i.e. short events like 4-8s long), or Segment for large epochs (epochs that are minutes long)	
% a.	Type: Segment
% 6.	IF YOU CHOSE SEGMENT:
% a.	Enter the duration of baseline to be used for calculating baseline/recovery, etc. (may be different from original baseline if you want only a portion of the baseline).
% i.	Baseline: 180
% b.	Same thing, but the duration of the recovery period used to analyze recovery
% i.	Recovery: 180
% c.	Enter the duration (s) of new segments if you want your big epochs to be resampled (e.g. if you have 5 minute epochs, but want 1 minute samples).
% i.	Duration: 60
% 
% 7.	IF YOU CHOSE TRIAL:
% Next, please think about your trial structure and determine when your trial starts, when your critical window for assessing the physiological response starts, and when your proximal baseline comparative window starts and ends (all in relation to start of trial, which will play a role in the trial extraction - how much of a pre- and post-trigger period you choose):
% this will be used to then calculate IBI accelearation and deceleration
% within that window (cw_ecg_calc)
% Examples:
% 1.	Baseline for 4s -> trial starts -> get physio from 1-4s after trial starts
% ?	base window = [-4 0]; trial starts = 4; response window = [1 4]
% 2.	Trial starts -> Baseline for 4s -> get physio from 1-4s after baseline
% ?	base window = [0 4]; trial starts = 0; response window = [5 8]
% 
% 
%  For EMG:
% 
% 1.	Enter next the precise channel name used for the EMG directory (can find by looking in emg\CHANNELNAME)
% a.	Channel name: Corr_filtered
% 2.	Enter the number of EMG channels (number of muscles that you will be examining)
% a.	Number: 2
% 3.	For each of those channels, enter the name of it (as noted in the .acq file)
% a.	Name 1: Corr_filtered
% b.	Name 2: Zyg_filtered
% 4.	Next, please think about your trial structure and determine when your trial starts, 
%       when your critical window for assessing the physiological response starts, 
%       and when your proximal baseline comparative window starts and ends 
%       (all in relation to start of trial, which will play a role in the trial extraction - 
%       how much of a pre- and post-trigger period you choose):
%
% Examples:
% 	Baseline for 4s -> trial starts -> get physio from 1-4s after trial starts
% ?	base window = [-4 0]; trial starts = 4; response window = [1 4]
% 	Trial starts -> Baseline for 4s -> get physio from 1-4s after baseline
% ?	base window = [0 4]; trial starts = 0; response window = [5 8]
% 
% 5.	Determine how you want to process each trial’s EMG data according
% to baseline or not
% a.	>> Proximal (each trial's proximal baseline is subtracted from its
% response)
%        >> z (each trial is standardized according to some pretask baseline
%        >> nobase (the average or max EMG is provided on each trial with no transformation)
% 6.    Determine whether you want to take the mean or max response on each
% trial
% b.	>> Mean – takes the average of the entire EMG response during the window
%   	>> Max – takes the max peak of the EMG response during the window
% 
% FOR SCL
% 1.	Enter next the precise channel name used for the EDA directory (can find by looking in eda\CHANNELNAME)
% a.	>> Channel name: SCL
% 2.	Enter the type of analysis desired: either Trial for trial by trial analyses (i.e. short events like 4-8s long), or Segment for large epochs (epochs that are minutes long)	
% a.	>> Type: Trial
% 3.	IF YOU CHOSE TRIAL:
% Next, please think about your trial structure and determine when your trial starts, when your critical window for assessing the physiological response starts, and when your proximal baseline comparative window starts and ends (all in relation to start of trial, which will play a role in the trial extraction - how much of a pre- and post-trigger period you choose):
% Examples:
% 	Baseline for 4s -> trial starts -> get physio from 1-4s after trial starts
% ?	base window = [-4 0]; trial starts = 4; response window = [1 4]
% 	Trial starts -> Baseline for 4s -> get physio from 1-4s after baseline
% ?	base window = [0 4]; trial starts = 0; response window = [5 8]
% 
% 4.	Next, enter the type of SCR response you would like to record:
% a.	>> ‘proximal’ – takes the average of the proximal baseline period and subtracts that from the response period. Then outputs the max response in that response window (or 0 if the max response is negative). 
% b.	>> ‘z’ – standardizes the responses in the response window of each trial to some overall baseline period.
% c.	>> ‘minmax’ – useful if no proximal baseline period. Finds the max response in the response window and subtracts from it the min response (that occurred before the max response).
% 5.	Enter whether you would like to transform the SCR, and then enter whether you’d like to either log transform it or sqrt it. 
% a.    >> 'sqrt' or 'log'
% 6     Duration of bins you want the SCL data to be reasampled into (s) or 1 if trial by trial
%
%       IF YOU CHOSE SEGMENT:
% 1.    Duration of bins you want the SCL data to be resampled into (s)
% e.g., 60


if nargin == 0
    chans = {'EMG' 'SCL' 'ECG'};
else
    chans = varargin;
end

numchans = length(chans);

PHYS.suffix = input('Enter a suffix for this analysis (e.g. _cleaned): ', 's');
PHYS.base = input('Enter number of seconds for baseline (or ''Trial'' if part of your timing file): ', 's'); %baseline time for EMG (s)
if strcmp(PHYS.base, 'Trial') == 1
    PHYS.baseTrig = input('Enter the trigger number for baseline in your timing file: '); %what trigger corresponds to baseline in your timing file
end

iter = 0;
while iter < numchans
    switch chans{iter+1}
        case 'EMG'
            
            PHYS.EMG.channame = input('Enter channel name for EMG (in ANSLAB directory): ', 's'); %channel name in directory

            PHYS.EMG.num = input('Enter number of EMG channels: ');
            for i = 1:PHYS.EMG.num
                fprintf(1, 'Enter the name of EMG signal %1d (ex. COR)', i);
                PHYS.EMG.names{i} = input(': ', 's');
            end
            
            fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
            fprintf(['\nNext, please think about your trial structure and determine\n']);
            fprintf(['\nwhen your trial starts, when your critical window for assessing\n']);
            fprintf(['\nthe physiological response starts, and when your proximal baseline\n']);
            fprintf(['\ncomparative window starts and ends (all in relation to start of trial.\n']);
            fprintf(['\nExamples:\n']);
            fprintf(['\nBaseline for 4s -> trial starts -> get physio from 1-4s after trial starts\n']);
            fprintf(['\nbase window = [-4 0]; trial starts = 4; response window = [1 4]\n']);
            fprintf(['\nTrial starts -> Baseline for 4s -> get physio from 1-4s after baseline\n']);
            fprintf(['\nbase window = [0 4]; trial starts = 0; response window = [5 8]\n']);
            fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
            
            PHYS.EMG.trialstart = input('Enter time that the trial starts - the time you will use for a reference for the baseline and response windows (s): '); % window for calculating average emg response
            
            PHYS.EMG.baseWindow = input('Enter window for processing proximal baseline EMG response relative to trial starting (e.g. [-4 0]; if appropriate - enter 0 if not): '); % window for calculating average emg response

            PHYS.EMG.responseWindow = input('Enter window for processing EMG response to stimuli relative to trial starting (e.g. [1 4]): '); % window for calculating average emg response

            PHYS.EMG.basetype = input('Enter analysis type to find EMG response (either proximal or z): ', 's');
            
            PHYS.EMG.meanmax = input('Use the mean EMG response or max EMG response?', 's');
            
            PHYS.EMG.binsize = input('What is the bin size for getting trial by trial data (i.e., not averaged across conditions; default = 1)?:');
            iter = iter + 1;
            
        case 'ECG'
            
            PHYS.ECG.channame = input('Enter channel name for ECG (in ANSLAB directory): ', 's'); %channel name in directory

            type = input('What type of ECG analysis? Enter - Trial - for trial by trial or - Segment - for long data phases (e.g. baseline, stress, etc.):', 's');
            
            
            
            if strcmp(type, 'Trial')
               fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
            fprintf(['\nNext, please think about your trial structure and determine\n']);
            fprintf(['\nwhen your trial starts, when your critical window for assessing\n']);
            fprintf(['\nthe physiological response starts, and when your proximal baseline\n']);
            fprintf(['\ncomparative window starts and ends (all in relation to start of trial.\n']);
            fprintf(['\nExamples:\n']);
            fprintf(['\nBaseline for 4s -> trial starts -> get physio from 1-4s after trial starts\n']);
            fprintf(['\nbase window = [-4 0]; trial starts = 4; response window = [1 4]\n']);
            fprintf(['\nTrial starts -> Baseline for 4s -> get physio from 1-4s after baseline\n']);
            fprintf(['\nbase window = [0 4]; trial starts = 0; response window = [5 8]\n']);
            fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
            
            PHYS.ECG.trialstart = input('Enter time that the trial starts - the time you will use for a reference for the baseline and response windows (s): '); 
            
            PHYS.ECG.baseWindow = input('Enter window for processing proximal baseline ECG response relative to trial starting (e.g. [-4 0]; if appropriate - enter [0 0] if not): ');

            PHYS.ECG.responseWindow = input('Enter window for processing ECG response to stimuli relative to trial starting (e.g. [1 4]): '); % window for calculating IBI accel/decel
 
            
            
                   elseif strcmp(type, 'Segment')
                PHYS.ECG.baserec = input('Duration of baseline period to be used ECG metrics (may or may not be same as regular baseline)?:' );
                
                PHYS.ECG.recwin = input('Duration of recovery period to be used for calculating ECG recovery metrics (may or may not be same as regular recovery)?:' );

                PHYS.ECG.binsize = input('Duration of new sub-segments for ECG data to be resampled into (in seconds; e.g. 60):');
            
            end
            iter = iter + 1;
            
        case 'SCL'
            
            PHYS.SCL.channame = input('Enter channel name for SCL (in ANSLAB directory): ', 's'); %channel name in directory
            
            type = input('What type of SCL analysis? Enter ''Trial'' for trial by trial or ''Segment'' for long data phases (e.g. baseline, stress, etc.):', 's');
             
            if strcmp(type, 'Trial')
                PHYS.SCL.trialstart = input('Enter time that the trial starts - the time you will use for a reference for the baseline and response windows (s): '); % window for calculating average emg response
            
                PHYS.SCL.baseWindow = input('Enter window for processing proximal baseline SCR response relative to trial starting (e.g. [-4 0]; if appropriate - enter [0 0] if not): '); % window for calculating average emg response

                PHYS.SCL.responseWindow = input('Enter window for processing SCR response to stimuli relative to trial starting (e.g. [1 4]): '); % window for calculating average emg response
 
                PHYS.SCL.scrtype = input('Enter analysis type to find SCRs (either proximal or z or minmax): ', 's');
            
                transform = input('Transform the SCR?:', 's');
                    if strcmp(transform, 'y') || strcmp(transform, 'yes')
                        PHYS.SCL.transform = input('Enter transformation method (log or sqrt):','s');
                    end
                   
                PHYS.SCL.binsize = input('Duration of new sub-segments for SCL data to be resampled into (e.g. 60 for segments, enter 1 for trial by trial):');
                        
            elseif strcmp(type, 'Segment')
            
                PHYS.SCL.binsize = input('Duration of new sub-segments for SCL data to be resampled into (e.g. 60):');
            end
            
            
            iter = iter + 1;
    end
end

save PHYS PHYS