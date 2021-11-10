% MATLAB Code for automated pre-processing and cleaning of EEG data from our paper "EEG-based clusters differentiate psychological distress, sleep quality and cognitive function in adolescents"

% Many thanks to my supervisor Paul Schwenn and to Makoto Miyakoshi (Makoto's preprocessing pipeline https://sccn.ucsd.edu/wiki/Makoto's_preprocessing_pipeline )



% 300920 v1.5 - updating for batch processing. Adding dipole fitting and IC
% rejection based on dipole residual variance and (inside brain) - copied
% from Paul S code "ti_eeg_clean"

% commented out [ALL EEG ...] lines at 61 and slightly below - not sure
% what they were doing...


%300720 v1.4 - updated labels of steps for write up
% 1.4 - step 3 line 80: adapted line noise frequencies and harmonics for Australia
% (50Hz,100,150,200)

% To run on a batch of files rather than just 1... set up lines 5-20

%####### Set directory
%cd('/Users/Current/Desktop/LABS_Data/PilotData_TI/Thompson_Institute_Pilot_EEG/EO_only')
cd('/Volumes/Owen_SSD/t1_AO/raw_bdfs_58')
rawDataFiles = dir('*.bdf');

outpath = '/Volumes/Owen_SSD/t1_AO/cleaned_58_pl_1_5'



% Channel location path, search for standard-10-5-cap385.elp in your eeglab plugins path e.g.
% C:\\Matlab\\eeglab\\plugins\\dipfit\\standard_BESA\\standard-10-5-cap385.elp
    chan_loc_path = '/Users/Current/Desktop/EEG_tools_Data/EEGLAB/eeglab2019_1/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp' 
    
    hdm_path = '/Users/Current/Desktop/EEG_tools_Data/EEGLAB/eeglab2019_1/plugins/dipfit/standard_BEM/standard_vol.mat'; % '[EEGLABroot]/eeglab/plugins/dipfit2.3/standard_BEM/elec/standard_1005.elc';
    bem_std_path = '/Users/Current/Desktop/EEG_tools_Data/EEGLAB/eeglab2019_1/plugins/dipfit/standard_BEM/standard_mri.mat'; % 'C:\\Work\\Matlab\\eeglab\\plugins\\dipfit\\standard_BEM\\standard_mri.mat'
   

 
    

% ------------ NOTES FROM CLEANING T1 Eyes Closed ----------------
% 300920 For some reason, MATLAB sees twice as many files in the directory, with
% the first half weirdly having the leading characters '._' in the
% filename. So here I am counting from n_files + 1 for this for loop (here
% 62 files for EC at baseline so going from 62)



 % Step 2: High-pass filter the data at 1-Hz. Note that EEGLAB uses pass-band edge, therefore 1/2 = 0.5 Hz.
      %230420 was getting an error saying I can't use 1650 as filter order
      %as minimum is 1690, recommended 3380. Moved to 3380 as per eeglab
      %recommendation
%    EEG = pop_eegfiltnew(EEG, 1, 0, 3380, 0, [], 0);
%Removing 6 channel(s)...
%eeg_insertbound(): 2 boundary (break) events added.
%eeg_insertbound(): event latencies recomputed and 389 events removed.
%eeg_checkset note: upper time limit (xmax) adjusted so (xmax-xmin)*srate+1 = number of frames
%Warning: converting all event types to strings
%eeg_checkset warning: 1/4 events had out-of-bounds latencies and were removed
%eeg_insertbound(): 2 boundary (break) events added.
%eeg_insertbound(): event latencies recomputed and 390 events removed.
%BUG 1971 WARNING: IF YOU ARE USING A SCRIPT WITTEN FOR A PREVIOUS VERSION OF EEGLAB (<2017)
%TO CALL THIS FUNCTION, BECAUSE YOU ARE REJECTING THE ONSET OF THE DATA, EVENTS MIGHT HAVE
%BEEN CORRUPTED. EVENT LATENCIES ARE NOW CORRECT (SEE https://sccn.ucsd.edu/wiki/EEGLAB_bug1971);
%Warning: Discrepency when recomputing event latency.
%Try to reproduce the problem and send us your dataset 
%> In eeg_eegrej (line 173) 
%Data break detected and taken into account for resampling
%resampling data 512.0000 Hz
%1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 
%resampling event latencies...
%resampling finished



% ------------% ------------% ------------% ------------% ------------% ------------




for subjID = 59:116 
    loadName = rawDataFiles(subjID).name;
    dataName = loadName(1:end-4);
 
 %####### Import data.
    EEG = pop_biosig(loadName);
    EEG.setname = dataName


% MAKE SURE DOUBLE PRECISION (64-bit number storage) IS SET IN EEGLAB: 
% File --> Preferences --> Memory options --> un-tick first box


% OLD single versionfullpath = '/Users/Current/Desktop/LABS_Data/PilotData_TI/Thompson_Institute_Pilot_EEG/1234/1234_EC.bdf'


% Get file parts and name
%    [filepath,name,ext] = fileparts(fullpath);
%    EEG.setname=name;

% Start EEGLAB
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG.etc.eeglabvers = '2019.1'; % this tracks which version of EEGLAB is being used, you may ignore it

% Import data
    EEG = pop_biosig(loadName, 'importannot','off','ref',[33 34] ,'refoptions',{'keepref' 'off'});
    %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',name,'gui','off'); 

% load channel locations
    EEG=pop_chanedit(EEG, 'lookup', chan_loc_path);
    EEG = eeg_checkset( EEG );

% Select EEG channels only
    EEG = pop_select( EEG, 'channel',{'Fp1' 'AF3' 'F7' 'F3' 'FC1' 'FC5' 'T7' 'C3' 'CP1' 'CP5' 'P7' 'P3' 'Pz' 'PO3' 'O1' 'Oz' 'O2' 'PO4' 'P4' 'P8' 'CP6' 'CP2' 'C4' 'T8' 'FC6' 'FC2' 'F4' 'F8' 'AF4' 'Fp2' 'Fz' 'Cz'});
    %[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[name ' EEG chs '],'gui','off'); 
    all_chanlocs = EEG.chanlocs;

  % Trim leading and trailing data corresponding to pre-start and trailing end
    % of session. See the EEG.event structure for more info.
    event_id_start = 1;
   event_id_stop = 2;
    event_idxs = [EEG.event.latency];  % Get all event indexes
    trim_idx = event_idxs(([EEG.event.type]==event_id_start | [EEG.event.type]==event_id_stop)); % start_event=1 stop_event=2
    % reject porition of continuous data in an EEGLAB dataset with
    % eeg_eegrej. Takes an array of regions tp reject each region defined
    % by start_idx and stop_idx to trim out.
    EEG = eeg_eegrej( EEG, [1 trim_idx(1); trim_idx(2) length(EEG.times)]);
    EEG = eeg_checkset( EEG );
    
    
 
    % Makoto...
    
    % Step 1: Downsample the data to 512 Hz (half the original at 1024) to
    % make ICA run faster
    EEG = pop_resample(EEG, 512);
    %EEG = pop_resample(EEG, 250, 0.8, 0.4); % For SIFT to ease antialiasing filter slope.
 
      % Step 2: High-pass filter the data at 1-Hz. Note that EEGLAB uses pass-band edge, therefore 1/2 = 0.5 Hz.
      %230420 was getting an error saying I can't use 1650 as filter order
      %as minimum is 1690, recommended 3380. Moved to 3380 as per eeglab
      %recommendation
    EEG = pop_eegfiltnew(EEG, 1, 0, 3380, 0, [], 0);
    
     % Step 3: Run cleanLineNoise (requires PREP Pipeline toolbox) The original CleanLine was replaced (25/04/2020).
    signal      = struct('data', EEG.data, 'srate', EEG.srate);
    lineNoiseIn = struct('lineNoiseMethod', 'clean', ...
                         'lineNoiseChannels', 1:EEG.nbchan,...
                         'Fs', EEG.srate, ...
                         'lineFrequencies', [50 100 150 200],...
                         'p', 0.01, ...
                         'fScanBandWidth', 2, ...
                         'taperBandWidth', 2, ...
                         'taperWindowSize', 4, ...
                         'taperWindowStep', 1, ...
                         'tau', 100, ...
                         'pad', 2, ...
                         'fPassBand', [0 EEG.srate/2], ...
                         'maximumIterations', 10);
    [clnOutput, lineNoiseOut] = cleanLineNoise(signal, lineNoiseIn);
    EEG.data = clnOutput.data;
 
    % Step 4: Apply clean_rawdata() to reject bad channels and correct continuous data using Artifact Subspace Reconstruction (ASR). Note this code specifies fixed amount of RAM (here 8GB, beware the large size) which disables adaptive RAM assignment to guarantee result reproducibility.
    originalEEG = EEG;
    EEG = clean_rawdata(EEG, 5, -1, 0.85, 4, 20, 0.25) %, 'availableRAM_GB', 8);  %% commented out 'availableRAM_GB' as it was causing an error "Too many arguments" for clean_rawdata
 
    
    
    % 2808 Check here how many channels have been rejected
    
    
    % Step 5: Interpolate all the removed channels
    EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');
 
    % Step 6: Re-reference the data to average
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1,:) = zeros(1, EEG.pnts);
    EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref(EEG, []);
    EEG = pop_select( EEG,'nochannel',{'initialReference'});
 
    % Step 7: Run AMICA using calculated data rank with 'pcakeep' option
    if isfield(EEG.etc, 'clean_channel_mask')
        dataRank = min([rank(double(EEG.data')) sum(EEG.etc.clean_channel_mask)]);
    else
        dataRank = rank(double(EEG.data'));
    end
    runamica15(EEG.data, 'num_chans', EEG.nbchan,...
        'outdir', ['/Volumes/Owen_SSD/t1_EC/cleaned_pl_1_5/amica_out_1_5'],...  % need to change this for batch processing
        'pcakeep', dataRank, 'num_models', 1,...
        'do_reject', 1, 'numrej', 15, 'rejsig', 3, 'rejint', 1);
    EEG.etc.amica  = loadmodout15(['/Volumes/Owen_SSD/t1_EC/cleaned_pl_1_5/amica_out_1_5/']); % need to change this for batch processing
    %npcs = size(EEG.etc.amica.W,1)
    EEG.etc.amica.S = EEG.etc.amica.S(1:EEG.etc.amica.num_pcs, :); % Weirdly, I saw size(S,1) be larger than rank. This process does not hurt anyway.
    EEG.icaweights = EEG.etc.amica.W;
    EEG.icasphere  = EEG.etc.amica.S;
    EEG = eeg_checkset(EEG, 'ica');
    
    

    
    
    % Step 8: From AMICA, dentify ICs for rejection using ICLabel scores
    % (keeping only components labeled >70% 'brain')
    EEG = iclabel(EEG, 'default'); % Original output from eegh has a bug: the second input ('default') is without ''--this is fixed here.
    brainIdx  = find(EEG.etc.ic_classification.ICLabel.classifications(:,1) >= 0.7); % > 70% brain == good ICs.
    
    
    % Step 8.2 - include Dipoles
    
    EEG = pop_dipfit_settings( EEG, 'hdmfile',hdm_path,'coordformat','MNI',...
        'mrifile',bem_std_path,'chanfile',chan_loc_path,...
        'coord_transform',[0.83215 -15.6287 2.4114 0.081214 0.00093739 -1.5732 1.1742 1.0601 1.1485] ,...
        'chansel',[1:32]);
    
    EEG = pop_multifit(EEG, 1:EEG.nbchan,'threshold', 100, 'dipplot','off','plotopt',{'normlen' 'on'});

    % Step 12: Search for and estimate symmetrically constrained bilateral dipoles
%     EEG = fitTwoDipoles(EEG, 'LRR', 35);

    % Save the dataset
    EEG = pop_saveset( EEG, 'filename', [dataName '_ica'], 'filepath', outpath);
    EEG_ica = EEG;
    
    % Perform IC rejection using ICLabel scores and r.v. from dipole fitting.
    EEG = iclabel(EEG, 'default'); % Original output from eegh has a bug: the second input ('default') is without ''--this is fixed here.
    brainIdx  = find(EEG.etc.ic_classification.ICLabel.classifications(:,1) >= 0.70); % > 70% brain == good ICs.
 
    % Perform IC rejection using residual variance of the IC scalp maps.
    rvList    = [EEG.dipfit.model.rv];
    goodRvIdx = find(rvList < 0.15)'; % < 15% residual variance == good ICs.
 
    % Perform IC rejection using inside brain criterion.
    load(EEG.dipfit.hdmfile); % This returns 'vol'.
    dipoleXyz = zeros(length(EEG.dipfit.model),3);
    for icIdx = 1:length(EEG.dipfit.model)
        dipoleXyz(icIdx,:) = EEG.dipfit.model(icIdx).posxyz(1,:);
    end
    depth = ft_sourcedepth(dipoleXyz, vol);
    depthThreshold = 1;
    insideBrainIdx = find(depth<=depthThreshold);
 
    % Take AND across the three criteria.
    goodIcIdx = intersect(brainIdx, goodRvIdx);
    goodIcIdx = intersect(goodIcIdx, insideBrainIdx);
 
    
    
    
     % Step 9. Perform IC rejection.
     
    % Perform IC rejection.
    EEG = pop_subcomp(EEG, goodIcIdx, 0, 1);
     
   
    % EEG = pop_subcomp(EEG, brainIdx, 0, 1);
    
   
 
    % Step 10. Post-process to update ICLabel data structure.
    EEG.etc.ic_classification.ICLabel.classifications = EEG.etc.ic_classification.ICLabel.classifications(brainIdx,:);
 
    % Step 11. Post-process to update EEG.icaact.
    EEG.icaact = [];
    EEG = eeg_checkset(EEG, 'ica');
    
    
    % NB at this stage skipping Makoto recommended steps of using dipole information for
    % IC rejection
    
    
    % Step 12: save file
    
    EEG = pop_saveset( EEG, 'filename', dataName, 'filepath', outpath);
    
    %NB 100820 - changed second setting in EEGLAB preferences to disable
    %saving separate .set and .fdt files so it is all saved as 1 .set file
    
    % add wavelet and multitaper here
    
    
    end %(end of for loop)
    