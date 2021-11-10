% MATLAB Code to apply multitaper analysis to get spectrograms from cleaned EEG files - from our paper "EEG-based clusters differentiate psychological distress, sleep quality and cognitive function in adolescents"


%v 1.3 061020 - trying out MT with 3 tapers (2hz res, 2sec windows) on all
% baseline, eyes closed, cleaned data

% see parameter advice at https://www.researchgate.net/post/How_do_you_choose_your_taper_value_for_multi-taper_spectral_analysis
% also here under "Procedure 2"... https://journals.physiology.org/doi/full/10.1152/physiol.00062.2015

% Multitaper pipeline
% Set directory
cd('/Volumes/Owen_SSD/t1_EO/cleaned_pl_1_5')
eeglab()
rawDataFiles = dir('*.set');


for subjID = 1:length(rawDataFiles) 
    loadName = rawDataFiles(subjID).name;
    dataName = loadName(1:end-4); %strip .set or .bdf off end of filename
 
% ####### Import data.
    EEG = pop_loadset(loadName);
    EEG.setname = dataName

    % W, t and p define the taper parameters for the multitaper analysis
    % W = bandwidth, T = duration of sampling window, p is an integer such that 2TW-p tapers are used

    W = 2 %2 hz bandwidth
    movingwin = [2 1] % 2 second windows, 1 second step size/overlap
    t_win = movingwin(1);  % window in s
    p = 1;
    chronux_params = struct();
    chronux_params.Fs = 512; % Sampling frequency - set to twice max frequency of interest
    chronux_params.tapers = [W t_win p]; %  W, T, p
    chronux_params.fpass = [0 256]; % frequency band to be used [minfreq maxfreq]
    chronux_params.err = [2 p];% Jackknife error bars
    
    EEG.data = transpose(EEG.data) %Chronux wants time in columns, channel in rows - EEG set is initially transpose of this (originally time by channel)


    [S,t,f,Serr] = mtspecgramc(EEG.data(:,:), movingwin, chronux_params)
    
    % Save multitaper spectrogram output to file
    file_for = sprintf('%s_mt_1_3_cl_1_5.mat',dataName) %%%%%%%!!!!!!!!!!@@@@@@@@@ CHange filename here
    save(file_for,'S','t','f','Serr') % change file path here
    
    % In future - can specify full filepath to save output to a particular directory
    % e.g. save('C:\My folder\filename','varname'
end


% Plot 32 spectrograms
%figure
%for r=1:32
%    axes(r) = subplot(4,8,r);
%    imagesc(t/60.,f,10*log10(S(:,:,r)'))
%    xlabel('Time [min]')
%    ylabel('Frequency [Hz]')
%    title("")
%    ylim([0 256])
%    set(gca,'YDir','normal')
%    colormap jet   
%end
%linkaxes(axes, 'xy')
