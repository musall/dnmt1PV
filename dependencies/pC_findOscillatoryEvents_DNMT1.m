function out = pC_findOscillatoryEvents_DNMT1_SM(recInfo, recLabels, optsIn)
%%
if ~exist('optsIn', 'var')
    optsIn = [];
end

checkTraces = true;

opts.showPlots = true;
opts.winSize = [0 10];
opts.verbose = true;
opts.reload = false;
opts.savePath = 'D:\RD10_project\';
opts.targDepth = 500;
opts.brainThresh = 1;
opts.eventThresh = 2;
opts.eventFreqRange = [1, 3];
opts.loadSC = false;
opts.scChan = nan;
opts = copyCommonFields(optsIn, opts); %check input opts and update fields


%% get data
recName = recInfo{strcmpi(recLabels, 'Folder')};
recPath = recInfo{strcmpi(recLabels, 'Path')};
probNr = recInfo{strcmpi(recLabels, 'Probe')};
% useChan = str2num(recInfo{strcmpi(recLabels, 'v1Channel')});

lfpDir = [fullfile(recPath, recName) filesep recName '_g0' filesep  recName '_g0_imec' probNr];
% lfpDir = strrep(lfpDir, 'naskampa.kampa-10g', 'naskampa');
lfMetaName = dir(fullfile(lfpDir,'*.lf.meta'));
lfMeta = pC_readSpikeGLXmeta(fullfile(lfMetaName.folder,lfMetaName.name));
lfFile = dir(fullfile(lfpDir,'*.lf.bin'));
lfFile = fullfile(lfpDir, lfFile.name);

savePath = [opts.savePath filesep recName filesep];
saveFile = fullfile(savePath, 'DNMT1_oscillations.mat');
if opts.loadSC
    saveFile = fullfile(savePath, 'DNMT1_oscillations_SC.mat');
end

fprintf('Current path: %s. ', recName);
if ~opts.reload && exist(saveFile, 'file')
    
    %% get data from local save file
    disp('Loading local data ... ');
    load(saveFile, 'out');
        
else
    
    if ~opts.loadSC
    %% get mean standard deviations of different parts of the recording to determine the brain surface.
    % check for bad channels results - tends to be more accurate if saline was applied
    badChansFile = fullfile(lfpDir, 'bad_channel_ids_list_fixed.npy');
    if ~exist(badChansFile, 'file')
        badChansFile = fullfile(lfpDir, 'bad_channel_ids_list.npy');
    end
    
    chanSteps = 1 : lfMeta.nChans-1; %spacing between channels when making spatial figure
    if exist(badChansFile, 'file')
        badChans = readNPY(badChansFile);
        if min(badChans) == 191
            badChans(badChans == 191) = []; %ignore the reference channel
        end
        
        % find brain surface based on last consecutive block of bad channels
        badChans = flipud(badChans+1);
        cIdx = min([find(diff(badChans) < -1, 1), length(badChans)]);
        cBrainIdx = true(1, lfMeta.nSites);
        if ~isempty(badChans)
            cBrainIdx(badChans(cIdx):end) = false;
        end
        cBrainIdx = fliplr(cBrainIdx);
        cBrainIdx = cBrainIdx(chanSteps);
        meanStd = cBrainIdx';
        
    else
        % if badchannels are not available use LFP method
        meanStd = medfilt1(pC_findBrainSurface_lfp(lfFile, lfMeta));
        meanStd = flipud(meanStd(chanSteps));
        meanStd = meanStd ./ std(meanStd);
        meanStd = meanStd - meanStd(1);
        cBrainIdx = meanStd' > opts.brainThresh; %contacts in the brain
    end
    
    depthRange = (1 : length(cBrainIdx)) *10;
    brainStart = find(cBrainIdx, 1); %cortical surface
    depthRange = depthRange - depthRange(brainStart);
    useChan = lfMeta.nChans -1 - find(depthRange >= opts.targDepth, 1); %get channel, based on opts.targDepth
    fprintf('Target depth = %i, Using channel: %i\n', opts.targDepth, useChan);
    
    else
        useChan = opts.scChan;
    end
    
    %get current channel
    lfpData = pC_extractAnalogChannel(lfFile, lfMeta.nChans, useChan, [], [], opts);
    lfpData = single(lfpData) .* lfMeta.uV_per_bit_lfp ./ 1000;
    
    
    %% spectral analysis
    fs = round(lfMeta.imSampRate); % in Hz
    lfpData(1:round(300*fs)) = 0; %exclude first 5 minutes because of noise in some recordings

    % Window size in seconds and convert it to samples
    windowSize_sec = 2; % seconds
    windowLength = windowSize_sec * fs; % in samples
    
    % Step size in seconds and convert it to samples
    stepSize_sec = 0.1; % seconds
    noverlap = windowLength - stepSize_sec * fs; % Number of overlapping samples
    
    %% Desired frequency resolution
    freqResolution = 0.1; % in Hz
    
    % Frequency range of interest
    freqMin = 0.5; % Hz
    freqMax = 10;  % Hz
    
    % Calculate nfft to achieve the desired frequency resolution
    nfft = fs / freqResolution; % nfft = 2500 / 0.5 = 5000
    
    % Compute the spectrogram
    [S, F, T] = spectrogram(lfpData, windowLength, noverlap, nfft, fs);
    
    % Find indices of frequencies between 0.5 Hz and 10 Hz
    freqIndices = (F >= freqMin) & (F <= freqMax);
    
    % Extract the desired frequency components
    S_desired = S(freqIndices, :);
    F_desired = F(freqIndices);
    
    %% Plot the spectrogram
    cData = 10*log10(abs(S_desired));
    medianPower = nanmedian(cData, 2);
    cData = cData - medianPower;
    cData(isinf(cData)) = 0;
    cData = smoothCol(cData, 2, 150, 'exp');
    
    if opts.showPlots
        figure;
        imagesc(T, F_desired, cData);
        axis xy; % Ensure the y-axis is in the correct direction
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title('Spectrogram of lfpData from 0.5 Hz to 10 Hz');
        colorbar;
        caxis([0 10]);
    end
    
    out.medianPower = medianPower;
    out.specFreq_all = F_desired;
    out.specPower_all = 10*log10(abs(S_desired));
    
    %% filter spectral power and find events
    Wn = 0.0005;               % Normalized cutoff frequency
    filterOrder = 500;         % Filter order
    
    % Design the FIR high-pass filter
    b = fir1(filterOrder, Wn, 'high', kaiser(filterOrder + 1, 3));
    
    % select frequency range to detect events. Has to be be within 0.5 and 10 Hz
    freqIdx = F_desired > opts.eventFreqRange(1) & F_desired < opts.eventFreqRange(2);
    
    % Step 3: Apply the High-Pass Filter
    eventTrace = nanmean(cData(freqIdx, :));
    eventTrace = detrendTrace(eventTrace, round(size(eventTrace,2) / 10));
    eventTrace = filtfilt(b, 1, double(eventTrace));
    eventTrace(1:400) = 0; %skip first part of recording
    eventTrace = eventTrace ./ std(eventTrace);
    eventTrace = eventTrace - prctile(eventTrace,30);
    eventTrace(1:400) = 0; %skip first part of recording
    eventTrace = smoothCol(eventTrace, 2, 20, 'exp');
    
    cEvents = pC_digitalToTimestamp(eventTrace > opts.eventThresh, 1 / mean(diff(T)));
    cEvents{1} = cEvents{1} + T(1) - mean(diff(T));
    cEvents{2} = cEvents{2} + T(1) - mean(diff(T));
    cEvents{3} = cEvents{3} + T(1) - mean(diff(T));
    
    if length(cEvents{2}) > length(cEvents{3})
        cEvents{2}(end) = [];
    end
    eventDurations = (cEvents{3} - cEvents{2});
    useEventsIdx = eventDurations > abs(diff(opts.winSize));
    out.eventDurations = eventDurations(useEventsIdx);
    
    useEvents = cEvents{2}(useEventsIdx);
    out.useEvents = useEvents;
    
%     disp(length(useEvents));
%     disp(sum(useEvents > spontOn & useEvents < spontOff))
    
    if opts.showPlots && ~isempty(useEvents)
        figure; 
        subplot(2,1,1); 
        imagesc(T, F_desired, cData);
%         imagesc(T, F_desired, smoothCol(cData, 2, 20, 'exp'));
        axis xy; % Ensure the y-axis is in the correct direction
        xlabel('Time (s)');
        ylabel('Frequency (Hz)');
        title('Spectrogram of lfpData from 0.5 Hz to 10 Hz');
%         colorbar;
        caxis([0 10]);
        xlim([T(1), T(end)]);
        
        subplot(2,1,2); 
        plot(T, eventTrace);
        nvline(useEvents, 'k');
        nvline(cEvents{3}(useEventsIdx), 'g');
        xlim([T(1), T(end)]);
        drawnow;
    end
    
    
    %%
    if ~isempty(useEvents)
        winSize = [opts.winSize]*fs;
        winSize = round(winSize(1) : winSize(2));
        org_events = round(useEvents .* fs);
        useIdx = org_events + winSize;
        useIdx = useIdx(:);
        
        testWins = lfpData(sort(useIdx));
        testWins = reshape(testWins, length(winSize), []);
        
        nrTrials = size(testWins,2);
    else
        nrTrials = 0;
        testWins = [];
    end
    
    S_all = cell(nrTrials, 1);    % For storing spectrogram data
    for iTrials = 1 : nrTrials
        
        trialData = testWins(:, iTrials);
        
        % Compute the spectrogram
        [S, F, ~] = spectrogram(trialData, windowLength, noverlap, nfft, fs);
        
        % Extract frequencies between freqMin and freqMax
        freqIndices = (F >= freqMin) & (F <= freqMax);
        S_desired1 = S(freqIndices, :);
        F_desired1 = F(freqIndices);
        
        % Store the spectrogram
        S_all{iTrials} = 10*log10(abs(S_desired1));
    end
    
    %     meanSpec = cat(3, S_all{:});
    % meanSpec = meanSpec - meanSpec(:,1,:);
    % meanSpec = meanSpec - median(cData,2);
    
    
    %%
    traceTime = opts.winSize(1) : 1/fs : opts.winSize(2);
    out.traceTime = traceTime;
    out.testWins = testWins;
    out.specTime = T;
    out.specFreq = F_desired;
    out.specPower = S_all;
    
    if opts.showPlots && checkTraces
        figure
    end
    
    freqIdx = F_desired > opts.eventFreqRange(1) & F_desired < opts.eventFreqRange(2);
    useFreq = F_desired(freqIdx);
    
    for x = 1 : length(useEvents)
        
        cData = S_all{x}(freqIdx,:) - medianPower(freqIdx);
        
        [a, timeIdx] = max(cData,[],2);
        [b, maxFreqIdx] = max(a,[],1);
        out.maxPower(x) = b;
        out.maxPowerFreq(x) = useFreq(maxFreqIdx);
        
        [out.peakMeanPower(x), meanFreqIdx] = max(nanmean(cData, 2));
        out.peakMeanFreq(x) = useFreq(meanFreqIdx);
        
        downsample_factor = 10;
        signal = double(testWins(:, x));
        signal = signal - prctile(signal, 75);
        Fs_down = fs / downsample_factor;
        signal_down = downsample(signal, downsample_factor);
        
        %         % Redefine Nyquist frequency
        %         Fn_down = Fs_down / 2;
        %
        %         % Redefine cutoff frequencies
        %         Wn_down = opts.eventFreqRange / Fn_down;
        %
        %         % Design filter with new parameters
        %
        %         [b_down, a_down] = butter(10, Wn_down, 'bandpass');
        
        % Filter the downsampled signal
        %         filtered_signal_down = filtfilt(b_down, a_down, signal_down');
        %         ds_traceTime = 0:1/Fs_down:5; % 2 seconds of data
        
        % Redefine Nyquist frequency
        Fn_down = Fs_down / 2;
        %
        % Redefine cutoff frequencies
        Wn_down = [0.5 20] / Fn_down;
        
        % Design filter with new parameters
        N_bp = 150;
        b_bp = fir1(N_bp, Wn_down, 'bandpass', hamming(N_bp+1));
        
        % Filter the downsampled signal
        filtered_signal_down = filtfilt(b_bp, 1, signal_down);
        ds_traceTime = opts.winSize(1) : 1/Fs_down : opts.winSize(2); % 2 seconds of data
        
%         filtered_signal_down = signal_down;
        
        if opts.loadSC
            peakProminence = 0.35;
        else
            peakProminence = 0.7;
        end
        [troughValues, troughLocs] = findpeaks(-filtered_signal_down, ds_traceTime, 'MinPeakProminence', peakProminence); % Invert the signal
        out.peakVals{x} = -troughValues; % Correct the sign for the trough heights
        
        if isempty(troughValues) || length(troughValues) <= opts.winSize(2) || out.peakMeanPower(x) < 3
            out.maxPower(x) = NaN;
            out.peakMeanPower(x) = NaN;
            out.eventDurations(x) = NaN;
            out.peakMeanFreq(x) = NaN;
            out.peakVals{x} = NaN;
        end
        
        if opts.showPlots && checkTraces
            % imagesc(T, F_desired, nanmean(meanSpec,3));
            subplot(1,2,1); cla;
%             plot(ds_traceTime, filtered_signal_down .* lfMeta.uV_per_bit_lfp ./ 1000);
            plot(ds_traceTime, filtered_signal_down);
            hold on;
%             plot(troughLocs, -troughValues .* lfMeta.uV_per_bit_lfp ./ 1000, 'ro', 'linewidth', 2);
            plot(troughLocs, -troughValues, 'ro', 'linewidth', 2);
            ylabel('signal voltage (mV)');
            axis square;
            xlabel('Time (s)');
            title('Filtered LFP with troughs');
            
            subplot(1,2,2);
            imagesc(T, useFreq, cData);
            nvline(T(timeIdx(maxFreqIdx)), 'r');
            nhline(useFreq(maxFreqIdx), 'r');
            axis xy; % Ensure the y-axis is in the correct direction
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            title(sprintf('Max power: %.2f - Mean peak power: %.2f', out.maxPower(x), out.peakMeanPower(x)));
            colorbar;
            axis square
            caxis([0 10])
            pause;
        end
    end
    
    %% save data locally
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end
    out.opts = opts;
    save(saveFile, 'out');
    
    % also save data in source lfp folder
    saveFile = fullfile(lfpDir, 'DNMT1_oscillations.mat');
    save(saveFile, 'out');
    
end
end