function [allSpecs, plotFreqs, oscRejCnt] = pC_checkVisualGamma_SM(recName, recPath, opts)

%% basic variables
stepSize = 100; %stepsize for channel analysis in um (should be a divider of 20).
GammaParadigmName = 'vStim_VisualPuffer'; %expected name of paradigm for gamma stimuli
GammaParadigmTrials = 640; %expected number of trials for gamma stimuli
plotSpectraChan = 5; %channel after which the power spectra is shown

%% set paths for where to find the data
lfpDir = [recPath filesep recName filesep recName '_g0' filesep  recName '_g0_imec' opts.imecNr];
nidaqDir = fileparts(lfpDir);
savePath = [opts.savePath filesep recName filesep];

if opts.useTaper
    saveFile = fullfile(savePath, ['lfpData_taper_imec' opts.imecNr '.mat']);
else
    saveFile = fullfile(savePath, ['lfpData_pwelch_imec' opts.imecNr '.mat']);
end

%% load data or compute new spectra
fprintf('Current path: %s. ', recName);
oscRejCnt = zeros(1,2);
if ~opts.reload && exist(saveFile, 'file')
    %% get data from local save file
    fprintf('Loading local data ... ');
    load(saveFile, 'allSpecs', 'plotFreqs', 'oscRejCnt');
    
else
    fprintf('Loading raw data ... ');
    
    %% get meta files and digital triggers
    lfMetaName = dir(fullfile(lfpDir,'*.lf.meta'));
    lfMeta = pC_readSpikeGLXmeta(fullfile(lfMetaName.folder,lfMetaName.name));
    
    niMetaName = dir(fullfile(nidaqDir,'*.nidq.meta'));
    niMeta = pC_readSpikeGLXmeta(fullfile(niMetaName.folder,niMetaName.name));
    trigDat = pC_extractDigitalChannel(nidaqDir, niMeta.nChans, niMeta.nChans, 'nidq'); %get digital channel from ni-daq. usually the lat channel
    
    %% load settings from stimulation software
    relevantVarNames = ["PulseCount", "PulseDur",  'PulseGap', 'OptoRamp', 'StimType', 'RedPowerOne', 'RedPowerTwo', 'BluePowerOne', 'BluePowerTwo', 'BlueOne', 'BlueTwo', 'RedOne', 'RedTwo', 'StimAngle'];
    [allStim, allStimData] = pC_getSessionData(lfpDir, relevantVarNames); %get settings for each trial from stimulation software
    
    % compare with trial triggers
    trialStart = trigDat{1}{2}{2};
    if length(trialStart) ~= size(allStim,2)
        trialStart(end) = []; %there is one recording with 1 trial too many
    end
    
    % load oscillation file and check if stimuli fall in oscillatory event
    oscFile = fullfile(lfpDir, 'DNMT1_oscillations.mat');
    load(oscFile, 'out');
    useIdx = find(~isnan(out.eventDurations));
    if ~isempty(useIdx)
        for x = 1 : length(useIdx)
            cEvent = out.useEvents(useIdx(x)); %onset of current
            cDur = out.eventDurations(useIdx(x)); %duration of current
            rejIdx = trialStart > cEvent & trialStart < (cEvent + cDur); %check for stimuli that fall into current event
            trialStart(rejIdx) = nan;
        end
    end
    oscRejCnt(1) = sum(isnan(trialStart)); %total amount of removed trials

    %% isolate StimData from last session. This is where full field visual stimuli were presented
    %find correct Gamma paradigm
    gammaRecIdx = [];
    for iRecs = 1 : length(allStimData)
        if strcmpi(allStimData{iRecs}.Paradigm, GammaParadigmName) && ...
                allStimData{iRecs}.nrTrials == GammaParadigmTrials
            
            gammaRecIdx = [gammaRecIdx iRecs];
        end
    end
    
    %check there is only one paradigm that matches criteria
    if length(gammaRecIdx) ~= 1
        error('Found no or too many gamma parardigms!!')
    end
    
    % do other things
    [trialCases, caseCnt, caseIdx] = checkStimData(allStimData{gammaRecIdx}, relevantVarNames);
    stimDur = min(trialCases(strcmpi(allStimData{gammaRecIdx}.VarNames, 'StimDuration'), :)); %duration of visual stimulus
    
    % remove optogenetic trials at this point
    optoIdx = strcmpi(relevantVarNames, 'BlueOne') | strcmpi(relevantVarNames, 'BlueTwo');
    useIdx = sum(trialCases(optoIdx, :),1) == 0; %remove laser cases
    trialCases = trialCases(:, useIdx);
    caseCnt = caseCnt(useIdx);
    nrTrials = sum(caseCnt); %all trials
    
    % find trials in larger structure
    gammaTrialIdx = false(1, size(allStim, 2));
    for i = 1 : size(trialCases, 2)
        nanIdx = isnan(trialCases(:,i));
        cIdx = sum(allStim(~nanIdx,:) == trialCases(~nanIdx,i)) == sum(~nanIdx);
        gammaTrialIdx(cIdx) = true; %find gamma trials in larger structure for all trials
    end
    
    if sum(gammaTrialIdx) ~= nrTrials
        error('Number of trials in larger structure doesnt match number of trials in the selected paradigm');
    end
    cTrials = trialStart(gammaTrialIdx); %gamma trials in the ephys recording
    if any(isnan(cTrials))
        fprintf('Removed %i/%i trials due to oscillatory events\n', sum(isnan(cTrials)), length(cTrials));
        oscRejCnt(2) = oscRejCnt(2) + sum(isnan(cTrials));
        cTrials(isnan(cTrials)) = [];
    end
    nrTrials = length(cTrials);
    
    %% get analog data
    lfFile = dir(fullfile(lfpDir,'*.lf.bin'));
    lfFile = fullfile(lfpDir, lfFile.name);
    winSize = round(1 * lfMeta.sRateHz); %size of data to be used for spectral analyis in datapoints
    
    % check if median average should be used
    if opts.removeAvg
        preStimAvg = single(pC_extractMedianAnalogChannel(lfFile, lfMeta.nChans, cTrials, -winSize:-1)); %get current channel
        postStimAvg = single(pC_extractMedianAnalogChannel(lfFile, lfMeta.nChans, cTrials, 1:winSize)); %get current channel
    end
    
    if opts.useTaper
        [~, freqs] = pmtm(rand(winSize,1), [], winSize,  lfMeta.sRateHz, 'DropLastTaper', false);
    else
        [~,freqs] = pwelch(rand(winSize,1), winSize/2, [], winSize, lfMeta.sRateHz, 'psd');
    end
    freqRange = freqs < 200; %define which frequencies to keep
    
    
    %% get mean standard deviations of different parts of the recording to determine the brain surface.
    badChansFile = fullfile(lfpDir, 'bad_channel_ids_list_fixed.npy');
    if ~exist(badChansFile, 'file')
        badChansFile = fullfile(lfpDir, 'bad_channel_ids_list.npy');
    end
    
    if exist(badChansFile, 'file')
        badChans = readNPY(badChansFile);
        badChans(badChans == 191) = []; %ignore the reference channel
        
        % find brain surface based on last consecutive block of bad channels
        badChans = flipud(badChans+1);
        cIdx = min([find(diff(badChans) < -1, 1), length(badChans)]);
        cBrainIdx = true(1, lfMeta.nSites);
        if ~isempty(badChans) && badChans(1) >= 384
            cBrainIdx(badChans(cIdx):end) = false;
        end
        cBrainIdx = fliplr(cBrainIdx);
        meanStd = cBrainIdx';
        
    else
        % if badchannels are not available use LFP method
        meanStd = medfilt1(pC_findBrainSurface_lfp(lfFile, lfMeta));
        meanStd = flipud(meanStd);
        meanStd = meanStd ./ std(meanStd);
        meanStd = meanStd - meanStd(1);
        cBrainIdx = smooth(meanStd)' > opts.brainThresh; %contacts in the brain
    end
    
    % get depth of individual contacts
    depthRange = (1 : length(cBrainIdx)/2) * 20;
    depthRange = sort(repmat(depthRange,1,2));
    brainStart = find(cBrainIdx, 1); %cortical surface
    depthRange = depthRange - depthRange(brainStart);
    
    depthIdx = depthRange >= 0 & depthRange < opts.brainRange;
    if sum(cBrainIdx) < sum(depthIdx)
        error('Not enough contacts in brain to cover requested depth');
    end
    cBrainIdx = cBrainIdx & depthIdx; %contacts at desired depth range
    depthRange = depthRange(cBrainIdx);
    
    %% load analog data and compute spectra    
    stepChanIdx = find(rem(depthRange, stepSize) < 20); %channel index for stepping through depth
    if sum(diff(stepChanIdx) == 1) < length(stepChanIdx)-1
        rejIdx = [false, diff(stepChanIdx) == 1]; %remove every second channel here since columns will get averaged.
        stepChanIdx(rejIdx) = [];
    else %when using all channels, remove every second contact here for column averaging
        if depthRange(2) == 0
            stepChanIdx = stepChanIdx(1:2:length(stepChanIdx));
        else
            stepChanIdx = stepChanIdx([1 2:2:length(stepChanIdx)]);
        end
    end
    
    useChans = length(cBrainIdx) - find(cBrainIdx); %channels within brain range. this should match 'depthRange'
    allSpecs = NaN(sum(freqRange), ceil(length(stepChanIdx)/2), 2);
    plotChecker = true;
    for iChans = 1 : length(stepChanIdx)
        
        % current channel
        cChan = useChans(stepChanIdx(iChans));
        
        clear preStim* postStim*
        preStim = single(pC_extractAnalogChannel(lfFile, lfMeta.nChans, cChan, cTrials, -winSize:-1)); %get current channel
        postStim = single(pC_extractAnalogChannel(lfFile, lfMeta.nChans, cChan, cTrials, 1:winSize)); %get current channel
        
        % average contacts from both columns together. only exception if first
        % contact has no second channel to be averaged with.
        if ~(iChans == 1 && depthRange(2) ~= 0)
            temp = single(pC_extractAnalogChannel(lfFile, lfMeta.nChans, cChan+1, cTrials, -winSize:-1)); %get current channel
            preStim = nanmean(cat(3, preStim, temp), 3);
            
            temp = single(pC_extractAnalogChannel(lfFile, lfMeta.nChans, cChan+1, cTrials, 1:winSize)); %get current channel
            postStim = nanmean(cat(3, postStim, temp), 3);
        end
        
        if opts.removeAvg
            preStim = preStim - preStimAvg;
            postStim = postStim - postStimAvg;
        end
        
        beforeSpec = zeros(sum(freqRange), nrTrials);
        afterSpec = zeros(sum(freqRange), nrTrials);
        for iTrials = 1 : nrTrials
            if opts.useTaper
                [temp, freqs] = pmtm(preStim(:,iTrials), [], winSize,  lfMeta.sRateHz, 'DropLastTaper', false, 'unity');
                beforeSpec(:, iTrials) = temp(freqRange);
                
                [temp, freqs] = pmtm(postStim(:,iTrials), [], winSize,  lfMeta.sRateHz, 'DropLastTaper', false, 'unity');
                afterSpec(:, iTrials) = temp(freqRange);
                
            else
                temp = pwelch(preStim(:,iTrials), winSize/2, [], winSize, lfMeta.sRateHz);
                beforeSpec(:, iTrials) = temp(freqRange);
                
                temp = pwelch(postStim(:,iTrials), winSize/2, [], winSize, lfMeta.sRateHz);
                afterSpec(:, iTrials) = temp(freqRange);
            end
            clear temp
        end
        allSpecs(:, iChans, 1) = nanmean(beforeSpec, 2);
        allSpecs(:, iChans, 2) = nanmean(afterSpec, 2);
        
        % show spectra comparison before and after
        plotFreqs = freqs(freqRange);
        if plotChecker && iChans > plotSpectraChan && opts.makePlot
            figure('renderer' ,'painters');
            hold on;
            cLine(1) = stdshade(beforeSpec', 0.5, 'g', plotFreqs);
            cLine(2) = stdshade(afterSpec', 0.5, 'k', plotFreqs);
            cLine(1).Parent.YScale = "log";
            %         ylim([10E-3 10E2]);
            legend(cLine, {'before' 'after'});
            title([recName '; depth: ' num2str(depthRange(iChans)) '; last channel is 0 here']);
            xlabel('Frequency(Hz)');
            plotChecker = false;
            xlim([0 200]); axis square
        end
    end
    
    % save results
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end
    save(saveFile, 'allSpecs', 'plotFreqs', 'oscRejCnt');
end

%%
if opts.makePlot
    h = figure('renderer' ,'painters');
    for x = 1 : 2
        subplot(1,2,x);
        
        if x == 1
            cMap = (log10(allSpecs(:,:,2)))';
        elseif x == 2
            cMap = ((allSpecs(:,:,2) - allSpecs(:,:,1)) ./ allSpecs(:,:,1))';
        end
        cImg = imagesc(cMap);
        ax = cImg.Parent;
        colormap(ax, colormap_blueblackred(256));
        maxFreq = plotFreqs(end);
        cRange = abs(prctile(cMap(:),97.5));
        caxis([-cRange cRange]);
        axis image;
        
        depthRange = (1:size(cMap,1)) ./ size(cMap,1) * opts.brainRange;
        cIdx = 1 : 5 : size(cMap,1);
        ax.YTick = cIdx;
        ax.YTickLabels = depthRange(cIdx);
        nhline(cIdx, 'w')
        ylabel('Depth (um)');
        
        useFreqs = 1 : ceil(20 / mean(diff(plotFreqs))) : length(plotFreqs);
        ax.XTick = useFreqs;
        ax.XTickLabels = round(plotFreqs(useFreqs),2);
        nvline(useFreqs, 'w')
        xlabel('Frequency (Hz)');
        title([recName ' - removeAvg = ' num2str(opts.removeAvg) ' - useTaper = ' num2str(opts.useTaper)]);
        axis square
        niceFigure
    end
end
fprintf('done\n');
