function [cluster, pulse, ramp, sensory, waveMetrics] = pC_passiveFRcomparisonSOM_SM(opts)
% function to extract optogenetic and sensory responses, as well as firing
% rate and some other features from Neuropixels recordings with PV-Cre and
% PV-DNMT1 mice.

savePath = [opts.localPath filesep opts.recName filesep];
saveFile = fullfile(savePath, [opts.recName '_sensoryStim.mat']);
if ~opts.redoAnalysis
    try
        load(saveFile, 'cluster', 'pulse', 'ramp', 'sensory', 'waveMetrics');
        
    catch
        opts.redoAnalysis = true;
        disp('Unable to load local save file. Loading raw data.');
    end
end

if opts.redoAnalysis
    %% set paths for where to find the data
    lfpDir = [fullfile(opts.recPath, opts.recName) filesep opts.recName '_g0' filesep  opts.recName '_g0_imec' opts.probNr];
    lfMetaName = dir(fullfile(lfpDir,'*.lf.meta'));
    lfMeta = pC_readSpikeGLXmeta(fullfile(lfMetaName.folder,lfMetaName.name));
    
    nidaqDir = fileparts(lfpDir);
    niFile = dir(fullfile(nidaqDir,'*.nidq.bin'));
    niMetaName = dir(fullfile(nidaqDir,'*.nidq.meta'));
    niMeta = readSpikeGLXmeta(fullfile(niMetaName.folder,niMetaName.name));
    
    % load data
    ksDir = fullfile(lfpDir, '\spikeinterface_KS2_5_output\sorter_output');
    [sp, trigDat, params] = pC_loadKSdir(ksDir, opts);
    stimStart = trigDat{1}{2}{2}; %start of stimuli in current session
    
    % load oscillation file and check if stimuli fall in oscillatory event
    oscFile = fullfile(lfpDir, 'DNMT1_oscillations.mat');
    if exist(oscFile, 'file')
        load(oscFile, 'out');
    else
        out.eventDurations = [];
    end
    useIdx = find(~isnan(out.eventDurations));
    if ~isempty(useIdx)
        for x = 1 : length(useIdx)
            cEvent = out.useEvents(useIdx(x)); %onset of current
            cDur = out.eventDurations(useIdx(x)); %duration of current
            rejIdx = stimStart > cEvent & stimStart < (cEvent + cDur); %check for stimuli that fall into current event
            stimStart(rejIdx) = nan;
        end
    end
    if any(isnan(stimStart))
        fprintf('Removed %i/%i from total trials due to oscillatory events\n', sum(isnan(stimStart)), length(stimStart));
    end
    cluster.oscRejCnt(1) = sum(isnan(stimStart)); %total amount of removed trials
    cluster.oscRejCnt(2) = 0; %total amount of removed trials in relevant conditions
    
    
    %get channelindex for brain surface
    badChansFile = fullfile(lfpDir, 'bad_channel_ids_list_fixed.npy');
    if ~exist(badChansFile, 'file')
        badChansFile = fullfile(lfpDir, 'bad_channel_ids_list.npy');
    end
    
    % find brain surface based on last consecutive block of bad channels
    badChans = readNPY(badChansFile);
    badChans(badChans == 191) = []; %ignore the reference channel
    badChans = flipud(badChans+1);
    cIdx = min([find(diff(badChans) < -1, 1), length(badChans)]);
    cBrainIdx = true(1, lfMeta.nSites);
    if ~isempty(badChans) && badChans(1) >= 384
        cBrainIdx(badChans(cIdx):end) = false;
    end
    cBrainIdx = fliplr(cBrainIdx);
    depthRange = (1 : length(cBrainIdx)/2) * 20;
    depthRange = sort(repmat(depthRange,1,2));
    ctxStart = depthRange(find(cBrainIdx, 1)); %estimated cortex surface
    
    %%  find depth and clusters in cortex
    sp.trueSpikeDepths = 3840 - sp.spikeDepths - ctxStart; %shold give correct depths (from top to bottom)
    clusterIDs = unique(sp.clu);
    cluster.clusterIDs = clusterIDs;
    cluster.clustDepth = nan(length(cluster.clusterIDs),1);
    Cnt = 0;
    for iClust = 1 : length(clusterIDs)
        Cnt = Cnt + 1;
        cDepth = round(nanmedian(sp.trueSpikeDepths(sp.clu == clusterIDs(iClust))));
        if cDepth > -1000 %keep cluster that are close to surface
            cluster.clustDepth(Cnt) = cDepth;
        end
    end
    
    %% get waveform metrics
    clustIdx = ismember(sp.cids, clusterIDs); %only used clusters   
    waveMetrics = sp.metrics(clustIdx,:);
        
%     opts.waveLoc =   fullfile(lfpDir, 'waveforms_kilosort2.5\waveforms');
%     opts.recomputeWM = false; % dont load raw waveforms by defailt. They are usually not located on the server because of space constraints
%     waveMetrics = pC_computeWaveMetrics_SM(sp, opts);
    
    %% get Opto tagging (5Hz)
    [allStim, allStimData] = pC_getSessionData(lfpDir, opts.relevantVarNames);
    useTrials = false(1, size(allStim,2));
    Cnt = 0;
    for iFiles = 1 : length(allStimData)
        [trialCases, caseCnt, caseIdx, diffVarNames] = checkStimData(allStimData{iFiles});
        optoRamp = trialCases(strcmpi(allStimData{iFiles}.VarNames, 'OptoPyramid'), :);
        OptoOne = trialCases(strcmpi(allStimData{iFiles}.VarNames, 'BlueOne'),:) > 0;
        OptoTwo = trialCases(strcmpi(allStimData{iFiles}.VarNames, 'BlueTwo'),:) > 0;
        optoOnly = trialCases(strcmpi(allStimData{iFiles}.VarNames, 'StimType'),:) == 8;
        pulseCount = trialCases(strcmpi(allStimData{iFiles}.VarNames, 'PulseCount'),:) == 5;
        
        % session without ramps
        cIdx = (OptoOne  | OptoTwo) & ~optoRamp & optoOnly & pulseCount; % PulseCount
        cTrials = sum(caseIdx(:, cIdx),2); %use these trials
        useTrials(Cnt + find(cTrials)) = true; %use these trials for analysis
        Cnt = Cnt + sum(caseCnt);
    end
    
    % get trials with 5Hz and maximal laser ampltiude
    cTarget = opts.optoLaser; %target laser
    cData = allStim(strcmpi(opts.relevantVarNames, cTarget), :); %power levels for use trials
    cIdx = useTrials & cData == max(allStim(strcmpi(opts.relevantVarNames, cTarget), useTrials)); %max power for blue laser 1
    
    % get exact times for each pulse in the 5Hz sequence
    cStims = stimStart(cIdx);
    if any(isnan(cStims))
        fprintf('Removed %i/%i trials due to oscillatory events\n', sum(isnan(cStims)), length(cStims));
        cluster.oscRejCnt(2) = cluster.oscRejCnt(2) + sum(isnan(cStims));
        cStims(isnan(cStims)) = [];
    end
    pulse.nrTrials = length(cStims);
    psthLength= int16((diff(opts.stimWindow)/opts.stimBin));
    pulse.seqPSTH = zeros(length(clusterIDs),psthLength);
    for iClusters = 1:length(clusterIDs)
        spT = sp.st(sp.clu == clusterIDs(iClusters));
        [pulse.seqPSTH(iClusters,:), pulse.seqBins, ~, ~, ~, ~] = psthAndBA(spT, cStims, opts.stimWindow, opts.stimBin);
    end
    disp([cTarget ', 5Hz stimulus PSTH - done'])
    
    % get exact times for each pulse in the 5Hz sequence and compute individual pulse PSTHs
    analogOffset = 100; %add some datapoints as baseline when checking analog signal
    analogDat = pC_extractAnalogChannel(fullfile(niFile.folder,niFile.name), niMeta.nChans, 5, cStims, -analogOffset: niMeta.sRateHz);
    pulseStarts = [];
    for x = 1 : length(cStims)
        db = diff([0; analogDat(:,x)]);
        db = zscore(single(db)) > 5;
        cPulses = find(diff(db) > 0) - analogOffset; %start of current pulses
        cPulses(find(diff(cPulses) < 100) + 1) = []; %remove pulses if time differences from previous event is too short
        pulseStarts = [pulseStarts; cStims(x) + (cPulses / niMeta.sRateHz)];
    end
    
    % get single pulse PSTH and compute AUCs
    psthLength= int16((diff(opts.pulseWindow)/opts.stimBin));
    pulse.PSTH = zeros(length(clusterIDs),psthLength);
    spikeAvg = nan(length(pulseStarts)*2, length(clusterIDs));
    for iClusters = 1:length(clusterIDs)
        spT = sp.st(sp.clu == clusterIDs(iClusters));
        [pulse.PSTH(iClusters,:), pulse.pulseBins, ~, ~, ~, binnedArray] = psthAndBA(spT, pulseStarts, opts.pulseWindow, opts.stimBin);
        spikeAvg(1:length(pulseStarts), iClusters) = mean(binnedArray(:, pulse.pulseBins < 0), 2);
        spikeAvg(length(pulseStarts)+1:end, iClusters) = mean(binnedArray(:, pulse.pulseBins > 0), 2);
    end
    disp([cTarget ', 5Hz pulse PSTH - done'])
    
    % get AUCs
    y = [ones(length(pulseStarts), 1); zeros(length(pulseStarts),1)];
    [pulse.AUCs, ~, pulseP] = colAUC(spikeAvg, y, 'abs',false);
    
    % get median spike times after each pulse
    cSpikeT = []; cSpikeC = [];
    for iPulses = 1 : length(pulseStarts)
        cIdx = (sp.st > pulseStarts(iPulses)+opts.pulseWindow(1)) & (sp.st < pulseStarts(iPulses)+opts.pulseWindow(2));
        cSpikeT = [cSpikeT; sp.st(cIdx) - pulseStarts(iPulses)];
        cSpikeC = [cSpikeC; sp.clu(cIdx)];
    end
    
    % compute median delay for response latency
    pulse.delays = zeros(length(clusterIDs), 1);
    for iClusters = 1:length(clusterIDs)
        pulse.delays(iClusters) = median(cSpikeT(cSpikeC == clusterIDs(iClusters)));
    end
    
    %% get Opto ramping
    useTrials = false(1, size(allStim,2));
    Cnt = 0;
    for iFiles = 1 : length(allStimData)
        [trialCases, caseCnt, caseIdx, diffVarNames] = checkStimData(allStimData{iFiles});
        optoRamp = trialCases(strcmpi(allStimData{iFiles}.VarNames, 'OptoPyramid'), :);
        OptoOne = trialCases(strcmpi(allStimData{iFiles}.VarNames, 'BlueOne'),:) > 0;
        OptoTwo = trialCases(strcmpi(allStimData{iFiles}.VarNames, 'BlueTwo'),:) > 0;
        optoOnly = trialCases(strcmpi(allStimData{iFiles}.VarNames, 'StimType'),:) == 8;
        pulseCount = trialCases(strcmpi(allStimData{iFiles}.VarNames, 'PulseCount'),:) == 2;
        
        % session without ramps
        cIdx = (OptoOne  | OptoTwo) & optoRamp & optoOnly & pulseCount; % PulseCount
        cTrials = sum(caseIdx(:, cIdx),2); %use these trials
        useTrials(Cnt + find(cTrials)) = true; %use these trials for analysis
        Cnt = Cnt + sum(caseCnt);
    end
    
    % get trials with maximal laser amplitude
    cData = allStim(strcmpi(opts.relevantVarNames, cTarget), :); %power levels for use trials
    cIdx = useTrials & cData == max(allStim(strcmpi(opts.relevantVarNames, cTarget), useTrials)); %max power for blue laser 1
    optoShift = unique(allStim(strcmpi(opts.relevantVarNames, 'OptoShift'), cIdx)); %start point of the opto ramp
    
    % get ramp PSTHs
    cStims = stimStart(cIdx);
    if any(isnan(cStims))
        fprintf('Removed %i/%i trials due to oscillatory events\n', sum(isnan(cStims)), length(cStims));
        cluster.oscRejCnt(2) = cluster.oscRejCnt(2) + sum(isnan(cStims));
        cStims(isnan(cStims)) = [];
    end
    ramp.nrTrials = length(cStims);
    psthLength = int16((diff(opts.stimWindow)/opts.stimBin));
    ramp.PSTH = zeros(length(clusterIDs),psthLength);
    spikeAvg = nan(length(cStims)*2, length(clusterIDs));
    for iClusters = 1:length(clusterIDs)
        spT = sp.st(sp.clu == clusterIDs(iClusters));
        [ramp.PSTH(iClusters,:), ramp.rampBins, ~, ~, ~, binnedArray] = psthAndBA(spT, cStims, opts.stimWindow, opts.stimBin);
        spikeAvg(1:length(cStims), iClusters) = mean(binnedArray(:, ramp.rampBins < optoShift & ramp.rampBins > optoShift - opts.rampRange), 2);
        spikeAvg(length(cStims)+1:end, iClusters) = mean(binnedArray(:, ramp.rampBins > optoShift & ramp.rampBins < optoShift + opts.rampRange), 2); %ramp starts 0.5 seconds before 0, so AUC is looking at time before time 0 to get ramping part
    end
    disp([cTarget ', Ramp stimulus PSTH - done']);
    
    % get AUCs
    y = [ones(length(cStims), 1); zeros(length(cStims),1)];
    [ramp.AUCs, ~, rampP] = colAUC(spikeAvg, y, 'abs',false);

    %% get sensory responses and firing rate
    useTrials = false(1, size(allStim,2));
    Cnt = 0;
    for iFiles = 1 : length(allStimData)
        [trialCases, caseCnt, caseIdx, diffVarNames] = checkStimData(allStimData{iFiles});
        optoRamp = trialCases(strcmpi(allStimData{iFiles}.VarNames, 'OptoPyramid'), :);
        visualOnly = trialCases(strcmpi(allStimData{iFiles}.VarNames, 'StimType'),:) == opts.stimType;
        pulseCount = trialCases(strcmpi(allStimData{iFiles}.VarNames, 'PulseCount'),:) == 2;
        
        % session without ramps
        cIdx = optoRamp & visualOnly & pulseCount; % PulseCount
        cTrials = sum(caseIdx(:, cIdx),2); %use these trials
        useTrials(Cnt + find(cTrials)) = true; %use these trials for analysis
        Cnt = Cnt + sum(caseCnt);
    end
    
    % get trials with opto ramp and visual stimulation
    cTarget = opts.optoLaser; %target laser
    cData = allStim(strcmpi(opts.relevantVarNames, cTarget), :); %power levels for use trials
    cIdx = useTrials & cData == max(allStim(strcmpi(opts.relevantVarNames, cTarget), useTrials)); %max power for blue laser 1
    optoShift = unique(allStim(strcmpi(opts.relevantVarNames, 'OptoShift'), cIdx)); %start point of the opto ramp
        
    % get sensory ramp PSTHs
    cStims = stimStart(cIdx);
    if any(isnan(cStims))
        fprintf('Removed %i/%i trials due to oscillatory events\n', sum(isnan(cStims)), length(cStims));
        cluster.oscRejCnt(2) = cluster.oscRejCnt(2) + sum(isnan(cStims));
        cStims(isnan(cStims)) = [];
    end
    sensory.ramp_nrTrials = length(cStims);
    psthLength = int16((diff(opts.stimWindow)/opts.stimBin));
    sensory.rampPSTH = zeros(length(clusterIDs),psthLength);
    spikeAvg = nan(length(cStims)*2, length(clusterIDs));
    for iClusters = 1:length(clusterIDs)
        spT = sp.st(sp.clu == clusterIDs(iClusters));
        [sensory.rampPSTH(iClusters,:), rampBins, ~, ~, ~, binnedArray] = psthAndBA(spT, cStims, opts.stimWindow, opts.stimBin);
        spikeAvg(1:length(cStims), iClusters) = mean(binnedArray(:, rampBins < optoShift & rampBins > optoShift - opts.rampRange), 2); %this is the time before ramp started
        spikeAvg(length(cStims)+1:end, iClusters) = mean(binnedArray(:, rampBins > 0 & rampBins < opts.rampRange), 2); %this looks at the time when sensory stimulus was presented (100 ms)
    end
    disp([cTarget ', Visual + Ramp stimulus PSTH - done']);
    sensory.rampBins = rampBins;
    
    % get AUCs
    y = [ones(length(cStims), 1); zeros(length(cStims),1)];
    [sensory.rampAUCs, ~, rampP] = colAUC(spikeAvg, y, 'abs',false);
    
    % get firing rate from baseline before sensory stimulation
    sensory.optoBaselineFR =  nanmean(sensory.rampPSTH(:, rampBins < optoShift), 2);
    sensory.optoRampFR =  nanmean(sensory.rampPSTH(:, rampBins > optoShift & rampBins < 0), 2);
    
    % get trials visual stimulation
    cData = allStim(strcmpi(opts.relevantVarNames, 'BluePowerOne'), :); %power levels for use trials
    cIdx = useTrials & cData == 0; %no laser trials
        
    % get visual stimulus PSTHs
    cStims = stimStart(cIdx);
    if any(isnan(cStims))
        fprintf('Removed %i/%i trials due to oscillatory events\n', sum(isnan(cStims)), length(cStims));
        cluster.oscRejCnt(2) = cluster.oscRejCnt(2) + sum(isnan(cStims));
        cStims(isnan(cStims)) = [];
    end
    sensory.nrTrials = length(cStims);
    psthLength = int16((diff(opts.stimWindow)/opts.stimBin));
    sensory.PSTH = zeros(length(clusterIDs),psthLength);
    spikeAvg = nan(length(cStims)*2, length(clusterIDs));
    for iClusters = 1:length(clusterIDs)
        spT = sp.st(sp.clu == clusterIDs(iClusters));
        [sensory.PSTH(iClusters,:), stimBins, ~, ~, ~, binnedArray] = psthAndBA(spT, cStims, opts.stimWindow, opts.stimBin);
        spikeAvg(1:length(cStims), iClusters) = mean(binnedArray(:, 1 : round(opts.visRange(2) / diff(stimBins(1:2)))), 2); %this is the time before ramp started
%         spikeAvg(1:length(cStims), iClusters) = mean(binnedArray(:, stimBins < optoShift & stimBins > optoShift - opts.visRange(2)), 2); %this is the time before ramp started
        spikeAvg(length(cStims)+1:end, iClusters) = mean(binnedArray(:, stimBins >  0 & stimBins < opts.visRange(2)), 2); %this looks at the time when sensory stimulus was presented (100 ms)
    end
    sensory.stimBins = stimBins;
    disp([cTarget ',Visual PSTH - done']);

    % get AUCs
    y = [ones(length(cStims), 1); zeros(length(cStims),1)];
    [sensory.AUCs, ~, sensory.AUC_pVals] = colAUC(spikeAvg, y, 'abs',false);
    
    % get firing rate from baseline before sensory stimulation
    sensory.baselineFR = nanmean(sensory.PSTH(:, stimBins < 0), 2);
    
    % return template waveforms from dominant channel
    cluster.tempWFs = zeros(size(sp.temps,2), size(clusterIDs,1));
    for iClusters = 1:length(clusterIDs)
        cIdx = sp.cids == clusterIDs(iClusters);
        if sum(cIdx) == 0
            cluster.clusterLabels{iClusters} = 'noise';
        else
            cData = squeeze(sp.temps(cIdx,:,:));
            [~,ind] = min(min(cData,[],1));
            cluster.tempWFs(:, iClusters) = cData(:,ind);
            
            %add decoder labels
            if iscell(sp.clusterLabels)
                cluster.clusterLabels{iClusters} = sp.clusterLabels{cIdx};
            elseif sp.clusterLabels(cIdx) == 0
                cluster.clusterLabels{iClusters} = 'noise';
            elseif sp.clusterLabels(cIdx) == 1
                cluster.clusterLabels{iClusters} = 'mua';
            elseif sp.clusterLabels(cIdx) == 2
                cluster.clusterLabels{iClusters} = 'sua';
            elseif sp.clusterLabels(cIdx) == 3
                cluster.clusterLabels{iClusters} = 'unassigned';
            end
        end
    end
    
    %% save local file
    save(saveFile, 'cluster', 'pulse', 'ramp', 'sensory', 'waveMetrics');
    
end
end