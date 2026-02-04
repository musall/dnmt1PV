function [pulse, cluster, ramp, waveMetrics] = pC_optoResponseSOM_SM(opts)
% function to extract optogenetic and sensory responses, as well as firing
% rate and some other features from Neuropixels recordings with PV-Cre and
% PV-DNMT1 mice.

savePath = [opts.localPath filesep opts.recName filesep opts.recName '_g0' filesep  opts.recName '_g0_imec' opts.probeNr];
savePath = strrep(savePath, '\\', '\');
saveFile = fullfile(savePath, [opts.recName '_optoRampResults.mat']);
if ~opts.redoAnalysis
    try
        load(saveFile, 'pulse', 'cluster', 'ramp', 'waveMetrics');
        
    catch
        opts.redoAnalysis = true;
        disp('Unable to load local save file. Loading raw data.');
    end
end

if opts.redoAnalysis
    %% set paths for where to find the data
    lfpDir = [fullfile(opts.recPath, opts.recName) filesep opts.recName '_g0' filesep  opts.recName '_g0_imec' opts.probeNr];
    lfMetaName = dir(fullfile(lfpDir,'*.lf.meta'));
    lfMeta = pC_readSpikeGLXmeta(fullfile(lfMetaName.folder,lfMetaName.name));
    
    nidaqDir = fileparts(lfpDir);
    niFile = dir(fullfile(nidaqDir,'*.nidq.bin'));
    niMetaName = dir(fullfile(nidaqDir,'*.nidq.meta'));
    niMeta = readSpikeGLXmeta(fullfile(niMetaName.folder,niMetaName.name));
    
    lfFile = dir(fullfile(lfpDir,'*.lf.bin'));
    lfFile = fullfile(lfpDir, lfFile.name);
    
    % load data
    ksDir = fullfile(lfpDir, '\spikeinterface_KS2_5_output\sorter_output');
    [sp, trigDat, params] = pC_loadKSdir(ksDir, opts);
    stimStart = trigDat{1}{2}{2}; %start of stimuli in current session
    
    %get channelindex for brain surface
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
        
    else
        % if badchannels are not available use LFP method
        meanStd = medfilt1(pC_findBrainSurface_lfp(lfFile, lfMeta));
        meanStd = flipud(meanStd(chanSteps));
        meanStd = meanStd ./ std(meanStd);
        meanStd = meanStd - meanStd(1);
        cBrainIdx = meanStd' > 1.5; %contacts in the brain
    end
    
    %     % find brain surface based on last consecutive block of bad channels
    %     badChans = readNPY(badChansFile);
    %     badChans(badChans == 191) = []; %ignore the reference channel
    %     badChans = flipud(badChans+1);
    %     cIdx = min([find(diff(badChans) < -1, 1), length(badChans)]);
    %     cBrainIdx = true(1, lfMeta.nSites);
    %     if ~isempty(badChans) && badChans(1) >= 384
    %         cBrainIdx(badChans(cIdx):end) = false;
    %     end
    %     cBrainIdx = fliplr(cBrainIdx);
    
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
    
    %% get waveform metrics from bombcell code
    %     clustIdx = ismember(sp.cids, clusterIDs); %only used clusters
    %     waveMetrics = sp.metrics(clustIdx,:);

    waveMetrics = pC_computeWaveMetricsFromTemplate_SM(sp, opts);
    cIdx = ismember(sp.cids, clusterIDs)';
    waveMetrics = selectBehaviorTrials(waveMetrics, cIdx, length(cIdx));
    
    [bcMetrics, bcLabel, forGUI] = pC_computeBombcellMetrics_fromTemplate_SM(opts);
    bcMetrics.bcLabel = bcLabel;
    
    for v = bcMetrics.Properties.VariableNames
        waveMetrics.(v{1}) = bcMetrics.(v{1})';
    end
    
    %% get Opto ramping
    laserPower = {'BluePowerOne', 'BluePowerTwo'};
    [allStim, allStimData] = pC_getSessionData(lfpDir, opts.relevantVarNames);
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
    
    ramp.PSTH = cell(1, length(laserPower));
    ramp.AUCs = cell(1, length(laserPower));
    for iLaser = 1 : length(laserPower)
        
        % get trials with maximal laser amplitude
        cData = allStim(strcmpi(opts.relevantVarNames, laserPower{iLaser}), :); %power levels for use trials
        cIdx = useTrials & cData == max(allStim(strcmpi(opts.relevantVarNames, laserPower{iLaser}), useTrials)); %max power for blue laser 1
        optoShift = unique(allStim(strcmpi(opts.relevantVarNames, 'OptoShift'), cIdx)); %start point of the opto ramp
        
        % get ramp PSTHs
        cStims = stimStart(cIdx);
        psthLength = int16((diff(opts.stimWindow)/opts.stimBin));
        ramp.PSTH{iLaser} = zeros(length(clusterIDs),psthLength);
        spikeAvg = nan(length(cStims)*2, length(clusterIDs));
        for iClusters = 1:length(clusterIDs)
            spT = sp.st(sp.clu == clusterIDs(iClusters));
            [ramp.PSTH{iLaser}(iClusters,:), ramp.rampBins, ~, ~, ~, binnedArray] = psthAndBA(spT, cStims, opts.stimWindow, opts.stimBin);
            spikeAvg(1:length(cStims), iClusters) = mean(binnedArray(:, ramp.rampBins < optoShift & ramp.rampBins > optoShift - opts.rampRange), 2);
            spikeAvg(length(cStims)+1:end, iClusters) = mean(binnedArray(:, ramp.rampBins > optoShift & ramp.rampBins < optoShift + opts.rampRange), 2); %ramp starts 0.5 seconds before 0, so AUC is looking at time before time 0 to get ramping part
        end
        disp([laserPower{iLaser} ', Ramp stimulus PSTH - done']);
        
        % get AUCs
        y = [ones(length(cStims), 1); zeros(length(cStims),1)];
        [ramp.AUCs{iLaser}, ~, rampP] = colAUC(spikeAvg, y, 'abs',false);
    end
    
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
    psthLength= int16((diff(opts.stimWindow)/opts.stimBin));
    pulse.seqPSTH = zeros(length(clusterIDs),psthLength);
    for iClusters = 1:length(clusterIDs)
        spT = sp.st(sp.clu == clusterIDs(iClusters));
        [pulse.seqPSTH(iClusters,:), pulse.seqBins, ~, ~, ~, ~] = psthAndBA(spT, cStims, opts.stimWindow, opts.stimBin);
    end
    disp([cTarget ', 5Hz stimulus PSTH - done'])
    
    % get exact times for each pulse in the 5Hz sequence and compute individual pulse PSTHs
    if strcmpi(opts.optoLaser, 'BluePowerOne')
        trigLine = 5;
    elseif strcmpi(opts.optoLaser, 'BluePowerTwo')
        trigLine = 6;
    end
    
    pulseStarts = [];
    pulseTrigStart = trigDat{1}{trigLine}{2}; %start of laser pulses
    for x = 1 : length(cStims)
        pulseStarts = [pulseStarts; pulseTrigStart(find(pulseTrigStart > cStims(x)-0.1, 1))];
    end
    
%     analogDat = pC_extractAnalogChannel(fullfile(niFile.folder,niFile.name), niMeta.nChans, 7, cStims, -analogOffset: niMeta.sRateHz);
%     pulseStarts = [];
%     for x = 1 : length(cStims)
%         db = diff([0; analogDat(:,x)]);
%         db = zscore(single(db)) > 5;
%         cPulses = find(diff(db) > 0) - analogOffset; %start of current pulses
%         cPulses(find(diff(cPulses) < 100) + 1) = []; %remove pulses if time differences from previous event is too short
%         pulseStarts = [pulseStarts; (cStims(x) + (cPulses / niMeta.sRateHz))];
%     end
    
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
        cIdx = sp.st > pulseStarts(iPulses) & sp.st < pulseStarts(iPulses)+opts.pulseWindow(2);
        cSpikeT = [cSpikeT; sp.st(cIdx) - pulseStarts(iPulses)];
        cSpikeC = [cSpikeC; sp.clu(cIdx)];
    end
    
    % compute median delay for response latency
    pulse.delays = zeros(length(clusterIDs), 1);
    for iClusters = 1:length(clusterIDs)
        pulse.delays(iClusters) = median(cSpikeT(cSpikeC == clusterIDs(iClusters)));
    end
    
    %% return template waveforms from dominant channel
%     cluster.tempWFs = forGUI.tempWv'; %use from bombcell
    cluster.tempWFs = zeros(size(sp.temps,2), size(clusterIDs,1));
    for iClusters = 1:length(clusterIDs)
        cIdx = sp.cids == clusterIDs(iClusters);
        if sum(cIdx) == 0
            cluster.clusterLabels{iClusters} = 'noise';
        else
            cData = squeeze(sp.temps(cIdx,:,:));
%             [~,ind] = min(min(cData,[],1));
            [~,ind] = max(max(cData,[],1)-min(cData,[],1)); %find max channel
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
    save(saveFile, 'pulse', 'cluster', 'ramp', 'waveMetrics');
    
end
end