function [csdResp, lfpResp, meanStd, normDepth, lfMeta, oscRejCnt] = ds_loadCSDresponse(recInfo, recLabels, opts)

%% basic variables
stimParadigmName = 'vStim_MultimodalStim'; %expected name of paradigm
stimParadigmTrials = 440; %expected number of trials for tactile stimuli

if isfield(opts, 'stimParadigmTrials')
    stimParadigmTrials = opts.stimParadigmTrials;
end

%% set paths for where to find the data
recName = recInfo{strcmpi(recLabels, 'Folder')};
recPath = recInfo{strcmpi(recLabels, 'Path')};
probNr = recInfo{strcmpi(recLabels, 'Probe')};

lfpDir = [recPath filesep recName filesep recName '_g0' filesep  recName '_g0_imec' probNr];
lfpDir = strrep(lfpDir, '.kampa-10g', ''); % dont use 10g network
nidaqDir = fileparts(lfpDir);

%% load data or compute new responses
savePath = [opts.savePath filesep recName filesep];
saveFile = fullfile(savePath, ['csdResponse_imec' probNr '.mat']);
fprintf('Current path: %s. ', recName);

if ~opts.reload && exist(saveFile, 'file')
    %% get data from local save file
    disp('Loading local data ... ');
    load(saveFile, 'csdResp', 'lfpResp', 'meanStd', 'normDepth', 'lfMeta');
    
else
    
    %% get meta files and digital triggers
    lfMetaName = dir(fullfile(lfpDir,'*.lf.meta'));
    lfMeta = pC_readSpikeGLXmeta(fullfile(lfMetaName.folder,lfMetaName.name));
    
    niMetaName = dir(fullfile(nidaqDir,'*.nidq.meta'));
    niMeta = pC_readSpikeGLXmeta(fullfile(niMetaName.folder,niMetaName.name));
    trigDat = pC_extractDigitalChannel(nidaqDir, niMeta.nChans, niMeta.nChans, 'nidq'); %get digital channel from ni-daq. usually the lat channel
    
    %% load settings from stimulation software
    relevantVarNames = ["PuffDur", "AudioAmp", "PulseCount", "PulseDur",  'PulseGap', 'OptoRamp', 'StimType', 'RedPowerOne', 'RedPowerTwo', 'BluePowerOne', 'BluePowerTwo', 'BlueOne', 'BlueTwo', 'RedOne', 'RedTwo', 'StimAngle'];
    [allStim, allStimData] = pC_getSessionData(lfpDir, relevantVarNames); %get settings for each trial from stimulation software
    
    % compare with trial triggers
    trialStart = trigDat{1}{2}{2};
    trialDiff = size(allStim,2) - length(trialStart);
    if trialDiff ~= 0
        error('StimData and digital Data have not the same number of trials');
        warning(['!!!Trial difference is ' num2str(trialDiff) '!!!']);
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
    if any(isnan(trialStart))
        fprintf('Removed %i/%i from total trials due to oscillatory events\n', sum(isnan(trialStart)), length(trialStart));
    end
    oscRejCnt(1) = sum(isnan(trialStart)); %total amount of removed trials
    oscRejCnt(2) = 0; %total amount of removed trials in relevant conditions
    
    chanSteps = 1 : opts.stepSize : lfMeta.nChans-1; %spacing between channels when making spatial figure
    chanDepth = lfMeta.nChans * 10 - 10 : -20 : 0; %spacing between contacts is 20um apart
    chanDepth = [chanDepth; chanDepth]; %two contacts next to each other at same depth
    chanDepth = chanDepth(:);
    
    %% find correct paradigm
    stimRecIdx = [];
    for iRecs = 1 : length(allStimData)
        if strcmpi(allStimData{iRecs}.Paradigm, stimParadigmName) && ...
                allStimData{iRecs}.nrTrials >= min(stimParadigmTrials) && allStimData{iRecs}.nrTrials <= max(stimParadigmTrials)
            
            stimRecIdx = [stimRecIdx iRecs];
        end
    end
    
    %check there is only one paradigm that matches criteria
    % if length(stimRecIdx) ~= 1
    %     error('Found no or too many gamma parardigms!!')
    % end
    stimRecIdx = stimRecIdx(1); %make sure there is only one paradigm to use
    
    % do other things
    [trialCases, caseCnt, caseIdx] = checkStimData(allStimData{stimRecIdx}, relevantVarNames);
    
    % check if optogenetic trials should be used or not
    optoIdx = strcmpi(relevantVarNames, 'BlueOne') | strcmpi(relevantVarNames, 'BlueTwo');
    stimTypeIdx = strcmpi(relevantVarNames, 'StimType'); %tactile only trials
    useIdx = sum(trialCases(optoIdx, :),1) == opts.optoGenetics & trialCases(stimTypeIdx, :) == opts.stimType; % use non opto sensory trials
    trialCases = trialCases(:, useIdx);
    caseCnt = caseCnt(useIdx);
    nrTrials = sum(caseCnt); %all trials
    
    % find trials in larger structure
    trialIdx = false(1, size(allStim, 2));
    for i = 1 : size(trialCases, 2)
        nanIdx = isnan(trialCases(:,i));
        cIdx = sum(allStim(~nanIdx,:) == trialCases(~nanIdx,i)) == sum(~nanIdx);
        trialIdx(cIdx) = true; %find trials in larger structure for all trials
    end
    
    if sum(trialIdx) ~= nrTrials
        warning('Number of trials in larger structure doesnt match number of trials in the selected paradigm');
        trialIdx = find(trialIdx);
        trialIdx = trialIdx(1:nrTrials);
    end
    cTrials = trialStart(trialIdx); %trials in the ephys recording
    cTrials(1:5) = nan;
    if any(isnan(cTrials))
        fprintf('Removed %i/%i trials due to oscillatory events\n', sum(isnan(cTrials)), length(cTrials));
        oscRejCnt(2) = oscRejCnt(2) + sum(isnan(cTrials));
        cTrials(isnan(cTrials)) = [];
    end
    
    %% get analog data
    lfFile = dir(fullfile(lfpDir,'*.lf.bin'));
    lfFile = fullfile(lfpDir, lfFile.name);
    baseSize = round(opts.baseDur * lfMeta.sRateHz);
    winSize = round(opts.postStim * lfMeta.sRateHz);
    
    % check if median average should be used
    if opts.removeAvg
        postStimAvg = single(pC_extractMedianAnalogChannel(lfFile, lfMeta.nChans, cTrials, -baseSize+1:winSize)); %get current channel
    end
    
    %% get mean standard deviations of different parts of the recording to determine the brain surface.
    % check for bad channels results - tends to be more accurate if saline was applied
    badChansFile = fullfile(lfpDir, 'bad_channel_ids_list_fixed.npy');
    if ~exist(badChansFile, 'file')
        badChansFile = fullfile(lfpDir, 'bad_channel_ids_list.npy');
    end
    
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
    
    depthRange = (1 : length(cBrainIdx)) *10*opts.stepSize;
    brainStart = find(cBrainIdx, 1); %cortical surface
    depthRange = depthRange - depthRange(brainStart);
    
    depthIdx = depthRange >= 0 & depthRange < opts.brainRange(2);
    if sum(cBrainIdx) < sum(depthIdx)
        error('Not enough contacts in brain to cover requested depth');
    end
    
    cBrainIdx = cBrainIdx & depthIdx; %contacts at desired depth range
    depthRange = depthRange(cBrainIdx);
    
    %% load analog data
    cData = zeros(sum(cBrainIdx), winSize+baseSize, 'single');
    Cnt = 0;
    for iChans = chanSteps(fliplr(cBrainIdx))
        Cnt = Cnt + 1;
        temp = pC_extractAnalogChannel(lfFile, lfMeta.nChans, iChans, cTrials, -baseSize+1:winSize, opts); %get current channel
        cData(Cnt, :) = nanmean(temp,2); %get trial average
    end
    lfpResp = flipud(cData .* lfMeta.uV_per_bit_lfp ./ 1000);
    
    %% compute CSD
    [csdResp, normDepth] = pC_SplineCSD(lfpResp, depthRange);
    csdResp = csdResp ./ 1000;
    
    % save data
    save(saveFile, 'csdResp', 'lfpResp', 'meanStd', 'normDepth', 'lfMeta');

end
