function [waveMetrics, opts] = pC_computeWaveMetricsFromTemplate_SM(sp, opts)

if ~isfield(opts, 'recWin') || isempty(opts.recWin)
    opts.recWin = round(sp.sample_rate / 1000 / 3); %compute recWin in samples to get 1/3 miliseconds
end

if ~isfield(opts, 'showSpike') || isempty(opts.showSpike)
    opts.showSpike = false; %flag to show spikes with duration
end

if opts.showSpike
    h1 = figure('renderer', 'painters');
end

%%
tic
nrClusts = size(sp.temps,1);
ptv = nan(1, nrClusts);
ptratio = nan(1, nrClusts);
hw = nan(1, nrClusts);
repolSlope = nan(1, nrClusts);
recSlope = nan(1, nrClusts);
spikeDur = nan(1, nrClusts);
firstDur = nan(1, nrClusts);
secondDur = nan(1, nrClusts);
WF = nan(size(sp.temps,2), nrClusts);

%% for bombcell code
cRec = [fullfile(opts.recPath, opts.recName) filesep opts.recName '_g0' filesep  opts.recName '_g0_imec' opts.probeNr];
ephysMetaDir = dir(fullfile(cRec, '*.ap.meta'));
ephysMetaDir = ephysMetaDir(1);
ephysKilosortPath = fullfile(cRec, 'spikeinterface_KS2_5_output', 'sorter_output');
savePath = [cRec 'qMetrics']; % where you want to save the quality metrics 

kilosortVersion = 2.5; % if using kilosort4, you need to have this value kilosertVersion=4. Otherwise it does not matter. 
gain_to_uV = NaN; % use this if you are not using spikeGLX or openEphys to record your data. this value, 
% when mulitplied by your raw data should convert it to  microvolts. 

[spikeTimes_samples, spikeClusters, templateWaveforms, templateAmplitudes, pcFeatures, ...
    pcFeatureIdx, channelPositions] = bc.load.loadEphysData(ephysKilosortPath, savePath);


param = bc.qm.qualityParamValues(ephysMetaDir, NaN, ephysKilosortPath);
meta = bc.dependencies.SGLX_readMeta.ReadMeta(ephysMetaDir.name, ephysMetaDir.folder);
[AP, ~, SY] = bc.dependencies.SGLX_readMeta.ChannelCountsIM(meta);
param.nChannels = AP + SY;
param.nSyncChannels = SY;

maxChannels = bc.qm.helpers.getWaveformMaxChannel(sp.temps);
waveformBaselineWindow = [param.waveformBaselineWindowStart, param.waveformBaselineWindowStop];
channelPositions = [sp.xcoords, sp.ycoords];

%%
for iClust = 1:nrClusts
    
    meanWave = squeeze(sp.temps(iClust,:,:));
    [~, targChan] = max(sum(abs(meanWave),1)); %find max channel
    
    cWave = squeeze(sp.temps(iClust,:,targChan))';
    WF(:, iClust) = cWave;
    minRange = round(size(cWave,1) / 5);
    spikeRange = minRange+1:size(cWave,1)-minRange; %define range for spike occurence
    
    [maxWave, maxWaveInd] = max(cWave(spikeRange));
    [minWave, minWaveInd] = min(cWave(spikeRange));
    
    maxWaveInd = maxWaveInd + minRange;
    minWaveInd = minWaveInd + minRange;
    
    %%
    ptv(iClust) = abs(maxWaveInd - minWaveInd) / sp.sample_rate * 10^6; %time difference in microseconds
    
    %% compute peak_to_valley ratio
    ptratio(iClust) = abs(maxWave / minWave);
    
    %% compute full-width half maximum for negative spike
    threshold = minWave * 0.5;
    
    beforeThresh = find(cWave(1:minWaveInd) < threshold); %find sample before threshold crossing
    if ~isempty(find(diff(beforeThresh) > 1,1))
        beforeThresh = beforeThresh(find(diff(beforeThresh) > 1,1)+1 : end);
    end
    if ~isempty(beforeThresh)
        beforeThresh = beforeThresh(1) - 1;
    end
    
    threshDiff = threshold - cWave(beforeThresh); %change to threshold
    sampleDiff = diff(cWave(beforeThresh:beforeThresh+1)); %change to next sample
    beforeThresh = beforeThresh + (threshDiff / sampleDiff); %exact time of threshold crossing before spike
    
    afterThresh = find(cWave(minWaveInd:end) < threshold); %find last sample before threshold crossing
    if ~isempty(find(diff(afterThresh) > 1,1))
        afterThresh = afterThresh(1 : find(diff(afterThresh) > 1,1));
    end
    if ~isempty(beforeThresh)
        afterThresh = afterThresh(end) + minWaveInd - 1; %find last sample before threshold crossing
    end
    if afterThresh == size(cWave,1); afterThresh = size(cWave,1) - 1; end %make sure this is not the last sample of the waveform
    
    threshDiff = threshold - cWave(afterThresh); %change to threshold
    sampleDiff = diff(cWave(afterThresh:afterThresh+1)); %change to next sample
    afterThresh = afterThresh + (threshDiff / sampleDiff); %exact time of threshold crossing after spike
    
    if ~isempty(afterThresh - beforeThresh)
        hw(iClust) = (afterThresh - beforeThresh) / sp.sample_rate * 1000; %time for negative spike
    else
        hw(iClust) = 0;
    end
    
    %%
    repolarizationTime = find(cWave(minWaveInd:end) > 0, 1, 'first') + minWaveInd - 1;
    a = fitlm(1:length(minWaveInd:repolarizationTime), cWave(minWaveInd:repolarizationTime));
    repolSlope(iClust) = a.Coefficients.Estimate(2);
    
    %% compute recovery slope
    %     Return the recovery slope of input waveforms. After repolarization,
    %     the neuron hyperpolarizes untill it peaks. The recovery slope is the
    %     slope of the actiopotential after the peak, returning to the baseline
    %     in dV/dT. The slope is computed within a user-defined window after
    %     the peak.
    
    [~, lateMaxWaveInd] = max(cWave(minWaveInd:end));
    lateMaxWaveInd = lateMaxWaveInd + minWaveInd - 1;
    
    ctime = lateMaxWaveInd : lateMaxWaveInd + opts.recWin - 1;
    if ctime(end) <= length(cWave)
        a = fitlm(1:length(ctime), cWave(ctime)); %get slope in defined time window after positive peak
        recSlope(iClust) = a.Coefficients.Estimate(2);
    end
    
    %%
%     [nPeaks, nTroughs, spatialDecaySlope, waveformBaseline, scndPeakToTroughRatio, mainPeakToTroughRatio, peak1ToPeak2Ratio,...
%         troughToPeak2Ratio, mainPeak_before_width, mainPeak_after_width, mainTrough_width, peakLocs, troughLocs, waveformDuration_peakTrough(thisUnit), ...
%         spatialDecayPoints, thisWaveform, spatialDecayPoints_loc, spatialDecayFit_1] = bc.qm.waveformShape(sp.temps, ...
%         thisUnit, maxChannels(thisUnit), param, channelPositions, waveformBaselineWindow);
    
    [nPeaks(iClust), nTroughs(iClust), spatialDecaySlope(iClust), waveformBaseline(iClust), scndPeakToTroughRatio(iClust), mainPeakToTroughRatio(iClust), peak1ToPeak2Ratio(iClust),...
        troughToPeak2Ratio(iClust), mainPeak_before_width(iClust), mainPeak_after_width(iClust), mainTrough_width(iClust), peakLocs{iClust}, troughLocs{iClust}, waveformDuration_peakTrough(iClust), ...
        spatialDecayPoints{iClust}, ~, spatialDecayPoints_loc{iClust}, spatialDecayFit_1(iClust)] = bc.qm.waveformShape(sp.temps, ...
        iClust, maxChannels(iClust), param, channelPositions, waveformBaselineWindow);
    
    %% compute spike duration
    threshold = std(cWave)/4;
    firstPeak = min([maxWaveInd, minWaveInd]);
    a = abs(smooth(cWave(1:firstPeak),10));
    a = a - prctile(a,1);
    spikeStart = find(a > threshold);
    if ~isempty(find(diff(spikeStart) > 1,1))
        spikeStart = spikeStart(find(diff(spikeStart) > 1,1)+1 : end);
    end
    if isempty(spikeStart)
        spikeStart = 1;
    else
        spikeStart = spikeStart(1);
    end
    
    lastPeak = max([maxWaveInd, minWaveInd]);
    a = abs(smooth(cWave(lastPeak:end), 10));
    a = a - prctile(a,1);
    spikeEnd = find(a < threshold);
    if ~isempty(find(diff(spikeEnd) > 1,1))
        spikeEnd = spikeEnd(1 : find(diff(spikeEnd) > 1,1));
    end
    if isempty(spikeEnd)
        spikeEnd = size(cWave,1);
    end
    spikeEnd = spikeEnd(1) - 1 + lastPeak; %find sample before threshold crossing
    spikeDur(iClust) = (spikeEnd - spikeStart) / sp.sample_rate * 10^6; %spike duration in us
    
    %% check for first and second spike duration
    %         % check for sign and invert if needed so first spike is negative
    %         if cWave(firstPeak) < 0
    %             currWave = smooth(cWave);
    %         else
    %             currWave = -smooth(cWave);
    %         end
    
    temp = cWave;
    temp(abs(temp) < std(temp) * 0.25) = 0;
    currWave = cWave;
    currWave(1:firstPeak) = temp(1:firstPeak);
    fNeg = find(currWave(1:minWaveInd)>=0,1, 'last')+1; %first negative value
    sNeg = (minWaveInd-2)+find(currWave(minWaveInd:end)>=0,1); %second negative value
    tNeg = sNeg+2+find(currWave(sNeg+3:end)<=0,1); %third negative value
    
    firstCross = (fNeg-1)+(currWave(fNeg-1)/abs(currWave(fNeg)-currWave(fNeg-1))); %first zero-line crossing
    secondCross = (sNeg)+(abs(currWave(sNeg))/abs(currWave(sNeg+1)-currWave(sNeg))); %second zero-line crossing
    thirdCross = (tNeg-1)+(currWave(tNeg-1)/abs(currWave(tNeg)-currWave(tNeg-1))); %third zero-line crossing
    
    if length([firstCross, secondCross, thirdCross]) == 3
        firstDur(iClust) = ((secondCross-firstCross)/sp.sample_rate)*10^6; %duration of the first spike in us
        secondDur(iClust) = ((thirdCross-secondCross)/sp.sample_rate)*10^6; %duration of the second spike in us
    end
    
    if rem(iClust, round(nrClusts/10)) == 0
        fprintf('Loaded %.0f/%.0f clusters (%.0f percent)\n',iClust, nrClusts, iClust/nrClusts*100)
    end
    
    %% show some results
    if opts.showSpike
        figure(h1); cla;
        title(sprintf('custer %d; firstDur: %gus\n', iClust, firstDur(iClust)));
        plot(cWave);
        %             nvline([spikeStart, spikeEnd]);
        nhline(0, 'k');
        nvline([firstCross, secondCross, thirdCross], 'k')
        pause;
    end
end

%% collect results in single structure and save to local file
waveMetrics.ptv = ptv;
waveMetrics.ptratio = ptratio;
waveMetrics.hw = hw;
waveMetrics.repolSlope = repolSlope;
waveMetrics.recSlope = recSlope;
waveMetrics.spikeDur = spikeDur;
waveMetrics.firstDur = firstDur;
waveMetrics.secondDur = secondDur;
waveMetrics.WF = WF;

waveMetrics.spatialDecaySlope = spatialDecaySlope;
waveMetrics.waveformBaseline = waveformBaseline;
waveMetrics.scndPeakToTroughRatio = scndPeakToTroughRatio;
waveMetrics.mainPeakToTroughRatio = mainPeakToTroughRatio;
waveMetrics.peak1ToPeak2Ratio = peak1ToPeak2Ratio;

waveMetrics.troughToPeak2Ratio = troughToPeak2Ratio;
waveMetrics.mainPeak_before_width = mainPeak_before_width;
waveMetrics.mainPeak_after_width = mainPeak_after_width;
waveMetrics.mainTrough_width = mainTrough_width;
waveMetrics.peakLocs = peakLocs;
waveMetrics.nPeaks = nPeaks;
waveMetrics.nTroughs = nTroughs;
waveMetrics.troughLocs = troughLocs;
waveMetrics.troughLocs = troughLocs;
waveMetrics.waveformDuration_peakTrough = waveformDuration_peakTrough;
waveMetrics.spatialDecayPoints = spatialDecayPoints;
waveMetrics.spatialDecayPoints_loc = spatialDecayPoints_loc;
waveMetrics.spatialDecayFit_1 = spatialDecayFit_1;
waveMetrics.maxChannels = maxChannels;
