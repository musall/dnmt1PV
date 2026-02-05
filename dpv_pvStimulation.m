%% basic variables
clear opts
opts.Location = {'S1'};
opts.makePlot = false;
opts.removeAvg = false;
opts.useTaper = true;
opts.groups = {'PV-Cre', 'PV-DNMT1'};
opts.brainRange = 3000;
opts.localPath = 'F:\DNMT1_project\';
opts.reload = false;
opts.gid = 'PassiveStimulation';
opts.stimType = 1;

opts.optoThresh = 0.05; %thresold for detection of excitatory responses to optogenetics
opts.excludeNoise = false;
opts.redoAnalysis = false;
opts.loadRaw = false;
opts.gaussLength = 10;
opts.relevantVarNames = ["PulseCount", "PulseDur",  'PulseGap', 'OptoRamp', 'OptoShift', ... 
        'StimType', 'RedPowerOne', 'RedPowerTwo', 'BluePowerOne', 'BluePowerTwo'];
opts.optoLaser = 'BluePowerTwo';
opts.stimBin = 0.001; %1ms bins
opts.stimWindow = [-1 1.5]; %psth range for stimulus sequence
opts.pulseWindow = [-0.005 0.005]; %psth range for single pulses
opts.rampRange = 0.5; %ramp range to compute AUCs
opts.visRange = [0.015 0.035]; %sensory range to compute AUCs

ctxRange = [0 1500];
scRange = [1800 2500];
ctrlColor = [212 212 212]./255;
KOcolor = [255 160 64]./255;
set(0,'defaultfigurecolor',[1 1 1])

%% run over groups
failedRecs = cell(1, length(opts.groups));
allInfo = cell(1, length(opts.groups));
allBaselineFR = cell(1, length(opts.groups));
allRecs = cell(1, length(opts.groups));
allLabels = cell(1, length(opts.groups));
allDepth = cell(1, length(opts.groups));
allOptoPulse = cell(1, length(opts.groups));
allOptoRamp = cell(1, length(opts.groups));
allOptoPulseDelay = cell(1, length(opts.groups));
allRampPSTH = cell(1, length(opts.groups));
allRampVisPSTH = cell(1, length(opts.groups));
allVisPSTH = cell(1, length(opts.groups));
allFirstDur = cell(1, length(opts.groups));
allSecondDur = cell(1, length(opts.groups));
allWFs = cell(1, length(opts.groups));
allVis = cell(1, length(opts.groups));
allVisP = cell(1, length(opts.groups));
allClustIDs = cell(1, length(opts.groups));
for iGroups = 1 : length(opts.groups)
    
    disp(['Loading data from group ' num2str(iGroups) ': ' opts.groups{iGroups}]);
    cPath = fullfile(opts.localPath, opts.groups{iGroups});
    recs = dir([cPath filesep '**' filesep '*npxResults.mat']);
    for iRecs = 1 : length(recs)
        
        clear cluster pulse ramp sensory waveMetrics
        saveFile = fullfile(recs(iRecs).folder, recs(iRecs).name);
        load(saveFile, 'cluster', 'pulse', 'ramp', 'sensory', 'waveMetrics');
                
        % collect data for analysis
        allRecs{iGroups} = [allRecs{iGroups}; ones(size(sensory.baselineFR)).*iRecs];
        allBaselineFR{iGroups} = [allBaselineFR{iGroups}; mean(cat(2, sensory.baselineFR, sensory.optoBaselineFR),2)];
        allLabels{iGroups} = [allLabels{iGroups}, cluster.clusterLabels];
        allClustIDs{iGroups} = [allClustIDs{iGroups}; cluster.clusterIDs];
        allDepth{iGroups} = [allDepth{iGroups}; cluster.clustDepth];
        allOptoPulse{iGroups} = [allOptoPulse{iGroups}; pulse.AUCs'];
        allOptoRamp{iGroups} = [allOptoRamp{iGroups}; ramp.AUCs'];
        allVis{iGroups} = [allVis{iGroups}; sensory.AUCs'];
        allVisP{iGroups} = [allVisP{iGroups}; sensory.AUC_pVals'];
        allOptoPulseDelay{iGroups} = [allOptoPulseDelay{iGroups}; pulse.delays];
        allRampPSTH{iGroups} = [allRampPSTH{iGroups}; ramp.PSTH];
        allVisPSTH{iGroups} = [allVisPSTH{iGroups}; sensory.PSTH];
        allRampVisPSTH{iGroups} = [allRampVisPSTH{iGroups}; sensory.PSTH];
        allFirstDur{iGroups} = [allFirstDur{iGroups}; waveMetrics.firstDur'];
        allSecondDur{iGroups} = [allSecondDur{iGroups}; waveMetrics.secondDur'];
        allWFs{iGroups} = [allWFs{iGroups}, waveMetrics.meanWaves];
        
    end
end


%% compare firing rates
frData = cell(1, length(opts.groups));
for iGroups = 1 : length(opts.groups)
    cData = allBaselineFR{iGroups}(strcmpi(allLabels{iGroups}, 'sua')');
    frData{iGroups} = cData(cData<prctile(cData, 99));
end

% make figure
figure('renderer', 'painters');
subplot(2,2,1);
h(1) = Violin(frData{1}, 0, 'ViolinColor', ctrlColor, 'EdgeColor', ctrlColor, 'Bandwidth', 3, 'BoxColor', [0 0 0]);
h(2) = Violin(frData{2}, 1, 'ViolinColor', KOcolor, 'EdgeColor', KOcolor, 'Bandwidth', 3, 'BoxColor', [0 0 0]);
nhline(0, 'k--')
ylim([-5.5 50.5])
axis square
ax = h(1).ViolinPlot.Parent;
ax.XTick = 0:1;
ax.XTickLabel = {'Control', 'DNMT1-KO'};
ylabel('Firing rate (Hz)');
title('Spont. firing rate') 
niceFigure

%% compare firing rate for PV cells
optoData = cell(1, length(opts.groups));
delayData = cell(1, length(opts.groups));
for iGroups = 1 : length(opts.groups)
%     cIdx = allOptoPulse{iGroups} > 0.6 &  allOptoPulseDelay{iGroups} > 0.001; %responds to pulse with at least 1ms delay
    pvIdx = allOptoRamp{iGroups} > 0.55 & (strcmpi(allLabels{iGroups}, 'sua')');
    cData = allBaselineFR{iGroups}(pvIdx);
    optoData{iGroups} = cData(cData<prctile(cData, 95));
end

% make figure
subplot(2,2,2);
h(1) = Violin(rmoutliers(optoData{1}), 0, 'ViolinColor', ctrlColor, 'EdgeColor', ctrlColor, 'Bandwidth', 2, 'BoxColor', [0 0 0]);
h(2) = Violin(rmoutliers(optoData{2}), 1, 'ViolinColor', KOcolor, 'EdgeColor', KOcolor, 'Bandwidth', 2, 'BoxColor', [0 0 0]);
nhline(0, 'k--')
ylim([-5.5 50.5])
axis square
ax = h(1).ViolinPlot.Parent;
ax.XTick = 0:1;
ax.XTickLabel = {'Control', 'DNMT1-KO'};
ylabel('Firing rate (Hz)');
title('Spont. firing rate - PV cells')
niceFigure

%% check optogenetic suppression
rampPSTH = cell(1, length(opts.groups));
for iGroups = 1 : length(opts.groups)
    pyrIdx = allOptoRamp{iGroups} < 0.45;
    cData = allRampPSTH{iGroups}(pyrIdx,:);
    cData = cData - mean(cData(:, 1:500),2); %baseline correction
    rampPSTH{iGroups} = cData;
end

% show optogenetic impact
subplot(2,2,3);
xlim([-0.5 2]);
nhline(0,'k--');
lines(1) = stdshade(rampPSTH{1}, 0.3, ctrlColor, ramp.rampBins+0.5, 20);
lines(2) = stdshade(rampPSTH{2}, 0.3, KOcolor, ramp.rampBins+0.5, 20);
axis square
legend(lines, 'Ctrl', 'DNMT-KO', 'Location', 'northwest');
title('Optogenetically induced suppression');
xlabel('time from laser onset (s)')
niceFigure;

% quantify impact as the average change in firing while ramp is on
rampIdx = ramp.rampBins>-0.5 & ramp.rampBins<(-0.5+opts.rampRange); %ramp ontime
% rampIdx = ramp.rampBins>-0.5 & ramp.rampBins<-0.3; %ramp ontime
for iGroups = 1 : length(opts.groups)
    cData = rmoutliers(nanmean(rampPSTH{iGroups}(:, rampIdx), 2), 'percentiles', [3 95]);
%     cData = (nanmean(rampPSTH{iGroups}(:, rampIdx), 2));
    optoImpact{iGroups} = cData;
end
% optoImpact{2} = nanmean(rampPSTH{2}(:, cIdx), 2);

subplot(2,2,4);
h(1) = Violin(optoImpact{1}, 0, 'ViolinColor', ctrlColor, 'EdgeColor', ctrlColor, 'Bandwidth', 1, 'BoxColor', [0 0 0]);
h(2) = Violin(optoImpact{2}, 1, 'ViolinColor', KOcolor, 'EdgeColor', KOcolor, 'Bandwidth', 1, 'BoxColor', [0 0 0]);
nhline(0, 'k--')
% ylim([-1 10])
ylim([-15 2])
axis square
ax = h(1).ViolinPlot.Parent;
ax.XTick = 0:1;
ax.XTickLabel = {'Control', 'DNMT1-KO'};
ylabel('FR change vs baseline (Hz)');
title('PV-mediated suppresion')
niceFigure