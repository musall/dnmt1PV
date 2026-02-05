%% define experiments to use
recLabels = {'Animal','ExperimentDate','Location','Expertise','Folder','Probe','Background','Path','Ramppower','useRec','hasGamma'};
recInfo{1} = {'PV2555','26/01/2022','V1','Naive','PV2555_20221026','0','PV-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'PV2555','28/01/2022','V1','Naive','PV2555_20221028','0','PV-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'PV2555','29/01/2022','V1','Naive','PV2555_20221029','0','PV-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'PV2556','26/01/2022','V1','Naive','PV2556_20221026','0','PV-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'PV2556','27/01/2022','V1','Naive','PV2556_20221027','0','PV-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'PV2556','28/01/2022','V1','Naive','PV2556_20221028','0','PV-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'PV2556','29/01/2022','V1','Naive','PV2556_20221029','0','PV-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control'};
recInfo{2} = {'PV2632','09/02/2023','V1','Naive','PV2632_20230209','0','PV-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','PV-DNMT1';'PV2632','10/02/2023','V1','Naive','PV2632_20230210','0','PV-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','PV-DNMT1';'PV2632','11/02/2023','V1','Naive','PV2632_20230211','0','PV-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','PV-DNMT1';'PV2633','24/01/2023','V1','Naive','PV2633_20230124','0','PV-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','PV-DNMT1';'PV2633','25/01/2023','V1','Naive','PV2633_20230125','0','PV-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','PV-DNMT1';'PV2633','26/01/2023','V1','Naive','PV2633_20230126','0','PV-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','PV-DNMT1';'PV2633','27/01/2023','V1','Naive','PV2633_20230127','0','PV-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','PV-DNMT1'};

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

%% visual responses - V1
visPSTH = cell(1, length(opts.groups));
visFWHM = cell(1, length(opts.groups));
visPeak = cell(1, length(opts.groups));
visPeakAmp = cell(1, length(opts.groups));
visMean = cell(1, length(opts.groups));
visAUC = cell(1, length(opts.groups));
visResponseCells = nan(1, length(opts.groups));
groupCells = nan(1, length(opts.groups));

for iGroups = 1 : length(opts.groups)
    cIdx = ~strcmpi(allLabels{iGroups}, 'noise')' & allDepth{iGroups} < 1500; %all non-noise clusters in cortex
    visIdx = allVisP{iGroups} < 0.05 & allVis{iGroups} > 0.5 & cIdx; %all responsive neurons
    cData = allVisPSTH{iGroups}(visIdx,:);
    cData = cData - mean(cData(:, sensory.rampBins<0),2); %baseline correction
    visPSTH{iGroups} = cData; %get PSTH
    visAUC{iGroups} = allVis{iGroups}(visIdx); %get AUC
        
    % count cells
    visResponseCells(iGroups) = sum(visIdx);
    groupCells(iGroups) = sum(cIdx);
    
    % check amplitude and temporal precision of sensory response
    respTime = sensory.rampBins > 0 & sensory.rampBins < 0.2;
    baseTime = sensory.rampBins < 0 & sensory.rampBins > -0.2;
    cTime = sensory.rampBins(respTime);
    for iClust = 1 : size(cData,1)
        if visAUC{iGroups}(iClust) > 0.6 %needs strong enough response to work
            
            visMean{iGroups}(iClust) = mean(cData(iClust, respTime)) - mean(cData(iClust, baseTime));

            % width of response
            cCluster = cData(iClust, respTime);
            cCluster(1) = 0;
            cCluster(end) = 0;
            cCluster = smooth(cCluster, 10);
            visFWHM{iGroups}(iClust) = fwhm(cTime, cCluster)*1000; %full width half maximum in miliseconds

            %time of response
            [visPeakAmp{iGroups}(iClust), temp] = max(cCluster);
            visPeak{iGroups}(iClust) = cTime(temp) * 1000; %time of response peak in miliseconds 
        else
            visFWHM{iGroups}(iClust) = nan;
            visPeak{iGroups}(iClust) = nan;
            visMean{iGroups}(iClust) = nan;
            visPeakAmp{iGroups}(iClust) = nan;
        end
    end
end
[~, fracP] = binomCompare(visResponseCells(1), groupCells(1), visResponseCells(2), groupCells(2));
visResponseFraction = visResponseCells ./ groupCells;

% show visual response trace
figure
xlim([-0.2 1]);
ylim([-5 25])
nhline(0,'k--');
nvline([0 0.5], 'k--');
lines(1) = stdshade(visPSTH{1}, 0.3, ctrlColor, ramp.rampBins, 20);
lines(2) = stdshade(visPSTH{2}, 0.3, KOcolor, ramp.rampBins, 20);
axis square
legend(lines, 'Ctrl', 'DNMT-KO', 'Location', 'northeast');
title('Visual response (V1)');
xlabel('time from stimulus (s)');
ylabel('firing rate (Hz)');
niceFigure;


%% show peak amplitude, delay and temporal precision of sensory response
figure
subplot(1,3,1); clear h cData;
cData{1} = rmoutliers(visPeak{1});
cData{2} = rmoutliers(visPeak{2});
h(1) = Violin(cData{1}, 0, 'ViolinColor', ctrlColor, 'EdgeColor', ctrlColor, 'BoxColor', [0 0 0]);
h(2) = Violin(cData{2}, 1, 'ViolinColor', KOcolor, 'EdgeColor', KOcolor, 'BoxColor', [0 0 0]);
ylim([35 225])
axis square
ax = h(1).ViolinPlot.Parent;
ax.XTick = 0:1;
ax.XTickLabel = {'Control', 'DNMT1-KO'};
ylabel('Time of visual response peak (ms)');
title('Visual response peak after stimulus (V1)');
niceFigure
%

% show visual peak response
subplot(1,3,2); clear h cData;
cData{1} = rmoutliers(visPeakAmp{1});
cData{2} = rmoutliers(visPeakAmp{2});
h(1) = Violin(cData{1}, 0, 'ViolinColor', ctrlColor, 'EdgeColor', ctrlColor, 'BoxColor', [0 0 0]);
h(2) = Violin(cData{2}, 1, 'ViolinColor', KOcolor, 'EdgeColor', KOcolor, 'BoxColor', [0 0 0]);
% nhline(0.5, 'k--')
% ylim([0.45 1])
axis square
ax = h(1).ViolinPlot.Parent;
ax.XTick = 0:1;
ax.XTickLabel = {'Control', 'DNMT1-KO'};
ylabel('Visual peak response [Hz]');
title('Visual peak response (V1)');
niceFigure

subplot(1,3,3); clear h cData;
cData{1} = rmoutliers(visFWHM{1});
cData{2} = rmoutliers(visFWHM{2});
h(1) = Violin(cData{1}, 0, 'ViolinColor', ctrlColor, 'EdgeColor', ctrlColor, 'BoxColor', [0 0 0]);
h(2) = Violin(cData{2}, 1, 'ViolinColor', KOcolor, 'EdgeColor', KOcolor, 'BoxColor', [0 0 0]);
ylim([0 150])
axis square
ax = h(1).ViolinPlot.Parent;
ax.XTick = 0:1;
ax.XTickLabel = {'Control', 'DNMT1-KO'};
ylabel('Full-width-half-maximum (ms)');
title('Visual response width (V1)');
niceFigure
