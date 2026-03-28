%% basic variables
clear opts
opts.makePlot = false;
opts.removeAvg = false;
opts.useTaper = true;
opts.groups = {'PV-Cre', 'PV-DNMT1'};
opts.brainRange = 3000;
opts.localPath = 'F:\DNMT1_project\';
opts.reload = false;
opts.gid = 'PassiveStimulation';
opts.stimType = 1;

opts.optoThresh = 0.05; % threshold for detection of excitatory responses to optogenetics
opts.excludeNoise = false;
opts.redoAnalysis = false;
opts.loadRaw = false;
opts.gaussLength = 10;
opts.relevantVarNames = ["PulseCount", "PulseDur",  'PulseGap', 'OptoRamp', 'OptoShift', ...
        'StimType', 'RedPowerOne', 'RedPowerTwo', 'BluePowerOne', 'BluePowerTwo'];
opts.optoLaser = 'BluePowerTwo';
opts.stimBin = 0.001; %1ms bins
opts.stimWindow = [-1 1.5]; % psth range for stimulus sequence
opts.pulseWindow = [-0.005 0.005]; % psth range for single pulses
opts.rampRange = 0.5; % ramp range to compute AUCs
opts.visRange = [0.015 0.035]; % sensory range to compute AUCs

ctxRange = [0 1500]; 
scRange = [1800 2500]; 
ctrlColor = [212 212 212]./255;
KOcolor = [255 160 64]./255;
set(0,'defaultfigurecolor',[1 1 1])

%% run over groups
nGroups = numel(opts.groups);

allBaselineFR = cell(1, nGroups);
allRecs = cell(1, nGroups);
allLabels = cell(1, nGroups);
allDepth = cell(1, nGroups);
allOptoPulse = cell(1, nGroups);
allOptoRamp = cell(1, nGroups);
allOptoPulseDelay = cell(1, nGroups);
allRampPSTH = cell(1, nGroups);
allRampVisPSTH = cell(1, nGroups);
allVisPSTH = cell(1, nGroups);
allFirstDur = cell(1, nGroups);
allSecondDur = cell(1, nGroups);
allWFs = cell(1, nGroups);
allVis = cell(1, nGroups);
allVisP = cell(1, nGroups);
allClustIDs = cell(1, nGroups);

% new
allAnimalID = cell(1, nGroups);
allRecID = cell(1, nGroups);

for iGroups = 1:nGroups

    disp(['Loading data from group ' num2str(iGroups) ': ' opts.groups{iGroups}]);
    cPath = fullfile(opts.localPath, opts.groups{iGroups});
    recs = dir([cPath filesep '**' filesep '*npxResults.mat']);

    for iRecs = 1:numel(recs)

        clear cluster pulse ramp sensory waveMetrics
        saveFile = fullfile(recs(iRecs).folder, recs(iRecs).name);
        load(saveFile, 'cluster', 'pulse', 'ramp', 'sensory', 'waveMetrics');

        nCells = numel(sensory.baselineFR);

        % expected path: ...\GROUP\animal_X\rec_Y\file.mat
        recFolder = recs(iRecs).folder;
        [animalFolder, recName] = fileparts(recFolder);
        [~, animalName] = fileparts(animalFolder);

        animalID = [opts.groups{iGroups} '_' animalName];
        recID = [opts.groups{iGroups} '_' animalName '_' recName];

        % collect data for analysis
        allRecs{iGroups} = [allRecs{iGroups}; repmat(iRecs, nCells, 1)];
        allAnimalID{iGroups} = [allAnimalID{iGroups}; repmat({animalID}, nCells, 1)];
        allRecID{iGroups} = [allRecID{iGroups}; repmat({recID}, nCells, 1)];

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
frData = cell(1, nGroups);
frAnimalID = cell(1, nGroups);
frRecID = cell(1, nGroups);

for iGroups = 1:nGroups
    suaIdx = strcmpi(allLabels{iGroups}, 'sua')';
    cData = allBaselineFR{iGroups}(suaIdx);
    cAnimal = allAnimalID{iGroups}(suaIdx);
    cRec = allRecID{iGroups}(suaIdx);

    keepIdx = ~isoutlier(cData, 'percentiles', [0 99]);
    frData{iGroups} = cData(keepIdx);
    frAnimalID{iGroups} = cAnimal(keepIdx);
    frRecID{iGroups} = cRec(keepIdx);
end

%% compare firing rate for PV cells
optoData = cell(1, nGroups);
optoAnimalID = cell(1, nGroups);
optoRecID = cell(1, nGroups);

for iGroups = 1:nGroups
    pvIdx = allOptoRamp{iGroups} > 0.55 & strcmpi(allLabels{iGroups}, 'sua')';
    cData = allBaselineFR{iGroups}(pvIdx);
    cAnimal = allAnimalID{iGroups}(pvIdx);
    cRec = allRecID{iGroups}(pvIdx);

    keepIdx = ~isoutlier(cData, 'percentiles', [0 99]);
    optoData{iGroups} = cData(keepIdx);
    optoAnimalID{iGroups} = cAnimal(keepIdx);
    optoRecID{iGroups} = cRec(keepIdx);
end

%% check optogenetic suppression
rampPSTH = cell(1, nGroups);
pyrAnimalID = cell(1, nGroups);
pyrRecID = cell(1, nGroups);

for iGroups = 1:nGroups
    pyrIdx = allOptoRamp{iGroups} < 0.45;
    cData = allRampPSTH{iGroups}(pyrIdx,:);
    cData = cData - mean(cData(:, 1:500),2); % baseline correction
    rampPSTH{iGroups} = cData;

    pyrAnimalID{iGroups} = allAnimalID{iGroups}(pyrIdx);
    pyrRecID{iGroups} = allRecID{iGroups}(pyrIdx);
end

% quantify impact as average change in firing while ramp is on
rampIdx = ramp.rampBins > -0.5 & ramp.rampBins < (-0.5 + opts.rampRange);

optoImpact = cell(1, nGroups);
impactAnimalID = cell(1, nGroups);
impactRecID = cell(1, nGroups);

for iGroups = 1:nGroups
    cData = nanmean(rampPSTH{iGroups}(:, rampIdx), 2);
    cAnimal = pyrAnimalID{iGroups};
    cRec = pyrRecID{iGroups};

    keepIdx = ~isoutlier(cData, 'percentiles', [5 95]);
    optoImpact{iGroups} = cData(keepIdx);
    impactAnimalID{iGroups} = cAnimal(keepIdx);
    impactRecID{iGroups} = cRec(keepIdx);
end

%% make figure
figure('renderer', 'painters');

% 1) spontaneous firing rate
statsFR = plotViolinWithLME( ...
    2, 2, 1, ...
    frData, frAnimalID, ...
    ctrlColor, KOcolor, ...
    'Firing rate (Hz)', ...
    'Spont. firing rate', ...
    [-5.5 50.5], ...
    3, true);

% 2) spontaneous firing rate - PV cells
statsPV = plotViolinWithLME( ...
    2, 2, 2, ...
    optoData, optoAnimalID, ...
    ctrlColor, KOcolor, ...
    'Firing rate (Hz)', ...
    'Spont. firing rate - PV cells', ...
    [-5.5 50.5], ...
    2, true);

% 3) trace plot
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

% 4) PV-mediated suppression
statsSupp = plotViolinWithLME( ...
    2, 2, 4, ...
    optoImpact, impactAnimalID, ...
    ctrlColor, KOcolor, ...
    'FR change vs baseline (Hz)', ...
    'PV-mediated suppression', ...
    [-15 2], ...
    1, true);

%% LME summary
statsSummary = table( ...
    {'Spont. firing rate'; ...
     'Spont. firing rate - PV cells'; ...
     'PV-mediated suppression'}, ...
    [statsFR.p; statsPV.p; statsSupp.p], ...
    [statsFR.t; statsPV.t; statsSupp.t], ...
    'VariableNames', {'Measure', 'pValue_LME', 'tStat_LME'});

disp(' ');
disp('========== LME summary ==========');
disp(statsSummary);


%% plotting function for violins
function stats = plotViolinWithLME(nRows, nCols, plotIdx, dataByGroup, animalIDByGroup, ctrlColor, KOcolor, yLabelTxt, titleTxt, yLimits, bw, drawZeroLine, transformType)

subplot(nRows, nCols, plotIdx); hold on;

plotData = cell(1,2);
plotAnimalID = cell(1,2);
plotRecID = cell(1,2);

for iG = 1:2
    keepIdx = ~isnan(dataByGroup{iG});
    plotData{iG} = dataByGroup{iG}(keepIdx);
    plotAnimalID{iG} = animalIDByGroup{iG}(keepIdx);
end

h(1) = Violin(plotData{1}, 0, ...
    'ViolinColor', ctrlColor, ...
    'EdgeColor', ctrlColor, ...
    'Bandwidth', bw, ...
    'BoxColor', [0 0 0]);

h(2) = Violin(plotData{2}, 1, ...
    'ViolinColor', KOcolor, ...
    'EdgeColor', KOcolor, ...
    'Bandwidth', bw, ...
    'BoxColor', [0 0 0]);

if drawZeroLine
    nhline(0, 'k--')
end

if ~isempty(yLimits)
    ylim(yLimits)
end

if ~exist('transformType', 'var')
    transformType = 'none';
end

axis square
ax = h(1).ViolinPlot.Parent;
ax.XTick = 0:1;
ax.XTickLabel = {'Control', 'DNMT1-KO'};
ylabel(yLabelTxt);

% build LME inputs
dataIn = [plotData{1}; plotData{2}];
% dataIn = log(dataIn - min(dataIn) + 1);
conditionID = [ones(numel(plotData{1}),1); 2*ones(numel(plotData{2}),1)];
animalVar = [plotAnimalID{1}; plotAnimalID{2}];

[pVal, tStat, fullmodel] = LME_compareMulti_SM(dataIn, conditionID, animalVar, 'normal', transformType);

title(sprintf('%s\nLME p = %.4g', titleTxt, pVal));
niceFigure

stats = struct;
stats.p = pVal;
stats.t = tStat;
stats.fullmodel = fullmodel;
stats.data = plotData;
stats.animalID = plotAnimalID;
stats.recID = plotRecID;
end
