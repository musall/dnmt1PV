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
opts.stimBin = 0.001; % 1 ms bins
opts.stimWindow = [-1 1.5]; % PSTH range for stimulus sequence
opts.pulseWindow = [-0.005 0.005]; % PSTH range for single pulses
opts.rampRange = 0.5; % ramp range to compute AUCs
opts.visRange = [0.015 0.035]; % sensory range to compute AUCs

ctxRange = [0 1500];
scRange = [1800 2500]; 
ctrlColor = [212 212 212]./255;
KOcolor = [255 160 64]./255;
set(0,'defaultfigurecolor',[1 1 1])

%% initialize containers
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

% new: animal and recording IDs per neuron
allAnimalID = cell(1, nGroups);
allRecID = cell(1, nGroups);

%% load data
for iGroups = 1:nGroups

    fprintf('Loading data from group %d: %s\n', iGroups, opts.groups{iGroups});
    cPath = fullfile(opts.localPath, opts.groups{iGroups});
    recs = dir([cPath filesep '**' filesep '*npxResults.mat']);

    for iRecs = 1:numel(recs)

        clear cluster pulse ramp sensory waveMetrics
        saveFile = fullfile(recs(iRecs).folder, recs(iRecs).name);
        load(saveFile, 'cluster', 'pulse', 'ramp', 'sensory', 'waveMetrics');

        nCells = numel(sensory.baselineFR);

        % path expected as: ...\GROUP\animal_X\rec_Y\file.mat
        recFolder = recs(iRecs).folder;
        [animalFolder, recName] = fileparts(recFolder);
        [~, animalName] = fileparts(animalFolder);

        % make IDs unique across whole dataset
        animalID = [opts.groups{iGroups} '_' animalName];
        recID = [opts.groups{iGroups} '_' animalName '_' recName];

        % collect data for analysis
        allRecs{iGroups} = [allRecs{iGroups}; repmat(iRecs, nCells, 1)];
        allAnimalID{iGroups} = [allAnimalID{iGroups}; repmat({animalID}, nCells, 1)];
        allRecID{iGroups} = [allRecID{iGroups}; repmat({recID}, nCells, 1)];

        allBaselineFR{iGroups} = [allBaselineFR{iGroups}; mean(cat(2, sensory.baselineFR, sensory.optoBaselineFR), 2)];
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
visPSTH = cell(1, nGroups);
visFWHM = cell(1, nGroups);
visPeak = cell(1, nGroups);
visPeakAmp = cell(1, nGroups);
visMean = cell(1, nGroups); 
visAUC = cell(1, nGroups);

% matching IDs for visually responsive neurons
visAnimalID = cell(1, nGroups);
visRecID = cell(1, nGroups);

visResponseCells = nan(1, nGroups);
groupCells = nan(1, nGroups);

for iGroups = 1:nGroups
    cIdx = ~strcmpi(allLabels{iGroups}, 'noise')' & allDepth{iGroups} < 1500; % non-noise cortical clusters
    visIdx = allVisP{iGroups} < 0.05 & allVis{iGroups} > 0.5 & cIdx; % responsive neurons

    cData = allVisPSTH{iGroups}(visIdx,:);
    cData = cData - mean(cData(:, sensory.rampBins < 0), 2); % baseline correction

    visPSTH{iGroups} = cData;
    visAUC{iGroups} = allVis{iGroups}(visIdx);
    visAnimalID{iGroups} = allAnimalID{iGroups}(visIdx);
    visRecID{iGroups} = allRecID{iGroups}(visIdx);

    visResponseCells(iGroups) = sum(visIdx);
    groupCells(iGroups) = sum(cIdx);

    % preallocate
    nResp = size(cData,1);
    visMean{iGroups} = nan(nResp,1);
    visFWHM{iGroups} = nan(nResp,1);
    visPeak{iGroups} = nan(nResp,1);
    visPeakAmp{iGroups} = nan(nResp,1);

    respTime = sensory.rampBins > 0 & sensory.rampBins < 0.2;
    baseTime = sensory.rampBins < 0 & sensory.rampBins > -0.2;
    cTime = sensory.rampBins(respTime);

    for iClust = 1:nResp
        if visAUC{iGroups}(iClust) > 0.6
            visMean{iGroups}(iClust) = mean(cData(iClust, respTime)) - mean(cData(iClust, baseTime));

            cCluster = cData(iClust, respTime);
            cCluster(1) = 0;
            cCluster(end) = 0;
            cCluster = smooth(cCluster, 10);

            visFWHM{iGroups}(iClust) = fwhm(cTime, cCluster) * 1000; % ms
            [visPeakAmp{iGroups}(iClust), temp] = max(cCluster);
            visPeak{iGroups}(iClust) = cTime(temp) * 1000; % ms
        end
    end
end

[~, fracP] = binomCompare(visResponseCells(1), groupCells(1), visResponseCells(2), groupCells(2)); 
visResponseFraction = visResponseCells ./ groupCells; 

%% show visual response trace
figure
xlim([-0.2 1]);
ylim([-5 25]);
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

statsPeak = plotViolinWithLME( ...
    1, 3, 1, ...
    visPeak, visAnimalID, visRecID, ...
    ctrlColor, KOcolor, ...
    'Time of visual response peak (ms)', ...
    'Visual response peak after stimulus (V1)', ...
    [35 225]);

statsAmp = plotViolinWithLME( ...
    1, 3, 2, ...
    visPeakAmp, visAnimalID, visRecID, ...
    ctrlColor, KOcolor, ...
    'Visual peak response [Hz]', ...
    'Visual peak response (V1)', ...
    []);

statsFWHM = plotViolinWithLME( ...
    1, 3, 3, ...
    visFWHM, visAnimalID, visRecID, ...
    ctrlColor, KOcolor, ...
    'Full-width-half-maximum (ms)', ...
    'Visual response width (V1)', ...
    [0 150]);

%% inspect LME model outputs
statsSummary = table( ...
    {'Visual response peak time'; ...
     'Visual peak amplitude'; ...
     'Visual response width (FWHM)'}, ...
    [statsPeak.p; statsAmp.p; statsFWHM.p], ...
    [statsPeak.t; statsAmp.t; statsFWHM.t], ...
    'VariableNames', {'Measure','pValue_LME','tStat_LME'});

disp(' ');
disp('===== Visual response LME summary =====');
disp(statsSummary);

%% plotting function for violins
function stats = plotViolinWithLME(nRows, nCols, plotIdx, dataByGroup, animalIDByGroup, recIDByGroup, ctrlColor, KOcolor, yLabelTxt, titleTxt, yLimits)

subplot(nRows, nCols, plotIdx); hold on;

% use same filtering for plotting and statistics
keepIdx = cell(1,2);
plotData = cell(1,2);
plotAnimalID = cell(1,2);
plotRecID = cell(1,2);

for iG = 1:2
    cVals = dataByGroup{iG};
    keepIdx{iG} = ~isnan(cVals) & ~isoutlier(cVals);

    plotData{iG} = cVals(keepIdx{iG});
    plotAnimalID{iG} = animalIDByGroup{iG}(keepIdx{iG});
    plotRecID{iG} = recIDByGroup{iG}(keepIdx{iG});
end

% violin plot
h(1) = Violin(plotData{1}, 0, ...
    'ViolinColor', ctrlColor, ...
    'EdgeColor', ctrlColor, ...
    'BoxColor', [0 0 0]);

h(2) = Violin(plotData{2}, 1, ...
    'ViolinColor', KOcolor, ...
    'EdgeColor', KOcolor, ...
    'BoxColor', [0 0 0]);

axis square
if ~isempty(yLimits)
    ylim(yLimits);
end

ax = h(1).ViolinPlot.Parent;
ax.XTick = 0:1;
ax.XTickLabel = {'Control', 'DNMT1-KO'};
ylabel(yLabelTxt);

% build inputs for LME
dataIn = [plotData{1}; plotData{2}];
conditionID = [ones(numel(plotData{1}),1); 2*ones(numel(plotData{2}),1)];
animalVar = [plotAnimalID{1}; plotAnimalID{2}];

[pVal, tStat, fullmodel, modelCompare] = LME_compareMulti_SM(dataIn, conditionID, animalVar);

title(sprintf('%s\nLME p = %.4g', titleTxt, pVal));
niceFigure

% return useful outputs
stats = struct;
stats.p = pVal;
stats.t = tStat;
stats.fullmodel = fullmodel;
stats.modelCompare = modelCompare;
stats.data = plotData;
stats.animalID = plotAnimalID;
stats.recID = plotRecID;
end

