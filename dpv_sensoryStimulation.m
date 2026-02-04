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
opts.localPath = 'F:\DNMT1_project\Ephys_data';
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
usedRecs = cell(1, length(opts.groups));
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
allRecNames = cell(1, length(opts.groups));
allClustIDs = cell(1, length(opts.groups));
allMouseID = cell(1, length(opts.groups));
allRecID = cell(1, length(opts.groups));
allOscRejCnt = cell(1, length(opts.groups));
sessionCnt = 0;
for iGroups = 1 : length(opts.groups)
    
    disp(['Loading data from group ' num2str(iGroups) ': ' opts.groups{iGroups}]);
    allInfo{iGroups} = recInfo{iGroups};
    
    recIdx = strcmpi(recLabels, 'Folder');
    pathIdx = strcmpi(recLabels, 'Path');
    allMice = unique(recInfo{iGroups}(:, strcmpi(recLabels, 'Animal'))); %mice for current group
    for iRecs = 1 : size(recInfo{iGroups}, 1)
        try
            
            % include current recordinf info to opts
            opts.recName = recInfo{iGroups}{iRecs,recIdx};
            opts.recPath = recInfo{iGroups}{iRecs,pathIdx};
            opts.probNr = recInfo{iGroups}{iRecs, strcmpi(recLabels, 'Probe')};
            cAnimal = recInfo{iGroups}{iRecs, strcmpi(recLabels, 'Animal')};
            
            % collect data for analysis
            [cluster, pulse, ramp, sensory] = pC_passiveFRcomparisonPV_SM(opts);

            allRecNames{iGroups} = [allRecNames{iGroups}; repmat({opts.recName}, size(sensory.baselineFR))];
            allRecs{iGroups} = [allRecs{iGroups}; ones(size(sensory.baselineFR)).*iRecs];
            allBaselineFR{iGroups} = [allBaselineFR{iGroups}; mean(cat(2, sensory.baselineFR, sensory.optoBaselineFR),2)];
            allLabels{iGroups} = [allLabels{iGroups}, cluster.clusterLabels];
            allClustIDs{iGroups} = [allClustIDs{iGroups}; cluster.clusterIDs];
            allDepth{iGroups} = [allDepth{iGroups}; cluster.clustDepth];
            allOscRejCnt{iGroups} = [allOscRejCnt{iGroups}; cluster.oscRejCnt];
            allOptoPulse{iGroups} = [allOptoPulse{iGroups}; pulse.AUCs'];
            allOptoRamp{iGroups} = [allOptoRamp{iGroups}; ramp.AUCs'];
            allVis{iGroups} = [allVis{iGroups}; sensory.AUCs'];
            allVisP{iGroups} = [allVisP{iGroups}; sensory.AUC_pVals'];
            allOptoPulseDelay{iGroups} = [allOptoPulseDelay{iGroups}; pulse.delays];
            allRampPSTH{iGroups} = [allRampPSTH{iGroups}; ramp.PSTH];
            allVisPSTH{iGroups} = [allVisPSTH{iGroups}; sensory.PSTH];
            allRampVisPSTH{iGroups} = [allRampVisPSTH{iGroups}; sensory.PSTH];
            usedRecs{iGroups} = [usedRecs{iGroups}; recInfo{iGroups}{iRecs,recIdx}];
            
            % keep track of animal ID and sesseion counter
            mouseIdx = repmat(find(strcmpi(allMice, cAnimal)), length(cluster.clusterIDs), 1);
            allMouseID{iGroups} = [allMouseID{iGroups}; mouseIdx];

            sessionCnt = sessionCnt + 1;
            sessionIdx = repmat(sessionCnt, length(cluster.clusterIDs), 1);
            allRecID{iGroups} = [allRecID{iGroups}; mouseIdx];

        catch
            fprintf('%s failed\n', recInfo{iGroups}{iRecs,recIdx});
            failedRecs{iGroups} = [failedRecs{iGroups}; recInfo{iGroups}{iRecs,recIdx}];
        end
    end
end

%% tactile responses - S1 - compare different depths, superficial vs deep
layerRange = [400, 2000];
visRange = [0, 0.025];
visPSTH = cell(length(opts.groups),2);
visFWHM = cell(length(opts.groups),2);
visPeak = cell(length(opts.groups),2);
visPeakAmp = cell(length(opts.groups),2);
visBaseline = cell(length(opts.groups),2);
visIDs = cell(2, length(opts.groups));
visAUC = cell(2, length(opts.groups));
visRespEarly = cell(2, length(opts.groups));
visRespEarlyTime = cell(2, length(opts.groups));
visRespLate = cell(2, length(opts.groups));
visResponseCells = nan(2, length(opts.groups));
groupCells = nan(2, length(opts.groups));
for iDepths = 1 : 2
    for iGroups = 1 : length(opts.groups)
        if iDepths == 1
            cIdx = allDepth{iGroups} < layerRange(1) & allDepth{iGroups} > -1000; %all non-noise clusters in cortex
        elseif iDepths == 2
            cIdx = allDepth{iGroups} > layerRange(1) & allDepth{iGroups} < layerRange(2); %all non-noise clusters in cortex
        end
        baselinePower = mean(allVisPSTH{iGroups}(~strcmpi(allLabels{iGroups}, 'noise')',sensory.rampBins<0 & sensory.rampBins > -0.1),2);
        visBaseline{iGroups,iDepths} = baselinePower; %get baseline
        visIdx = allVisP{iGroups} < 0.08 & allVis{iGroups} > 0.5 & cIdx; %all responsive neurons
        visIDs{iGroups,iDepths} = allMouseID{iGroups}(visIdx);
        cData = allVisPSTH{iGroups}(visIdx,:);
        cData = cData - mean(cData(:, sensory.rampBins<0 & sensory.rampBins > -0.1),2); %baseline correction
        visPSTH{iGroups,iDepths} = cData; %get PSTH
        visAUC{iGroups,iDepths} = allVis{iGroups}(visIdx); %get AUC
        visRespEarly{iGroups,iDepths} = nanmean(cData(:, sensory.rampBins > visRange(1) & sensory.rampBins < visRange(2)),2); %get early response
        [~, visRespEarlyTime{iGroups,iDepths}] = max(smoothCol(cData(:, sensory.rampBins > 0 & sensory.rampBins < 0.05), 2),[],2); %get early response
        useIdx = visRespEarly{iGroups,iDepths} > 0;
        visRespEarly{iGroups,iDepths} = visRespEarly{iGroups,iDepths}(useIdx);
        visRespEarlyTime{iGroups,iDepths} = visRespEarlyTime{iGroups,iDepths}(useIdx) * mean(diff(sensory.rampBins)) * 1000; %convert to miliseconds
        visRespLate{iGroups,iDepths} = nanmean(cData(useIdx, sensory.rampBins > 0.05 & sensory.rampBins < 0.2),2); %get late response
        visIDs{iGroups,iDepths} = visIDs{iGroups,iDepths}(useIdx);
        
        % count cells
        visResponseCells(iGroups,iDepths) = sum(visIdx);
        groupCells(iGroups,iDepths) = sum(cIdx);
        
        % check amplitude and temporal precision of sensory response
        respTime = sensory.rampBins > 0 & sensory.rampBins < 0.05;
        cTime = sensory.rampBins(respTime);
        for iClust = 1 : size(cData,1)
            if visAUC{iGroups,iDepths}(iClust) > 0.6 %needs strong enough response to work
                % width of response
                cCluster = cData(iClust, respTime);
                cCluster(1) = 0;
                cCluster(end) = 0;
                cCluster = smooth(cCluster, 10);
                visFWHM{iGroups,iDepths}(iClust) = fwhm(cTime, cCluster)*1000; %full width half maximum in miliseconds
                
                %time of response
                [visPeakAmp{iGroups,iDepths}(iClust),temp] = max(cCluster);
                visPeak{iGroups,iDepths}(iClust) = cTime(temp) * 1000; %time of response peak in miliseconds
            else
                visFWHM{iGroups,iDepths}(iClust) = nan;
                visPeak{iGroups,iDepths}(iClust) = nan;
                visPeakAmp{iGroups,iDepths}(iClust) = nan;
            end
        end
        visFWHM{iGroups,iDepths}(visFWHM{iGroups,iDepths}>50) = NaN;
    end
end

% show response trace
figure
t = tiledlayout(2, 3);
t.Padding = 'compact'; % 'compact' or 'none'
cTitles = {'Supragranular', 'Infragranular'};
clear h
for iDepth = 1 : 2
    % show traces
    nexttile;
    xlim([-0.05 0.4]);
    nhline(0,'k--');
    nvline([0 0.5], 'k--');
    cData = visPSTH(:, iDepth);
    lines(1) = stdshade(cData{1}, 0.3, ctrlColor, ramp.rampBins, 5);
    lines(2) = stdshade(cData{2}, 0.3, KOcolor, ramp.rampBins+0.005, 5);
    axis square
    legend(lines, 'Ctrl', 'DNMT-KO', 'Location', 'northeast');
    title(['Tactile response - ' cTitles{iDepth}]);
    xlabel('time from stimulus (s)');
    ylabel('firing rate (Hz)');
    ylim([-2 35]);
    niceFigure;
    
    %show early visual response
    nexttile;
    cData = visRespEarly(:, iDepth);
    h(1) = Violin(rmoutliers(cData{1}, 'percentiles', [2 98]), 0, 'ViolinColor', ctrlColor, 'EdgeColor', ctrlColor, 'BoxColor', [0 0 0]);
    h(2) = Violin(rmoutliers(cData{2}, 'percentiles', [2 98]), 1, 'ViolinColor', KOcolor, 'EdgeColor', KOcolor, 'BoxColor', [0 0 0]);
    xlim([-0.5 1.5]);
    nhline(0, 'k--');
    ylim([-2 70]);
    axis square
    ax = h(1).ViolinPlot.Parent;
    ax.XTick = 0:1;
    ax.XTickLabel = {'Control', 'DNMT1-KO'};
    ylabel('firing rate (hz)');
    
    %LME test
    mouseIdx = visIDs(:,iDepth);
    mouseIdx{2} = mouseIdx{2} + max(unique(mouseIdx{1}));
    mouseIdx = cat(1, mouseIdx{:});
    actData = cat(1, cData{:});
    conditionID = cat(1, ones(length(cData{1}),1), ones(length(cData{2}),1)*2);
    [pVal, cTStat] = LME_compare(actData, conditionID, mouseIdx);
    
    title(sprintf('%s: Early tactile response (0-30ms; S1)\n Ctrl = %.2f hz; DNMT1-KO = %.2f hz\n LME pVal = %d,  t-stat = %.6f' , ...
        cTitles{iDepth}, median(cData{1}),  median(cData{2}), pVal, cTStat));
    niceFigure
    
    disp('==================')
    disp([cTitles{iDepth} ': Early spiking response:'])
    fprintf('%s: Peak response %.2f %c %.2fHz\n'  , opts.groups{1}, mean(cData{1}), char(177), sem(cData{1}));
    fprintf('%s: Peak response %.2f %c %.2fHz\n'  , opts.groups{2}, mean(cData{2}), char(177), sem(cData{2}));
    fprintf('LME test: pVal = %.6f, tStat = %.2f\n', pVal, cTStat);
    fprintf('Nr of neural clusters %s: %i; %s: %i:\n', opts.groups{1}, visResponseCells(1,1), opts.groups{2},  visResponseCells(2,1));

    %show late visual response
    nexttile;
    cData = visRespLate(:, iDepth);    
    h(1) = Violin(rmoutliers(cData{1}, 'percentiles', [5 95]), 0, 'ViolinColor', ctrlColor, 'EdgeColor', ctrlColor, 'BoxColor', [0 0 0]);
    h(2) = Violin(rmoutliers(cData{2}, 'percentiles', [5 95]), 1, 'ViolinColor', KOcolor, 'EdgeColor', KOcolor, 'BoxColor', [0 0 0]);
    xlim([-0.5 1.5]);
    nhline(0, 'k--');
    ylim([-2 27.5]);
    axis square
    ax = h(1).ViolinPlot.Parent;
    ax.XTick = 0:1;
    ax.XTickLabel = {'Control', 'DNMT1-KO'};
    ylabel('firing rate (hz)');
    
    %LME test
    actData = cat(1, cData{:});
    [pVal, cTStat] = LME_compare(actData, conditionID, mouseIdx);
    
    title(sprintf('%s: Late tactile response (50-100ms; S1)\n Ctrl = %.2f hz; DNMT1-KO = %.2f hz\n LME pVal = %d,  t-stat = %.6f' , ...
        cTitles{iDepth}, median(cData{1}),  median(cData{2}), pVal, cTStat));
    
    niceFigure
    
    disp('==================')
    disp([cTitles{iDepth} ': Late spiking response:'])
    fprintf('%s: Peak response %.2f %c %.2fHz\n'  , opts.groups{1}, mean(cData{1}), char(177), sem(cData{1}));
    fprintf('%s: Peak response %.2f %c %.2fHz\n'  , opts.groups{2}, mean(cData{2}), char(177), sem(cData{2}));
    fprintf('LME test: pVal = %.6f, tStat = %.2f\n', pVal, cTStat);
    disp('==================')
    
end

%% show early and late quantification a bit differently in the same figure
% figure
% 
% %show early visual response
% Cnt = 0;
% violDist = 0.75;
% violIdx = [];
% for iDepth = 1 : 2
%     cData = visRespEarly(:, iDepth);
%     cData{1} = rmoutliers(cData{1}, 'percentiles', [2 95]);
%     cData{2} = rmoutliers(cData{2}, 'percentiles', [2 95]);
%     h(1) = Violin(cData{1}, Cnt, 'ViolinColor', ctrlColor, 'EdgeColor', ctrlColor, 'BoxColor', [0 0 0]);
%     h(2) = Violin(cData{2}, violDist + Cnt, 'ViolinColor', KOcolor, 'EdgeColor', KOcolor, 'BoxColor', [0 0 0]);
%     violIdx = [violIdx, Cnt, Cnt + violDist];
%     Cnt = Cnt + 1.75;
%     fprintf('Supragranular clusters; %s: %i; %s: %i:\n', opts.groups{1}, length(cData{1}), opts.groups{2},  length(cData{2}));
% end
% 
% % show late visual response
% Cnt = Cnt + 1;
% for iDepth = 1 : 2
%     cData = visRespLate(:, iDepth);
%     cData{1} = rmoutliers(cData{1}, 'percentiles', [2 95]);
%     cData{2} = rmoutliers(cData{2}, 'percentiles', [2 95]);
%     h(1) = Violin(cData{1}, Cnt, 'ViolinColor', ctrlColor, 'EdgeColor', ctrlColor, 'BoxColor', [0 0 0]);
%     h(2) = Violin(cData{2}, violDist + Cnt, 'ViolinColor', KOcolor, 'EdgeColor', KOcolor, 'BoxColor', [0 0 0]);
%     violIdx = [violIdx, Cnt, Cnt + violDist];
%     Cnt = Cnt + 1.75;
%     fprintf('Infragranular clusters; %s: %i; %s: %i:\n', opts.groups{1}, length(cData{1}), opts.groups{2},  length(cData{2}));
% end
% 
% ax = h(1).ViolinPlot.Parent;
% ax.XTick = violIdx;
% ax.XTickLabel = {'ctrl', 'KO'};
% ylabel('firing rate [hz]');
% ylim([-5 45]);
% nhline(0, 'k--');
% xlim([-0.75 Cnt-0.25]);
% niceFigure

%% check spontaneous firing rate versus depth
layerRange = [400, 1000];

figure; hold on;
clear h
spontFiring = cell(length(opts.groups), 2);
spontFiringIDs = cell(length(opts.groups), 2);
for iGroups = 1 : length(opts.groups)
    
    basePower = mean(allVisPSTH{iGroups}(strcmpi(allLabels{iGroups}, 'sua')',sensory.rampBins < 0 & sensory.rampBins > -1),2);
    baseDepth = allDepth{iGroups}(strcmpi(allLabels{iGroups}, 'sua'));
    cMouseID = allMouseID{iGroups}(strcmpi(allLabels{iGroups}, 'sua'));

    useIdx = baseDepth > 0 & baseDepth < layerRange(2);
    baseDepth = baseDepth(useIdx);
    basePower = basePower(useIdx);
    baseAnimalID = cMouseID(useIdx);
    
    % keep activity for statistics
    spontFiring{iGroups, 1} = basePower(baseDepth < layerRange(1));
    spontFiring{iGroups, 2} = basePower(baseDepth > layerRange(1) & baseDepth <= layerRange(2));

    % keep mouse IDs for statistics
    spontFiringIDs{iGroups, 1} = baseAnimalID(baseDepth < layerRange(1));
    spontFiringIDs{iGroups, 2} = baseAnimalID(baseDepth > layerRange(1) & baseDepth <= layerRange(2));
        
    % show error bar plot across depth
    baseDepth = baseDepth / 1000;
    binEdges = 0 : 0.2 : 1; % Binning depths
    binCenters = binEdges(1:end-1) + diff(binEdges)/2; % Calculate bin centers
    
    % Calculate mean and SEM for each bin
    [~, ~, binIdx] = histcounts(baseDepth, binEdges); % Bin indices
    meanActivity = arrayfun(@(i) median(basePower(binIdx == i & basePower > 0), 'omitnan'), 1:numel(binCenters));
    semActivity = arrayfun(@(i) std(basePower(binIdx == i & basePower > 0), 'omitnan') / sqrt(sum(binIdx == i)), 1:numel(binCenters));
    meanActivity(isnan(meanActivity)) = 0;
    semActivity(isnan(semActivity)) = 0;
    
    % Create the error bar plot
    if iGroups == 1
        useColor = ctrlColor;
    elseif iGroups == 2
        useColor = KOcolor;
    end
    h(iGroups) = errorbar(binCenters, meanActivity, semActivity, '-o', 'Color', useColor, 'LineWidth', 1.5, 'MarkerSize', 6);
    xlabel('depth [mm]');
    ylabel('mean firing rate [hz]');
    
end
xlim([0, 1]);
nhline(0, 'k--');
grid on
niceFigure;
h(1).LineWidth = 3;
h(1).MarkerSize = 8;
h(2).LineWidth = 3;
h(2).MarkerSize = 8;
legend(h, opts.groups);
axis square
view([90, 90]);  % This flips x and y axes

disp('==================')
useData = spontFiring(:,1);
disp('Superficial spont firing:')
fprintf('%s: Firing rate %.2f %c %.2fHz\n'  , opts.groups{1}, median(useData{1}(useData{1}>0)), char(177), sem(useData{1}(useData{1}>0)));
fprintf('%s: Peak response %.2f %c %.2fHz\n'  , opts.groups{2}, median(useData{2}(useData{2}>0)), char(177), sem(useData{2}(useData{2}>0)));
% fprintf('pVal ranksum test: %f\n', ranksum(useData{1}(useData{1}>0), useData{2}(useData{2}>0)))
fprintf('Nr of neural clusters %s: %i; %s: %i:\n', opts.groups{1}, length(useData{1}), opts.groups{2},  length(useData{2}));

% LME test
mouseIdx = spontFiringIDs(:,1);
mouseIdx{2} = mouseIdx{2} + max(unique(mouseIdx{1}));
mouseIdx = cat(1, mouseIdx{:});
actData = cat(1, useData{:});
conditionID = cat(1, ones(length(useData{1}),1), ones(length(useData{2}),1)*2); 

% do test and give feedback
[pVal, cTStat] = LME_compare(actData, conditionID, mouseIdx);
fprintf('LME test: pVal = %.6f, tStat = %.2f\n', pVal, cTStat);
disp('==================')

useData = spontFiring(:,2);
disp('Infragranular spont firing:')
fprintf('%s: Firing rate %.2f %c %.2fHz\n'  , opts.groups{1}, median(useData{1}(useData{1}>0)), char(177), sem(useData{1}(useData{1}>0)));
fprintf('%s: Peak response %.2f %c %.2fHz\n'  , opts.groups{2}, median(useData{2}(useData{2}>0)), char(177), sem(useData{2}(useData{2}>0)));
% fprintf('pVal ranksum test: %f\n', ranksum(useData{1}(useData{2}>0), useData{2}(useData{2}>0)))
fprintf('Nr of neural clusters %s: %i; %s: %i:\n', opts.groups{1}, length(useData{1}), opts.groups{2},  length(useData{2}));

% LME test
mouseIdx = spontFiringIDs(:,2);
mouseIdx{2} = mouseIdx{2} + max(unique(mouseIdx{1}));
mouseIdx = cat(1, mouseIdx{:});
actData = cat(1, useData{:});
conditionID = cat(1, ones(length(useData{1}),1), ones(length(useData{2}),1)*2); 

% do test and give feedback
[pVal(2), cTStat(2)] = LME_compare(actData, conditionID, mouseIdx);
fprintf('LME test: pVal = %.6f, tStat = %.2f\n', pVal(2), cTStat(2));
disp('==================')

title(sprintf('Spontaneous firing vs depth\n pVal_s_u_p = %.6f, tStat_s_u_p = %.2f\n pVal_i_n_f_r_a = %.6f, tStat_i_n_f_r_a = %.2f' ...
    , pVal(1), cTStat(1), pVal(2), cTStat(2)));