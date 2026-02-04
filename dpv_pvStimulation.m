%% define experiments to use
recLabels = {'Animal','ExperimentDate','Location','Expertise','Folder','Probe','Background','Path','Ramppower','useRec','hasGamma'};
recInfo{1} = {'PV2555','26/01/2022','V1','Naive','PV2555_20221026','0','PV-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'PV2555','28/01/2022','V1','Naive','PV2555_20221028','0','PV-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'PV2555','29/01/2022','V1','Naive','PV2555_20221029','0','PV-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'PV2556','26/01/2022','V1','Naive','PV2556_20221026','0','PV-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'PV2556','27/01/2022','V1','Naive','PV2556_20221027','0','PV-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'PV2556','28/01/2022','V1','Naive','PV2556_20221028','0','PV-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'PV2556','29/01/2022','V1','Naive','PV2556_20221029','0','PV-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control'};
recInfo{2} = {'PV2632','09/02/2023','V1','Naive','PV2632_20230209','0','PV-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','PV-DNMT1';'PV2632','10/02/2023','V1','Naive','PV2632_20230210','0','PV-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','PV-DNMT1';'PV2632','11/02/2023','V1','Naive','PV2632_20230211','0','PV-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','PV-DNMT1';'PV2633','24/01/2023','V1','Naive','PV2633_20230124','0','PV-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','PV-DNMT1';'PV2633','25/01/2023','V1','Naive','PV2633_20230125','0','PV-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','PV-DNMT1';'PV2633','26/01/2023','V1','Naive','PV2633_20230126','0','PV-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','PV-DNMT1';'PV2633','27/01/2023','V1','Naive','PV2633_20230127','0','PV-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','PV-DNMT1'};

%% basic variables
clear opts
opts.Location = {'V1', 'S1'};
opts.groups = {'PV-Cre', 'PV-DNMT1'};
opts.cortexRange = 1000;
opts.localPath = 'D:\DNMT1_project\';
opts.gid = 'PassiveStimulation';
opts.redoAnalysis = false;
opts.loadRaw = false;
opts.relevantVarNames = ["PulseCount", "PulseDur",  'PulseGap', 'OptoRamp', 'OptoShift', ... 
        'StimType', 'RedPowerOne', 'RedPowerTwo', 'BluePowerOne', 'BluePowerTwo'];
opts.stimBin = 0.001; %1ms bins
opts.stimWindow = [-1 1.5]; %psth range for stimulus sequence
opts.pulseWindow = [-0.005 0.005]; %psth range for single pulses
opts.rampRange = 0.5; %ramp range to compute AUCs
opts.showSpike = false;

ctxRange = [0 1500];
ctrlColor = [212 212 212]./255;
KOcolor = [255 160 64]./255;
set(0,'defaultfigurecolor',[1 1 1])
opts.repoSavePath = 'F:\DNMT1_project\';

%% run over groups
usedRecs = cell(1, length(opts.groups));
failedRecs = cell(1, length(opts.groups));
allClustIDs = cell(1, length(opts.groups));
allMouseID = cell(1, length(opts.groups));
allRecID = cell(1, length(opts.groups));
allRecNames = cell(1, length(opts.groups));
allDepth = cell(1, length(opts.groups));
allOptoRamp = cell(2, length(opts.groups));
allRampPSTH = cell(2, length(opts.groups));
allOptoPulse = cell(1, length(opts.groups));
allWF = cell(1, length(opts.groups));
allPulsePSTH = cell(2, length(opts.groups));
allPulseSeqPSTH = cell(1, length(opts.groups));
allPulseDelay = cell(2, length(opts.groups));
allLabels = cell(1, length(opts.groups));
allMetrics = cell(1, length(opts.groups));
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
            opts.probeNr = recInfo{iGroups}{iRecs, strcmpi(recLabels, 'Probe')};
            cAnimal = recInfo{iGroups}{iRecs, strcmpi(recLabels, 'Animal')};
            if strcmpi(opts.probeNr, '0')
                opts.optoLaser = 'BluePowerOne';
            elseif strcmpi(opts.probeNr, '1')
                opts.optoLaser = 'BluePowerTwo';
            end
            
            % get data from current recording
            disp(opts.recName);
            
            % collect data for analysis
            [pulse, cluster, ramp, waveMetrics] = pC_optoResponsePV_SM(opts);

            if size(cluster.tempWFs,2) ~= length(pulse.AUCs)
                break
            end
            
            allRecNames{iGroups} = [allRecNames{iGroups}; repmat({opts.recName}, length(ramp.AUCs))];
            allLabels{iGroups} = [allLabels{iGroups}, cluster.clusterLabels];
            allDepth{iGroups} = [allDepth{iGroups}; cluster.clustDepth];
            allMetrics{iGroups} = appendBehavior(allMetrics{iGroups},waveMetrics);
            
            allOptoRamp{iGroups,1} = [allOptoRamp{iGroups,1}; ramp.AUCs{1}'];
            allOptoRamp{iGroups,2} = [allOptoRamp{iGroups,2}; ramp.AUCs{2}'];
            allRampPSTH{iGroups,1} = [allRampPSTH{iGroups,1}; ramp.PSTH{1}];
            allRampPSTH{iGroups,2} = [allRampPSTH{iGroups,2}; ramp.PSTH{2}];
            
            allOptoPulse{iGroups} = [allOptoPulse{iGroups}; pulse.AUCs'];
            allPulsePSTH{iGroups} = [allPulsePSTH{iGroups}; pulse.PSTH];
            allPulseSeqPSTH{iGroups} = [allPulseSeqPSTH{iGroups}; pulse.seqPSTH];
            allPulseDelay{iGroups} = [allPulseDelay{iGroups}; pulse.delays];
            allWF{iGroups} = [allWF{iGroups}; cluster.tempWFs'];
            allMetrics{iGroups} = appendBehavior(allMetrics{iGroups},waveMetrics);
            usedRecs{iGroups} = [usedRecs{iGroups}; fullfile(opts.recPath, opts.recName) ];
            
            % keep track of animal ID and sesseion counter
            mouseIdx = repmat(find(strcmpi(allMice, cAnimal)), length(pulse.AUCs), 1);
            allMouseID{iGroups} = [allMouseID{iGroups}; mouseIdx];

            sessionCnt = sessionCnt + 1;
            sessionIdx = repmat(sessionCnt, length(cluster.clusterIDs), 1);
            allRecID{iGroups} = [allRecID{iGroups}; mouseIdx];

        catch
            fprintf('%s failed\n', recInfo{iGroups}{iRecs,recIdx});
            lfpDir = [fullfile(opts.recPath, opts.recName) filesep opts.recName '_g0' filesep  opts.recName '_g0_imec' opts.probeNr];
            apFile = dir([lfpDir filesep '*ap.meta']);
            apFile = strrep(fullfile(lfpDir, apFile.name), 'ap.meta', 'ap.bin');
            failedRecs{iGroups} = [failedRecs{iGroups}; apFile];
        end
    end
end

%% show optogenetic suppression
figure; clear h rampOut meanOut
minThresh = 2;
optoData = cell(1, length(opts.groups));
rampPSTH = cell(1, length(opts.groups));
for iGroups = 1 : length(opts.groups)
    
    cIdx = strcmpi(allLabels{iGroups}, 'sua')' & allDepth{iGroups} > 0 & allDepth{iGroups} < 2000; %all sua clusters in cortex
    for x = 1 : 2
        cData = allRampPSTH{iGroups,x};
        cData = cData - mean(cData(:, 1:500),2); %baseline correction
        pyrIdx = cIdx & nanmean(cData(:, ramp.rampBins>-0.5 & ramp.rampBins < 0), 2) < -minThresh & nanmean(cData(:, ramp.rampBins>0.5 & ramp.rampBins < 1), 2) < -minThresh;
        rampOut{x} = cData(pyrIdx,:);
        meanOut{x} = mean(cData(pyrIdx, ramp.rampBins>-0.5 & ramp.rampBins<1), 2);
    end
    
    rampPSTH{iGroups} = cat(1, rampOut{:});
    optoData{iGroups} = cat(1, meanOut{:});
    optoData{iGroups} = rmoutliers(optoData{iGroups}(optoData{iGroups}<0), 'percentile', [5 95]);
    
end

subplot(1,2,1);
xlim([-0.5 2]);
nhline(0,'k--');
lines(2) = stdshade(rampPSTH{2}, 0.3, KOcolor, ramp.rampBins+0.5, 100);
lines(1) = stdshade(rampPSTH{1}, 0.3, ctrlColor, ramp.rampBins+0.5, 100);
axis square
legend(lines, 'Ctrl', 'DNMT-KO', 'Location', 'northwest');
title('Optogenetically induced suppression');
xlabel('time from laser onset (s)')
ylabel('Firing rate (Hz)');
ylim([-8 4])
niceFigure;

subplot(1,2,2);
h(1) = Violin((optoData{1}), 0, 'ViolinColor', ctrlColor, 'EdgeColor', ctrlColor, 'BoxColor', [0 0 0]);
h(2) = Violin((optoData{2}), 1, 'ViolinColor', KOcolor, 'EdgeColor', KOcolor, 'BoxColor', [0 0 0]);
nhline(0, 'k--');
ylim([-14 2]);
axis square
ax = h(1).ViolinPlot.Parent;
ax.XTick = 0:1;
ax.XTickLabel = {'Control', 'DNMT1-KO'};
ylabel('Firing rate (Hz)');
title(sprintf('Mean optogenetically-induced suppression\n Ctrl FR = %.2f Hz; DNMT1-KO FR = %.2f Hz; pVal = %d', ...
    median(optoData{1}),  median(optoData{2}), ranksum(optoData{1}, optoData{2})))
niceFigure

disp('==================')
disp('Optogenetic suppression:')
fprintf('%s: Peak response %.2f %c %.2fHz\n'  , opts.groups{1}, mean(optoData{1}), char(177), sem(optoData{1}));
fprintf('%s: Peak response %.2f %c %.2fHz\n'  , opts.groups{2}, mean(optoData{2}), char(177), sem(optoData{2}));
fprintf('pVal ranksum test: %f\n', ranksum(optoData{1}, optoData{2}))
fprintf('Nr of neural clusters %s: %i; %s: %i:\n', opts.groups{1}, length(optoData{1}), opts.groups{2},  length(optoData{2}));
disp('==================')
