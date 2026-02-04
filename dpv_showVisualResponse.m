%% define experiments to use
recLabels = {'Animal','ExperimentDate','Location','Expertise','Folder','Probe','Background','Path','Ramppower','useRec','hasGamma'};
recInfo{1} = {'SOM2554','26/01/2022','S1','Naive','SOM2554_20221026','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2554','27/01/2022','S1','Naive','SOM2554_20221027','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2554','28/01/2022','S1','Naive','SOM2554_20221028','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2554','29/01/2022','S1','Naive','SOM2554_20221029','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2563','26/01/2022','S1','Naive','SOM2563_20221026','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2563','27/01/2022','S1','Naive','SOM2563_20221027','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2563','28/01/2022','S1','Naive','SOM2563_20221028','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2563','29/01/2022','S1','Naive','SOM2563_20221029','1','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1'};
recInfo{2} = {'SOM2626','08/02/2023','S1','Naive','SOM2626_20230208','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2626','09/02/2023','S1','Naive','SOM2626_20230209','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2626','10/02/2023','S1','Naive','SOM2626_20230210','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2626','11/02/2023','S1','Naive','SOM2626_20230211','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2627','24/01/2023','S1','Naive','SOM2627_20230124','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2627','25/01/2023','S1','Naive','SOM2627_20230125','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2627','26/01/2023','S1','Naive','SOM2627_20230126','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1';'SOM2627','27/01/2023','S1','Naive','SOM2627_20230127','1','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1'};

%% basic variables
opts.Location = {'S1'};
opts.makePlot = false;
opts.removeAvg = false;
opts.useTaper = false;
opts.groups = {'SOM-Cre', 'SOM-DNMT1'};
opts.baseDur = 0.5;
opts.postStim = 1;
opts.stepSize = 4;
opts.stimType = 3; %1 for vision, 2 for audio, 3 for tactile
opts.brainRange = [-200, 1500]; %depth range for plotting
opts.verbose = false; %flag to supress some of the text outputs
opts.brainThresh = 1;
opts.groupColors = {([212, 212, 212]./255), ([255, 160, 64]./255), [0, 0, 1], [1, 0, 0]};
opts.optoGenetics = false; %flag to isolate sensory responses during optogenetics
opts.layerRange = [500, 1000]; %range for computing CSD traces
opts.layerDepths = 25 : 75 : 1000; %depths of different layers in microns
ctrlColor = [212 212 212]./255;
KOcolor = [255 160 64]./255;
opts.savePath = 'D:\DNMT1_project\';
opts.reload = false;

%% run over groups
nrGroups = length(opts.groups);
usedRecs = cell(1, length(opts.groups));
failedRecs = cell(1, length(opts.groups));
allCSD = cell(1, length(opts.groups));
allLFP = cell(1, length(opts.groups));
allStimStd = cell(1, length(opts.groups));
for iGroups = 1 : nrGroups
    
    recIdx = strcmpi(recLabels, 'Folder');
    pathIdx = strcmpi(recLabels, 'Path');
    allInfo{iGroups} = recInfo{iGroups};
    
    for iRecs = 1 : size(recInfo{iGroups}, 1)
        try
            
            [csdResp, lfpResp, meanStd, normDepth, lfMeta] = ds_loadCSDresponse(recInfo{iGroups}(iRecs,:), recLabels, opts);
            recName = recInfo{iGroups}{iRecs,strcmpi(recLabels, 'Folder')};
            recPath = recInfo{iGroups}{iRecs,strcmpi(recLabels, 'Path')};

            % keep info in larger array
            if isempty(allCSD{iGroups})
                allCSD{iGroups} = {csdResp};
                allLFP{iGroups} = {lfpResp};
                allStimStd{iGroups} = {meanStd};
                usedRecs{iGroups} = {fullfile(recInfo{iGroups}{iRecs,pathIdx}, recInfo{iGroups}{iRecs,recIdx})};
            else
                allCSD{iGroups} = [allCSD{iGroups}; csdResp];
                allLFP{iGroups} = [allLFP{iGroups}; lfpResp];
                allStimStd{iGroups} = [allStimStd{iGroups}; meanStd];
                usedRecs{iGroups} = [usedRecs{iGroups}; {fullfile(recInfo{iGroups}{iRecs,pathIdx}, recInfo{iGroups}{iRecs,recIdx})}];
            end
        catch ME
            disp(['Aborted: ' ME.message]);
            failedRecs{iGroups} = {failedRecs{iGroups}; recInfo{iGroups}{iRecs,recIdx}};
        end
    end
end

%% make brain surface figure
% figure
% for iGroups = 1 : nrGroups
%     % show alignment to brain surface to check for inconsistencies
%     subplot(1,nrGroups,iGroups); hold on;
%     cStd = cat(2,allStimStd{iGroups}{:});
%     plot(cStd)
%     title([opts.groups{iGroups} '; Brain surface for each recording']);
%     axis square
%     nhline(opts.brainThresh);
% end

%% combined figures - LFP and CSD for each group
normStim = cell(1, length(opts.groups));
cRange = 1.5*1E6; %for tactile
lfpRange = 0.7; %for vision
stimLabel = 'Tactile';
traceRange = [-7 2]; %for CSD

if opts.stimType == 1
    lfpRange = 0.7;
    csdRange = 6*1E5;
elseif opts.stimType == 3
    lfpRange = 0.7;
    csdRange = 1*1E6;
end

clear lines
depthLabel = {'Supragranular', 'Infragranular'};
respLabels = {'LFP absolute', 'CSD absolute', 'CSD rescale'};
h1 = figure('renderer' ,'painters');
h2 = figure('renderer' ,'painters');
xRange = [(opts.baseDur - 0.05), (opts.baseDur + 0.1)] .*lfMeta.sRateHz;
baseDur = round(opts.baseDur .*lfMeta.sRateHz);
nrLayers = length(opts.layerDepths);
peakResp = cell(nrGroups,2);
peakTime = cell(nrGroups,2);
depthPeaks = cell(nrGroups,nrLayers);

for iGroups = 1 : nrGroups
    for x = 1 : 2
        
        figure(h1);
        subplot(2, nrGroups, (x-1)*nrGroups + iGroups)
        if x == 1
            mergeStim = cat(3,allLFP{iGroups}{:});
            mergeStim = mergeStim - nanmean(mergeStim(:, 1:baseDur, :), 2);
            cRange = lfpRange;
            depthRange = (1 : opts.stepSize*10 : size(mergeStim,1)*opts.stepSize*10)-1;
            
        elseif x == 2
            mergeStim = cat(3,allCSD{iGroups}{:});
            mergeStim = mergeStim - nanmean(mergeStim(:, 1:baseDur, :), 2);
%             cRange = csdRange;
            cRange = prctile(abs(mergeStim(:)),99.5);
            depthRange = round(normDepth * 1E6);
            
        elseif x == 3
            mergeStim = cat(3,allCSD{iGroups}{:});
            mergeStim = mergeStim - nanmean(mergeStim(:, 1:baseDur, :), 2);
            cRange = prctile(abs(mergeStim(:)),99.7);
            depthRange = round(normDepth * 1E6);
        end
        
        % show image
        cImg = imagesc(nanmean(mergeStim, 3));
        ax = cImg.Parent;
        ax.XTick = 1 : lfMeta.sRateHz*0.05 : size(mergeStim,2);
        ax.XTickLabel = (((0 : lfMeta.sRateHz*0.05 : size(mergeStim,2)-1)./lfMeta.sRateHz) - opts.baseDur)*1000;
        xlim(xRange);
        caxis([-cRange cRange]);
        colormap(ax, parula(256));
        colorbar
        niceFigure;
        
        title([opts.Location{1} ' - ' opts.groups{iGroups} '; ' stimLabel ' stimulation; ' respLabels{x} '; optogenetics: ' num2str(opts.optoGenetics)]);
        xlabel('time(ms)')
        nvline(opts.baseDur*lfMeta.sRateHz, 'w', 'linewidth', 2);
        
        depthSteps = round(200 / mean(diff(depthRange)));
        ax.YTick = 1 : depthSteps : size(mergeStim,1);
        ax.YTickLabels = depthRange(1 : depthSteps : length(depthRange));
        [~, ctxStart] = min(abs(depthRange));
        nhline(ctxStart-0.5, 'w', 'linewidth', 2);
        ylabel('Depth (um)');
        ylim([0 find(depthRange > 1000, 1)]);
        axis square;
        drawnow;
        
        % plot traces from different depths
        if x == 2
                        
            % get depth profile
            for iSteps = 1 : nrLayers
                if iSteps == 1
                    cIdx = depthRange < opts.layerDepths(iSteps);
                else
                    cIdx = depthRange > opts.layerDepths(iSteps-1) & depthRange < opts.layerDepths(iSteps);
                end
                
                cData = squeeze(nanmean(mergeStim(cIdx, :, :), 1));
                t = (0 : size(cData,1)-1)/ lfMeta.sRateHz - opts.baseDur;
                stimOn = find(abs(t) == min(abs(t)));
                cData = cData - cData(stimOn, :);
                cIdx = t>0 & t < 0.1;
                depthPeaks{iGroups,iSteps} = min(cData(cIdx,:));
            end

            % get traces from different depths
            mergeTrace = cell(1,2);
            cIdx = depthRange < opts.layerRange(1);
            mergeTrace{1} = squeeze(nanmean(mergeStim(cIdx, :, :), 1));
            cIdx = depthRange > opts.layerRange(1) & depthRange < opts.layerRange(2);
            mergeTrace{2} = squeeze(nanmean(mergeStim(cIdx, :, :), 1));
            
            figure(h2)
            % show image again
            subplot(2,2, iGroups);
            cImg = imagesc(nanmean(mergeStim, 3));
            ax = cImg.Parent;
            ax.XTick = 1 : lfMeta.sRateHz*0.05 : size(mergeStim,2);
            ax.XTickLabel = (((0 : lfMeta.sRateHz*0.05 : size(mergeStim,2)-1)./lfMeta.sRateHz) - opts.baseDur)*1000;
            xlim(xRange);
            colorRange = prctile(abs(mergeStim(:)),99.5);
            caxis([-colorRange colorRange]);
            colormap(ax, parula(256));
            colorbar
            niceFigure;
            
            title([opts.Location{1} ' - ' opts.groups{iGroups} '; ' stimLabel ' stimulation; ' respLabels{x} '; optogenetics: ' num2str(opts.optoGenetics)]);
            xlabel('time(ms)')
            nvline(opts.baseDur*lfMeta.sRateHz, 'w', 'linewidth', 2);
            
            depthSteps = round(200 / mean(diff(depthRange)));
            ax.YTick = 1 : depthSteps : size(mergeStim,1);
            ax.YTickLabels = depthRange(1 : depthSteps : length(depthRange));
            [~, ctxStart] = min(abs(depthRange));
            nhline(ctxStart-0.5, 'w', 'linewidth', 2);
            ylabel('Depth (um)');
            ylim([0 find(depthRange > 1000, 1)]);
            axis square;
            drawnow;
            
            for iDepth = 1 : 2
                subplot(2,2,iDepth+2)
                t = (0 : size(mergeTrace{iDepth},1)-1)./ lfMeta.sRateHz - opts.baseDur;
                xlim([-0.05 0.25]); axis square;
                ylim(traceRange)
                nvline(0, '--k');
                nhline(0, '--k');
                mergeTrace{iDepth} = mergeTrace{iDepth} - mergeTrace{iDepth}(abs(t) == min(abs(t)),:);
                lines(iGroups) = stdshade(mergeTrace{iDepth}', 0.5, opts.groupColors{iGroups}, t);
                title([opts.Location{1} '; ' stimLabel ' stimulation; ' respLabels{x} '; ' depthLabel{iDepth}]);
                ylabel('LFP deflection (uV)');
                xlabel('time (s)');
                niceFigure;
                cIdx = t>0 & t < 0.1;
                [peakResp{iGroups, iDepth}, peakTime{iGroups, iDepth}] = min(mergeTrace{iDepth}(cIdx,:));
            end

        end
    end
end
legend(lines, opts.groups)

%% show some statistics
for iDepth = 1 : 2
disp('====================')
disp(depthLabel(iDepth));
disp('CSD peak response:')
fprintf('%s: Peak response %.2f %c %.2fuV\n'  , opts.groups{1}, mean(peakResp{1, iDepth}), char(177), sem(peakResp{1, iDepth}));
fprintf('%s: Peak response %.2f %c %.2fuV\n'  , opts.groups{2}, mean(peakResp{2, iDepth}), char(177), sem(peakResp{2, iDepth}));
fprintf('pVal ranksum test: %f\n', ranksum(peakResp{1, iDepth}, peakResp{2, iDepth}))

disp('====================')
zeroTime = find(t > 0, 1);
disp('CSD peak time:')
fprintf('%s: Peak time %.2f %c %.2fms\n'  , opts.groups{1}, mean(t(peakTime{1, iDepth}+zeroTime)*1000), char(177), sem(t(peakTime{1, iDepth}+zeroTime)*1000));
fprintf('%s: Peak time %.2f %c %.2fms\n'  , opts.groups{2}, mean(t(peakTime{2, iDepth}+zeroTime)*1000), char(177), sem(t(peakTime{2, iDepth}+zeroTime)*1000));
fprintf('pVal ranksum test: %f\n', ranksum(peakTime{1}, peakTime{2}))
disp('====================')
end

%% more CSD analysis
figure; hold on;
cData = cellfun(@mean,depthPeaks);
cError = cellfun(@sem,depthPeaks);

clear cLines
plotDephts = opts.layerDepths - opts.layerDepths(1);
cLines(1) = errorshade(plotDephts, cData(1,:), cError(1,:), cError(1,:), ctrlColor, 0.2);
cLines(2) = errorshade(plotDephts, cData(2,:), cError(2,:), cError(2,:), KOcolor, 0.2);

axis square;
title('Mean response difference');
view(270,270);
legend(opts.groups, 'location', 'northwest');
niceFigure
for x = 1:length(cLines)
    cLines(x).LineWidth = 4;
end
ylabel('current source density [uV/mm^2]');
xlabel('depth [mm]');
grid on

%% more CSD analysis
layerLabels = cell(1, length(opts.layerDepths));
layerLabels{1} = ['0 - ' num2str(opts.layerDepths(1))];

% Generate labels for each layer
for i = 1:length(opts.layerDepths)
    layerLabels{i} = sprintf('%i', opts.layerDepths(i));
end

figure; hold on
clear cLine
cLine(1) = errorbar(cellfun(@mean,depthPeaks(1,:))', cellfun(@sem, depthPeaks(1,:))', '-o', 'linewidth', 2, 'color', ctrlColor);
cLine(2) = errorbar(cellfun(@mean,depthPeaks(2,:))', cellfun(@sem, depthPeaks(2,:))', '-o', 'linewidth', 2, 'color', KOcolor);
xlim([0, nrLayers+1]);
axis square;
title('Peak CSD response');
view(270,270);
ax = gca;
ax.XTick = 1 : 2 : nrLayers;
ax.XTickLabel = layerLabels(1 : 2 : nrLayers);
legend(opts.groups, 'location', 'northwest');
ax.YTick = -10 : 2 : 0;
niceFigure
for x = 1:length(cLine)
    cLine(x).LineWidth = 6;
    cLine(x).MarkerSize = 12;
end
ylabel('current source density [uV/mm^2]');
xlabel('depth [um]');
ylim([-10 0]);
grid on;