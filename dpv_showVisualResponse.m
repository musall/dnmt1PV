%% basic variables
opts.Location = {'S1'};
opts.makePlot = false;
opts.removeAvg = false;
opts.useTaper = false;
opts.groups = {'PV-Cre', 'PV-DNMT1'};
opts.baseDur = 0.5;
opts.postStim = 1;
opts.stepSize = 4;
opts.stimType = 1;

opts.brainRange = [-200, 1500]; %depth range for plotting
opts.verbose = false; %flag to supress some of the text outputs
opts.brainThresh = 1;
opts.groupColors = {([212, 212, 212]./255), ([255, 160, 64]./255), [0, 0, 1], [1, 0, 0]};
opts.layerRange = [500, 1000]; %range for computing CSD traces
opts.layerDepths = 25 : 75 : 1000; %depths of different layers in microns
ctrlColor = [212 212 212]./255;
KOcolor = [255 160 64]./255;
opts.savePath = 'F:\DNMT1_project\';
opts.reload = false;
opts.imecNr = '0';

%% run over groups
nrGroups = length(opts.groups);
usedRecs = cell(1, length(opts.groups));
failedRecs = cell(1, length(opts.groups));
allCSD = cell(1, length(opts.groups));
allLFP = cell(1, length(opts.groups));
allStimStd = cell(1, length(opts.groups));
for iGroups = 1 : nrGroups
    
    cPath = fullfile(opts.savePath, opts.groups{iGroups});
    recs = dir([cPath filesep '**' filesep 'csdResponse_imec' opts.imecNr '.mat']);
    
    allCSD{iGroups} = cell(1, length(recs));
    allLFP{iGroups} = cell(1, length(recs));
    allStimStd{iGroups} = cell(1, length(recs));
    for iRecs = 1 : length(recs)
        
        saveFile = fullfile(recs(iRecs).folder, recs(iRecs).name);
        cFile = load(saveFile, 'csdResp', 'lfpResp', 'meanStd', 'normDepth', 'lfMeta');
        
        allCSD{iGroups}{iRecs} = cFile.csdResp;
        allLFP{iGroups}{iRecs} = cFile.lfpResp;
        allStimStd{iGroups}{iRecs} = cFile.meanStd;

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
cRange = 1.5*1E6;
lfpRange = 0.7;
stimLabel = 'Visual';
traceRange = [-3 1]; %for CSD
lfMeta = cFile.lfMeta;

if opts.stimType == 1
    lfpRange = 0.4;
    csdRange = 2.5;
elseif opts.stimType == 3
    lfpRange = 150;
    csdRange = 1.3*1E6;
end

clear lines
depthLabel = {'Supragranular', 'Infragranular'};
respLabels = {'LFP absolute', 'CSD absolute', 'CSD rescale'};
h1 = figure('renderer' ,'painters');
h2 = figure('renderer' ,'painters');
xRange = [(opts.baseDur - 0.05), (opts.baseDur + 0.2)] .*lfMeta.sRateHz;
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
            depthRange = round(cFile.normDepth * 1E6);
            
        elseif x == 3
            mergeStim = cat(3,allCSD{iGroups}{:});
            mergeStim = mergeStim - nanmean(mergeStim(:, 1:baseDur, :), 2);
            cRange = prctile(abs(mergeStim(:)),99.7);
            depthRange = round(cFile.normDepth * 1E6);
        end
        
        % show image
        cImg = imagesc(nanmean(mergeStim, 3));
        ax = cImg.Parent;
        ax.XTick = 1 : lfMeta.sRateHz*0.05 : size(mergeStim,2);
        ax.XTickLabel = (((0 : lfMeta.sRateHz*0.05 : size(mergeStim,2)-1)./lfMeta.sRateHz) - opts.baseDur)*1000;
        xlim(xRange);
        caxis([-cRange cRange]);
        colormap(ax, parula(256));
        niceFigure;
        
        title([opts.Location{1} ' - ' opts.groups{iGroups} '; ' stimLabel ' stimulation; ' respLabels{x}]);
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
            niceFigure;
            
            title([opts.Location{1} ' - ' opts.groups{iGroups} '; ' stimLabel ' stimulation; ' respLabels{x}]);
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
                ylabel('CSD deflection');
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