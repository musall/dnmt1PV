%% basic variables
clear opts
opts.Location = {'V1'};
opts.makePlot = false;
opts.removeAvg = false;
opts.useTaper = true;
opts.groups = {'PV-Cre' , 'PV-DNMT1'};
opts.brainRange = 3000;
opts.hasGamma = '1';
opts.useRec = '1';
opts.savePath = 'F:\DNMT1_project\';
opts.reload = false;
opts.gid = 'PassiveStimulation';
opts.brainThresh = 2;
opts.imecNr = '0';

ctxRange = [0 1000];
scRange = [1800 2500];

%% run over groups
usedRecs = cell(1, length(opts.groups));
allSpecs = cell(1, length(opts.groups));
allFreqs = cell(1, length(opts.groups));
allInfo = cell(1, length(opts.groups));
allRejCnt = cell(1, length(opts.groups));
for iGroups = 1 : length(opts.groups)
    
    cPath = fullfile(opts.savePath, opts.groups{iGroups});
    recs = dir([cPath filesep '**' filesep 'lfpData_*imec' opts.imecNr '.mat']);
    allSpecs{iGroups} = cell(1, length(recs));
    allFreqs{iGroups} = cell(1, length(recs));
    allRejCnt{iGroups} = cell(1, length(recs));
    for iRecs = 1 : length(recs)
        
        saveFile = fullfile(recs(iRecs).folder, recs(iRecs).name);
        cFile = load(saveFile, 'allSpecs', 'plotFreqs', 'oscRejCnt');
        
        % keep info in larger array
        allSpecs{iGroups}{iRecs} = cFile.allSpecs;
        allFreqs{iGroups}{iRecs} = cFile.plotFreqs;
        allRejCnt{iGroups}{iRecs} = cFile.oscRejCnt;
        
    end
end


%% full range across the probes
lfpSize = size(allSpecs{1}{1},2);
depthRange = (1:lfpSize) ./ lfpSize * opts.brainRange;
ctxIdx =  depthRange >= ctxRange(1) & depthRange <= ctxRange(2);
scIdx =  depthRange >= scRange(1) & depthRange <= scRange(2);

h = figure('renderer' ,'painters');
Cnt = 0;
for iGroups = 1 : length(allSpecs)
    
    mergeSpecs = nanmean(cat(4,allSpecs{iGroups}{:}), 4);
    plotFreqs = nanmean(cat(4,allFreqs{iGroups}{1}), 4);

    for x = 1 : 2
        Cnt = Cnt + 1;
        subplot(length(allSpecs), 2, Cnt);
        
        if x == 1
            cMap = (log10(mergeSpecs(:,:,2)))';
        elseif x == 2
            cMap = ((mergeSpecs(:,:,2) - mergeSpecs(:,:,1)) ./ mergeSpecs(:,:,1))';
        end
        cImg = imagesc(cMap);
        ax = cImg.Parent;
        colormap(ax, colormap_blueblackred(256));
        maxFreq = plotFreqs(end);
        cRange = abs(prctile(cMap(:),97.5));
        caxis([-cRange cRange]);
        axis image;
        
        depthlines = [find(ctxIdx,1), find(ctxIdx,1,'last'), find(scIdx,1), find(scIdx,1,'last')];
        ax.YTick = depthlines;
        ax.YTickLabels = depthRange(depthlines);
        nhline(depthlines, 'w--')
        ylabel('Depth (um)');
        
        useFreqs = 1 : ceil(20 / mean(diff(plotFreqs))) : length(plotFreqs);
        ax.XTick = useFreqs;
        ax.XTickLabels = round(plotFreqs(useFreqs),2);
        nvline(useFreqs, 'w--')
        xlabel('Frequency (Hz)');
        xlim([0 find(plotFreqs > 100, 1)])
        title([opts.groups{iGroups} ' - removeAvg = ' num2str(opts.removeAvg) ' - useTaper = ' num2str(opts.useTaper)]);
        axis square
        niceFigure
    end
end

%% show change in oscillatory power
depthRange = (1:size(mergeSpecs,2)) ./ size(mergeSpecs,2) * opts.brainRange;
ctxIdx =  depthRange >= ctxRange(1) & depthRange <= ctxRange(2);
scIdx =  depthRange >= scRange(1) & depthRange <= scRange(2);

h3 = figure('name' , 'Traces');
useColors = {[212 212 212]./255, [255 160 64]./255};

locLabels = 'cortex';
Cnt = 0;
clear cLines
for iGroups = 1 : length(allSpecs)
    mergeSpecs = cat(4,allSpecs{iGroups}{:});
    plotFreqs = nanmean(cat(4,allFreqs{iGroups}{1}), 4);
    
    for x = 1 : 2
        if x == 1
            cData = (log10(mergeSpecs(:,:,2,:)));
            cRange = 0.65;
        elseif x == 2
            cData = ((mergeSpecs(:,:,2,:) - mergeSpecs(:,:,1,:)) ./ mergeSpecs(:,:,1,:));
            cRange = 0.65;
        end
        cData = squeeze(cData);
        
        % make image for cortex or SC only
        Cnt = Cnt + 1;
        figure(h3);
        subplot(2,2,iGroups+2);
        cIdx = ctxIdx;
        
        if x == 2
            cImg = imagesc(nanmean(cData(:,cIdx, :), 3)');
            ax = cImg.Parent;
            colormap(ax, colormap_blueblackred(256));
            maxFreq = plotFreqs(end);
            %             cRange = abs(prctile(cData(:),95));
            
            caxis([-cRange cRange]);
            axis square;
            stepIdx = 1 : 2 : sum(cIdx);
            ax.YTick = stepIdx;
            cDepth = depthRange(cIdx);
            ax.YTickLabels = cDepth(stepIdx)/1000;
            nhline(stepIdx, 'w--')
            ylabel('Depth [mm]');
            
            useFreqs = 1 : ceil(20 / mean(diff(plotFreqs))) : length(plotFreqs);
            ax.XTick = useFreqs;
            ax.XTickLabels = round(plotFreqs(useFreqs),2);
            nvline(useFreqs, 'w--')
            xlabel('Frequency [Hz]');
            xlim([0 find(plotFreqs > 100, 1)])
            title([opts.groups{iGroups} ' - removeAvg = ' num2str(opts.removeAvg) ' - useTaper = ' num2str(opts.useTaper)]);
            axis square
            niceFigure
        end
        
        % get traces
        layerIdx =  depthRange >= 500 & depthRange <= 900;
        cTrace = squeeze(nanmean(cData(:,layerIdx, :), 2)).*100;
        
        % show traces
        figure(h3);
        if x == 1
            subplot(2,2,1);
        elseif x == 2
            subplot(2,2,2);
        end
        
        hold on;
        cLines(iGroups) = stdshade((cTrace)', 0.1, useColors{iGroups}, plotFreqs);
        axis square
        niceFigure
        xlabel('frequency (Hz)');
        title(locLabels);
        xlim([0 115]);
        
        if x == 1
            ylabel('power spectral density [uV^2 / Hz]');
        elseif x == 2
            ylabel('stimulus-induced difference [%]');
            nhline(0, 'k--');
        end
        if iGroups == 2
            legend(cLines, opts.groups);
        end
    end
end