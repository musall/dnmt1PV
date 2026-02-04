%% define experiments to use
recLabels = {'Animal','ExperimentDate','Location','Expertise','Folder','Probe','Background','Path','Ramppower','useRec','hasGamma'};
recInfo{1} = {'SOM2554','26/01/2022','V1','Naive','SOM2554_20221026','0','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'SOM2554','28/01/2022','V1','Naive','SOM2554_20221028','0','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'SOM2554','29/01/2022','V1','Naive','SOM2554_20221029','0','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'SOM2563','26/01/2022','V1','Naive','SOM2563_20221026','0','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'SOM2563','27/01/2022','V1','Naive','SOM2563_20221027','0','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'SOM2563','28/01/2022','V1','Naive','SOM2563_20221028','0','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control';'SOM2563','29/01/2022','V1','Naive','SOM2563_20221029','0','SOM-Cre','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','Control'};
recInfo{2} = {'SOM2626','09/02/2023','V1','Naive','SOM2626_20230209','0','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','SOM-DNMT1';'SOM2626','10/02/2023','V1','Naive','SOM2626_20230210','0','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','SOM-DNMT1';'SOM2626','11/02/2023','V1','Naive','SOM2626_20230211','0','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','SOM-DNMT1';'SOM2627','24/01/2023','V1','Naive','SOM2627_20230124','0','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','SOM-DNMT1';'SOM2627','25/01/2023','V1','Naive','SOM2627_20230125','0','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','SOM-DNMT1';'SOM2627','26/01/2023','V1','Naive','SOM2627_20230126','0','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','SOM-DNMT1';'SOM2627','27/01/2023','V1','Naive','SOM2627_20230127','0','SOM-DNMT1','\\naskampa.kampa-10g\lts\invivo_ephys\Neuropixels\','high','1','1','SOM-DNMT1'};

%% basic variables
clear opts
opts.Location = {'V1'};
opts.makePlot = false;
opts.removeAvg = false;
opts.useTaper = true;
opts.groups = {'SOM-Cre' , 'SOM-DNMT1'};
opts.brainRange = 3000;
opts.hasGamma = '1';
opts.useRec = '1';
opts.savePath = 'F:\DNMT1_project\Ephys_Data\';
opts.reload = false;
opts.gid = 'PassiveStimulation';
opts.brainThresh = 2;

if strcmpi(opts.Location, 'V1')
    opts.imecNr = '0';
elseif strcmpi(opts.Location, 'S1')
    opts.imecNr = '1';
end
ctxRange = [0 1000];
scRange = [1800 2500];

%% run over groups
usedRecs = cell(1, length(opts.groups));
failedRecs = cell(1, length(opts.groups));
allSpecs = cell(1, length(opts.groups));
allFreqs = cell(1, length(opts.groups));
allInfo = cell(1, length(opts.groups));
for iGroups = 1 : length(opts.groups)
    
    allInfo{iGroups} = recInfo{iGroups};
    
    recIdx = strcmpi(recLabels, 'Folder');
    pathIdx = strcmpi(recLabels, 'Path');
    for iRecs = 1 : size(recInfo{iGroups}, 1)
        try
            
            [cSpecs, cFreqs, oscRejCnt] = pC_checkVisualGamma_SM(recInfo{iGroups}{iRecs,recIdx}, recInfo{iGroups}{iRecs,pathIdx}, opts); drawnow;
            disp(oscRejCnt)
            
            % keep info in larger array
            if isempty(allSpecs{iGroups})
                allSpecs{iGroups} = {cSpecs};
                allFreqs{iGroups} = {cFreqs};
                usedRecs{iGroups} = {recInfo{iGroups}{iRecs,recIdx}};
                allRejCnt{iGroups} = {oscRejCnt};
            else
                allSpecs{iGroups} = [allSpecs{iGroups}; cSpecs];
                allFreqs{iGroups} = [allFreqs{iGroups}; cFreqs];
                allRejCnt{iGroups} = [allRejCnt{iGroups}; oscRejCnt];
                usedRecs{iGroups} = [usedRecs{iGroups}; recInfo{iGroups}{iRecs,recIdx}];
            end
        catch
            failedRecs{iGroups} = [failedRecs{iGroups}; recInfo{iGroups}{iRecs,recIdx}];
        end
    end
end


%% work in progress
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
        nhline(depthlines, 'w--');
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

%% show plots for specific depth range
gammaChange = cell(length(allSpecs), 2, 2);
depthRange = (1:size(mergeSpecs,2)) ./ size(mergeSpecs,2) * opts.brainRange;
ctxIdx =  depthRange >= ctxRange(1) & depthRange <= ctxRange(2);
scIdx =  depthRange >= scRange(1) & depthRange <= scRange(2);

% cortex and SC figures
h3 = figure('name' , 'Traces');
useColors = {[212 212 212]./255, [255 160 64]./255};
gammaRng = [60 70];

locLabels = 'cortex';
Cnt = 0;
clear cLines
for iGroups = 1 : length(allSpecs)
    mergeSpecs = cat(4,allSpecs{iGroups}{:});
    plotFreqs = nanmean(cat(4,allFreqs{iGroups}{1}), 4);
    
    for x = 1 : 2
        if x == 1
            cData = (log10(mergeSpecs(:,:,2,:)));
            cRange = 1;
        elseif x == 2
            cData = ((mergeSpecs(:,:,2,:) - mergeSpecs(:,:,1,:)) ./ mergeSpecs(:,:,1,:));
            cRange = 1;
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
        layerIdx =  depthRange >= 400 & depthRange <= 700;
        cTrace = squeeze(nanmean(cData(:,layerIdx, :), 2)).*100;
        
        % for statistics
        freqIdx = plotFreqs > gammaRng(1) & plotFreqs < gammaRng(2);
        gammaChange{iGroups, 1, x} = mean(cTrace(freqIdx, :), 1);
        
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

%% give statistics for gamma power
disp('====================')
x = 1; %1 for absolute, 2 for stimulus-induced difference
fprintf('Absolute gamma change cortex.');
disp('Absolute gamma power difference between groups');
for iGroups = 1 : 2
    fprintf('%s: Gamma power %.2f %c %.2f uV^2\n'  , opts.groups{iGroups}, mean(gammaChange{iGroups,1,x}), char(177), sem(gammaChange{iGroups,1,x}));
end
fprintf('pVal ranksum test: %f\n', ranksum(gammaChange{1,1,x}, gammaChange{2,1,x}))

disp('====================')
x = 2; %1 for absolute, 2 for stimulus-induced difference
disp('Relative gamma change cortex.');
disp('Stimulus-induced gamma change between groups');
for iGroups = 1 : 2
    fprintf('%s: Gamma power %.2f %c %.2f percent\n'  , opts.groups{iGroups}, mean(gammaChange{iGroups,1,x}), char(177), sem(gammaChange{iGroups,1,x}));
end
fprintf('pVal ranksum test: %f\n', ranksum(gammaChange{1,1,x}, gammaChange{2,1,x}))
disp('====================')
