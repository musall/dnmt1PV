function [qMetric, newUnitType, forGUI, opts] = pC_computeBombcellMetrics_fromTemplate_SM(opts)
% function to compute quality metrics and units with bombcell code
% based on the gettingStarted.mlx file from the bombcell repository - SM 20250501

cRec = [fullfile(opts.recPath, opts.recName) filesep opts.recName '_g0' filesep  opts.recName '_g0_imec' opts.probeNr];
ephysMetaDir = dir(fullfile(cRec, '*.ap.meta'));
ephysMetaDir = ephysMetaDir(1);
ephysRawFile = "NaN"; % path to your raw .bin or .dat data
ephysKilosortPath = fullfile(cRec, 'spikeinterface_KS2_5_output', 'sorter_output');
savePath = [cRec filesep 'bc_qMetrics']; % where you want to save the quality metrics 

kilosortVersion = 2.5; % if using kilosort4, you need to have this value kilosertVersion=4. Otherwise it does not matter. 
gain_to_uV = NaN; % use this if you are not using spikeGLX or openEphys to record your data. this value, when mulitplied by your raw data should convert it to  microvolts. 

% load data
[spikeTimes_samples, spikeClusters, templateWaveforms, templateAmplitudes, pcFeatures, ...
    pcFeatureIdx, channelPositions] = bc.load.loadEphysData(ephysKilosortPath, savePath);

% get params
param = bc.qm.qualityParamValues(ephysMetaDir, ephysRawFile, ephysKilosortPath, gain_to_uV, kilosortVersion);
meta = bc.dependencies.SGLX_readMeta.ReadMeta(ephysMetaDir.name, ephysMetaDir.folder);
[AP, ~, SY] = bc.dependencies.SGLX_readMeta.ChannelCountsIM(meta);
param.nChannels = AP + SY;
param.nSyncChannels = SY;

%% compute quality metrics
[qMetric, unitType] = bc.qm.runAllQualityMetrics(param, spikeTimes_samples, spikeClusters, ...
    templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, savePath);

forGUI = []; %should usually exist but just to make sure
if ~exist('forGUI', 'var') || ~isempty(dir([savePath, filesep, 'templates.qualityMetricDetailsforGUI.mat']))
    load([savePath, filesep, 'templates.qualityMetricDetailsforGUI.mat'], 'forGUI')
end

%% optional: Inspect
% After running quality metrics, espacially the first few times, it's a good idea to inspect your data and the quality metrics using the built-in GUI. Use your keyboard to navigate the GUI: 
% left/right arrow : toggle between units 
% u  : brings up a input dialog to enter the unit you want to go to
% g  : go to next good unit 
% m : go to next multi-unit 
% n  : go to next noise unit
% a  : go to next non-somatic unit ("a" is for axonal)
% up/down arrow : toggle between time chunks in the raw data

% bc.load.loadMetricsForGUI;
% 
% unitQualityGuiHandle = bc.viz.unitQualityGUI_synced(memMapData, ephysData, qMetric, forGUI, rawWaveforms, ...
%     param, probeLocation, unitType, loadRawTraces);

%% to adjust thresholds, edit qualityParamValues and run again
% param = bc.qm.qualityParamValues(ephysMetaDir, ephysRawFile, ephysKilosortPath, gain_to_uV, kilosortVersion);
% unitType = bc.qm.getQualityUnitType(param, qMetric, savePath);
% bc.qm.plotGlobalQualityMetric(qMetric, param, unitType, uniqueTemplates, forGUI.tempWv);

%% save as tsv file
% goodUnits = unitType == 1;
% muaUnits = unitType == 2;
% noiseUnits = unitType == 0;
% nonSomaticUnits = unitType == 3; 

% adjust to sortingview convention
newUnitType = repmat("noise", size(unitType,1), 1);
newUnitType(unitType == 2) = "MUA";
newUnitType(unitType == 1) = "SUA";
newUnitType(unitType == 3) = "nonSoma";

% (for use with another language: output a .tsv file of labels. You can then simply load this)
label_table = array2table([qMetric.phy_clusterID, newUnitType], 'VariableNames', {'cluster_id', 'group'});
writetable(label_table,[savePath filesep 'bc_unit_labels.tsv'],'FileType', 'text','Delimiter','\t');  
    
    