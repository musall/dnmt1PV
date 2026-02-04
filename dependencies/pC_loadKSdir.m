function [sp, trigDat, params] = pC_loadKSdir(ksDir, varargin)
% function to load spike data from kilosort.
% the resulting struct sp contains all information from the spikesorting.
% sp.cids are the cluster ids which are used to identify different
% clusters. Their labels from the decoder or human labels are given in 
% sp.clusterLabels. 

%% check params for specific settings
if ~isempty(varargin)
    params = varargin{1};
else
    params = [];
end

if ~isfield(params, 'localPath')
    params.localPath = [];
end

if ~isfield(params, 'loadRaw')
    params.loadRaw = false;
end

if ~isfield(params, 'excludeNoise')
    params.excludeNoise = false;
end

if ~isfield(params, 'loadPCs')
    params.loadPCs = true;
end

if ~isfield(params, 'niCorrection')
    params.niCorrection = true;
end

%make sure the path does not end with path separator
if ksDir(end) == filesep
    ksDir(end) = [];
end

% make savepath
lfpDir = fileparts(fileparts(ksDir));
[nidaqDir, lfpFolder] = fileparts(lfpDir);
[recFolder, recName] = fileparts(nidaqDir);
[~, baseFolder] = fileparts(recFolder);
saveName = [recName '_npxData.mat']; %filename for saved data

if isempty(params.localPath)
    params.saveFile = [lfpDir filesep saveName];
else
    params.saveFile = [params.localPath filesep baseFolder filesep recName filesep lfpFolder filesep saveName];
end

%% load data
checkFile = exist(params.saveFile, 'file') > 0;

if ~checkFile || params.loadRaw
    %% load raw spike data
    sp = pC_loadParamsPy(fullfile(ksDir, 'params.py'));
    ss = readNPY(fullfile(ksDir, 'spike_times.npy'));
    
    % get sampling rate from meta file
    apMetaName = dir(fullfile(lfpDir,'*.ap.meta'));
    apMeta = pC_readSpikeGLXmeta(fullfile(apMetaName.folder,apMetaName.name));
    st = double(ss)/apMeta.sRateHz;
    
    spikeTemplates = readNPY(fullfile(ksDir, 'spike_templates.npy')); % note: zero-indexed
    
    if exist(fullfile(ksDir, 'spike_clusters.npy'), 'file')
        clu = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
    else
        clu = spikeTemplates;
    end

    tempScalingAmps = readNPY(fullfile(ksDir, 'amplitudes.npy'));
    if params.loadPCs
        pcFeat = readNPY(fullfile(ksDir,'pc_features.npy')); % nSpikes x nFeatures x nLocalChannels
        pcFeatInd = readNPY(fullfile(ksDir,'pc_feature_ind.npy')); % nTemplates x nLocalChannels
    else
        pcFeat = [];
        pcFeatInd = [];
    end
    
    cgsFile = '';
    if exist(fullfile(ksDir, 'decoder_output_dataframe.csv'), 'file')
        cgsFile = fullfile(ksDir, 'decoder_output_dataframe.csv');
    elseif exist(fullfile(lfpDir,'metrics', 'decoder_output_dataframe.csv'), 'file') %check alternative location if not directly in folder
        cgsFile = fullfile(lfpDir,'metrics', 'decoder_output_dataframe.csv');
    end
    
    if exist(fullfile(ksDir, 'sorting-curation.tsv'), 'file')
        cgsFile = fullfile(ksDir, 'sorting-curation.tsv');
    end

    % check for decoder output
    metrics = [];
    if ~isempty(cgsFile)

        % check for manual labels from sortingview
        [~, cFile] = fileparts(cgsFile); %get filename for selected curation file
        labelFile = 'sorting-curation.tsv'; %sortingview filename
        
        cids = [];
        sp.clusterFile = [];
        if strcmpi(cFile, labelFile)
            disp('Found manual sortingview labels for recording. Using them to provide cluster labels.')
            [cids, clusterLabels] = readClusterGroupsTSV(cgsFile, 'sortingview');
            sp.clusterFile = labelFile;
        end
        
        % check decoder output and use if no manual labels are present.
        % Also keep quality metrics.
        labelFile = 'decoder_output_dataframe'; %decoder filename
        if strcmpi(cFile, labelFile)
            
            % get decoder output
            metrics = readtable(cgsFile);
            if isempty(cids)
                disp('Found no manual labels for recording. Using decoder output to identify noise clusters.')
                clusterLabels = zeros(1, length(metrics.decoder_label), 'single'); %non-mua or sua clusters should be noise
                clusterLabels(cellfun(@(x)strcmpi(x,'mua'),metrics.decoder_label')) = 1; %mua cluster
                clusterLabels(cellfun(@(x)strcmpi(x,'sua'),metrics.decoder_label')) = 2; %sua cluster
                cids = metrics.cluster_id;
                sp.clusterFile = labelFile;
            end
            
            % set unlabeled clusters to 0
            useIdx = ismember(metrics.cluster_id, cids)';
            tempLabels = clusterLabels;
            clusterLabels = zeros(1, length(metrics.cluster_id), 'single'); %non-mua or sua clusters should be noise
            clusterLabels(useIdx) = tempLabels;
            cids = metrics.cluster_id;
        end

%         % when removing noise, just remove them from the cluster id
%         if params.excludeNoise
%             disp('Removing noise clusters based on noise labels');
%             useIdx = clusterLabels ~= 0;
%             cids = cids(useIdx);
%             clusterLabels = clusterLabels(useIdx);
%         end
% 
%         % only use cluster ids for which we have labels
%         useIdx = ismember(clu, cids);
%         spikeTemplates = spikeTemplates(useIdx);
%         tempScalingAmps = tempScalingAmps(useIdx);
%         clu = clu(useIdx);
%         if params.loadPCs
%             pcFeat = pcFeat(useIdx, :,:);
%         end
%         metrics = metrics(ismember(metrics.cluster_id, cids),:);
        
    else
        disp('!!! No decoder output or manual labels found. Clusters are unclassified. !!!')
        clu = spikeTemplates;
        amps = readtable(fullfile(ksDir, 'cluster_Amplitude.tsv'), 'FileType', 'text', 'Delimiter', '\t');
        cids = amps.cluster_id;
        clusterLabels = 3*ones(size(cids));
    end

    coords = readNPY(fullfile(ksDir, 'channel_positions.npy'));
    ycoords = coords(:,2); xcoords = coords(:,1);
    temps = readNPY(fullfile(ksDir, 'templates.npy'));
    winv = readNPY(fullfile(ksDir, 'whitening_mat_inv.npy'));

    % compute spikedepths
    spikeFeatInd = pcFeatInd(spikeTemplates+1,:);
    spikeFeatYcoords = ycoords(spikeFeatInd+1); % 2D matrix of size #spikes x 12
    redPcFeat = squeeze(pcFeat(:,1,:)); % take first PC only
    redPcFeat(redPcFeat<0) = 0; % some entries are negative, but we don't really want to push the CoM away from there.
    spikeDepths = sum(spikeFeatYcoords.*redPcFeat.^2,2)./sum(redPcFeat.^2,2); % center of mass is sum(coords.*features)/sum(features)

    %% make correction against NI DAQ
    metaName = dir(fullfile(nidaqDir,'*.nidq.meta'));
    niMeta = pC_readSpikeGLXmeta(fullfile(metaName.folder,metaName.name));
    trigDat = pC_extractDigitalChannel(nidaqDir, niMeta.nChans, niMeta.nChans, 'nidq'); %get digital channel from ni-daq. usually the lat channel

    if params.niCorrection
        lfpSyncDat = pC_extractSyncFromAnalog(lfpDir, 385, 385, 'lf'); %extract digital signal from analog trace

        % make correction
        [~,b] = makeCorrection(trigDat{1}{1}{2}, lfpSyncDat{2}, false); % adjust time base between spike and event times
        st = applyCorrection(st, b);
        sp.b = b; %save correction so its easy to check if alignmend was done
    else
        sp.b = NaN;
    end

    %% create struct and save data
    sp.st = st;
    sp.spikeTemplates = spikeTemplates;
    sp.clu = clu;
    sp.tempScalingAmps = tempScalingAmps;
    sp.cids = cids;
    sp.xcoords = xcoords;
    sp.ycoords = ycoords;
    sp.winv = winv;
    sp.pcFeat = pcFeat;
    sp.spikeDepths = spikeDepths;
    sp.clusterLabels = clusterLabels; %these are either manual or decoder labels. Will use the same name here for later processing.
    sp.metrics = metrics;
    sp.temps = temps(sp.cids+1, :, :);
    sp.pcFeatInd = pcFeatInd(sp.cids+1, :);

    if ~exist(fileparts(params.saveFile), 'dir')
        mkdir(fileparts(params.saveFile))
    end
    save(params.saveFile, 'sp', 'trigDat', 'params','-v7.3'); % save raw data

else
    % load data
    disp('Loading local file');
    load(params.saveFile, 'sp', 'trigDat', 'params'); % load data from file
    params.saveFile = [params.localPath filesep baseFolder filesep recName filesep lfpFolder filesep saveName];
    params.savePath =  fileparts(params.saveFile);

end

