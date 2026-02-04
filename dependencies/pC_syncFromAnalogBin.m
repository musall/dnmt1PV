function syncDat = pC_syncFromAnalogBin(filename, numChans, syncChanIndex, loadRaw, trigThresh)
% usage: syncDat = pC_syncFromAnalogBin(filename, numChans, syncChanIndex, loadRaw, trigThresh)
% extract digital sync channel from a binary file, like the LFP.bin file from spikeGLX. 
% extraChanIndices are 1-indexed

if ~exist('loadRaw', 'var') || isempty(loadRaw)
    loadRaw = false; %use this if triggers should be loaded from raw file
end

if ~exist('trigThresh', 'var') || isempty(loadRaw)
    trigThresh = 1; %threshold to detect logical events in standard deviations
end

[folder,fn] = fileparts(filename);
syncFname =  fullfile(folder, [fn '_trig' num2str(syncChanIndex) '.mat']);

% check if sync file exist and get data from there if not empty
checker = false;
if exist(syncFname, 'file')
    syncSize = dir(syncFname);
    
    if syncSize.bytes > 0
        checker = true;
    end
end

if checker && ~loadRaw
    % get syncData from dedicated file
    load(syncFname, 'syncDat');
    
else
    
    % get analog data
    syncDat = pC_extractAnalogChannel(filename, numChans, syncChanIndex);
    
    % translate digital channel to events
    syncDat = syncDat > std(single(syncDat)) * trigThresh;
    metaName = strrep(filename, '.bin', '.meta');
    metaVars = readSpikeGLXmeta(metaName);
    syncDat = pC_digitalToTimestamp(syncDat, metaVars.sRateHz); %get sync times from digital channel

    save(syncFname, 'syncDat'); %save for later
end
