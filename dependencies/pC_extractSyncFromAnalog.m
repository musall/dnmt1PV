function syncDat = pC_extractSyncFromAnalog(folder, numChans, syncChanIndex, datFile, loadRaw, trigThresh)
% usage: syncDat = pC_extractSyncFromAnalog(folder, numChans, syncChanIndex, datFile, loadRaw, trigThresh)
% extract digital sync signal from analog trace in a int16 binary file.
% This can be used to get the sync channel from the 'lf' binary file from spikeGLX.
% To save time digital traces are saved in a sync file and will be loaded
% from there if the file exists. Use loadRaw = true to force loading from
% binary data file. trigThresh is the threshold with which the analog trace
% is converted to digital (default is 1).
% extraChanIndices are 1-indexed

if ~exist('datFile', 'var') || isempty(datFile)
    datFile = 'lf'; %load LFP sync by default
end

if ~exist('trigThresh', 'var') || isempty(loadRaw)
    trigThresh = 1; %threshold to detect logical events in standard deviations
end

if ~exist('loadRaw', 'var') || isempty(loadRaw)
    loadRaw = false; %load sync file by default
end

dataFiles = dir(fullfile(folder,['*.' datFile '.bin']));
syncDat = cell(length(dataFiles), 3);
for d = 1:length(dataFiles)
    filename = fullfile(folder, dataFiles(d).name);
    syncDat(d,:) = pC_syncFromAnalogBin(filename, numChans, syncChanIndex, loadRaw, trigThresh);
end