function analogDat = pC_extractMedianAnalogChannel(filename, numChans, eventTimes, winIdx)
% usage: analogDat = pC_extractMedianAnalogChannel(filename, numChans, eventTimes, winIdx)
% function to extract analog channel from binary data file.
% fileName is the name of the binary file, numChans says how many channels
% there are in total, chanIdx identifies the channel of interest, eventTimes
% are the timpoints around which data should be returned and winIdx is a
% timeIndex that defines how much data is returned around each event in
% eventTimes.

maxReadSize = floor(1e9 / numChans); %maximum number of bytes to be read in one step

if ~exist('eventTimes', 'var')
    eventTimes = []; %no events, get the whole trace
end

if ~exist('winIdx', 'var')
    winIdx = []; %no events, get the whole trace
end

% get datafile and nr of samples
fprintf(1,' loading median over all analog channels %s\n', filename);
d = dir(filename);
nSamp = d.bytes/2/numChans;
if isempty(eventTimes)
    analogDat = zeros(1, nSamp, 'int16');
else
    analogDat = zeros(length(winIdx), length(eventTimes), 'int16');
end

% skip over the first samples of the other channels
fid = fopen(filename, 'r');

if isempty(eventTimes)
    % whole trace
    nBatch = floor(nSamp/maxReadSize);
    for b = 1:nBatch
        tempDat = fread(fid, [1, maxReadSize*numChans], 'int16=>int16'); % load data from all channels
        tempDat = reshape(tempDat, numChans, maxReadSize);
        analogDat(1, (b-1)*maxReadSize+1:b*maxReadSize) = median(tempDat(1:end-1,:), 1);
    end
    tempDat = fread(fid, [1, Inf], 'int16=>int16'); % load data from all channels
    tempDat = reshape(tempDat, numChans, []);
    analogDat(1, nBatch*maxReadSize+1:end) = median(tempDat(1:end-1,:), 1);
    
else
    metaName = strrep(filename, '.bin', '.meta');
    metaVars = readSpikeGLXmeta(metaName);
    
    % get specific times
    Cnt = 0;
    for x = 1 : length(eventTimes)
        
        cOn = eventTimes(x) * metaVars.sRateHz + winIdx(1) ; %onset for current window
        
        fseek(fid, (numChans*2)* round(cOn - Cnt), 'cof'); %move pointer to current window onset
        tempDat = fread(fid, [1, length(winIdx)*numChans], 'int16=>int16'); %read data according to window size
        tempDat = reshape(tempDat, numChans, length(winIdx));
        analogDat(:,x) = median(tempDat(1:end-1,:), 1);
        
        Cnt = cOn + length(winIdx);
        
    end
end
fclose(fid);


