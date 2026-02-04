function [allStim, allStimData] = pC_getSessionData(lfpDir, relevantVarNames)
% function to get all settings files from the 'vStim' stimulation software 
% and create a single array for all trials. The output array allStim contains 
% all settings for the relevant variables, defined by 'relevantVarNames'.
% allStimdata also gives out all raw StimData variables from all files.
% example:
% relevantVarNames = ["PulseCount", "PulseDur",  'PulseGap', 'OptoRamp', 'StimType', 'RedPowerOne', 'RedPowerTwo', 'BluePowerOne', 'BluePowerTwo'];

settingsDir = fileparts(fileparts(lfpDir));
settingsFiles = dir([settingsDir filesep '*settings.mat']);
allStim = [];
allStimData = cell(1, length(settingsFiles));
for iFiles = 1  : length(settingsFiles)

    load([settingsDir filesep settingsFiles(iFiles).name])
    disp(settingsFiles(iFiles).name)
    StimData.nrTrials = length(StimData.TimeStamps);
    
    if StimData.nrTrials < size(StimData.VarVals, 2)
       warning('Not enough exectued trials in StimData !!'); 
       warning(sprintf('Found %.0f trials instead of preset %.0f trials', StimData.nrTrials, size(StimData.VarVals, 2))); 
    end
    
    StimData.VarVals = StimData.VarVals(:,1:StimData.nrTrials); %only include completed trials
    allStimData{iFiles} = StimData;

    
    cStim = NaN(length(relevantVarNames), StimData.nrTrials);
     for i = 1 : length(relevantVarNames)
        if sum(ismember(StimData.VarNames,relevantVarNames(i)))>0
        cStim(i,:) =  StimData.VarVals(ismember(StimData.VarNames,relevantVarNames(i)),1:StimData.nrTrials);
        else
        disp('VarName '+relevantVarNames(i)+' does not exist')
        cStim(i,:) = NaN;
        end
     end
    
    if isempty(allStim)
        allStim = cStim;
    else
        allStim = [allStim, cStim];
    end 
end