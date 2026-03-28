function [pVal_group, tStat_group, fullmodel, modelCompare] = ...
    LME_compare_groups(dataIn1, dataIn2, animalID1, animalID2)

% dataIn1   : recordings from group 1
% dataIn2   : recordings from group 2
% animalID1 : animal ID for each entry in dataIn1
% animalID2 : animal ID for each entry in dataIn2

% convert to column vectors
dataIn1   = double(dataIn1(:));
dataIn2   = double(dataIn2(:));
animalID1 = animalID1(:);
animalID2 = animalID2(:);

% sanity checks
assert(length(dataIn1) == length(animalID1), ...
    'dataIn1 and animalID1 must have same length');
assert(length(dataIn2) == length(animalID2), ...
    'dataIn2 and animalID2 must have same length');

% make animal IDs unique across groups
animalID1 = categorical("G1_A" + string(animalID1));
animalID2 = categorical("G2_A" + string(animalID2));

% combine all observations
cData   = [dataIn1; dataIn2];
group   = categorical([repmat("G1", length(dataIn1), 1); ...
                       repmat("G2", length(dataIn2), 1)]);
animal  = [animalID1; animalID2];

% build table
tbl = table(cData, group, animal, ...
    'VariableNames', {'cData','group','animal'});

% fit model
fullmodel = fitlme(tbl, 'cData ~ group + (1|animal)');

% extract group effect
coefNames = fullmodel.CoefficientNames;
idx = find(strcmp(coefNames, 'group_G2') | strcmp(coefNames, 'group_G1'), 1, 'last');

% safer alternative: just use row 2 if only intercept + one group term
if isempty(idx)
    idx = 2;
end

pVal_group = fullmodel.Coefficients.pValue(idx);
tStat_group = fullmodel.Coefficients.tStat(idx);

if nargout > 3
    nullmodel = fitlme(tbl, 'cData ~ 1 + (1|animal)');
    modelCompare = compare(nullmodel, fullmodel, 'CheckNesting', true);
end
end
