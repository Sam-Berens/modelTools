function [conTable] = runHCons(mdl,conTable)
% runHCons.m
% Apply hypothesis contrasts to a fitted model and compute test statistics.
% Sam Berens (s.berens@sussex.ac.uk)
% 05/04/2025
%
% conTable = runHCons(mdl, conTable)
%
% DESCRIPTION:
%   This function processes a table of contrast matrices (one per
%   experimental effect) and applies them to the fitted model to compute
%   corresponding test statistics and p-values. For each contrast, it
%   determines whether a t-test (for single-degree-of-freedom contrasts) or
%   an F-test is appropriate, and augments the table with the test type,
%   test statistic, and p-value.
%
% INPUTS:
%   mdl      - A fitted model object containing coefficients, covariance,
%              and degrees of freedom.
%   conTable - A table where each row contains a contrast matrix (H) for a
%              specific effect. The row names should denote the effect
%              names.
%
% OUTPUT:
%   conTable - The input table augmented with the following columns:
%                testType  - A string describing the test (e.g., 't(df)' or
%                            'F(df1,df2)').
%                testStat  - The computed t-value or F statistic.
%                pValue    - The p-value of the hypothesis test.
%
% USAGE EXAMPLE:
%   conTable = runHCons(mdl, conTable);
%
% NOTES:
%   - This function assumes that each contrast in conTable is represented 
%     as a single variable.
%   - For contrasts with one degree of freedom, the test statistic is 
%     converted to a t-value.
% 
conName = conTable.Row;
testType = cell(size(conName));
testStat = nan(size(conName));
pValue = nan(size(conName));
Hs = conTable.Variables;
if size(Hs,2) > 1
    error(['The conTable input must have a single variable ',...
        'corresponding to the desired contrasts.']);
end
for iCon = 1:numel(conName)
    [pValue(iCon,1),fStat,df1,df2,rowSign] = ...
        runFCon(Hs{iCon},mdl);
    if df1 == 1
        testType{iCon,1} = sprintf('t(%i)',df2);
        testStat(iCon,1) = sqrt(fStat) * rowSign;
    else
        testType{iCon,1} = sprintf('F(%i,%i)',df1,df2);
        testStat(iCon,1) = fStat;
    end
end
conTable = [conTable,table(testType,testStat,pValue)];
return