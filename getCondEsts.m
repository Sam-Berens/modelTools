function [condXs] = getCondEsts(mdl,condXs)
% getCondEsts.m
% Compute condition estimates and 95% CIs for experimental conditions.
% Sam Berens (s.berens@sussex.ac.uk)
% 05/04/2025
%
% condXs = getCondEsts(mdl, condXs)
%
% DESCRIPTION:
%   This function calculates the predicted response and corresponding 95%
%   confidence interval for each experimental condition. It iterates over
%   each condition weight vector provided in the input table, computes the
%   predicted response using getXEUL, and then constructs a table with the
%   estimates and confidence intervals.
%
% INPUTS:
%   mdl    - A fitted model containing coefficient estimates, covariance
%            matrix, and degrees of freedom.
%   condXs - A table with row names indicating condition labels and a 
%            column ('X') containing the condition weight vectors.
%
% OUTPUT:
%   condXs - An updated table with two new columns:
%            est  - The predicted response for each condition.
%            Ci95 - A two-column matrix with the lower and upper bounds of 
%                   the 95% confidence interval.
%
% USAGE EXAMPLE:
%   condEsts = getCondEsts(mdl, condXs);
%
% NOTES:
%   - The function leverages getXEUL to handle interaction terms and apply
%     the model's inverse link function if necessary.
% 
condName = condXs.Row;
est = nan(size(condName));
Ci95 = nan(size(condName,1),2);
Xs = num2cell(condXs.X);
for iX = 1:size(Xs,1)
    cX = Xs(iX,:);
    [~,est(iX),Ci95(iX,1),Ci95(iX,2)] = getXEUL(mdl,cX{:});
end
condXs = table(est,Ci95,'RowNames',condName);
return