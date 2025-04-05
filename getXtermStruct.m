function [xTermStruct] = getXtermStruct(coefName)
% getXtermStruct.m
% Parse coefficient names to determine interaction structure.
% Sam Berens (s.berens@sussex.ac.uk)
% 05/04/2025
% 
% xTermStruct = getXtermStruct(coefName) processes a cell array of
% coefficient names to identify which coefficients correspond to main
% effects and which represent interactions.
% 
% INPUT:
%   coefName - A cell array of strings representing the names of 
%              coefficients. Interaction terms are indicated by
%              colon-separated names.
% 
% OUTPUT:
%   xTermStruct - A cell array where each element corresponds to a
%                 coefficient. For a main effect, the element is its own
%                 index. For an interaction term, it contains the indices
%                 of the main effects that form the interaction.
% 
% USAGE EXAMPLE:
%   coefName = {'(Intercept)', 'A', 'B', 'A:B'};
%   xTermStruct = getXtermStruct(coefName);
% 
isXterm = cellfun(@(s)contains(s,':'),coefName);
xTermStruct = cell(size(isXterm));
for iInter = 1:1:numel(xTermStruct)
    if isXterm(iInter)
        [~,xTermStruct{iInter}] = ...
            ismember(split(coefName{iInter},':'),coefName(~isXterm));
    else
        xTermStruct{iInter} = iInter;
    end
end
return