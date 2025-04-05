function [X] = addXterms(X,xTermStruct)
% addXterms.m
% Augment a design matrix with interaction terms.
% Sam Berens (s.berens@sussex.ac.uk)
% 05/04/2025
% 
% X = addXterms(X, xTermStruct) calculates interaction terms for the design
% matrix X based on the provided structure that defines which main effects
% are combined.
% 
% INPUTS:
%   X           - A design matrix with columns corresponding to main 
%                 effects.
%   xTermStruct - A cell array where each element specifies the indices of 
%                 main effects that contribute to a term. For a main
%                 effect, the element is simply its own column index; for 
%                 interactions, it is an array of indices.
% 
% OUTPUT:
%   X - The updated design matrix that includes both the original main 
%       effect columns and additional columns for the interaction terms, 
%       computed as the product of the corresponding main effect values.
% 
% USAGE EXAMPLE:
%   X = addXterms(X, xTermStruct);
% 
k = numel(xTermStruct);
nCols2Add = k - size(X,2);
X = [X,nan(size(X,1),nCols2Add)];
for iK = 1:k
    sX = X(:,xTermStruct{iK});
    X(:,iK) = prod(sX,2);
end
return