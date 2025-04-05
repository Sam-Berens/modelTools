function [pValue,fStat,df1,df2,rowSign] = runFCon(H,glme)
% runFCon.m
% Compute F contrast test statistics for a given contrast matrix.
% Sam Berens (s.berens@sussex.ac.uk)
% 05/04/2025
%
% [pValue, fStat, df1, df2, rowSign] = runFCon(H, glme)
%
% DESCRIPTION:
%   This function performs an F-test for a specified contrast matrix by
%   calculating the F statistic, the degrees of freedom for the contrast
%   (df1) and the model error (df2), and the associated p-value.
%   Additionally, it returns the sign of the contrast effect (i.e., the
%   sign of H*b) to indicate the direction of the effect.
%
% INPUTS:
%   H    - A contrast matrix specifying a linear combination of model 
%          coefficients.
%   glme - A fitted generalized linear model (or similar) with fields:
%          Coefficients, CoefficientCovariance, and DFE (degrees of freedom
%          for error).
%
% OUTPUTS:
%   pValue  - The p-value corresponding to the F-test.
%   fStat   - The computed F statistic.
%   df1     - The degrees of freedom associated with the contrast (i.e., 
%             the rank of H).
%   df2     - The degrees of freedom for the model error.
%   rowSign - The sign of the contrast effect, computed as sign(H*b).
%
% USAGE EXAMPLE:
%   [p, f, df1, df2, sign] = runFCon(H, glme);
%
% NOTES:
%   - Ensure that the contrast matrix H is properly specified and that glme
%     contains all necessary fields.
% 
b = glme.Coefficients.Estimate;
C = glme.CoefficientCovariance;
df1 = rank(H);
df2 = glme.DFE;
fStat = (H*b)'*(inv(H*C*H'))*(H*b) / df1;
pValue = 1-fcdf(fStat,df1,df2);
rowSign = sign(H*b);
return