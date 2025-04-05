function [X,est,low,upp] = getXEUL(mdl,varargin)
% getXEUL.m
% Compute predicted responses and 95% CIs from a design matrix.
% Sam Berens (s.berens@sussex.ac.uk)
% 05/04/2025
%
% [X, est, low, upp] = getXEUL(mdl, X1 ,X2, ..., Xn)
%
% DESCRIPTION:
%   This function builds a design matrix from the supplied predictor values
%   (either full model terms or only simple effects) and augments it with
%   interaction terms via addXterms. It then computes the estimated
%   response values and their 95% confidence intervals using the model's
%   coefficients and covariance matrix. If the model specifies a link
%   function, the inverse link is applied to the outputs.
%
% INPUTS:
%   mdl       - A fitted statistical model (e.g., LinearModel) containing 
%               coefficient estimates, covariance, degrees of freedom, and, 
%               optionally, a link function.
%   predictor - A variable number of inputs representing predictor values.
%               Provide either all model terms or only the simple
%               (non-interaction) effects.
%
% OUTPUTS:
%   X   - The design matrix constructed from the inputs, including added 
%         interaction terms.
%   est - The estimated predicted responses (after applying the inverse
%         link function if applicable).
%   low - The lower bounds of the 95% confidence intervals.
%   upp - The upper bounds of the 95% confidence intervals.
%
% USAGE EXAMPLE:
%   [X, est, low, upp] = getXEUL(mdl, in1, in2, in3, ...);
%
% NOTES:
%   - An error is raised if the number or dimensions of the provided
%     predictors are not compatible with the model.
%
%% Extract some things
d = mdl.NumCoefficients;
dIn = numel(varargin);
coefName = mdl.CoefficientNames';
numSimple = sum(cellfun(@(s)~contains(s,':'),coefName));

%% Check the inputs
if (dIn < d) && (dIn ~= numSimple)
    error(['You must provide either all the model terms ',...
        'or only the simple effects.']);
end

%% Start filling up X
n = max(cellfun(@(v)numel(v),varargin));
X = zeros(n,d);
for ii = 1:dIn
    if numel(varargin{ii}) == 1
        X(:,ii) = varargin{ii};
    else
        if numel(varargin{ii}) == n
            X(:,ii) = varargin{ii}(:);
        else
            error('Incompatible inputs dimensions.');
        end
    end
end

%% Add any remaining interaction terms
xTermStruct = getXtermStruct(coefName);
X = addXterms(X,xTermStruct);

%% Maths time
b = mdl.Coefficients.Estimate;
C = mdl.CoefficientCovariance;
df = mdl.DFE;
est = X*b;
Ci95 = sqrt(diag(X*C*X'))*[tinv(0.025,df),tinv(0.975,df)] + est;

if isprop(mdl,'Link')
    invLink = @(l) mdl.Link.Inverse(l);
else
    invLink = @(l) l;
end
est = invLink(est);
Ci95 = invLink(Ci95);

%% Reshape the outputs
est = reshape(est,[n,1]);
low = reshape(Ci95(:,1),[n,1]);
upp = reshape(Ci95(:,2),[n,1]);
return