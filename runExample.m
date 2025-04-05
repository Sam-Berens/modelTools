function [] = runExample()
% runExample.m
% Execute an example analysis pipeline using the included function.
% Sam Berens (s.berens@sussex.ac.uk)
% 05/04/2025
%
% runExample()
%
% DESCRIPTION:
%   This function provides a complete, integrated example of the included
%   functions. It performs the following steps:
%     1. Generates an example model and associated factor structure using a
%        helper function.
%     2. Computes ANOVA contrast vectors and condition weight matrices with
%        getTreatWeights.
%     3. Applies hypothesis contrasts via runHCons and computes condition 
%        estimates with getCondEsts.
%     4. Plots the average effect of a linear predictor (Complexity) over a
%        specified domain.
%     5. Produces an HTML report summarising all the results using 
%        reportResults.
%
% INPUTS:
%   None.
%
% OUTPUT:
%   None. The function displays a figure and generates an HTML report.
%
% USAGE EXAMPLE:
%   runExample();
%
% NOTES:
%   - A synthetic dataset is generated within this function to simulate 
%     model fitting.
%   - This example illustrates how the individual components integrate into
%     a complete analysis pipeline.
% 
%% Produce an example model
[mdl,factorStruct] = getModel();
coefName = mdl.CoefficientNames';

%% Get the ANOVA contrast vectors and condition weights
[anovaHs,condXs,linXs] = getTreatWeights(coefName,factorStruct);

%% Run ANOVA contrasts and get the condition estimates
anovaHs = runHCons(mdl,anovaHs);
condEsts = getCondEsts(mdl,condXs);

%% Plot the average effect of Complexity as a linear predictor
domain = (-4:0.1:4)';
X = repmat(mean(linXs.Complexity(2:2:end,:).X(:,1:6),1),numel(domain),1);
X(isnan(X)) = domain;
in = num2cell(X,1);
[~,est,low,upp] = getXEUL(mdl,in{:});

Complexity = figure;
plot(domain,est,'DisplayName','est');
hold on;
plot(domain,low,'DisplayName','low');
plot(domain,upp,'DisplayName','upp');
legend;
xlabel('Complexity');
ylabel('y');
hold off;

%% Report the results
reportResults(mdl,anovaHs,condEsts,Complexity);
return

function [mdl,factorStruct] = getModel()
formula = 'y ~ 1 + Session*ActiveDrug*(Complex+Complexity)';
factorStruct(1).name   = 'Session';
factorStruct(1).coding = 'bin';
factorStruct(1).coefs  = [2,3];
factorStruct(2).name   = 'Treatment';
factorStruct(2).coding = 'bin';
factorStruct(2).coefs  = 4;
factorStruct(3).name   = 'Trialtype';
factorStruct(3).coding = 'bin';
factorStruct(3).coefs  = 5;
factorStruct(4).name   = 'Complexity';
factorStruct(4).coding = 'lin';
factorStruct(4).coefs  = 6;
factorStruct(4).centre = 0;
n = 12^3;
Session = randi([0,2],n,1);
Session01 = double(Session==1);
Session02 = double(Session==2);
Session = categorical(Session);
ActiveDrug = randi([0,1],n,1);
Complex = randi([0,1],n,1);
Complexity = nan(size(Complex));
for iD = 1:size(Complex)
    if Complex(iD)
        Complexity(iD) = randn(1,1);
    end
end
nanzscore = @(v) (v-nanmean(v))./nanstd(v);
Complexity = nanzscore(Complexity);
Complexity(isnan(Complexity)) = 0;
y = randn(n,1);
DataTable = table(...
    Session,Session01,Session02,ActiveDrug,Complex,Complexity,y);
mdl = fitlm(DataTable,formula);
return