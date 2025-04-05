function [] = reportResults(varargin)
% reportResults.m
% Generate an HTML report summarising model results and figures.
% Sam Berens (s.berens@sussex.ac.uk)
% 05/04/2025
%
% reportResults(arg1, arg2, ..., argN)
%
% DESCRIPTION:
%   This function compiles various outputs—such as model details, contrast
%   tables, condition estimates, and figures—into a formatted HTML file.
%   It writes a temporary text file capturing all non-figure inputs and 
%   saves figure inputs as PNG images. The final HTML report, which uses a 
%   black background with green courier text, is then automatically opened 
%   in a web browser.
%
% INPUTS:
%   varargin - A variable number of inputs, which may include:
%              - Textual or numeric data (e.g., model summaries, tables).
%              - MATLAB figure objects; these are saved as images and 
%                embedded in the report.
%
% OUTPUT:
%   - No direct outputs. The function creates an HTML file (with a
%     timestamped filename) that displays the compiled results and opens it
%     for viewing.
%
% USAGE EXAMPLE:
%   reportResults(mdl, anovaHs, condEsts, figureHandle);
%
% NOTES:
%   - Non-figure inputs are printed directly, while figures are processed 
%     and saved into a 'Figs' directory.
%   - The function uses MATLAB’s diary functionality to capture console 
%     output.
% 
%% Record the inputs in a temporary text file
nInputs = nargin;
format short;
diary('temp.txt');
diary on;
for iInputs = 1:1:nInputs
    if ~isa(varargin{iInputs},'matlab.ui.Figure')
        fprintf('<h1>%s</h1>%c',inputname(iInputs),10)
        disp(varargin{iInputs});
        fprintf('%c%c%c%c%c',10,10,10,10,10)
    else
        if ~exist('Figs','dir')
            mkdir('Figs');
        end
        print(varargin{iInputs},sprintf('.%sFigs%s%s.png',filesep,...
            filesep,inputname(iInputs)),'-dpng');
        fprintf('<h1>%s</h1>%c',inputname(iInputs),10)
        fprintf('<img src=".%sFigs%s%s.png">%c%c','/','/',...
            inputname(iInputs),10,10)
        close(varargin{iInputs});
    end
end
diary off;

%% Select data from the text so it may be inserted into an HTML body
FileId_In = fopen('temp.txt');
TextData = fread(FileId_In);
fclose all;
Str1 = double('diary on;')';
Str2 = double('diary off;')';

ToCopy = ones(size(TextData));
for iTextData = (size(TextData,1)-size(Str1,1)+1):-1:1
    ToMatch = TextData(iTextData:(iTextData+size(Str1,1)-1),1);
    Match = sum(double(ToMatch == Str1),1) == size(Str1,1);
    if Match
        ToCopy(iTextData:(iTextData+size(Str1,1)-1),1) = zeros(size(Str1));
    end
end
for iTextData = (size(TextData,1)-size(Str2,1)+1):-1:1
    ToMatch = TextData(iTextData:(iTextData+size(Str2,1)-1),1);
    Match = sum(double(ToMatch == Str2),1) == size(Str2,1);
    if Match
        ToCopy(iTextData:(iTextData+size(Str2,1)-1),1) = zeros(size(Str2));
    end
end
ToCopy = logical(ToCopy);
Html_Body = TextData(ToCopy);
Html_Body = Html_Body(1:end,:);

%% Construct the rest of the HTML file
Html_Head = double(sprintf(...
    ['<!DOCTYPE html>%c<html>%c<body style="background-color:#000000;',...
    'color:#00ff40;font-family:courier;white-space:pre;">%c'],10,10,10))';
Htlm_End = double(sprintf('%c</html>',10))';
Html_All = [Html_Head;Html_Body;Htlm_End];

%% Delete the temporary file
delete('temp.txt');

%% Write the output file
FileName_Out = sprintf('ResultsDisp_%04d%02d%02d%02d%02d%i.html',...
    round(clock));
FileID_Out = fopen(FileName_Out,'w');
fwrite(FileID_Out,Html_All);
fclose all;
if ispc
    winopen(FileName_Out);
elseif ismac
    system(['open ', FileName_Out]);
else % Assume Linux or Unix
    system(['xdg-open ', FileName_Out]);
end
return