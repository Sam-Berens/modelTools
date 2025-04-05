function [anovaHs,condXs,linXs] = getTreatWeights(coefName,factorStruct)
% getTreatWeights.m
% Compute contrast matrices & condition weights for treatment-coded models.
% Sam Berens (s.berens@sussex.ac.uk)
% 05/04/2025
%
% [anovaHs, condXs, linXs] = getTreatWeights(coefName, factorStruct)
% 
% This function computes contrast matrices (H) for ANOVA analyses and
% condition weight matrices (X) for a statistical model based on treatment
% effects. It takes as inputs a cell array of coefficient names and a
% factor structure that describes the predictors (factors) in the model,
% including their coding schemes and associated coefficients.
% 
% INPUTS:
%   coefName     - A cell array of strings representing the names of the 
%                  fixed effects (main effects and interaction terms).
%                  Interaction terms must be specified using colons (':').
% 
%   factorStruct - A structure array where each element defines a factor
%                  with the following fields:
%                      .name   : A string with the factor's name.
%                      .coding : A string indicating the coding type, 
%                                either 'bin' for binary or 'lin' for
%                                linear.
%                      .coefs  : A numeric vector with indices 
%                                corresponding to the factor's coefficients
%                                in coefName.
%                      .centre : (Required for 'lin' factors only)
%                                The centre value used for scaling
%                                continuous predictors.
% 
% OUTPUTS:
%   anovaHs - A table where each row contains a contrast matrix (H) for a 
%             main effect or interaction term. The row names indicate the 
%             effect (e.g., 'Session:Treatment').
% 
%   condXs  - A table of condition weight vectors computed from a reduced 
%             design matrix that encodes the experimental conditions based
%             on the factors.
% 
%   linXs   - A structure with fields corresponding to each linearly coded 
%             factor. Each field contains a table of condition weight
%             vectors where continuous predictors are coded with NaNs.
% 
% DESCRIPTION:
%   The function first validates the inputs and sets default testing 
%   parameters if no inputs are provided. It then:
%     - Constructs a design matrix (X) by enumerating all possible
%       conditions based on the levels of each factor.
%     - Computes the structure of main effects and interactions from the 
%       coefficient names.
%     - Derives G matrices for main effects that, when multiplied by X,
%       yield the contrast matrices for each effect.
%     - Scales and adjusts the design matrix for continuous (linear)
%       factors based on their centre values.
%     - Produces two sets of condition weight matrices: one (condXs) for
%       the binary-coded conditions and one (linXs) specific to continuous
%       predictors.
% 
% EXAMPLE:
%   % Define coefficient names (including interaction terms with ':')
%   coefName = {'(Intercept)'; 'Session01'; 'Session02'; 'Activedrug';
%               'Complex'; 'Complexity'; 'Session01:Activedrug';
%               'Session02:Activedrug'; 'Session01:Complex';
%               'Session02:Complex'; 'Activedrug:Complex';
%               'Session01:Complexity'; 'Session02:Complexity';
%               'Activedrug:Complexity'; ...
%               'Session01:Activedrug:Complex';
%               'Session02:Activedrug:Complex';
%               'Session01:Activedrug:Complexity';
%               'Session02:Activedrug:Complexity'};
% 
%   % Define factor structure
%   factorStruct(1).name   = 'Session';
%   factorStruct(1).coding = 'bin';
%   factorStruct(1).coefs  = [2, 3];
% 
%   factorStruct(2).name   = 'Treatment';
%   factorStruct(2).coding = 'bin';
%   factorStruct(2).coefs  = 4;
% 
%   factorStruct(3).name   = 'Trialtype';
%   factorStruct(3).coding = 'bin';
%   factorStruct(3).coefs  = 5;
% 
%   factorStruct(4).name   = 'Complexity';
%   factorStruct(4).coding = 'lin';
%   factorStruct(4).coefs  = 6;
%   factorStruct(4).centre = 0;
% 
%   % Compute contrast and condition weight matrices
%   [anovaHs, condXs, linXs] = getTreatWeights(coefName, factorStruct);
% 
% NOTES:
%   - If no inputs are provided, the function runs in a test mode using a
%     default example model.
%   - The function supports only binary coding for factors with more than
%     two levels.
% 
%% Check the inputs
if nargin == 0
    %% You must be testing. Here is an example model to play with:
    % 'y ~ 1 + (Session01+Session02)*Activedrug*(Complex+Complexity)'
    coefName = {...
        '(Intercept)';                     % 01
        'Session01';                       % 02
        'Session02';                       % 03
        'Activedrug';                      % 04
        'Complex';                         % 05
        'Complexity';                      % 06
        'Session01:Activedrug';            % 07
        'Session02:Activedrug';            % 08
        'Session01:Complex';               % 09
        'Session02:Complex';               % 10
        'Activedrug:Complex';              % 11
        'Session01:Complexity';            % 12
        'Session02:Complexity';            % 13
        'Activedrug:Complexity';           % 14
        'Session01:Activedrug:Complex';    % 15
        'Session02:Activedrug:Complex';    % 16
        'Session01:Activedrug:Complexity'; % 17
        'Session02:Activedrug:Complexity'; % 18
        };
    factorStruct = struct();
    factorStruct(1,1).name = 'Session';
    factorStruct(1,1).coding = 'bin';
    factorStruct(1,1).coefs = [2,3];
    factorStruct(2,1).name = 'Treatment';
    factorStruct(2,1).coding = 'bin';
    factorStruct(2,1).coefs = 4;
    factorStruct(3,1).name = 'Trialtype';
    factorStruct(3,1).coding = 'bin';
    factorStruct(3,1).coefs = 5;
    factorStruct(4,1).name = 'Complexity';
    factorStruct(4,1).coding = 'lin';
    factorStruct(4,1).coefs = 6;
    factorStruct(4,1).centre = 0;
    
elseif nargin == 1
    error(['Two inputs are required:%c',...
        ' - a cell array of coefficient names;%c',...
        ' - a valid factor structure;%c'...
        'Type "edit %s" for more info.'],10,10,10,mfilename);
elseif nargin > 2
    error('Too many input arguments.');
end

%% Check the validity of the factorStuct input
for iFS = 1:numel(factorStruct)
    if ~(strcmpi(factorStruct(iFS).coding,'bin') || ...
            strcmpi(factorStruct(iFS).coding,'lin'))
        error('Unrecognised coding for factor %i.',iFS);
    end
    if (numel(factorStruct(iFS).coefs) > 1) && ...
            strcmpi(factorStruct(iFS).coding,'lin')
        error(['This function only supports binary coding for factors ',...
            'with more than 2 levels.']);
    end
    if strcmpi(factorStruct(iFS).coding,'lin') && ...
            ~isfield(factorStruct(iFS),'centre')
        error('Mean values must be specified for linearly coded factors.');
    end
end

%% Set some important variables
% Total number of fixed effects
k = numel(coefName);

% Number of factors
f = numel(factorStruct);

% Number of levels per factor
l = cellfun(@(v)numel(v)+1,{factorStruct.coefs});

% Get the xTermStruct (the structure of the interactions)
xTermStruct = getXtermStruct(coefName);

%% Construct a matrix (Xidx) that enumerates all possible conditions ...
% ... and specifies the level of each factor within each condition.
% Xidx will have a row for each condition and a column for each factor.
Xidx = getXidx(l,f);

%% Construct a reduced design matrix (X) ...
% ... that encodes parameter weights for all the conditions listed in Xidx
X = getX(factorStruct,k,xTermStruct,Xidx);

%% Construct a set of G matrices for all the main effects (G_main)
% G matrices are right multiplied by X (a design matrix) to produce
% parameter contrast matrices (H) by explicitly contrasting conditions
% coded within the rows of X. As such, G matrices have dimensions [r,m],
% where r is the rank of the parameter contrast matrix (H) and m is the
% number of conditions encoded across the rows of X.
G_main = getG_main(factorStruct,X);

%%
%getHs(factorStruct,)
linFacs = find(strcmp({factorStruct.coding},'lin'));
iH = 0;
Hname = cell((2^f)-1,1);
Hmat = cell((2^f)-1,1);
for iCol = 1:f % Looping through the cols of Pascal's triangle
    Fx = nchoosek(1:f,iCol);
    for iFx = 1:size(Fx,1)
        %% Extract the terms (factor indices) ...
        % ... in the current effect (row of Fx) being considered
        t = Fx(iFx,:);
        fxName = strjoin({factorStruct(t).name},':');
        %% Below is how we compute the rank of this effect
        % However, we don't actually need to do this explicitly ...
        % ... we get it for free.
        % rank = prod(cellfun(@(v)numel(v),{factorStruct(t).coefs}'));
        %% Construct a contrast matrix that operates on rows of X (gg)
        % gg matrices for main effects have been computed above. gg
        % matrices for interaction terms are computing by performing an
        % element-wise multiplication with every row of the gg matrices for
        % each constituent main effect, and then stacking the resulting
        % rows on top of each other.
        for iT = 1:numel(t)
            if iT == 1
                gg = G_main{t(iT)};
            else
                gg = crossGs(gg,G_main{t(iT)});
            end
        end
        %% Figure out which coefficients need to be scaled
        facsToScale = linFacs(~ismember(linFacs,t));
        [coefsToScale,centres] = getCoefsToScale(...
            coefName,factorStruct,facsToScale);
        %% Scale columns of the design matrix (SX) according to the above
        SX = X;
        if ~isempty(coefsToScale)
            SX(:,coefsToScale) = SX(:,coefsToScale).*centres;
        end
        %% Compute the final contrast and add it to the list
        H = gg*SX;
        if sum(H.^2,'all') > 0
            iH = iH + 1;
            H = H ./ min(nonzeros(abs(H)),[],'all'); % Scale the contrast
            Hname{iH} = fxName;
            Hmat{iH} = H;
        end
    end
end

%% Remove all empty contrasts and construct an table of contrasts
sel = cellfun(@(v)~isempty(v),Hmat);
Hname = Hname(sel);
Hmat = Hmat(sel);
anovaHs = cell2table(Hmat,'VariableNames',{'H'},'RowNames',Hname);

%% Construct tables of condition weights
% First, trim X to remove rows that have zeros for each linear predictor
% and remove all columns that correspond to interactions (we will add these
% back in later).
rowsToSel = all(Xidx(:,linFacs)==1,2);
X = X(rowsToSel,cellfun(@(v)numel(v)==1,xTermStruct)');

% Get the names for each weighting vector
Xnames = cellfun(@(c1,c2)...
    sprintf('%s%i',c1,c2),...
    repmat({factorStruct.name},size(Xidx,1),1),...
    num2cell(Xidx),'UniformOutput',false);
Xnames = Xnames(rowsToSel,~ismember((1:f),linFacs));
Xnames = cellstr(join(string(Xnames),':',2));

% Replace all the ones in each column of X corresponding to a linear
% predictor with the centre for that predictor.
for iLinFac = 1:numel(linFacs)
    iCol = factorStruct(linFacs(iLinFac)).coefs;
    centre = factorStruct(linFacs(iLinFac)).centre;
    X(:,iCol) = centre;
end

% Now we prepare a new set of conditions weights that have complex values
% in the place of each linear predictor. These values can be replaced with
% values of the continuous predictor to give continuous model estimates.
% Complex values are used here as they obey standard rules of arithmetic
% and so can be used to compute interaction terms. However, we later
% replace these complex values with nans (after the interactions have been
% calculated).
linXs = cell(numel(linFacs),1);
for iLinFac = 1:numel(linFacs)
    iCol = factorStruct(linFacs(iLinFac)).coefs;
    linXs{iLinFac} = X;
    linXs{iLinFac}(:,iCol) = 1i;
end

% Add the interactions back in
X = addXterms(X,xTermStruct);
for iLinFac = 1:numel(linFacs)
    linXs{iLinFac} = addXterms(linXs{iLinFac},xTermStruct);
end

% Add the weights to a table
condXs = cell2table(num2cell(X,2),...
    'VariableNames',{'X'},'RowNames',Xnames);
for iLinFac = 1:numel(linFacs)
    % Convert the complex values to nans
    linXs{iLinFac}(linXs{iLinFac}==1i) = NaN;
    linXs{iLinFac} = cell2table(num2cell(linXs{iLinFac},2),...
        'VariableNames',{'X'},'RowNames',Xnames);
end
linXs = cell2struct(linXs,{factorStruct(linFacs).name}',1);
return

function [Xid] = getXidx(l,f)
% getXidx
% Enumerate all possible experimental conditions.
% Sam Berens (s.berens@sussex.ac.uk)
% 05/04/2025
% 
% Xid = getXidx(l, f) returns a matrix where each row represents a unique
% combination of factor levels for the experiment.
% 
% INPUTS:
%   l - A vector of length f where each element represents the number of 
%       levels for a corresponding factor (typically, number of 
%       coefficients + 1).
%   f - The number of factors.
% 
% OUTPUT:
%   Xid - An m-by-f matrix (with m = prod(l)) enumerating all possible 
%         combinations of factor levels. Factor levels are zero-indexed.
% 
% USAGE EXAMPLE:
%   Xid = getXidx([3 2 4], 3);
% 
m = prod(l); % Provisionally, the number of conditions
Xid = nan(m,f);
scale = flip(cumprod(flip([l(2:end),1])));
for iM = 1:m
    c = iM-1;
    for ig = 1:f
        n = floor(c/scale(ig));
        Xid(iM,ig) = n;
        c = c - n*scale(ig);
    end
end
return

function [X] = getX(factorStruct,k,xTermStruct,Xidx)
% getX
% Construct a design matrix encoding condition parameter weights.
% Sam Berens (s.berens@sussex.ac.uk)
% 05/04/2025
% 
% X = getX(factorStruct, k, xTermStruct, Xidx) builds a design matrix that
% represents the experimental conditions based on the provided factor
% structure and the enumerated condition indices.
% 
% INPUTS:
%   factorStruct - A structure array defining each factor (with fields such
%                  as name, coding, and coefs).
%   k            - Total number of coefficients (typically, length of
%                  coefName).
%   xTermStruct  - A cell array that describes the interaction structure
%                  from the coefficient names.
%   Xidx         - A matrix where each row specifies a unique combination 
%                  of factor levels (as produced by getXidx).
% 
% OUTPUT:
%   X - A design matrix where each row corresponds to an experimental
%       condition and each column corresponds to a parameter weight. The
%       matrix includes both main effects and interaction terms (added via
%       addXterms).
% 
% NOTE:
%   Factor levels are assumed to be encoded as 0 (baseline) and increasing
%   integers for higher levels.
% 
f = numel(factorStruct);
m = size(Xidx,1); % Provisionally, the number of conditions
X = nan(m,k);
X(:,1) = 1;
for iX = 1:size(X,1)
    for iFac = 1:f
        level = Xidx(iX,iFac);
        ffx = factorStruct(iFac).coefs;
        for iffx = 1:numel(ffx)
            if level==0
                X(iX,ffx) = 0;
            elseif level==iffx
                X(iX,ffx(iffx)) = 1;
            else
                X(iX,ffx(iffx)) = 0;
            end
        end
    end
end
X = addXterms(X,xTermStruct);
return

function [G_main] = getG_main(factorStruct,X)
% getG_main
% Generate G matrices for main effects.
% Sam Berens (s.berens@sussex.ac.uk)
% 05/04/2025
% 
% G_main = getG_main(factorStruct, X) computes a set of G matrices for each
% main effect defined in factorStruct. These matrices are used to derive
% parameter contrasts from the design matrix X.
% 
% INPUTS:
%   factorStruct - A structure array where each element defines a factor,
%                  including its coefficient indices.
%   X            - The design matrix encoding experimental conditions.
% 
% OUTPUT:
%   G_main - A cell array where each element is a matrix corresponding to a
%            factor. Each matrix has a number of rows equal to the number 
%            of coefficients for that factor and columns equal to the 
%            number of conditions.
% 
% NOTE:
%   For each condition, if all entries corresponding to a factor's
%   coefficients are zero, the matrix entry is set to -1. If the condition
%   matches a specific coefficient, the entry is set to 1; otherwise, it is
%   0.
% 
f = numel(factorStruct);
G_main = cell(f,1);
for iFac = 1:f
    coefs = factorStruct(iFac).coefs;
    nRows = numel(coefs);
    G_current = nan(nRows,size(X,1));
    for iRow = 1:nRows
        for iX = 1:size(X,1)
            row = X(iX,:);
            if sum(row(coefs')) == 0
                G_current(iRow,iX) = -1;
            elseif row(coefs(iRow)) == 1
                G_current(iRow,iX) = 1;
            else
                G_current(iRow,iX) = 0;
            end
        end
    end
    G_main{iFac} = G_current;
end
return

function [coefsToScale,means] = getCoefsToScale(...
    coefName,factorStruct,facsToScale)
% getCoefsToScale
% Determine coefficients and means for scaling linear predictors.
% Sam Berens (s.berens@sussex.ac.uk)
% 05/04/2025
% 
% [coefsToScale, means] = getCoefsToScale(coefName, factorStruct,
% facsToScale) identifies the coefficients corresponding to linearly coded
% factors that require scaling and retrieves their associated centre
% values.
% 
% INPUTS:
%   coefName     - A cell array of strings representing coefficient names.
%   factorStruct - A structure array defining the factors, including their 
%                  coefficient indices and centre (mean) values.
%   facsToScale  - Indices of the factors (within factorStruct) for which 
%                  scaling is needed.
% 
% OUTPUTS:
%   coefsToScale - A row vector of indices corresponding to the 
%                  coefficients that need to be scaled.
%   means        - A row vector of centre values corresponding to these 
%                  coefficients.
% 
% NOTE:
%   This function uses string matching to determine which coefficients to
%   scale and sorts the resulting indices in ascending order.
% 
simpleCoeffsToScale = [factorStruct(facsToScale).coefs]; % Simple FX only
if ~isempty(simpleCoeffsToScale)
    simpleMeans = [factorStruct(facsToScale).centre]; % Same size as above
else
    simpleMeans = [];
end
coefsToScale = [];
means = [];
for ii = 1:numel(simpleCoeffsToScale)
    terms = find(contains(coefName,coefName(simpleCoeffsToScale(ii))));
    coefsToScale = [coefsToScale;terms]; %#ok<AGROW>
    means = [means;ones(size(terms)).*simpleMeans(ii)]; %#ok<AGROW>
end
[~,iSort] = sort(coefsToScale);
coefsToScale = coefsToScale(iSort)'; % Make a row vector
means = means(iSort)'; % Make a row vector
return

function [C] = crossGs(A,B)
% crossGs
% Compute a column-wise outer product of two G matrices.
% Sam Berens (s.berens@sussex.ac.uk)
% 05/04/2025
% 
% C = crossGs(A, B) computes, for each column, the outer product between
% the rows of matrices A and B. In other words, for each column index k,
% every element from column k in A is multiplied by every element from
% column k in B. The resulting products are then reshaped into a
% two-dimensional matrix.
% 
% This operation is used to construct interaction contrasts from main
% effect G matrices. It is not a simple element-wise (Hadamard) product;
% rather, it performs a column-wise outer product and flattens the result.
%
% INPUTS:
%   A - A matrix of size [n, k] representing the first set of contrast
%       weights.
%   B - A matrix of size [m, k] representing the second set of
%       contrast weights.
% 
% OUTPUT:
%   C - A matrix of size [n*m, k] where each column k is formed by
%       flattening the outer product of the k-th columns of A and B.
% 
% NOTE:
%   This function reshapes A and B to enable the outer product computation 
%   for each column, then converts the resulting three-dimensional array
%   into a two-dimensional matrix.
% 
[n,k] = size(A);
[m,~] = size(B);
C = reshape(A,[n,1,k]).*reshape(B,[1,m,k]);
C = reshape(C,[n*m,k]);
return