function [X2,y2,condVec2] = preprocess(X,y,XNorm,yNorm,selectionAbs,selectionStd,otherSelVar,condVec)

% PREPROCESS will standardize sampling arrays and responses, and select
% trials according to given criteria.
%
% X is the raw sampling array. The first dimension should be the trials,
% and following dimensions should be sampling dimensions.
%
% y is a 1 x nTrials vector of responses (e.g. -1 and 1 for accuracy).
%
% XNorm is a boolean indicating whether or not to standardize predictors
% within each trial (default: true). Highly recommended if a variable
% number of bubbles across trials is used.
%
% yNorm is a boolean indicating whether or not to standardize responses
% (default: true). If there is no standardization, responses should still
% sum to 0.
%
% selectionAbs is a 2-element vector indicating respectively below and
% above which values of y, if any, to reject trials. Optional.
%
% selectionStd is a scalar or 2-element vector indicating, in standard
% deviations, the criterion to use to reject extreme negative and positive
% values of y, within each condition. If selectionStd is a scalar, the 
% function will use the same absolute value, but with different signs, for 
% the lower and upper bounds. Optional.
%
% otherSelVar is a structure or a cell array of structures. Each structure
% allows to reject values of y based on values of another variable and some
% criterion. Each structure should contain a field 'y' containing the
% 1 x Ntrl vector of variables, a field 'selAbs' containing the 2-element
% vector of absolute criteria for lower and upper bounds, and a field
% 'selStd' containing a scalar or 2-element vector of criteria in standard
% deviations. Optional.
%
% condVec is a 1 x nTrials vector indicating which condition is associated
% with each trial. Optional.
%
% X2 is the preprocessed sampling array.
%
% y2 is the preprocessed vector of responses.
%
% condVec2 is the 1 x nTrials vector of conditions with rejected trials 
% removed.
%
%
% Copyright (c) 2017 Laurent Caplette
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.%
%
%
% TO DO
% Make it possible to have multidimensional response variables
% Add possibility of transforming response variable (median, quantile)
% Add possibility of padding X, given prop, if natural padding is lost

% Defaults
if isempty(XNorm), XNorm = true; end
if isempty(yNorm), yNorm = true; end
if isempty(selectionAbs), selectionAbs = [NaN NaN]; end
if isempty(selectionStd), selectionStd = [NaN NaN]; end
if isempty(condVec), condVec = ones(1,length(y)); end
    
% Extract selection structures
if isstruct(otherSelVar) 
    selVar{1} = otherSelVar; 
elseif iscell(otherSelVar)
    selVar = otherSelVar; 
elseif ~isempty(otherSelVar)
    error('otherSelVar must be a structure or cell array of structures');
end
selVar{length(otherSelVar)+1}.y = y;
selVar{length(otherSelVar)+1}.selAbs = selectionAbs;
selVar{length(otherSelVar)+1}.selStd = selectionStd;

% Input checks
if ~isvector(y), error('y must be a vector'); end
if size(X,1)~=length(y), error('X and y must have the same nb of trials'); end
if ndims(X)<2, error('X must be a N-dimensional array, N>=2'); end
if XNorm~=0 && XNorm~=1, error('XNorm must be 0 or 1'); end
if yNorm~=0 && yNorm~=1, error('yNorm must be 0 or 1'); end
expectedNames = {'selAbs'; 'selStd'; 'y'};
for ii = 1:length(selVar)
    v = selVar{ii};
    if numel(v.selAbs)~=2, error('selectionAbs must be a 2 element vector'); end
    if any(~strcmpi(sort(fieldnames(v)),expectedNames)), error('Structure fields are not ok.'); end
    if isscalar(v.selStd), selVar{ii}.selStd = [-abs(v.selStd) abs(v.selStd)]; end
end

% Parameters
nCond = length(unique(condVec));
nPred = size(X);
nPred = nPred(2:end);
X = X(:,:); % vectorize sampling dimensions of X

% Rejection of extreme trials
extTrials = find(sum(X,2)'==0); % find trials with no bubbles in case there were any
for ii = 1:length(selVar)
    v = selVar{ii};
    if ~isnan(v.selAbs(1)), extTrials = [extTrials, find(v.y<v.selAbs(1))]; end
    if ~isnan(v.selAbs(2)), extTrials = [extTrials, find(v.y>v.selAbs(2))]; end
    for cond = 1:nCond % for relative criteria, find extremes within each condition
        y_cond = v.y(condVec==cond);
        if ~isnan(v.selStd(1)), extTrials = [extTrials, find(zscore(y_cond)<v.selStd(1) | zscore(y_cond)>v.selStd(2))]; end
    end
end
selected = true(1,length(y));
selected(extTrials) = false;
X2 = X(selected,:);
y2 = y(selected);
condVec2 = condVec(selected);

% Standardizations
if yNorm,    
    for cond = 1:nCond, y2(condVec2==cond) = zscore(y2(condVec2==cond)); end
end
if XNorm, X2 = zscore(X2,1,2); end
X2 = reshape(X2, [length(y2) nPred]); 





