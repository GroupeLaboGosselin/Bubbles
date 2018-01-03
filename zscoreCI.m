function ZCI = zscoreCI(SCI,conds,varargin)

% ZSCORECI z-scores a (smooth) classification image or array of
% classification images according to one of seven methods. This must be
% performed before using the RFT method to compute statistics.
%
%
% SCI is the (smooth) classification image, or array of classification
% images (if there are many conditions).
%
% CONDS is a boolean indicating whether there are multiple conditions or
% not (as the last dimension of SCI). Optional.
%
% METHOD #1: Z-SCORE WITH ITSElF.
% ZCI = zscoreCI(SCI,conds)
% No other arguments required.
%
% METHOD #2: Z-SCORE WITH ASSUMED NO-SIGNAL REGION.
% ZCI = zscoreCI(SCI,conds,noiseMask)
% NOISEMASK is a boolean array the same size as the classification image.
% Ones indicate the no-signal region in the SCI classification image to use
% to compute noise mean and standard deviation.
%
% METHOD #3: Z-SCORE WITH PREDEFINED MEAN(S).
% ZCI = zscoreCI(SCI,conds,mu)
% MU is the predefined mean(s). One for each condition.
% Standard deviation(s) will be computed from SCI.
%
% METHOD #4: Z-SCORE WITH PREDEFINED MEAN(S) AND STANDARD DEVIATION(S).
% ZCI = zscoreCI(SCI,conds,mu,sigma)
% MU is the predefined mean(s). One for each condition.
% SIGMA is the predefined standard deviation(s). One for each condition.
%
% METHODS #5-6: Z-SCORE USING PERMUTATIONS.
% ZCI = zscoreCI(SCI,conds,permSCI,permMethod)
% PERMSCI is the array of permutation classification images (permutations
% as the first dimension).
% PERMMETHOD is a string indicating the z-scoring method. Can be 'across'
% if all predictors are averaged together or 'within' if z-scoring should
% be performed within each predictor.
%
% METHOD #7: Z-SCORE USING ANALYTICAL FORMULA.
% ZCI = zscoreCI(SCI,conds,nCorrect,nIncorrect,stdNoise,filter)
% NCORRECT is the number of correct trials.
% NINCORRECT is the number of incorrect trials.
% STDNOISE is the standard deviation of the white gaussian noise.
% FILTER is the Gaussian filter used to perform smoothing.
% For an explanation, see:
% 	Chauvin, Worsley, Schyns, Arguin & Gosselin (2005). Accurate
% 	statistical tests for smooth classification images. Journal of Vision,
% 	5, 659-667.
%
%
% Copyright (c) 2017, Laurent Caplette, Alan Chauvin & Frédéric Gosselin
% Code written by Laurent Caplette
% Some parts are based on code by Alan Chauvin & Frédéric Gosselin
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
% along with this program.  If not, see <https://www.gnu.org/licenses/>.


% Default
if ~exist('conds','var') || isempty(conds), conds = 0; end
if ~isscalar(conds), error('Conds must be 0 or 1'); end

% Parameters
nPred = size(SCI);
nPred = nPred(nPred~=1); % including conditions
if conds
    nCond = nPred(end);
    nElemNC = prod(nPred(1:end-1)); % excluding conditions
end

switch nargin
 
    case {1,2} % Option #1: Z-score using only CI

        if conds
            SCI = reshape(SCI, [nElemNC nCond]);
            ZCI = zeros(nElemNC,nCond);
            for cond = 1:nCond, ZCI(:,cond) = (SCI(:,cond)-mean(SCI(:,cond)))./std(SCI(:,cond)); end
            ZCI = reshape(ZCI, nPred);
        else
            ZCI = (SCI-mean(SCI(:)))./std(SCI(:));
        end
        
    case 3
        
        if ~isscalar(varargin{1}) && all(size(varargin{1})==size(SCI)) % Option #2: Z-score using no-signal region
            
            noiseMask = varargin{1};
            if conds
                SCI = reshape(SCI, [nElemNC nCond]);
                ZCI = zeros(nElemNC,nCond);
                noise = cell(1,nCond);
                for cond = 1:nCond
                    noise{cond} = noiseMask.*SCI(:,cond);
                    ZCI(:,cond) = (SCI(:,cond)-mean(noise{cond}(:)))./std(noise{cond}(:));
                end
                ZCI = reshape(ZCI, nPred);
            else
                noise = noiseMask.*SCI;
                ZCI = (SCI-mean(noise(:)))./std(noise(:));
            end
            
        elseif isvector(varargin{1}) && numel(varargin{1})<numel(SCI) % Option #3: Z-score using CI and given mean
        
            mu = varargin{1};
            if conds
                SCI = reshape(SCI, [nElemNC nCond]);
                ZCI = zeros(nElemNC,nCond);
                for cond = 1:nCond, ZCI(:,cond) = (SCI(:,cond)-mu(cond))./std(SCI(:,cond)); end
                ZCI = reshape(ZCI, nPred);
            else
                ZCI = (SCI-mu)./std(SCI(:));
            end
            
        else 
            
            error('Third argument must be a mask of the same size as SCI or a vector of mean(s)')
            
        end
        
    case 4
        
        if ~ischar(varargin{2}) % Option #4: Z-score using CI and given mean and std
            
            mu = varargin{1};
            sigma = varargin{2};
            if conds
                SCI = reshape(SCI, [nElemNC nCond]);
                ZCI = zeros(nElemNC,nCond);
                for cond = 1:nCond, ZCI(:,cond) = (SCI(:,cond)-mu(cond))./sigma(cond); end
                ZCI = reshape(ZCI, nPred);
            else
                ZCI = (SCI-mu)./sigma;
            end
            
        elseif strcmpi(varargin{2},'across') % Option #5: Z-score using permutations, across all elements
            
            permSCI = varargin{1};
            nPerm = size(permSCI,1);
            if conds
                SCI = reshape(SCI, [nElemNC nCond]);
                permSCI = reshape(permSCI, [nPerm*nElemNC nCond]);
                ZCI = zeros(nElemNC,nCond);
                for cond = 1:nCond, ZCI(:,cond) = (SCI(:,cond)-mean(permSCI(:,cond)))./std(permSCI(:,cond)); end
                ZCI = reshape(ZCI, nPred);
            else
                ZCI = (SCI-mean(permSCI(:)))./std(permSCI(:));
            end
            
        elseif strcmpi(varargin{2},'within') % Option #6: Z-score using permutations within each element
            
            permSCI = varargin{1};
            if size(permSCI,1) < 100, error('Not enough permutations to z-score within each element'); end
            ZCI = (SCI-squeeze(mean(permSCI)))./squeeze(std(permSCI)); % works with or without conds
                 
        else 
            
            error('Fourth argument must either be a vector of standard deviation(s) or the ''across'' or ''within'' string')
            
        end
        
    case 6 % Option #7: Z-score using analytical formula
        
        nc = varargin{1};
        ni = varargin{2};
        stdNoise = varargin{3};
        filter = varargin{4};
        n = ni + nc;
        nVarX = n*stdNoise^2*sum(filter(:).^2);
        nVarY = nc*ni/n;
        ZCI = sqrt(n)*SCI/sqrt(nVarX)/sqrt(nVarY);
              
    otherwise
        
        error('Incorrect number of arguments')
        
end

    