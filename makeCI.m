function [CI, permCI] = makeCI(X2,y2,nPerm,condVec)

% MAKECI computes a weighted sum (cross-product; equivalent to
% least-squares multiple linear regression when X is random, which it
% should be) between the independent variables across trials and the
% dependent variable. Will also compute weighted sums with the trial
% ordering randomly permuted to obtain a null distribution. Can handle the 
% presence of conditions.
%
% X2 is the (preprocessed) sampling array. The first dimension should be
% the trials and the following dimensions represent the sampling
% dimensions. Preferably, X should contain only bubble centers. These can
% be normalized or not, in which case X can be a logical array (this will
% quicken the analysis).
%
% Y2 is the (preprocessed) 1 x nTrials vector of responses. Responses 
% cannot be multidimensional for now.
%
% NPERM is the desired numper of permutations. Optional.
%
% CONDVEC is a 1 x nTrials vector representing the condition associated to
% each trial. Optional.
%
% CI is the resulting array of classification images representing which
% predictors/samples correlate with responses for each condition (condition
% being the last dimension).
%
% PERMCI is the array of permutation classification images (permutations 
% are the first dimension and condition is the last). Optional.
%
%
% Laurent Caplette, 2017
%
%
% TO DO
% Add possibility of multidimensional responses 
% Add possibility of using sampling matrices as shown (with cut-off)
% Add possibility of using permutations with resampling (bootstraps)
% Add possibility of using other functions (e.g. regression, mutual info)

% Input checks
if ~isvector(y2), error('y2 must be a vector'); end
if ndims(X2)<2, error('X2 must be a N-dimensional array, N>=2'); end
if size(y2,1)>size(y2,2), y2 = y2'; end % y in row form
if nPerm<0 || nPerm>100000, error('Invalid number of permutations.'); end
if isempty(condVec), condVec = ones(size(y2)); end

% Parameters
nPred = size(X2);
nPred = nPred(2:end);
nCond = length(unique(condVec));

X2 = X2(:,:); % vectorize sampling dimensions of X
nElem = size(X2,2);

% Compute weighted sums
CI = zeros(nElem,nCond); 
if nPerm~=0, permCI = zeros(nPerm,nElem,nCond); end % if 1 cond, dimension will disappear
for cond = 1:nCond
    y2_cond = y2(condVec==cond);
    X2_cond = X2(condVec==cond,:);
    CI(:,cond) = y2_cond*X2_cond; % weighted sum
    if nPerm~=0
        for perm = 1:nPerm
            index = randperm(length(y2_cond));
            permCI(perm,:,cond) = y2_cond(index)*X2_cond; % weighted sum
        end
    end 
end

% Reshape classification images
CI = reshape(CI,[nPred nCond]);
if nPerm~=0, permCI = reshape(permCI, [nPerm nPred nCond]); end
