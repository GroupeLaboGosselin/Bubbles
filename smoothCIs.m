function [SCI, permSCI] = smoothCIs(CI,permCI,sigma,padding,conds)

% SMOOTHCIS smooths (observed and permutation) classification images of
% arbitrary dimensionality, handling non-smooth dimensions and presence of
% conditions.
%
%
% CI is a non-smoothed observed classification image of any size (1D or 
% more). If there are more than one condition, these should be the last
% dimension.
%
% PERMCI is an array of permutation classification images (2D or more, the
% first dimension being the permutations; must be 1D more than the observed
% classification image). If there are more than one condition, these should 
% be the last dimension. Optional.
%
% SIGMA is a vector of standard deviations (in pixels) for the smoothing
% kernel in each dimension (excluding conditions). Standard deviations can 
% be 0 if the corresponding dimension is non-smooth. If the standard
% deviation is the same across dimensions, sigma can be a scalar even for
% >1D classification images. No sigma should be put for the conditions
% dimension if there is one.
%
% PADDING is a string indicating the padding method to be applied. Padding 
% can be 'none' (default), 'circular', 'replicate', or 'symmetric' (see 
% padarray function documentation). Padding can take a lot of time, so it 
% is better to have padded CIs beforehand (possibly obtained with padded 
% Xs); this will also give better results. If there is no padding, it is
% highly recommended to add some, or else edge artifacts will occur.
%
% CONDS indicates the presence of conditions (as the last dimension of the
% classification images). Can be 0 (default) or 1.
%
% SCI is the smoothed observed classification image.
%
% PERMSCI is an array of smoothed permutation classification images.
%
%
% (c) 2017 Laurent Caplette & Frédéric Gosselin
% Code written by Laurent Caplette, 2017
% Some parts are based on code by Frédéric Gosselin


% Defaults
if ~exist('padding','var') || isempty(padding), padding = 'none'; end
if ~exist('conds','var') || isempty(conds), conds = 0; end

% Parameters
nPred = size(CI); % nb of predictors in each dim, including conditions
nPred = nPred(nPred~=1); % singleton dimensions are not a real dimension
nDims = numel(nPred); % nb of dims, including conditions
nDimsNC = nDims - conds; % nb of dims, excluding conditions
nPerm = size(permCI,1);

% Input checks and transformations
expectedPadding = {'none', 'circular', 'replicate', 'symmetric'};
if ~any(strcmpi(padding, expectedPadding)), error('Invalid padding method.'); end
if conds~=0 && conds~=1, error('conds must be 0 or 1'); end
if nDims==1 && isrow(CI) 
    CI = CI'; % 1-D vector in column form
    transposeFlag = true;
end 
if numel(sigma)==1, sigma = repmat(sigma,[1 nDimsNC]); end
if numel(sigma)~=nDimsNC, error('Sigma size must be 1 or nb of dimensions'); end
if all(sigma==0)
    SCI = CI;
    permSCI = permCI;
    return
end
if conds, sigma = [sigma 0]; end % not smooth across condition dimension

% Determine kernel size
sigma2 = sigma(sigma~=0); % sigma without the zeros
kern = zeros(1,numel(sigma2));
for k = 1:numel(sigma2), kern(k) = ceil(6*sigma2(k)); end % size of 6*sigma

% Pad CIs
if strcmpi(padding,'none')
    PCI = CI;
    permPCI = permCI;
else % pad only dimensions to be smoothed
    if nDims==1, % necessarily smooth
        PCI = padarray(CI, [nPred 0], padding);
    else
        PCI = padarray(CI, nPred.*(sigma~=0), padding); 
    end
    permPCI = padarray(permCI, [0 nPred.*(sigma~=0)], padding);
end
nPredP = size(PCI);
nPredP = nPredP(nPredP~=1);

% Create mesh and bubble
inputs = cell(1,numel(sigma2));
for ii = 1:numel(sigma2), inputs{ii} = (0:kern(ii))-kern(ii)/2; end
outputs = cell(1,numel(sigma2));
[outputs{:}] = ndgrid(inputs{:}); % inputs & outputs as comma-separated lists for any nDims
bubble = 0;
for ii = 1:numel(sigma2), bubble = bubble - (outputs{ii}.^2/sigma2(ii).^2); end
bubble = exp(bubble);
bubble = bubble./max(bubble(:));

% Order (smooth, then non-smooth) and collapse non-smooth in 1 dimension
if nDims~=1
    order = [find(sigma~=0) find(sigma==0)]; % smooth, then non-smooth
    PCI = permute(PCI, order);
    permPCI = permute(permPCI, [1 order+1]); % perms, then smooth, then non-smooth
    PCI = reshape(PCI, [nPredP(sigma~=0) prod(nPredP(sigma==0))]);
    permPCI = reshape(permPCI, [nPerm nPredP(sigma~=0) prod(nPredP(sigma==0))]);
end

% Convolve CIs
SCI = zeros(prod(nPredP(sigma~=0)),prod(nPredP(sigma==0)));
permSCI = zeros(nPerm,prod(nPredP(sigma~=0)),prod(nPredP(sigma==0)));
S.subs = cell(1,numel(sigma2)+1);
for ii = 1:numel(sigma2), S.subs{ii} = ':'; end % subscripts for smooth dimensions in PCI
S.type = '()';
S2.subs = cell(1,numel(sigma2)+2);
for ii = 2:numel(sigma2)+1, S2.subs{ii} = ':'; end % subscripts for smooth dimensions in permPCI
S2.type = '()';
for pred = 1:prod(nPredP(sigma==0)) % convolution for each non-smooth predictor
    S.subs{numel(sigma2)+1} = pred; % subscript for non-smooth dimension in PCI
    S2.subs{numel(sigma2)+2} = pred; % subscript for non-smooth dimension in permPCI
    SCI(:,pred) = vect(convn(squeeze(subsref(PCI,S)), bubble, 'same'));
    for perm = 1:nPerm
        S2.subs{1} = perm; % subscript for perm dimension
        permSCI(perm,:,pred) = vect(convn(squeeze(subsref(permPCI,S2)), bubble, 'same'));
    end
end

% Unwrap and place back in original order
if nDims~=1
    SCI = reshape(SCI, [nPredP(sigma~=0) nPredP(sigma==0)]);
    permSCI = reshape(permSCI, [nPerm nPredP(sigma~=0) nPredP(sigma==0)]);
    [~, reverseOrder] = sort(order); 
    SCI = permute(SCI, reverseOrder);
    permSCI = permute(permSCI, [1 reverseOrder+1]);
end

% Unpad CIs
if ~strcmpi(padding,'none')
    if nDims==1
        SCI = SCI(nPred+1:2*nPred);
        if transposeFlag, SCI = SCI'; end
        if ~isempty(permCI), permSCI = permSCI(:,nPred+1:2*nPred); end
    else
        nPad = nPred.*(sigma~=0);
        S.type = '()';
        S.subs = cell(1,nDims);
        for dim = 1:nDims, S.subs{dim} = nPad(dim)+1:nPad(dim)+nPred(dim); end
        SCI = subsref(SCI,S);
        if ~isempty(permCI), 
            S.subs = cell(1,nDims+1);
            S.subs{1} = ':';
            for dim = 2:nDims+1, S.subs{dim} = nPad(dim-1)+1:nPad(dim-1)+nPred(dim-1); end
            permSCI = subsref(permSCI,S);
        end
    end
end
