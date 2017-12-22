function [sigPix, sigClus] = statCI_max(SCI,permSCI,smoothDims,tails,alpha_pix,clus_type,alpha_clus,thresh_clus,conds,c,mask)

% STATCI_MAX finds thresholds and indices of significant pixels in a
% classification image, accorging to a pixel- or cluster-level correction
% base on the maximum statistic method (Holmes, 1996). Calls FIND_PERM_CLUS
% and FIND_SIG_CLUS which do most of the job.
%
%
% SCI is the observed smooth classification image.
%
% PERMSCI is the array of smooth permutation classification images.
%
% SMOOTHDIMS is a 1 x nDim vector indicating which vectors are smooth
% (either with a 1 or with the size of the standard deviation).
%
% TAILS is a string indicating to which tail of the null distribution the
% observed statistics should be compared. Can be 'left', 'right', or
% 'both' (default).
%
% ALPHA_PIX is the p-value at which to apply the pixel-level correction.
% 
% CLUS_TYPE is the type of cluster test to be applied. This can be 'size',
% in which case significant clusters will be chosen only based on their
% size (nb of elements), or 'sum', in which case significant clusters will
% be chosen based on their total summed value.
%
% ALPHA_CLUS is the p-value at which to apply the cluster-level correction.
%
% THRESH_CLUS in the primary cluster threshold. The cluster test requires
% the choice of an arbitrary primary threshold: the test will then
% determine what is the required size of a cluster (or summed value) above
% this primary threshold for it to be significant. This can be a p-value or
% a raw regression coefficient (see smooth classification image).
%
% CONDS indicates the presence of conditions (as the last dimension of the
% classification images). Can be either 0 (no conditions; default) or 1.
%
% C specifies the connectivity to be used in the computation of the
% clusters (see bwlabel function documentation). Can be a scalar or an
% array as returned by the function conndef. Default is
% conndef(sum(smoothDims),'maximal).
%
% MASK is an array indicating on which pixels to perform the computations.
%
% SIGPIX is a structure indicating the threshold of the pixel test and the
% indices of the significant pixels.
%
% SIGCLUS is a structure indicating the threshold of the cluster test and
% the indices of the significant pixels.
%
%
% (c) 2017 Laurent Caplette


% Defaults
if ~exist('c','var') || isempty(c), c = conndef(sum(smoothDims),'maximal'); end % maximal connectivity as default
if ~exist('conds','var') || isempty(conds), conds = 0; end
if isempty(tails), tails = 'both'; end

% Parameters
nPerm = size(permSCI,1);
nPred = size(SCI);
nPred = nPred(nPred~=1);
nDims = numel(nPred);
nDimsNC = nDims - conds; % excluding conditions

% Input checks and transformations
if exist('mask','var') && ~isempty(mask) % CIs will contain NaNs where mask contains NaNs
    SCI = SCI.*mask; 
    permSCI = permSCI.*repmat(shiftdim(mask,-1),[nPerm 1]);
end
if conds~=0 && conds~=1, error('conds must be 0 or 1'); end
if numel(smoothDims)~=nDimsNC, error('smoothDims size must be nb of dimensions'); end
if conds, smoothDims = [smoothDims 0]; end % not smooth across conditions
smoothDims = logical(smoothDims); % if sigmas are put, will transform them in 1s
if nPerm<100, error('Too few permutations to use maximum statistic method'); end
if isscalar(thresh_clus) && abs(thresh_clus)<1
    if ~strcmpi(tails,'both') && nPerm<inv(thresh_clus)
        error('Not enough permutations for the chosen primary p-value')
    elseif strcmpi(tails,'both') && nPerm<inv(thresh_clus/2)
        error('Not enough permutations for the chosen primary p-value')
    end
end
expectedTails = {'left','right','both'};
if ~any(strcmpi(tails,expectedTails)), error('tails must be ''left'', ''right'' or ''both'''); end
if alpha_pix<=0 || alpha_pix>=1, error('alpha_pix must be between 0 and 1 exclusively'); end
if alpha_clus<=0 || alpha_clus>=1, error('alpha_clus must be between 0 and 1 exclusively'); end
expectedTypes = {'size','sum'};
if ~any(strcmpi(clus_type,expectedTypes)), error('clus_type must be ''size'' or ''sum'''); end
if strcmpi(tails,'both') && isscalar(thresh_clus) 
    if abs(thresh_clus)>=1
        thresh_clus = [-abs(thresh_clus) abs(thresh_clus)];
    else 
        thresh_clus = repmat(thresh_clus, [1 2]);
    end
end
if sum(smoothDims)==1 && c~=8, error('Connectivity for 1D smooth vectors must be 8'); end % vectors are treated as matrices
iptcheckconn(c,'statCI_max','c',10) % check if c is valid (for nDims>1), using existing matlab fct

% Transform p-value in coefficient if p-value is given as primary threshold, based on average data
if nargout>1 && all(abs(thresh_clus)<1)
    if strcmpi(tails,'left')
        thresh_clus = mean(quantile(permSCI(:,:),thresh_clus));
    elseif strcmpi(tails,'right')
        thresh_clus = mean(quantile(permSCI(:,:),1-thresh_clus));
    else % both
        thresh_clus(1) = mean(quantile(permSCI(:,:),thresh_clus(1)/2));
        thresh_clus(2) = mean(quantile(permSCI(:,:),1-(thresh_clus(2)/2)));
    end
end

% Pixel test
if strcmpi(tails,'left')
    sigPix.thresh = quantile(min(permSCI(:,:),[],2),alpha_pix);
    sigPix.ind = find(SCI<sigPix.thresh);
elseif strcmpi(tails,'right')
    sigPix.thresh = quantile(max(permSCI(:,:),[],2),1-alpha_pix);
    sigPix.ind = find(SCI>sigPix.thresh);
else % both
    sigPix.thresh(1) = quantile(min(permSCI(:,:),[],2),alpha_pix/2);
    sigPix.thresh(2) = quantile(max(permSCI(:,:),[],2),1-(alpha_pix/2));
    sigPix.ind{1} = find(SCI<sigPix.thresh(1));
    sigPix.ind{2} = find(SCI>sigPix.thresh(2));
end

% Cluster test
if nargout>1
    if strcmpi(tails,'left')
        clusVal = find_perm_clus(-permSCI, smoothDims, clus_type, abs(thresh_clus), c);
        sigClus.thresh = quantile(clusVal, 1-alpha_clus);
        sigClus.ind = find_sig_clus(-SCI, smoothDims, clus_type, abs(thresh_clus), sigClus.thresh, c);
        sigClus.thresh = -sigClus.thresh;
    elseif strcmpi(tails,'right')
        clusVal = find_perm_clus(permSCI, smoothDims, clus_type, abs(thresh_clus), c);
        sigClus.thresh = quantile(clusVal, 1-alpha_clus);
        sigClus.ind = find_sig_clus(SCI, smoothDims, clus_type, abs(thresh_clus), sigClus.thresh, c);
    else % both
        clusVal = find_perm_clus(-permSCI, smoothDims, clus_type, abs(thresh_clus(1)), c);
        sigClus.thresh(1) = quantile(clusVal, 1-(alpha_clus/2));
        sigClus.ind{1} = find_sig_clus(-SCI, smoothDims, clus_type, abs(thresh_clus(1)), sigClus.thresh(1), c);
        clusVal = find_perm_clus(permSCI, smoothDims, clus_type, abs(thresh_clus(2)), c);
        sigClus.thresh(2) = quantile(clusVal, 1-(alpha_clus/2));
        sigClus.ind{2} = find_sig_clus(SCI, smoothDims, clus_type, abs(thresh_clus(2)), sigClus.thresh(2), c);
        if strcmpi(clus_type,'sum'), sigClus.thresh = [-1 1] .* sigClus.thresh; end
    end
end
