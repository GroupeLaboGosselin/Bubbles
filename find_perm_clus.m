function clusVal = find_perm_clus(permSCI, smoothDims, clus_type, thresh_clus, c)

% FIND_PERM_CLUS finds max cluster values for each permutation. See
% STATCI_MAX function documentation for details.
%
% Laurent Caplette, 2017

% Parameters
nPerm = size(permSCI,1);
nPred = size(permSCI);
nPred = nPred(2:end);
nPredS = nPred(smoothDims); % nb of predictors in each smooth dimension
nPredNS = nPred(~smoothDims); % nb of predictors in each non-smooth dimension
nElemNS = prod(nPredNS); % nb of predictors in all non-smooth dimensions
if numel(nPredS)==1, nPredS = [nPredS 1]; end

% Order and reshape permSCI
order = [find(smoothDims) find(~smoothDims)]; 
permSCI = permute(permSCI, [1 order+1]); % perm, then smooth, then non-smooth
permSCI = reshape(permSCI, [nPerm*prod(nPredS) nElemNS]); % collapse perm & smooth in 1 dim, and non-smooth in 1 dim

tvar = permSCI>thresh_clus; % thresholded array

% Find max cluster value for each permutation and predictor in non-smooth dimensions
clusVal = zeros(nElemNS,nPerm);
for pred = 1:nElemNS
    temp = reshape(tvar(:,pred), [nPerm nPredS]);
    permSCI_temp = reshape(permSCI(:,pred), [nPerm nPredS]);
    for perm = 1:nPerm
        temp2 = reshape(temp(perm,:), nPredS);
        [L,n] = bwlabeln(temp2,c);
        if n==1 % clusVal will stay 0 if n==0
            if strcmpi(clus_type,'size')
                clusVal(pred,perm) = sum(vect(L));
            else % sum
                clusVal(pred,perm) = sum(vect(permSCI_temp(perm,L==1))); 
            end
        elseif n>1
            val = zeros(1,n);
            if strcmpi(clus_type,'size')
                for ii = 1:n, val(ii) = sum(vect(L==ii)); end
            else % sum
                for ii = 1:n, val(ii) = sum(vect(permSCI_temp(perm,L==ii))); end
            end
            clusVal(pred,perm) = max(val);
        end
    end
end
if nElemNS~=1, clusVal = max(clusVal); end % 1 x nPerm vector
