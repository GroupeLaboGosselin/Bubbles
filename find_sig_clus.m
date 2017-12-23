function ind = find_sig_clus(SCI, smoothDims, clus_type, thresh_clus, clusValThresh, c)

% FIND_SIG_CLUS finds if classification image contains any significant 
% clusters, and returns the linear indices of the significant indices. See
% STATCI_MAX function documentation for details.
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
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% Parameters
nPred = size(SCI);
nPred = nPred(nPred~=1);
nPredS = nPred(smoothDims);
nPredNS = nPred(~smoothDims);
nElemNS = prod(nPredNS);
if numel(nPredS)==1, nPredS = [nPredS 1]; end
if numel(nPredNS)==1, nPredNS = [nPredNS 1]; end
if ~exist('c','var') || isempty(c), c = conndef(sum(smoothDims),'maximal'); end % maximal connectivity as default

% Order and reshape
order = [find(smoothDims) find(~smoothDims)];
SCI = permute(SCI, order); % smooth, then non-smooth
SCI = reshape(SCI, [prod(nPredS) nElemNS]); % collapse smooth and non-smooth in 1 dim each

tvar = SCI>thresh_clus; % thresholded array

% Find significant clusters
SCI_sig = false(prod(nPredS),nElemNS);
for pred = 1:nElemNS
    temp = reshape(tvar(:,pred), nPredS);
    SCI_temp = reshape(SCI(:,pred), nPredS);
    [L,n] = bwlabeln(temp,c);
    cval = zeros(1,n);
    if strcmpi(clus_type,'size'),
        for ii = 1:n, cval(ii) = sum(vect(L==ii)); end
    else % sum
        for ii = 1:n, cval(ii) = sum(vect(SCI_temp(L==ii))); end
    end
    sigClusIdx = find(cval>clusValThresh);
    SCI_sig_temp = false(size(SCI_temp));
    for ii = 1:length(sigClusIdx)
        SCI_sig_temp(L==sigClusIdx(ii)) = true;
    end
    SCI_sig(:,pred) = SCI_sig_temp(:);
end

% Put back in original order and find indices of significant pixels
SCI_sig = squeeze(reshape(SCI_sig, [nPredS nPredNS]));
[~, reverseOrder] = sort(order); 
SCI_sig = permute(SCI_sig, reverseOrder);
ind = find(SCI_sig);
