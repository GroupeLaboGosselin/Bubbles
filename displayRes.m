function [SCIpix, SCIclus] = displayRes(CIsize,sigPix,sigClus,displayOpt,mask)

% DISPLAYRES creates images illustrating which pixels are significant at 
% the pixel and cluster levels, and displays the results.
%
%
% CISIZE is a vector representing the size of the classification image.
%
% SIGPIX is a structure obtained after a pixel test in STATCI_MAX.
%
% SIGCLUS is a structure obtained after a cluster test in STATCI_MAX.
% Optional.
%
% DISPLAYOPT is a boolean indicating whether to display figures or not.
% Note than if the classification image contains 4 dimensions or more,
% figures will not be displayed irrespective of this argument. (Default: 1)
%
% MASK is a boolean array indicating for which pixels to show the
% significance value. Optional.
%
% SCIPIX is an image indicating which pixels are significant at the pixel 
% level. Pixels can be 0 (not significant), 1 (significant at the upper
% tail) or -1 (significant at the lower tail).
%
% SCICLUS is an image indicating which pixels are significant at the 
% cluster level. Pixels can be 0 (not significant), 1 (significant at the 
% upper tail) or -1 (significant at the lower tail).
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
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% Parameters
nDims = sum(CIsize~=1);

% Input checks
if isempty(displayOpt), displayOpt = 1; end
if displayOpt~=0 && displayOpt~=1, error('displayOpt must be 0 or 1'); end
displayOpt = logical(displayOpt);
if nDims>3 && displayOpt, warning('Figures will not be displayed. Too many dimensions.'); end
if ~exist('mask','var'), mask = true(CIsize); end

% Pixel-level significance image
SCIpix = zeros(1,prod(CIsize));
if iscell(sigPix.ind)
    SCIpix(sigPix.ind{1}) = -1;
    SCIpix(sigPix.ind{2}) = 1;
else
    SCIpix(sigPix.ind) = 1;
end
SCIpix = reshape(SCIpix, CIsize);

% Cluster-level significance image
SCIclus = zeros(1,prod(CIsize));
if iscell(sigClus.ind)
    SCIclus(sigClus.ind{1}) = -1;
    SCIclus(sigClus.ind{2}) = 1;
else
    SCIclus(sigClus.ind) = 1;
end
SCIclus = reshape(SCIclus, CIsize);

% Display figures
if nDims==2 && displayOpt
    figure, imagesc(SCIpix.*mask), colorbar, title('Pixel test')
    figure, imagesc(SCIclus.*mask), colorbar, title('Cluster test')
elseif nDims==3 && displayOpt
    SCIpix = permute(SCIpix, [find(CIsize~=min(CIsize)) find(CIsize==min(CIsize))]);
    SCIclus = permute(SCIclus, [find(CIsize~=min(CIsize)) find(CIsize==min(CIsize))]);
    for pred = 1:size(SCIpix,3)
        figure, imagesc(SCIpix(:,:,pred).*mask(:,:,pred)), colorbar, title(sprintf('Pixel test, Predictor #%i of 3rd dimension', pred))
        figure, imagesc(SCIclus(:,:,pred).*mask(:,:,pred)), colorbar, title(sprintf('Cluster test, Predictor #%i of 3rd dimension', pred))
    end
    [~, reverseOrder] = sort([find(CIsize~=min(CIsize)) find(CIsize==min(CIsize))]);
    SCIpix = permute(SCIpix, reverseOrder);
    SCIclus = permute(SCIclus, reverseOrder);
end
