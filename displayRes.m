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
% level. If significance is tested at one tail, pixels can be 0 or 1; if
% significance is tested at both tails, pixels can be 0, -1 (significant at
% the lower tail, or 1 (significant at the upper tail).
%
% SCICLUS is an image indicating which pixels are significant at the 
% cluster level. If significance is tested at one tail, pixels can be 0 or 
% 1; if significance is tested at both tails, pixels can be 0, -1 
% (significant at the lower tail, or 1 (significant at the upper tail).
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
if ~exist('displayOpt','var') || isempty(displayOpt), displayOpt = 1; end
if displayOpt~=0 && displayOpt~=1, error('displayOpt must be 0 or 1'); end
displayOpt = logical(displayOpt);
if nDims>3 && displayOpt, warning('Figures will not be displayed. Too many dimensions.'); end
if ~exist('mask','var'), mask = true(CIsize); end

% Pixel-level significance image
SCIpix = [];
if exist('sigPix','var') && ~isempty(sigPix)
    SCIpix = zeros(1,prod(CIsize));
    if iscell(sigPix.ind)
        SCIpix(sigPix.ind{1}) = -1;
        SCIpix(sigPix.ind{2}) = 1;
    else
        SCIpix(sigPix.ind) = 1;
    end
    SCIpix = reshape(SCIpix, CIsize);
end

% Cluster-level significance image
SCIclus = [];
if exist('sigClus','var') && ~isempty(sigClus)
    SCIclus = zeros(1,prod(CIsize));
    if iscell(sigClus.ind)
        SCIclus(sigClus.ind{1}) = -1;
        SCIclus(sigClus.ind{2}) = 1;
    else
        SCIclus(sigClus.ind) = 1;
    end
    SCIclus = reshape(SCIclus, CIsize);
end

% Display figures
if nDims==2 && displayOpt
    if ~isempty(SCIpix), figure, imagesc(SCIpix.*mask), colorbar, title('Pixel test'), end
    if ~isempty(SCIclus), figure, imagesc(SCIclus.*mask), colorbar, title('Cluster test'), end
elseif nDims==3 && displayOpt
    order = [find(CIsize~=min(CIsize)) find(CIsize==min(CIsize))];
    [~, reverseOrder] = sort([find(CIsize~=min(CIsize)) find(CIsize==min(CIsize))]);
    if ~isempty(SCIpix)
        SCIpix = permute(SCIpix, order);
        for pred = 1:size(SCIpix,3)
            figure, imagesc(SCIpix(:,:,pred).*mask(:,:,pred)), colorbar, title(sprintf('Pixel test, Predictor #%i of 3rd dimension', pred))
        end
        SCIpix = permute(SCIpix, reverseOrder);
    end
    if ~isempty(SCIclus)
        for pred = 1:size(SCIclus,3)
            figure, imagesc(SCIclus(:,:,pred).*mask(:,:,pred)), colorbar, title(sprintf('Cluster test, Predictor #%i of 3rd dimension', pred))
        end
        SCIclus = permute(SCIclus, order);
        SCIclus = permute(SCIclus, reverseOrder);
    end
end
