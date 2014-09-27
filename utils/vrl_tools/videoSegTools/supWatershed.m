% Superpixels with the watershed transform
function L = supWatershed(I)
tic,
if( nargin < 1 )
    clear; close all; clc;
    %I = imread('peppers.png');
    % I = 255-imread('simpSquare.jpg'); I = cat(3, I, I, I);
    I = imread('frame17.jpg');
end
    
gradVar = .01;
I(:,:,1) = medfilt2(I(:,:,1), [5 5]); I(:,:,2) = medfilt2(I(:,:,2), [5 5]); I(:,:,3) = medfilt2(I(:,:,3), [5 5]);
I = im2double(I);

% Filter and Extract Gradient
If = imfilter( I, fspecial('gauss', [7 7], 1.2), 'conv', 'same' );
[Ix1 Iy1] = gradient(If(:,:,1)); 
[Ix2 Iy2] = gradient(If(:,:,2));
[Ix3 Iy3] = gradient(If(:,:,3));
Ix = max( cat(3, Ix1, Ix2, Ix3), [], 3 );
Iy = max( cat(3, Iy1, Iy2, Iy3), [], 3 );
gradMag = sqrt( Ix.^2 + Iy.^2 );

% Nonlinearity on Gradient
edgeInd = 1 - exp( -gradMag ./ gradVar );

% Threshold the Gradient
edgeInd( edgeInd > .9 ) = 1; edgeInd( edgeInd <=.9 ) = 0;

% % % % % Remove Border Artifact

% % % % edgeInd(1:2, :) = 0; edgeInd(:, 1:2) = 0;
% % % % edgeInd(end-1:end, :) = 0; edgeInd(:, end-1:end) = 0;

% Clean the Edge Map
CC = bwconncomp(edgeInd);
numPixels = cellfun(@numel,CC.PixelIdxList);
idx = find( numPixels < 300 );
for connIter = 1:numel(idx)
    edgeInd(CC.PixelIdxList{idx(connIter)}) = 0;
end


% Perform Watershed on Distance Transform
distMap = -bwdist( edgeInd );
initDistMap = distMap;
rmin = imregionalmin(distMap);
minMap = imdilate(rmin,strel('disk',2));
%figure; imagesc(minMap);
CC = bwconncomp(minMap);
for connIter = 1:CC.NumObjects
    for pixIter = 1:numel( CC.PixelIdxList{connIter} )
        distMap(CC.PixelIdxList{connIter}(pixIter)) = -10^6 + pixIter;
    end
end


L = watershed(distMap,8);
L = imdilate(L, strel('disk', 2) );
supVis = label2rgb(L);

if( nargin < 1 )
figure; imshow(I);
figure; imagesc(distMap);
figure; imagesc( supVis );
return;
end
toc