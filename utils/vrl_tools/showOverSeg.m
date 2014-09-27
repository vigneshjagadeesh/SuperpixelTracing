function showOverSeg( img, sp, figHand, supStruct)

if( nargin == 0 )
    img = imread('C:\researchCode\dataBase\sample_tracks\stack7\VikingFrame_2.bmp');
    img = adapthisteq( imresize(img, [512 512]) );
    sp = vgg_segment_ms( repmat(img, [1 1 3]), 5, 3, 100 );
    figHand = 1;
end

meanImg = zeros( size(img, 1), size(img, 2) );

[Ix Iy] = gradient(double(sp));
gradMag = sqrt( Ix.^2 + Iy.^2 );
gradMag = gradMag > 0;
% gradMag = imerode( gradMag, strel('disk',2) );
for supIter = 1:numel(unique(sp))
    pixId = sp == supIter;
    meanVal = mean( img(pixId) );
    meanImg( pixId ) = meanVal;
end

newMask = meanImg;
newMask( gradMag ) = 255;
newMask(:,:,2) = meanImg;
newMask(:,:,3) = meanImg;
figure(figHand);
imshow( uint8(newMask) );


if( nargin > 3 ) 
    numLabels = numel(supStruct(1).labelCost)-1;
    %% Display all the Superpixel Costs
        for supIter = 1:numel(unique(sp))
            pixId = sp == supIter;
            for labIter = 1:numLabels % This is the number of labels+1 for background
                labImg{labIter}( pixId ) = supStruct(supIter).labelCost(labIter);
            end
        end
        for labIter = 1:numLabels-1
            subplot( ceil((numLabels-1)/2), 2, labIter );
            imagesc( reshape(labImg{labIter}(:), [size(img,1) size(img,2)] ) );
        end
    
end