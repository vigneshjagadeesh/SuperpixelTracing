function [CCbg numBg] = modelBg(bgMask, erodeFlag)

%% Model Background

masker = bgMask;
masker = imerode(masker, strel('disk',7));
CCbg = bwconncomp(masker);
ccSize = cellfun(@numel,CCbg.PixelIdxList);
[ccSize sortedInds] = sort( ccSize, 'descend');
%CCbg.PixelIdxList(ccSize<50)=[];
%CCbg.NumObjects = CCbg.NumObjects - sum(ccSize<50);
for iter = 5:CCbg.NumObjects
    masker( CCbg.PixelIdxList{sortedInds(iter)} ) = 0;
end
if( erodeFlag )
    masker = imerode( masker, strel('disk', 4) );
end
CCbg = bwconncomp(masker);
figure(100); imagesc( masker ); title('Background Mask'); pause; close(100);
numBg = CCbg.NumObjects;
