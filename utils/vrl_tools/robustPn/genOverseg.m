clear;
close all;
clc;

sliceCtr = 1;

for stackIter = 0:11
    try
        I = adapthisteq( imread(['C:\Users\vignesh\Desktop\RetinaProject\splitMerge\splitMerge0\sm0_' num2str( (stackIter) ) '.bmp']) );
        I = imresize(I, [512 512]);
        preCompOverSeg{1}(:,:,sliceCtr) = vgg_segment_ms( repmat(I, [1 1 3]), 5, 3, 100);
        sliceCtr = sliceCtr + 1;
    catch me
        continue;
    end
end
save('C:\Users\vignesh\Desktop\RetinaProject\splitMerge\splitMerge0\msOverSeg53.mat', 'preCompOverSeg');