function paintResults(I, annoOne, labelColors)
imshow( I ); 
hold on;
if( size(annoOne,3) == 1 )
    annoOne(1:2,:)       = 0;
    annoOne(end-1:end,:) = 0;
    annoOne(:,1:2)       = 0;
    annoOne(:,end-1:end) = 0;
    for labIter = 1:size(labelColors,1)
        annoTemp(:,:,labIter) = ( annoOne == labIter);
    end
    annoOne = annoTemp;
end

for iter = 1:size(annoOne,3)
    if( sum(annoOne(:,:,iter)) == 0 )
        continue;
    end
    cc = contour(annoOne(:,:,iter), [0.5, 0.5]);
    cc2 = get_contours(cc);
    for ccIter = 1:numel(cc2)
        fill(cc2{ccIter}(1, :), cc2{ccIter}(2, :), labelColors(iter,:), 'FaceAlpha', 0.4);
    end
end
hold off;

axis tight;