clear; close all; clc; addpath(genpath('../../'));

if( ismac )
    inpDir = '/Users/vigneshjagadeesh/Dropbox/BEL_testing_histeq/belOutputs/';
    addpath( genpath( '/research/matlab_bgl/' ) );
else
    inpDir = 'C:\researchCode\dataBase\emChallenge\514_ipl_00_219\belOutputs\';
end

mapping = getmapping(8, 'u2');
supInfo = struct( 'mapping', mapping );
nodeOffset = 0;
nodeSpacer = [];
SPrev = zeros(512, 512);
noSlices = 5;
for iter = 1:noSlices
    fileNum = sprintf('%.2d', iter);
    fileName = [inpDir 'I000' fileNum '.tif'];
    fileNameb = [inpDir 'I000' fileNum '_prob.tif'];
    Ib = imread( fileNameb );
    Ib = 1-im2double( Ib );
    
    I = imread( fileName );
    IFull(:,:,iter) = I;
    Is = imfilter( double(I), fspecial('gauss', [11 11], 3) );
    [Ix Iy] = gradient(Is);
    gradMag = sqrt(Ix.^2 + Iy.^2);
    edgeMap = exp( - gradMag / 10 );
    
    %edgeMap = imfill(edgeMap, 'holes');
    
    Ib = medfilt2( Ib, [3 3] );
    
    synMap = edgeMap .* Ib;
    bwMap = 1-adaptivethresh(synMap);
    bwMap = imdilate(bwMap, ones(5,5));
    CC = bwconncomp(bwMap);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    idx = find( numPixels < 250 );
    for connIter = idx, bwMap(CC.PixelIdxList{connIter}) = 0; end
    %bwMap = imopen(bwMap, strel('disk', 2) );
    
    distMap = bwdist( bwMap );
    maxPts = distMap > imdilate(distMap, [1 1 1; 1 0 1; 1 1 1]);
    %maxPts = imextendedmax(distMap, 30);
    distMap = max( distMap(:) ) - distMap;
    % distMap( maxPts ) = -Inf;
    segRes = watershed( bwMap );
    segRes = imdilate( segRes, ones(3,3) );
    [op S supInfo] = shuffleSupIndFast(segRes, cat(3, I, I, I), SPrev, [], nodeOffset, supInfo, iter);
    SPrev = segRes;
    figure(iter); imshow( uint8( op ) ) ;pause(.5);
    supFull(:,:,iter) = segRes + nodeOffset;
    nodeOffset = nodeOffset + numel( unique( segRes ) )
    nodeSpacer = [nodeSpacer; nodeOffset];
end

nodeOffset = numel(supInfo);
adjMatDiv = zeros(nodeOffset);
adjMatCov = zeros(nodeOffset);
divg = zeros(nodeOffset);
covrage = zeros(nodeOffset);
for supIter = 1:nodeOffset
    for neIter = 1:numel( supInfo(supIter).neighbors )
        currN = supInfo(supIter).neighbors(neIter);
        areaRatio = abs( supInfo(supIter).area / supInfo(currN).area - 1 );
        interA = numel(intersect(supInfo(supIter).currPixelIdx, supInfo(currN).currPixelIdx));
        unionA = numel(union    (supInfo(supIter).currPixelIdx, supInfo(currN).currPixelIdx));
        if( supInfo(supIter).neiflagger(neIter) )
            temp1 = sum( min( supInfo(supIter).colorHist(:), supInfo(currN).colorHist(:))) ;
            temp2 = sum( min( supInfo(supIter).lbpHist(:),   supInfo(currN).lbpHist(:)  )) ;
            divg( supIter, currN ) = 0;      covrage( supIter, currN ) = 0;
            divg( currN, supIter ) = temp1 + temp2 ; covrage( currN, supIter ) = interA / unionA;
        else
            %(1 - sum( min( supInfo(supIter).colorHist(:), supInfo(currN).colorHist(:))));
            divg( supIter, currN ) = 0;      covrage( supIter, currN ) = 0;
            divg( currN, supIter ) = 0;      covrage(currN, supIter) = 0;
        end
    end
end
sigmaDiv = median(divg(divg~=0));  sigmaCov = median(covrage(covrage~=0));
adjMatDiv(divg~=0)    = exp(-divg(divg~=0)/sigmaDiv);
adjMatCov(covrage~=0) = exp(-covrage(covrage~=0)/sigmaCov);
adjMat = adjMatDiv + adjMatCov;
adjMat = sparse( adjMat );

%adjMat = sparse( divg );

[D P] = floyd_warshall_all_sp(adjMat);
D = D + 10^6 * eye( size(D) );
allSp = cell(10,10);
visMask = zeros(size(I,1), size(I,2), noSlices+2);
tarCtr = 1;
excepCount = 0;
firstSup = supFull(:,:,1); figure(102); imshow(IFull(:,:,1)); [XQ,YQ] = ginput; queryIndices = unique(firstSup(sub2ind(size(firstSup),ceil(YQ), ceil(XQ))));
firstSup( supInfo(queryIndices(1)).currPixelIdx ) = -100;   figure(103); imagesc(firstSup);
hyperedgeList = []; 
for queryIter = 1:nodeSpacer(1)  %1:numel(queryIndices) %
    tarIter = queryIter;         %queryIndices( queryIter );
    currCtr = 1;
    pickLastSlice = [10^6*ones(1, nodeSpacer(end-1)) ones(1,nodeSpacer(end)-nodeSpacer(end-1)) ];
    [minVal minInd] = min( pickLastSlice .* D(tarIter,:) );
    for currNode = nodeSpacer(end-1)+1:nodeSpacer(end)
        predNode = currNode;
        currPath = [currNode];
        while( predNode ~= tarIter )
            predNode1 = P(tarIter, predNode)
            if( predNode1 == predNode )
                excepCount = excepCount + 1;
                display('caught exception');
                break;
            else
                predNode = predNode1;
            end
            currPath = [currPath; predNode];
        end
        allSp{currCtr, tarCtr} = currPath;
        if( ~isempty(currPath) && numel(currPath) == noSlices )
            if(currNode == minInd)
                hyperedgeList = [hyperedgeList currPath];
                mask = zeros(size(I,1), size(I,2), noSlices);
                for sliceIter = 1:numel(currPath)
                    tempMask = zeros( size(I,1), size(I,2));
                    tempMask( supInfo(currPath(sliceIter)).currPixelIdx ) = 255;
                    mask(:,:,noSlices-sliceIter+1) = tempMask;
                    figure(99); subplot(3,2,noSlices-sliceIter+1); imshow(tempMask); title([num2str(noSlices-sliceIter+1)]);
                end
                visMask(:,:,2:end-1) = mask;
                hyperedgeMask{queryIter} = mask;
                
                figure(100); axis([1 size(I,2) 1 size(I,1) 1 noSlices]); title('SHORTEST PATH');
                % p = patch(isosurface(visMask));
                % set(p, 'FaceColor', [rand rand rand], 'EdgeColor', 'none');  % pause; 
            end
        end
        currCtr = currCtr+1;
    end
    tarCtr = tarCtr + 1;
end
excepCount
for hIter = 1:numel( hyperedgeMask )
    currMask = hyperedgeMask{hIter};
    save(['currHyperedges/topDownHyperEdges' num2str(hIter) '.mat'], 'currMask' );
end
    return;


%% Reconstruct Labelled Video
uniqueInds = unique(infLabels);
colorMapper = rand( numel(uniqueInds), 3 );
dispImg = zeros( noRows, noCols, noSlices );
for supIter = 1:numel( infLabels )
    currInd = supInfo(supIter).pixelIdx;
    % for pixIter = 1:size(currInd,1), reconVideo( currInd(pixIter,1), currInd(pixIter,2), currInd(pixIter,3) ) = infLabels( supIter ); end
    reconVideo( noRows*noCols*(currInd(:,3)-1) + noRows*(currInd(:,2)-1) + currInd(:,1) ) = infLabels(supIter);
    redChannel( noRows*noCols*(currInd(:,3)-1) + noRows*(currInd(:,2)-1) + currInd(:,1) ) = 255.*colorMapper( find( uniqueInds == infLabels(supIter)) ,1);
    grnChannel( noRows*noCols*(currInd(:,3)-1) + noRows*(currInd(:,2)-1) + currInd(:,1) ) = 255.*colorMapper( find( uniqueInds == infLabels(supIter)) ,2);
    bluChannel( noRows*noCols*(currInd(:,3)-1) + noRows*(currInd(:,2)-1) + currInd(:,1) ) = 255.*colorMapper( find( uniqueInds == infLabels(supIter)) ,3);
end

for num_image = 1:noSlices
    figure(5); imshow( uint8( IFull(:,:,num_image) ) ); axis off;
    figure(6); imshow( uint8( cat(3, redChannel(:,:,num_image), grnChannel(:,:,num_image), bluChannel(:,:,num_image) ) ) ); pause(.1); axis('image');axis off;
    pause;
end

uniqueLabs = unique(reconVideo(:));
IVis = zeros( size(reconVideo,1), size(reconVideo,2), size(reconVideo,3)+2);
IVis(:,:,2:end-1) = reconVideo;
close all; figure(100); axis([ 1 size(IVis,2) 1 size(IVis,1) 1 size(IVis,3) ]);
for labIter = 3:numel(uniqueLabs)-3
    visMask = (IVis == labIter);
    %isosurface( visMask); pause;
    p = patch(isosurface(visMask));
    set(p, 'FaceColor', [rand rand rand], 'EdgeColor', 'none');   pause(0.01);
end