function [A y adjMat numFg numBg] = constructAdjHKNN(supInfo, supSeg, I, MovAnnoFrameOne)
%% Construct A
[noRows noCols ns] = size( supSeg{1} );
noSlices = numel( supSeg );
nodeOffset = numel(supInfo);
hypConst = 0.001;
divCol = zeros(nodeOffset); adjMatCol = zeros(nodeOffset);
divTex = zeros(nodeOffset); adjMatTex = zeros(nodeOffset);
covrage = zeros(nodeOffset); adjMatCov = zeros(nodeOffset);
for supIter = 1:nodeOffset
    currSlice = supInfo(supIter).pixelIdx(1,3);
    for neIter = 1:numel( supInfo(supIter).neighbors )
        currN = supInfo(supIter).neighbors(neIter);
        neigSlice = supInfo(currN).pixelIdx(1,3);
        num = numel(intersect(supInfo(supIter).currPixelIdx, supInfo(currN).currPixelIdx));
        %den = numel(union(supInfo(supIter).currPixelIdx, supInfo(currN).currPixelIdx));
        den = min( numel(supInfo(supIter).currPixelIdx), numel(supInfo(currN).currPixelIdx) );
        covrage( supIter, currN ) = 2 * ( 1 - (num / den) ); % If high then no overlap since intersection is very low
        
        % divTex( supIter, currN ) = (1 - sum( min( supInfo(supIter).lbpHist(:), supInfo(currN).lbpHist(:) ) ));
        % divCol( supIter, currN ) = (1 - sum( min( supInfo(supIter).colorHist(:), supInfo(currN).colorHist(:)))) ;
        if( currSlice == neigSlice )
            divCol( supIter, currN ) = kld(supInfo(supIter).colorHist(:),supInfo(currN).colorHist(:),'sym');
            divTex( supIter, currN ) = kld(supInfo(supIter).lbpHist(:),supInfo(currN).lbpHist(:),'sym');
        else
            %%%%%%%%%%%%% REMOVE IS NEED BE covrage( supIter, currN )
            divCol( supIter, currN ) = kld(supInfo(supIter).colorHist(:),supInfo(currN).colorHist(:),'sym')*covrage( supIter, currN );
            divTex( supIter, currN ) = kld(supInfo(supIter).lbpHist(:),supInfo(currN).lbpHist(:),'sym')*covrage( supIter, currN );
        end
    end
end
% sigmaCol = 0.01; sigmaTex = 0.01; sigmaCov = 0.01;
adjMatCol(divCol~=0)    = exp(-divCol(divCol~=0));
adjMatTex(divTex~=0)    = exp(-divTex(divTex~=0));
adjMatCov(covrage~=0) = exp(-covrage(covrage~=0));
adjMat = adjMatCol .* adjMatTex;
adjMat = max( adjMat, adjMat');

%% Construct Hyperedges Using K Nearest Neighbors
numNode = size( adjMat, 1);
KNN = 5;
A = zeros(numNode,numNode);
colSpacer = linspace(25,200,KNN);
for iNode = 1 : numNode
    masker = zeros(noRows, noCols, noSlices);
    pixelMarks = supInfo(iNode).pixelIdx; twodPts = sub2ind([noRows noCols], pixelMarks(:,1), pixelMarks(:,2) ); masker( noRows*noCols*(pixelMarks(1,3)-1) + twodPts) = 255;
    [currNodeDistances,currNodeIndices] = sort(adjMat(iNode,:),'descend');
    currNodeDistances = currNodeDistances(1:KNN);
    currNodeIndices   = currNodeIndices(1:KNN);
    currNodeDistances = currNodeDistances( currNodeDistances ~= 0 );
    currNodeIndices = currNodeIndices( currNodeDistances ~= 0 );
    for neIter = 1:numel(currNodeIndices)
        pixelMarks = supInfo(currNodeIndices(neIter)).pixelIdx; twodPts = sub2ind([noRows noCols], pixelMarks(:,1), pixelMarks(:,2) ); masker( noRows*noCols*(pixelMarks(1,3)-1) + twodPts) = colSpacer(neIter);
    end
    A(iNode,currNodeIndices) = currNodeDistances;
    %for sliceIter = 1:noSlices
    %    subplot(ceil(noSlices/2), 2, sliceIter); imagesc( masker(:,:,sliceIter) );
    %end
    %pause(2);
end

%% Construct Hyperedges using Top Down Information for Longer Range Interactions across the z-direction
hPath = '/Users/vignesh/Desktop/hypClean/hypergraphClean/data/currHyperedges/';
hFiles = dir( [hPath 'topDown*'] );
fullSupSeg = supSeg{1}; nodeOffset = max( fullSupSeg(:) );
for sliceIter = 2:numel(supSeg), fullSupSeg = cat(3, fullSupSeg, supSeg{sliceIter}+nodeOffset ); nodeOffset = max( fullSupSeg(:) ); end  
%topDownHyperEdges = zeros(numel(hFiles), numNode);
for hIter = 1:numel( hFiles )
    % For each mat file find superpixels with intersection
    load([hPath hFiles(hIter).name]); % loaded file with this name
    hyperNodes = unique( fullSupSeg( currMask~= 0 ) );
    hfind = find( hyperNodes < max(supSeg{1}(:)) ); hfind = hfind(1);
    hyperNodes(hyperNodes==hfind) = [];% Add to existing hyperedges
    for nodeIter = 1:numel(hyperNodes)
        currN = hyperNodes( nodeIter );
        colDist = kld(supInfo(hfind).colorHist(:),supInfo(currN).colorHist(:),'sym');
        texDist = kld(supInfo(hfind).lbpHist(:),supInfo(currN).lbpHist(:),'sym');
        A( hfind, hyperNodes( nodeIter) ) = hypConst * exp( -colDist ) * exp( -texDist );
        %topDownHyperEdges(hIter, hyperNodes( nodeIter) ) = exp( -colDist ) * exp( -texDist );
    end
end
%A = [A; topDownHyperEdges];
A = max(A, A');

currSeg = supSeg{1};
if( 0  == 1 )
    %% Unsupervised Initialization of y 1
    numInitSuperpixels = numel( unique( currSeg ) );
    if( 1 == 0 )
        supPicker = ceil ( linspace(20, numInitSuperpixels-20, 100 ) );
        Img = zeros( size(currSeg,1), size(currSeg,2) );
        for objIter = 1:numel(supPicker)
            Img( currSeg == supPicker(objIter) ) = objIter;
        end
    else
        [supPicker tempo Img] = naive2Initialize( I, 100 );
    end
    uniqueIDs = unique(Img); uniqueIDs( uniqueIDs==0 ) = [];
    numFg = numel(unique(Img))-1;
    %[CCbg numBg] = modelusBg(Img==0, 1);
    numBg = 0;
else
    %% Heuristic Initialization of y
    %Img = imread('I_label_watershed.png');
    Img = MovAnnoFrameOne;
    uniqueIDs = unique(Img); uniqueIDs( uniqueIDs==0 ) = [];
    numFg = numel(unique(Img))-1;
    [CCbg numBg] = modelBg(Img==0, 0);
    
end


y = zeros(numNode,numFg+numBg);
for labelIter = 1 : numFg
    overlapSuperpixels = currSeg(Img== uniqueIDs(labelIter) );
    y(overlapSuperpixels,labelIter) = 1;
end
for labelIter = (numFg+1):(numFg+numBg)
    overlapSuperpixels = currSeg(CCbg.PixelIdxList{labelIter-numFg});
    y(overlapSuperpixels,labelIter) = 1;
end