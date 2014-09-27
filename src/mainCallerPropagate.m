clear; close all; clc;
codeDir = '/Users/vjagadeesh/Desktop/cleanCodes/hyperDiffClean/SuperpixelTracing/';
addpath( genpath( [codeDir 'data'] ) );
addpath( genpath( [codeDir 'utils/'] ) );
%addpath(          [codeDir 'src/annotateNValidate/'] );

% Initialize Parameters ... 
iter = 0;
noSlices = 5;
VALIDATE_RESULTS = false;
SUP = 1; 
load('genSuperpixels');
IMGSIZE = 512;
noRows = IMGSIZE; 
noCols = IMGSIZE;
mapping = getmapping(8, 'u2');
supInfo = struct( 'mapping', mapping );
dataSet = 1;

% Read Input Data
if( dataSet == 1 )
    MovObj = VideoReader('EMchallenge_down_8.avi');
    MovAnnoFrameOne = imread('I_label_watershed.png');
elseif( dataSet == 2 );
    MovObj = VideoReader('EM_Challenge_top.avi');
    MovAnnoFrameOne = imread('I_miccaiOldInit.png');
end
I =  read(MovObj,1);
%inpDir = 'C:\researchCode\dataBase\emChallenge\514_ipl_00_219\';
%I =  imread([inpDir '514_ipl_' sprintf('%.2d',1) '_' num2str(1+125) '.tif']);
if( size(I,3) == 1), I = cat(3, I, I, I); end
I = imresize(I, [512 512]);
Iprev = I;
reconVideo = zeros( noRows, noCols, noSlices);
redChannel = reconVideo; grnChannel = reconVideo; bluChannel = reconVideo;
SPrev = zeros( size( I ) );
tic,
nodeOffset = 0;

% Graph Construction and Generation of Oversegmentations
tic
COMPUTESEG = false;
if( COMPUTESEG )
    for iter = 1:noSlices
        I =  read(MovObj,iter);%imread([inpDir '514_ipl_' sprintf('%.2d',iter+1) '_' num2str(iter+125+1) '.tif']);
        I =  imresize(I, [IMGSIZE IMGSIZE]) ;
        I = histeq(I(:,:,1));
        I1 = I;
        I = medfilt2( I, [3 3] );
        % I = anisodiff(I, 10, 20, .25, 1);
        IFull(:,:,iter) = I;
        figure(400); imshow( uint8(I) );
        if( size(I,3) == 1), I = cat(3, I, I, I); end;
        if( size(I1,3) == 1), I1 = cat(3, I1, I1, I1); end
        %% The list of superpixel generators ...
        if( SUP == 1)
            currSup{iter} = vgg_segment_gb( uint8(I), .8, 100, 100);% was 25 for surfer % was 200 200  ...   let this be set to a constant value, or a value dependent on the stroke
        elseif( SUP == 2 )
            temp = uint8(I);
            %[to currSup] = edison_wrapper(temp,@RGB2Luv, 'SpatialBandWidth', 9, 'RangeBandwidth', 5, 'MinimumRegionArea', 500);
            currSup{iter} = vgg_segment_ms( temp, 8, 5, 100); % was 25 before .. constant or value dependent on the stroke ...
            currSup{iter} = double( currSup{iter} );
        elseif( SUP == 3 )
            currSup{iter} = slicSuperpixels( uint8(I), 400, 20);
        end
        [supVis supSeg{iter} supInfo] = shuffleSupIndFast( currSup{iter}, I1, SPrev, Iprev, nodeOffset, supInfo, iter ); % Superpixel Statistics
        figure(401); imshow( uint8( supVis ) );
        nodeOffset = nodeOffset + numel( unique( supSeg{iter} ) )
        SPrev = supSeg{iter};
        Iprev = I;
    end
    %save('genSuperpixels.mat', 'currSup', 'supInfo', 'supSeg', 'IFull');
    toc
else
    load('genSuperpixels.mat');
end
[A y adjMat numFg numBg] = BU_constructAdjHKNN(supInfo, supSeg, read(MovObj,1), MovAnnoFrameOne);

numBg
numFg

%% Actual Cut
f = partitionHypergraph(A, y);
%[f compInv] = runLGC(y, adjMat, 0);

[~,infLabels] = max(f,[],2);

%% Reconstruct Labelled Video
infLabels( infLabels > numFg ) = numFg + 1;
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
toc
if( VALIDATE_RESULTS )
    FMeas = zeros( numFg,1 );
    for targetIter = 1:numFg
        if( mod( targetIter, 10) == 0 ), display([num2str(targetIter)]); end
        currTar = reconVideo == targetIter;
        if( dataSet == 1 )
            [masker FMeas(targetIter)] = retrieveAnnotation(targetIter+1, currTar, 0);
        elseif( dataSet == 2 )
            [masker FMeas(targetIter)] = retrieveMiccai11Anno(targetIter, currTar, 0);
        end
    end
    plot( 1:numel(FMeas), FMeas, 'r' ); title(['Median F Measure is ' num2str( median(FMeas) ) ]);
    for iter = 1:noSlices
        if( dataSet == 1 )
            RAND_idx(iter) = retrieveRANDindex(reconVideo(:,:,iter), iter, MovObj);
        end
    end
end

Img = imread('I00000_label.tif');
for num_image = 1:noSlices
    figure(5*num_image); imshow( uint8( IFull(:,:,num_image) ) ); axis off;
    if( num_image == 1 )
        wShed = watershed( Img );
        hold on; contour( wShed, [0 0], 'r' ); hold off;
    end
    figure(6*num_image); imshow( uint8( cat(3, redChannel(:,:,num_image), grnChannel(:,:,num_image), bluChannel(:,:,num_image) ) ) ); pause(.1); axis('image');axis off;
    pause;
end

uniqueLabs = unique(reconVideo(:));
IVis = zeros( size(reconVideo,1), size(reconVideo,2), size(reconVideo,3)+2);
IVis(:,:,2:end-1) = reconVideo;

close all;
labelColors = rand( numFg, 3 );

vidObj = VideoWriter('tracedEM.avi'); vidObj.FrameRate = 1; open(vidObj);
fig = figure(1);
for sliceIter = 1:noSlices
    I = MovObj.read(sliceIter);
    currAnno = reconVideo(:,:,sliceIter); currAnno( currAnno == 95 ) = 0;
    paintResults(imresize(I,[512 512]), currAnno, labelColors);
    pause(1);
    print ('-f1', '-dpng', ['sampleFig' num2str(sliceIter) '.png']);
    currFrame = getframe;
    writeVideo(vidObj,currFrame);
end
close(vidObj);

vidObj = VideoWriter('tracedEM1.avi'); vidObj.FrameRate = 1; open(vidObj);
for iter = 1:noSlices
I = imread(['sampleFig' num2str(iter) '.png']);
writeVideo(vidObj, I);
end
close(vidObj);
return;
figure(100); axis([ 1 size(IVis,2) 1 size(IVis,1) 1 size(IVis,3) ]);
for labIter = 3:numel(uniqueLabs)-3
    visMask = (IVis == labIter);
    %isosurface( visMask); pause;
    p = patch(isosurface(visMask));
    set(p, 'FaceColor', [rand rand rand], 'EdgeColor', 'none');   pause;
end