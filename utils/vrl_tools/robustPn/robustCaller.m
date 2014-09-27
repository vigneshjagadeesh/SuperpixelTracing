clear; close all; clc;
addpath('C:\researchCode\codebase\vrl_tools\');
% H = adapthisteq( imread('C:\Users\vignesh\Desktop\RetinaProject\splitMerge\splitMerge0\sm0_0.bmp') );
H = ( imread('C:\researchCode\dataBase\multiCamera\epfl\terrace1\ter2\ter2_510.png') );

[noRows noCols noSlices] = size(H);
numObjects = 2;
objColor = [1 0 0; 0 1 0; 0 0 1];
summer = zeros(noRows, noCols);
% % for objIter = 1:numObjects
% %     userMask{objIter} = roipoly(H); 
% %     summer = summer + userMask{objIter};
% % end
% % userMask{objIter+1} = summer == 0 ;
% % save( 'initMask.mat', 'userMask');
load('initMask.mat');

%% Higher Order Potentials
msp = [ 3 2; 6 4];
gsp = [.8 50; .8 15; .8 100;.8 150];
hocCost = 100;
nlSigma = 2.2;
intCost =  2500;
Qhoc = .2;
for stackIter = 1:50
        % I = adapthisteq( imread( ['C:\Users\vignesh\Desktop\RetinaProject\splitMerge\splitMerge0\sm0_' num2str(stackIter) '.bmp'] ) );
        I = ( imread(['C:\researchCode\dataBase\multiCamera\epfl\terrace1\ter2\ter2_' num2str( (2*stackIter) +511) '.png']) );  
        if( size(I,3) ==  3)
            grayI = rgb2gray(I);
        else
            grayI = I;
        end
        [Vx Vy] = openCvFlow((int32([noRows noCols])), [0.5, 1, 5, 5, 5, 1.5, 0], double(rgb2gray(H)), double(rgb2gray(I)) );
        Vx = reshape(Vx, [noRows noCols]);         Vy = reshape(Vy, [noRows noCols]);
        optFlowMag = sqrt(Vx.^2 + Vy.^2);
        
        
        for supIter = 1:size(msp,1)
            if( size(I,3) == 1)
                segRes{supIter} = vgg_segment_ms( repmat(I, [1 1 3]), msp(supIter,1), msp(supIter,2), 40);
                % segRes{supIter} = vgg_segment_gb( repmat(I, [1 1 3]), gsp(supIter,1), gsp(supIter,2), 100);
            else
                segRes{supIter} = vgg_segment_ms( I, msp(supIter,1), msp(supIter,2), 40);
                % segRes{supIter} = vgg_segment_gb( I, gsp(supIter,1), gsp(supIter,2), 15);
            end
            figure(supIter);
            imagesc( segRes{supIter} );
        end
tbeta = 25; tv = 4.5; talpha = .8; tp = 1;
%% Unary Potential
    unary = [];
    for objIter = 1:numObjects+1
        Ifg{objIter} = find(userMask{objIter} == 1);
        if( size(I, 3) == 1)
            fgHist{objIter} = vrl_grayhist( H, Ifg{objIter}, 8);
            fgProject{objIter} = -log(vrl_grayhistbp( I, fgHist{objIter} )+0.05).*0 ;
        else
            fgHist{objIter} = vrl_colorhist( H, Ifg{objIter}, 8);
            if( objIter == numObjects + 1 )
                fgProject{objIter} = -log( vrl_colorhistbp( I, fgHist{objIter} ) + 0.05 )  ;
            else
                fgProject{objIter} = -log(vrl_colorhistbp( I, fgHist{objIter} ) .* (optFlowMag > 3) .* exp(-bwdist(userMask{objIter})/15)  +0.05)  ;
            end
        end
        unary = [unary fgProject{objIter}(:)];
    end
    unary = unary';
%% Interaction Potential
imSmooth = imfilter(double(I), fspecial('gauss', [15 15], 3));
im_test_smooth = double(reshape(I, noRows*noCols, noSlices));
[edges, no_nlinks] = vrl_Construct2DLattice([noRows noCols],1);
[weights,w_dist] = vrl_edgeweight(edges, im_test_smooth', [noRows noCols], nlSigma ); 
edges = edges';
vrl_visualize_interaction_map(imSmooth, nlSigma, 2);
% PnWrapper( [2 100 1000 10], [1 2 3 4], [1 2 3 4], [1 2 3 4]);

ctr = 1;
for supIter = 1:numel(segRes)
    allLabels = unique(segRes{supIter});
    for iter = 1:numel( unique(segRes{supIter}) )
        pixelIDs{ctr} = find( segRes{supIter} == allLabels(iter) );
        supStruct(ctr).numSp = numel( pixelIDs{iter} );
        supStruct(ctr).pixList = pixelIDs{iter}-1;
        
        
        currPixels = grayI( pixelIDs{iter} );
        varEstimate = std( double( currPixels ) );
        totalOverlap = 0;
        for objIter = 1:numObjects
            numOverlap = sum( userMask{objIter}( pixelIDs{iter} ) );
            totalOverlap = totalOverlap + numOverlap;
            supStruct(ctr).labelCost(objIter) = (1/varEstimate) .* ( 1 + hocCost.* ( 1-(numOverlap/numel(currPixels)) ) );            
        end
        supStruct(ctr).labelCost(numObjects+1) = (1/varEstimate) .* ( 1 + hocCost.* ( (totalOverlap/numel(currPixels)) ) );
        supStruct(ctr).labelCost(numObjects+2) = (1/varEstimate) .* ( 1 + hocCost );
        supStruct(ctr).labelCost(numObjects+3) = Qhoc;

        % % % % %                 currPixels = fgProject( pixelIDs{iter} );
        % % % % %                 varEstimate =  sum( (currPixels - mean(currPixels)).^2 ) / numel(currPixels) ;
        % % % % %                 G = exp( - tbeta * varEstimate);
        % % % % %                 supStruct(ctr).labelCost = [0 0];
        % % % % %                 supStruct(ctr).gMax = numel(currPixels)^talpha * (tp + tv * G );
        ctr = ctr + 1;
    end
end

segFg = PnWrapper( int32([numObjects+1, noRows*noCols, numel(weights)]), ...
                                            double( unary(:) ),  double( (edges(:)-1)' ), double(  intCost * ( weights)' ), ...
                                            supStruct);
segFg = reshape( segFg, [size(I,1) size(I,2)] );
figure(100); imshow(I);                      
for objIter = 1:numObjects+1
    userMask{objIter} = segFg == objIter-1;
    if(objIter ~= (numObjects+1) )
        hold on; contour( userMask{objIter}-0.5 , [0 0], 'Color', objColor(objIter, :));   hold off;
    end
end
H = I;
end