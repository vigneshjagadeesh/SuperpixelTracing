function [op S supInfo] = shuffleSupIndFast(S, I, SPrev, Iprev, nodeOffset, supInfo, currSlice, varargin)
USE_RGB = true;
if( nargin == 7 )
    batchIter = currSlice;
else
    batchIter = varargin{1};
end
if( USE_RGB )
    no_bins = 8;%16
    bin_width = 256 / no_bins;    
    % figure(402); imshow(flowCol); pause(.1);
    %% LAB Histogram Parameters
else
    [L1 A1 B1] = RGB2Lab(I);
    I = cat(3, L1, A1, B1);
    LABBINS = 32;
    BL = 100 / LABBINS;
    BA = 255 / LABBINS;
    CL = linspace(BL, 100, LABBINS );
    CA = linspace(-128+BA, 127, LABBINS );
end

% Texture
    lbpCode = padarray( lbp(rgb2gray(I), 1, 8, supInfo(1).mapping, 'i'), [1 1], 'replicate', 'both' );
    
    % Optical Flow .5 15 5  10 7  1.5
    try,
    Vx = zeros(size(I,1), size(I,2));  Vy = zeros(size(I,1), size(I,2));
        %[Vx Vy] = vrl_opencvFlow(double([size(I,1) size(I,2) .5 15 7  10 7  1.5]), rgb2gray(uint8(Iprev)), rgb2gray(uint8(I)));
    catch me
        display('Flow not computed');
        Vx = zeros(size(I,1), size(I,2));
        Vy = zeros(size(I,1), size(I,2));
    end
    Vx( abs(Vx) > 40 ) = 0; Vy( abs(Vy) > 40 ) = 0;
    Vx = medfilt2( Vx, [5 5] ); Vy = medfilt2( Vy, [5 5] );
    Im = (Vx+Vy) > 0;
    [X Y] = meshgrid( 1:size(I,2), 1:size(I,1) );
    
    Xwarp = round(X - Vx); Ywarp = round(Y - Vy);
    Xwarp( Xwarp < 1 ) = 1; Xwarp( Xwarp > size(I,2) ) = size(I,2); Ywarp( Ywarp < 1 ) = 1; Ywarp( Ywarp > size(I,1) ) = size(I,1);
    % warpedIm = vrl_interp2_wrapper(double( rgb2gray(Iprev) ), Xwarp, Ywarp);
    % figure(401); imshow(warpedIm, []); title('Warped Image');
    supWarpIndices = SPrev( sub2ind( [size(I,1) size(I,2)], Ywarp(:), Xwarp(:) ) );
    % flowCol = flowToColor( cat(3, Vx, Vy) );

display(['In Slice Number' num2str(currSlice) ])
I = double(I);
RED = I(:,:,1);
GREEN = I(:,:,2);
BLUE  = I(:,:,3);



prevNodeOffset = nodeOffset - numel( unique(SPrev) );
% Preprocess
Spre = zeros( size(S) );
uniqueInd = unique(S); % Dont want to deal with zeros
newId = 1;
for iter = 1:numel(uniqueInd)
    currMask = S == uniqueInd(iter);
    CC = bwconncomp(currMask);
    for ccIter = 1:CC.NumObjects
        Spre(CC.PixelIdxList{ccIter}) = newId; newId = newId + 1;
    end
end
S = Spre; uniqueInd = unique(S);

SNew = zeros( size(S) );
for iter = 1:numel(uniqueInd)
    SNew( S==uniqueInd(iter) ) = iter;
end
supVis1 = zeros(size(I,1), size(I,2)); supVis2 = zeros(size(I,1), size(I,2)); supVis3 = zeros(size(I,1), size(I,2));
colorCode = rand( numel(uniqueInd), 3).*255;
S = SNew;
%% Neighbors Computation
R = [S(:,2:end) S(:,end)];
L = [S(:,1)     S(:,1:end-1)];
B = [S(2:end,:); S(end,:)];
T = [S(1,:);     S(1:end-1,:)];
NW = S; NW(2:end, 2:end) = S(1:end-1,1:end-1);
NE = S; NE(2:end, 1:end-1) = S(1:end-1, 2:end);
SW = S; SW(1:end-1, 2:end) = S(2:end, 1:end-1);
SE = S; SE(1:end-1, 1:end-1) = S(2:end, 2:end);
Neighbors = [R(:) L(:) B(:) T(:) NW(:) NE(:) SW(:) SE(:)];

%% Superpixel Specific Statistics
for iter = 1:numel(uniqueInd)
    %%% You are now inside a particular superpixel in I(t)
    pickInd = (S==iter);
    [y x]   = find( pickInd);
    supInfo(iter+nodeOffset).center = [mean(x(:)) mean(y(:)) 1];
    %%% Detect Neighbors
    spatNeighbors = unique( Neighbors( pickInd, : ) ) ;
    spatNeighbors( spatNeighbors == iter ) = [];
    spatNeighbors = spatNeighbors + nodeOffset;
    timeNeighbors = [];
    timeNeighborsFlow = [];
    if( batchIter > 1 )
        timeNeighbors = unique( SPrev( pickInd ) )+prevNodeOffset; % across time
        timeNeighborsFlow = unique( supWarpIndices(pickInd) ) +prevNodeOffset ;
        timeNeighbors = unique([timeNeighbors; timeNeighborsFlow]);
    end
    
    try,
        supInfo(iter+nodeOffset).neighbors = [spatNeighbors(:); timeNeighbors ];
        supInfo(iter+nodeOffset).neiflagger = [zeros(numel(spatNeighbors(:)),1); ones(numel(timeNeighbors),1) ];
    catch me,
        me
        spatNeighbors
        timeNeighbors
    end
    
    %%% Construct Features
    supVis1(pickInd) = colorCode(iter,1);
    supVis2(pickInd) = colorCode(iter,2);
    supVis3(pickInd) = colorCode(iter,3);
    supInfo(iter+nodeOffset).area          = sum( pickInd(:) );
    supInfo(iter+nodeOffset).pixelIdx      = [y x currSlice*ones( numel(y) , 1)];
    supInfo(iter+nodeOffset).currPixelIdx  = find( pickInd );
    
    %%%%% LBP Extraction
    supInfo(iter+nodeOffset).lbpHist       = hist(lbpCode(pickInd),0:((supInfo(1).mapping.num)-1)); 
    supInfo(iter+nodeOffset).lbpHist       = supInfo(iter+nodeOffset).lbpHist./sum(supInfo(iter+nodeOffset).lbpHist);
    %supInfo(iter+nodeOffset).hoof          = vrl_hog(Vx, Vy, Im, pickInd, 16);
    
    if( USE_RGB )
% % %         h            = double( ceil( [RED(pickInd) GREEN(pickInd) BLUE(pickInd)] / bin_width ) );
% % %         h(h==0)      = 1;
% % %         color_hist    = zeros(no_bins, no_bins, no_bins);
% % %         try
% % %             for row_iter = 1:size(h,1)
% % %                 color_hist(h(row_iter,1),h(row_iter,2),h(row_iter,3)) = color_hist(h(row_iter,1),h(row_iter,2),h(row_iter,3)) + 1;
% % %             end
% % %         catch me
% % %             color_hist
% % %         end
        singChannel  = I(:,:,1);
        color_hist = histc(singChannel(pickInd), linspace(0, 255, no_bins+1));
        color_hist(end-1) = color_hist(end-1) + color_hist(end);
        color_hist(end) = [];
         color_hist = color_hist ./ sum( color_hist(:) );
         supInfo(iter+nodeOffset).colorHist = color_hist;
        
    else
        LH = hist(RED(pickInd), CL(1:end-1) );
        LA = hist(GREEN(pickInd), CA(1:end-1) );
        LB = hist(BLUE(pickInd), CA(1:end-1) );
        
        color_hist = [LH(:) ./ sum(LH(:)); LA(:) ./ sum(LA(:)); LB(:) ./ sum(LB(:))];
        
        supInfo(iter+nodeOffset).colorHist     = color_hist;
    end
end
op=cat(3, supVis1, supVis2, supVis3);
