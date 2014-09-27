function [op S supInfo] = shuffleSupInd(S, I, SPrev, SPrevPrev, nodeOffset, supInfo, currSlice)
DEBUG = false;

I = im2double(I);
prevNodeOffset = nodeOffset - numel( unique(SPrev) );
prevprevNodeOffset = prevNodeOffset - numel( unique(SPrevPrev) );
uniqueInd = unique(S); % Dont want to deal with zeros
S_out = zeros(size(S));
SNew = zeros( size(S) );
for iter = 1:numel(uniqueInd)    
    SNew( S==uniqueInd(iter) ) = iter;
end
S = SNew; clear SNew;
if(DEBUG), clc; end

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
    currNode = iter+nodeOffset;
    pickInd  = (S==iter);
    [y x]    = find( pickInd);
    S_out(pickInd) = mean( I(pickInd) ); %randomInd(iter, 1);
%     S1(pickInd) = mean( R(pickInd) ); %randomInd(iter, 1);
%     S2(pickInd) = mean( G(pickInd) ); %randomInd(iter, 2);
%     S3(pickInd) = mean( B(pickInd) ); %randomInd(iter, 3);    
    %%% Detect Neighbors
    spatNeighbors = unique( Neighbors( pickInd, : ) ) + nodeOffset;
    spatNeighbors = setdiff(spatNeighbors,currNode);
    timeNeighbors = [];
    
    if currSlice > 2
        timeNeighbors = unique( SPrevPrev( pickInd ) ) + prevprevNodeOffset;                    % across time t and t-2
        timeNeighbors = cat(1, timeNeighbors, unique( SPrev( pickInd ) ) + prevNodeOffset );    % across time t and t-1
    elseif( currSlice == 2 )
%     if currSlice > 1
        timeNeighbors = unique( SPrev( pickInd ) ) + prevNodeOffset; % across time t and t-1
%     end
        if( DEBUG ), timeNeighbors,  end
    end

    supInfo(currNode).neighbors = [spatNeighbors(:); timeNeighbors];
    
    %%% Construct Features
    sprintf('');
    supInfo(currNode).area          = sum( pickInd(:) );
    supInfo(currNode).pixelIdx      = [y x currSlice*ones( numel(y) , 1)];
    supInfo(currNode).center        = mean(supInfo(currNode).pixelIdx,1);
    supInfo(currNode).currPixelIdx  = find( pickInd ); 
    supInfo(currNode).colorHist     = vrl_grayhist(I, pickInd, 16);
    % supInfo(currNode).colorHist     = vrl_colorhist(I, pickInd, 8);
    % supInfo(currNode).hog           = vrl_hog(Ix, Iy, I_roi, 8);
    % supInfo(currNode).supShape      = vrl_shapeMoments(pickInd);
end

if(DEBUG), figure; subplot(211); imagesc(SPrev); subplot(212); imagesc(S); end
% op = cat(3, S1,S2,S3);
op = S_out;
end