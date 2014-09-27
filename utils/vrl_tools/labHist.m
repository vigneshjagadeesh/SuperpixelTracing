function fullHist = labHist(I, bins, pickInd)

%% I is the input image in the LAB color space
%% bins has the number of bins in the L, A and B color channels

L = I(:,:,1); A = I(:,:,2); B = I(:,:,3);

BL = 100 / bins(1);
BA = 255 / bins(2);
CL = linspace(BL, 100, bins(1) );
CA = linspace(-128+BA, 127, bins(2) );

LH = hist(L(pickInd), CL(1:end-1) );
LA = hist(A(pickInd), CA(1:end-1) );
LB = hist(B(pickInd), CA(1:end-1) );

LH = LH ./ sum(LH(:));
LA = LA ./ sum(LA(:));
LB = LB ./ sum(LB(:));

fullHist = [LH(:); LA(:); LB(:)];