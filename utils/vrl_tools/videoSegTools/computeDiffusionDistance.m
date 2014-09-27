function [diffDist d] = computeDiffusionDistance(v1, v2, sigVal, CHOICE)

% Compute Distances between two vectors v1 and v2 of same size and the CHOICE of norm
% may be 'ONE', 'TWO', 'HIST_INT'

if( numel(v1(:)) ~= numel(v2(:)) )
    error('Dimensions of input must match');
end

diffVal = v1(:) - v2(:);

if( strcmp(CHOICE, 'ONE' ) )
    d = norm(diffVal,1);
elseif( strcmp(CHOICE, 'TWO' ) )
    d = norm(diffVal,2);
elseif( strcmp(CHOICE, 'HIST_INT' ) )
    d = sum( min( v1(:), v2(:) ) );
elseif( strcmp(CHOICE, 'MATRIX_EUC' ) )
    aa=sum(v1.*v1,1); bb=sum(v2.*v2,1); ab=v1'*v2; 
    d = sqrt(abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab));
    sigVal = sigVal / 10;
else
    error('Unrecognized Distance');
end
diffDist = exp( -d ./ sigVal );