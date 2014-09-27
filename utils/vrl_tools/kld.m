function dist = kld(p,q,type)

%%
% Kullback–Leibler divergence (also information divergence, information gain, relative entropy, or KLIC).
%
% example: dist = kld( p , q , 'sym' );
% "p" and "q" must be 1-dimensional vectors and must sum to unity.

%%

if nargin < 3
    error('At least 3 input arguments required. Ex. dist = kld(p,q,''sym'')');
end
if ~(any(size(p,1)==1 || size(p,2)==1 || size(q,1)==1 || size(q,2)==1)) 
    error('First 2 inputs should be 1-dim vector!');
end
if size(p) ~= size (q)
    error('The two vector must have the same size!');
end
if any( abs(sum(p)- 1) > 0.00001 ) || any( abs(sum(q)- 1) > 0.00001 )
    error('KL distance can only be used for probability density vector, i.e. must sum to 1.');
end

disttypes = {'sym','asym'};
disttype = disttypes{ strcmpi( type, disttypes ) };
dist = 0;

switch disttype
case 'sym'
    temp = ( p > 0 & q > 0);
    dist = sum( p(temp) .* log( p(temp) ./ q(temp) ) );
    dist = sum( q(temp) .* log( q(temp) ./ p(temp) ) ) + dist ;
case 'asym'
    temp = ( p > 0 & q > 0);
    dist = sum( p(temp) .* log( p(temp) ./ q(temp) ) );
end

end