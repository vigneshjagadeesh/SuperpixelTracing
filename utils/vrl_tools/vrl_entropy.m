function entropyVal = vrl_entropy(p)

% Compute entropy of a probability distribution 
if( sum(p(:)) ~= 1 )
    p = p(:)./sum(p(:));
end
entropyVal = 0;
p(p<=0) = eps;
if( p < -.05 )
    error('Probability is too negative');
end

for iter = 1:numel(p)
    entropyVal = entropyVal - ( p(iter) * log2(p(iter)) );
end