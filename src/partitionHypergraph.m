function f = partitionHypergraph(A, y)
SOFT = 1; % 1: soft(probabilistic), 2: hard
gamma = 0.1;
if( SOFT )
    H = A;
else
    H = A > 0;
end
diag_w  = sum(A,1);
numNode = size(A,1);
W       = diag(diag_w);
diag_de = sum(H,2);
diag_de(diag_de==0) = eps;
diag_dv = H * diag_w';
diag_dv(diag_dv==0) = eps;

De_n1 = diag(1./diag_de);
Dv_n12 = diag(1./sqrt(diag_dv));
Theta = Dv_n12 * H * W * De_n1 * H' * Dv_n12;

delta_coeff = eye(numNode) - gamma * Theta;

f = delta_coeff \ y;