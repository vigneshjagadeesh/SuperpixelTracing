function queryVecs = wVC(AdjMat)
% minimum weight vertex cover for active learning

if ( nargin < 1 )
    % % % %     clear; close all; clc;
    % % % %
    % % % % f = [9 1 2 8 4]';
    % % % % fE = [1 1 2 3 4]';
    % % % % tE = [3 5 4 4 5]';
    % % % %
    % % % % numNodes = 5;
    % % % % numConstraints = numel(fE);
    % % % %
    % % % % A = zeros( numConstraints, numNodes );
    % % % %
    % % % % for conIter = 1:numConstraints
    % % % %     A( conIter, fE(conIter) ) = 1;
    % % % %     A( conIter, tE(conIter) ) = 1;
    % % % % end
    % % % %
    % % % %
    % % % % b = ones( numConstraints, 1 );
    % % % % lb = zeros( numNodes, 1 );
    % % % % ub = ones( numNodes, 1 );
    % % % %
    % % % % X = linprog( f, -A, -b, [], [], lb, ub);
    AdjMat = [0    0.7   0    0.3   0   0;
        0.7  0     0.2  0     0.5 0;
        0    0.2   0    0     0   .9;
        0.3  0     0    0     0   .1;
        0    0.5   0    0     0   0;
        0    0     0.9  0.1   0   0
        ];
end

%% Generate some basic statistics
upA = triu(AdjMat);
f = 1./sum(AdjMat,2);
numNodes = size(AdjMat,2);

%% Find all the edges
[fE tE] = find(upA);

numConstraints = numel(fE);

Ae = zeros( numConstraints, numNodes );

for conIter = 1:numConstraints
    Ae( conIter, fE(conIter) ) = 1;
    Ae( conIter, tE(conIter) ) = 1;
end

be = ones( numConstraints, 1 );
lb = zeros( numNodes, 1 );
ub = ones( numNodes, 1 );

X = linprog( f, -Ae, -be, [], [], lb, ub);
sol = round(X);

queryVecs = find( sol );

[val idx] = sort( f (queryVecs) );

queryVecs = queryVecs( idx );