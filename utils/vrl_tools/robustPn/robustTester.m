%% TestHoc
clear;
close all;
clc;

noRows = 4;
noCols = 4;
unary = [.1 .2 .5 .1 .5 .1 .5 .3 .5 .5 .7 .1 .6 .6 .8 .1;
        .2 .1 .1 .3 .1 .2 .6 .2 .1 .6 .7 .3 .3 .5 .6 .1]';
    
edges = [1 2; 2 3; 3 4; 5 6; 6 7; 7 8;9 10; 10 11; 11 12; 13 14; 14 15; 15 16; 1 5; 2 6; 3 7; 4 8; 5 9; 6 10; 7 11; 8 12; 9 13; 10 14; 11 15; 12 16]';
weights = [1 1 1 1 1 1 1 1 1 1 1 1; 2 2 2 2 2 2 2 2 2 2 2 11];

supStruct(1).num = 8;
supStruct(1).ID = [1 2 3 4 5 6 7 8]-1;
supStruct(2).num = 8;
supStruct(2).ID = [9 10 11 12 13 14 15 16]-1;

segRes = PnWrapper( int32([2, noRows*noCols, numel(weights), 4]), ...
                                            unary(:),  (edges(:)-1)', (100 * weights)', ...
                                            supStruct);