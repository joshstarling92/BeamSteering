close all
clear all
clc
tic

% S.h1 = figure;
% S.h2 = figure;
% S.h3 = figure;
% S.h4 = figure;
% S.h5 = figure;
% S.h6 = figure;
% S.h7 = figure;
% S.h8 = figure;
% S.h9 = figure;
% S.h10 = figure;


for i = 1:10
    B = randn(100*i^1.3);
    C = B.^2;
    %name = S.(sprintf('h%d',i));
    figure
    plot(C)
end
time = toc