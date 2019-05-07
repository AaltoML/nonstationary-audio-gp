% DemoGT created by Nick clark~ nick@ihr.mrc.ac.uk ~to show off some of the
% functionality of the GammaTone Filter Kit.   Run this and make sure the
% files decribed in the README pdf are in the working directory or your
% path. 14/06/07--> date UK style ;) ENJOY!


close all; clear all; clc;

sr = 16e3; % sampling rate
x = zeros(1, sr*25e-3); %create a 25ms input
x([1 100 200 300]) = 1; %make a click train

[forward,feedback,fc,ERB,B] = GammaToneMake(sr,50,200,3000,'Moore'); %make filterbank
y = GammaToneApply(x,forward,feedback); %filter into individual channels

figure; subplot(1,2,1); BMMplot(y, fc, sr, [1 10 20 30 40 50]); %Label a selection of centre frequencies
subplot(1,2,2); BMMplot(y([1 10 20 30 40 50],:), fc, sr, []); % dont plot the centre frequencies on the ordinate and plot a selection of centre frequencies.