function [varc,varx,mux,Info] = LearnMarginalParams(y,Opts)     
    
% function [varc,varx,mux,Info] = LearnMarginalParams(y,Opts)
%  
% Learns the parameters from the marginal data using a
% coarse grid search followed by a conjugate gradient
% optimisation routine

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Down sample
NumRetain = Opts.NumRetain; % number points to retain from the total

T = length(y);

DS = floor(T/NumRetain);

sigy = max(sqrt(var(y)),1e-12);

yDS = y(1:DS:NumRetain*DS);

yDS = yDS/sigy;
Opts.varcRng = Opts.varcRng/(sigy^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid Search
%disp('Grid Search for marginal parameters') % WW
[varc,varx,mux,InfoGrid] = GridMargParams(yDS,Opts);

Info.varcGrid = varc*sigy^2;
Info.varxGrid = varx;
Info.muxGrid = mux;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conjugate Gradients
%disp('Conjugate Gradient search for marginal parameters') % WW

% Super weak prior on varc to stop varc running away to
% really large values and mux reducing to compensate
Opts.logvarcMn = -log(sigy^2);
Opts.logvarcVar = 10;

[varc,varx,mux,InfoCG] = ConjGradMargParams(yDS,varc, varx,mux,Opts);

varc = varc*sigy^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
Info.Grid = InfoGrid;
Info.CG = InfoCG;