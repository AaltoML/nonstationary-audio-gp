%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard Options for Learning Marginal Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Options for both grid and conjugate gradient searches
Opts.NumXs = 250;% Number of xs to integrate over
Opts.zUpper = 5; % Upper limit of x
Opts.zLower = -5; % Lower limit of x
Opts.ChSz = 1000; % Length of chunk
Opts.NumRetain = 500; % Just use 500 samples to learn off

% Options for grid search
Opts.varcRng = [0.1,10];  
Opts.varxRng = [0.1,10];  
Opts.muxRng = [-3,3];  
Opts.NumGridPoints = 5;

% Options for conjugate gradients
Opts.NumItsMM = 20;

