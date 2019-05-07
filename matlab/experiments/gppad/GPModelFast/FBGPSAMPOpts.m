%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard options for FBGPMAP.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options for Inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Opts.tol = 5; % extra number of length-scales of latents
              % to add on to avoid edge effects

Opts.ITS =   [  [500,1000];...
      ones(9,1)*[0,50]]; % number of iterations of ESSMCMC

Opts.OptsHMCMC.epsilon = 0.0025; % step size in the dynamics
Opts.OptsHMCMC.NUpdateEps = 5; % number of samples to draw before updating epsilon 
Opts.OptsHMCMC.Tau = 200;       % number of leap-frog updates per proposal 
Opts.OptsHMCMC.MinAccRate = 0.75; % Minimum tolerable acceptance rate
Opts.OptsHMCMC.MaxAccRate = 0.95; % Maximum tolerable acceptance rate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options for Learning
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
Opts.NumIts = 20;

