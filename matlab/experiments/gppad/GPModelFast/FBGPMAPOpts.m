%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard options for FBGPMAP.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options for Inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Opts.tol = 5; % extra number of length-scales of latents
              % to add on to avoid edge effects
Opts.MinLength = 8; % Down sample to this length but no
                    % more in the cascade

if ~isfield(Opts,'NumIts')
  Opts.NumIts = 2000; % total number of iterations for
                      % MAP inference of the envelopes
end

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

% Number of conjugate gradient steps for marginal
% parameter estimation
Opts.NumItsMM = 20;

