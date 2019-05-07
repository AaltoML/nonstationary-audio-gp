%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONS FOR GRADIENT ASCENT

Opts.NumIts = 2000; % total number of iterations of
                          % gradient based optimisation

Opts.nev = 100; % Number of eigenvalues to compute

Opts.tol = 8; % extra number of length-scales of latents
              % to add on to avoid edge effects

Opts.MinLength = 8; % Down sample to this length but no
                    % more during inference

Opts.TruncPoint = 8; % Truncation point for
                     % eigenvalues/vectors

Opts.NumItsLL = 30;  % Number of iterations of online
                     % gradient ascent

Opts.loglenLR0 = log(1.2); % initial step size in the
                           % online gradient ascent.

Opts.m = 0.5; % momentum variable in the gradient ascent
