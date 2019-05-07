%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard options for learning the time-scales using GPPPAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Loading standard options for learning the time-scales using GPPAD')

Opts.NumIts = 2000; % total number of iterations of
                    % gradient based optimisation

Opts.nev = 100; % Number of eigenvalues to compute

Opts.tol = 5; % extra number of length-scales of latents
              % to add on to avoid edge effects

Opts.MinLength = 8; % Down sample to this length but no
                    % more during inference

Opts.TruncPoint = 8; % Truncation point for eigenvalues/vectors

Opts.NumItsBisect = 10; % Max number of iterations of the
                        % bisection algorithm. 

Opts.LenTol = 0.03; % terminate learning of length-scale
                    % when located within this fraction
                    % (e.g. 0.03 = within 3% of the optimum)  


Opts.Overlap = 5; % Overlap of chunks for the error-bar computation
