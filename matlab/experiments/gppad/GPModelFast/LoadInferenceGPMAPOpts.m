%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options for Inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(Opts,'tol')
  Opts.tol = 5; % extra number of length-scales of latents
                % to add on to avoid edge effects
end

if ~isfield(Opts,'MinLength')
  Opts.MinLength = 8; % Down sample to this length but no
                      % more in the cascade
end

if ~isfield(Opts,'NumIts')
  Opts.NumIts = 2000; % total number of iterations for
                      % MAP inference of the envelopes
end
