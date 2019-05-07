function [x,Info] = HMCMCGPFastV2(x,y,Params,Opts)

  % This function RANDOMISES the step size, epsilon, as
  % suggested by Iain, but it appears to make things _much_
  % less efficient. 
  
% function [x,Info] = HMCMCGPFastV2(x,y,Params,Dims,Opts)
%
% carriers out HMCMC Sampling for the FAST version of GPPAD
% returning the long version of the latent variables
% 
%
% INPUTS
% x = latent variables over which to sample [1,Tx]
% y = observations [T,1]
% Params = structure of parameters  
% Opts = a structure of options containing:
%         epsilon = size of the leap-frog step
%         N = number of samples to return
%
% OUTPUTS
% x = samples from the distribution [1,Tx]
% Info = structure of information about the sampling run including
%  acc = number of samples accepted
%  epsilon = final value of the step size parameter
  
T = length(y);
Tx = length(x);

Dims.T = T; Dims.Tx = Tx;

fftCov = GetFFTCovFast(Params.len,Tx);

mean_step_size = Opts.mean_step_size;
Tau = Opts.Tau;
N = Opts.N;
NUpdateEps = Opts.NUpdateEps;

% Get the initial objective and derivatives
[E,g] = GetGPObjFast(x,y,fftCov,Params,Dims);

xSamp = x; acc = 0; accLoc = 0;

for samp=1:N   % loop over samples
    p = randn(Tx,1);   % initialise momentum
    H = p'*p/2 + E;     % evaluate the Hamiltonian
    
    % Draw epsilon from an exponential distribution
    epsilon = - mean_step_size * log(1-rand);
    
    xnew = x; gnew = g;
    
    for tau = 1:Tau     % make tau leap-frog steps
        p = p - epsilon*gnew/2;         % make a half step in p
        xnew = xnew + epsilon*p;        % make a step in x
        
        [Enew,gnew] = GetGPObjFast(xnew,y,fftCov,Params,Dims);

        p = p - epsilon*gnew/2;         % make a half step in p
    end
    
    Hnew = p'*p/2+Enew; % find the new value of H
    dH = Hnew - H;      % decide whether to accept
    
    if dH < 0 
        accept = 1;
    elseif rand< exp(-dH)
        accept = 1;
    else
        accept = 0;
    end
    
    if accept==1
        g = gnew; x = xnew; E = Enew;
        acc = acc+1;
        accLoc = accLoc+1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if mod(samp,NUpdateEps)==0
      if accLoc/NUpdateEps<Opts.MinAccRate 
        mean_step_size = mean_step_size/2;
      elseif accLoc/NUpdateEps>Opts.MaxAccRate 
        mean_step_size = mean_step_size*2;
      else
      end
      accLoc = 0;
    end
    
    str1 = ['Fraction complete: ',num2str(round(samp/N*100)/100)];
    str2 = [' Fraction accepted: ',num2str(round(acc/samp*100)/100)];
    str3 = [' Mean step size: ',num2str(mean_step_size)];
    str4 = [' Epsilon: ',num2str(epsilon)];

    
    if samp<N
      fprintf(['\r',str1,str2,str3,str4,'     ']);
    else
      fprintf(['\r',str1,str2,str3,str4,'      \n\n']);
    end
    
end

Info.acc = acc;
Info.mean_step_size = mean_step_size;
