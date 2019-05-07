function [x,Info] = HMCMCGPMarg2(x,y,Params,Opts)

% function [x,Info] = HMCMCGPMarg2(x,y,Params,Opts)
%
% carriers out HMCMC for GPPAD on the MARGINAL DATA just
% for the latents and NOT FOR THE PARAMETERS (see
% HMCMCGPMarg2 for this)

T = length(y);

epsilon = Opts.epsilon;
Tau = Opts.Tau;
N = Opts.N;
NUpdateEps = Opts.NUpdateEps;

% Get the initial objective and derivatives
[E,g] = GetGPMargObj(x,y,Params);

xSamp = x; acc = 0; accLoc = 0;

for samp=1:N   % loop over samples
    p = randn(T,1);   % initialise momentum
    H = p'*p/2 + E;     % evaluate the Hamiltonian
    
    xnew = x; gnew = g;
    
    for tau = 1:Tau     % make tau leap-frog steps
        p = p - epsilon*gnew/2;         % make a half step in p
        xnew = xnew + epsilon*p;        % make a step in x
        
        [Enew,gnew] = GetGPMargObj(xnew,y,Params);

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
      if accLoc/NUpdateEps<0.2
        epsilon = epsilon/2;
      elseif accLoc/NUpdateEps>0.8
        epsilon = epsilon*2;
      else
      end
      accLoc = 0;
    end
    
    str1 = ['Fraction complete: ',num2str(round(samp/N*100)/100)];
    str2 = [' Fraction accepted: ',num2str(round(acc/samp*100)/100)];
    str3 = [' Epsilon: ',num2str(epsilon)];
    
    if samp<N
      fprintf(['\r',str1,str2,str3,'     ']);
    else
      fprintf(['\r',str1,str2,str3,'      \n\n']);
    end
    
end

Info.acc = acc;
Info.epsilon = epsilon;
