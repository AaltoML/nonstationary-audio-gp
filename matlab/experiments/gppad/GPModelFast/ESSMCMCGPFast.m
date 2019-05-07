function [x,Info] = ESSMCMCGPFast(x,y,Params,Opts)

% function [x,Info] = ESSMCMCGPFast(x,y,Params,Opts)
%
% carriers out Elliptical Slice Sampling for GPPAD
% returning the long form of the transformed envelopes
% 
%
% INPUTS
% x = variables over which to sample
% variables = additional variables to be passed to the gradient and
% finction evaluation routines [gradE(x,variables) and findE(x,variables) ] 
% epsilon = size of the leap-frog step
% no_samps = number of samples to return
%
% OUTPUTS
% XS = samples from the distribution
% acc = number of samples accepted

Tx = length(x);

% Covert from x to zero mean, unit variance, z   
[varx,len,mux,varc,vary] = UnpackParamsGP(Params); 
T = length(y);
zOld = (x-mux)/sqrt(varx);

% Compute the covariance of the GP
fftCov = GetFFTCovFast(Params.len,Tx);

% Start drawing N samples from the markov chain
N = Opts.N;
NumRejects= 0;

for samp=1:N   % loop over samples

  % Choose the ellipse
  zPert = SampleGP(fftCov);
  logy = log(rand)+GetGPLike(zOld,y,Params);

  % Draw an initial proposal
  theta = rand*2*pi;
  thetaLim = [theta-2*pi,theta];
  
  % slice sample around the ellipse until accepted
  accept = 0;
  
  while accept==0
    
    zNew = zOld*cos(theta)+zPert*sin(theta);
    
    if GetGPLike(zNew,y,Params)>logy
      % accepted
      zOld = zNew;
      accept = 1;
    else
      % rejected: reduce the bracket
      NumRejects = NumRejects+1;
      if theta<0
        thetaLim(1) = theta;
        theta = thetaLim(1)+rand*diff(thetaLim);
      else
        thetaLim(2) = theta;
        theta = thetaLim(1)+rand*diff(thetaLim);
      end
    end
    
  end

  % display information
  str1 = ['Fraction complete: ',num2str(round(samp/N*100)/100)];
  if samp<N
    fprintf(['\r',str1,'     ']);
  else
    fprintf(['\r',str1,'      \n\n']);
  end
  
end

Info.NumRejectsPerSample = NumRejects/N;
x = zOld*sqrt(varx)+mux;
