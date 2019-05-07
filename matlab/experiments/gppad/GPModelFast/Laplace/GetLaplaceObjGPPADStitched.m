function [Objs,Lens] = GetLaplaceObjGPPADStitched(LenRng, ...
                                              NumLens,y,Params,Opts)
  
% function [Objs,Lens] =
%     GetLaplaceObjGPPADStitched(LenRng,NumLens,y,Params,Opts)
%
% Evaluates the stitched together objective function
% (Laplace's approximation to the likelihood of the model,
% log p(y|\theta)). The objective has to be 'stitched
% together' because only local regions can be evaluated
% at the same time, due to constraints arising from the
% eigenvalue truncation.
%
% The function can be very slow for long signals, but may be
% useful especially if you believe the signal contains
% modulation on multiple different time-scales. These will
% show up as multiple optima in Objs.
%
% If you just want to learn a time-scale,
% GetLenGradAscent.m might be preferable.
%
% INPUTS
% LenRng = First and last time-scale to be evaluated [1,2] 
% NumLens = number of time-scales to evaluate. 
%           If NumLens = [], then the method figures out
%           the smallest value of NumLens that can cover
%           the range.
% y = Signal [T,1]
% Params = structure of parameters including
%     varc = carrier variance
%     varx = transformed envelope variance
%     mu = transformed envelope mean
%     vary = observation noise (scalar or [T,1])    
% Opts = structure of options including
%     NumIts = number of iterations in the inference step  
%     tol = amount of padding (missing) data (in time-scales)  
%     nev = number of eigenvalues to find (around 100 is a
%           good number)  
%     TruncPoint = Number of standard deviations to retain 
%     MinLength = small length-scale to allow (controls
%                 downsampling rates)
% OUTPUTS
% Objs = evaluation of objective function [NumLens,1]
% Lens = time-scales at which the objective function is
%        evaluated at [NumLens,1]


alpha = 3;

% If the user does not specify a number of length-scales
% then use the minimum number to cover the range.
if isempty(NumLens)
  NumLens = ceil((log(LenRng(2))-log(LenRng(1)))/log(alpha))+1;
  disp(['Number of time-scales to evalulate objective at: ',num2str(NumLens)])
end

logLens = linspace(log(LenRng(1)),log(LenRng(2)),NumLens);
Lens = exp(logLens);

if logLens(2)-logLens(1)>log(alpha)
  disp('Error in GetLaplaceObjGPPADStitched')
  disp('Separation between points is too great:')
  disp('Reduce LenRng or increase NumLens')
  return;
end

NumPerSeg = sum(logLens<logLens(1)+log(alpha));

NumSeg = ceil((NumLens-1)/(NumPerSeg-1));

disp(['Number of time-scales per segment is: ',num2str(NumPerSeg)])

Objs = repmat(NaN,[NumLens,1]);

for s=1:NumSeg
  
  % Information
  str1 = 'Stitched Objective computation progress: ';
  str2 = [num2str(s),'/',num2str(NumSeg)];
  disp([str1,str2])
  
  % Figuring out the range of length-scales to test
  sStart = 1+(s-1)*(NumPerSeg-1);
  sEnd = min([sStart-1+NumPerSeg,NumLens]);
  ind = sStart:sEnd;
  curLens = Lens(ind);

  % Evaluating the objectives for this section
  [Ls,curObj,InfoLap] = ...
      GetLaplaceObjGPPADMulti(curLens,y,Params,Opts);    

  % Stitching the objectives together 
  if s==1 
    % For the first set of objectives to be computed no
    % stitching is required.
    Objs(ind) = curObj;
  else
    % For subsequent objectives, subtracting off an offset is
    % used to stitch the objectives back together.
    delta = curObj(1)-Objs(ind(1));
    Objs(ind(2:end)) = curObj(2:end)-delta;
  end

  % clf;
  % subplot(2,1,1)
  % hold on
  % plot(Lens,Objs)
  
  % subplot(2,1,2)
  % plot(curLens,curObj,'r')
  % keyboard

  
end

