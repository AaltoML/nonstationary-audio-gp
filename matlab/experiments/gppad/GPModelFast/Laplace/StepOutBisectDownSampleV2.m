function [Len,Obj,Info] = StepOutBisectDownSampleV2(Len,y,Params,Opts)
  
  % function [Len,Obj,Info] = StepOutBisectDownSampleV2(Len,y,Params,Opts)
  %
  % Approximate Maximum Likelihood learning of the
  % time-scale in GPPAD using Laplace's approximation. A
  % bunch of (hopefully) intelligent heuristics are used to
  % complete this optimisation that have been tweaked using
  % many different stimuli.
  % 
  % The basic algorithm evaluates the Laplace objective at
  % three length-scales specified in Len. Then, if the
  % optimum lies in this range (i.e. between Len(1) and
  % Len(3)) the method bisects to reduce the range. If it
  % lies outside of this range, then the range is extended
  % by shifting Len on a log-scale. Some noise is added to
  % the shift to stop cycles back and forth betweem two
  % sets of identical length scales. This can happen
  % because the signals are split into chunks for
  % evaluating the objective, and the chunk length is set
  % using Len(1). Therefore, each Len(1) implicitly
  % defines a slightly different objective.
  %  
  % Down-sampling is used to speed up the method. The
  % down-sampling rate is initially set by the smallest
  % value in Len and by the minimum allowed length-scale
  % (Opts.MinLength). The down-sampling rate is reduced as
  % more bisections occur. So, the computation time per
  % iteration steadily rises, and the accuracy increases.
  %
  %  
  % INPUTS
  % Len = 3 length scales hopefull chosen so that the
  %       optimum lies somewhere between them [3,1]
  % y = signal [T,1]
  % Params = Structure of parameters including
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
  %     Overlap = amount to overlap segments (in time-scales)  
  %     MinLength = small length-scale to allow (controls
  %                 downsampling rates)
  %     NumItsBisect = maximum number of bisections to
  %                    allow
  %     LenTol = determine length-scale to within this
  %              percentage error (e.g. LenTol = 0.03
  %              means learn it to within 3%)
  %
  % OUTPUTS
  % Len = learned length-scales centre one is the best
  %       guess at the approximate ML length-scale
  %       [3,1]. Other two who the error in the location
  %       of the ML estimate - not the uncertainty
  %       arising from having finite data.
  % Obj = the objective at the three length-scales in Len [3,1]
  % Info = structure of information about the learning including  
  %        LenHist = the history of learned length scales
  %        ObjHist = the history of associated objectives  

  if length(Len)==1;
    LMin = 0;
    LMax = inf;
    Len = Len*[2/3,1,3/2];
  elseif length(Len)==2;
    LMin = Len(1);
    LMax = Len(2);
    Len = (Len(1)+Len(2))/2*[2/3,1,3/2];
    Len(1) = max([Len(1),LMin]);
    Len(3) = min([Len(3),LMax]);
  else
  end
    
  Len = sort(Len); % in case the length-scales are
                   % incorrectly ordered 
    
  NumBisect = 0;  % zero when criteria is reached
    
  LastMaxLoc = inf; % used to remember the last max location
  
  fprintf('\n\nSTEP OUT, BISECT AND DOWNSAMPLE ALGORITHM')
  
  numSteps = 1; % number of steps

  MinLength = Opts.MinLength; % minimum length-scale

  MaxSep = 1.5; % top and bottom lengths must be no more
                % than this fraction of the middle value
  MaxLength = 100; % maximum length-scale to use,
                   % downsample above this
  
  MaxSteps = 50; % maximum number of steps 
  
  LenHist = []; % history of length scales
  ObjHist = []; % history of objective evaluations
  DSHist = [];  % history of down sampling rates
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  while Opts.NumItsBisect>NumBisect % loop while number
                                     % of bisections
                                     % less than specified 
    numSteps = numSteps+1; % update number of steps  
        
    % Down sample to make things fast (if accuracy tolerance has
    % not been reached)
    if  Opts.LenTol>(1-Len(1)/Len(2)) & Opts.LenTol>(Len(3)/Len(2)-1)
      DS = floor(max([1,Len/MaxLength]));
    else
      DS = floor(max([Len/MinLength,1,Len/MaxLength]));
    end
    
    LenDS = Len/DS;
    yDS = y(1:DS:end);
  
    % Get the objective for the current three points
    [LenDS,Obj,InfoLap] = GetLaplaceObjGPPADMulti(LenDS,yDS, ...
                                                  Params,Opts);    
    
    % Update the histories
    LenHist = [LenHist;Len];
    ObjHist = [ObjHist;Obj];
    DSHist = [DSHist;DS];

    % Display useful information
    fprintf(['\nSOBDS iteration number: ',num2str(numSteps)])   
    fprintf(['\nSOBDS time-scales: ',num2str(Len)])   
    fprintf(['\nSOBDS objectives: ',num2str(Obj)])
    fprintf(['\nSOBDS down-sample rate: ',num2str(DS)])
    str = num2str(max([1-Len(1)/Len(2),Len(3)/Len(2)-1])*100);
    fprintf(['\nSOBDS proximity: ',str,'%%'])
    fprintf(['\nSOBDS bisections: ',num2str(NumBisect),'/'])
    fprintf([num2str(Opts.NumItsBisect),'\n'])

    % Exit if tolerance has been reached
    if Opts.LenTol>(1-Len(1)/Len(2)) & ...
          Opts.LenTol>(Len(3)/Len(2)-1) & InfoLap.MaxLoc ==0
      Info.NumBisect = NumBisect;
      break;
    end

    % Exit if max number of steps has been reached
    if numSteps>MaxSteps
      Info.NumBisect = NumBisect;
      break;
    end

    
    % Figure out whether the optimum is still at the
    % flankers or whether it is in the middle
    if InfoLap.MaxLoc ==0 % optimum in the middle - we bisect

      NumBisect=NumBisect+1; % Increase count of bisections
      %MinLength = sqrt(2)*MinLength; % Down sampling rate reduced
      MinLength = 1.4*MinLength; % Down sampling rate reduced
      
      % Decide whether to bisect to the right or to the
      % left of the middle point
      
      l1 = log(Len(2))-log(Len(1));
      l2 = log(Len(3))-log(Len(2));
      
      if l1>l2
        flag = -1;      
      elseif l1<l2
        flag = 1;
      elseif rand<1/2
        flag = -1;
      else
        flag = 1;
      end  
      
      if flag ==-1 % left bisection
        Len = [Len(1),(Len(1)+Len(2))/2,Len(2)];
      else % right bisection
        Len = [Len(2),(Len(2)+Len(3))/2,Len(3)];
      end
      
    elseif InfoLap.MaxLoc==-1 % optimum at the left flanker
           
      if Len(1)>LMin
        % If we are not up against the bottom of the range
        if LastMaxLoc==0
          % if the last iteration was a bisection then keep the
          % separation between points roughly constant...
          fac = Len(2)/Len(1); % mean factor to step out by
        elseif LastMaxLoc==1
          % if the last iteration was on the right then
          % decrease the separation
          fac = Len(2)/Len(1); % mean factor to step out by
          beta = -1/2;
          fac = min([1+(fac-1)*(1+beta),MaxSep]);
        else
          % if the last iteration was also on the left
          % flanker, then increase the separation 
          fac = Len(2)/Len(1); % mean factor to step out by
          beta = 1/2;
          fac = min([1+(fac-1)*(1+beta),MaxSep]);
        end
        
        Len = [Len(1)/fac,Len(1),Len(2)]; % step out in the new
                                          % direction
                                          % add noise to stop cycles
        Len(1) = Len(1)+(Len(2)-Len(1))/10*(rand-1/2); 
        Len(3) = Len(3)+(Len(3)-Len(2))/10*(rand-1/2);

        Len(1) = max([Len(1),LMin]);
                
      else
        % If we are up against the bottom of the range
        
        Len = [Len(1),(Len(2)+Len(1))/2,Len(2)]
      end
      
    elseif InfoLap.MaxLoc==1 % optimum at the right flanker
      
      if Len(3)<LMax
        % If we're not up against the top of the range
        if LastMaxLoc==0
          % if the last iteration was a bisection then keep the
          % separation between points roughly constant...
          fac = Len(3)/Len(2);
        elseif LastMaxLoc==-1
          % if the last iteration was on the left
          % flanker, then decrease the separation 
          fac = Len(3)/Len(2);
          beta = -1/2;
          fac = min([1+(fac-1)*(1+beta),2]);
        else
          % if the last iteration was also on the right
          % flanker, then increase the separation 
          fac = Len(3)/Len(2);
          beta = 1/2;
          fac = min([1+(fac-1)*(1+beta),2]);
        end
      
        Len = [Len(2),Len(3),Len(3)*fac];
        
        % add noise to stop cycles
        Len(1) = Len(1)+(Len(2)-Len(1))/10*(rand-1/2);
        Len(3) = Len(3)+(Len(3)-Len(2))/10*(rand-1/2);
      
        Len(3) = min([Len(3),LMax]);
                
      else
        % If we are up against the top of the range
        Len = [Len(2),(Len(3)+Len(2))/2,Len(3)];
      end
    else % error if we are here
      disp('Error in StepOutLength')
      return;
    end

    LastMaxLoc = InfoLap.MaxLoc;
   
  end
    
  % return various bits of information
  Info.LenHist = LenHist;
  Info.ObjHist = ObjHist;
  Info.DSHist = DSHist;
  Info.numSteps = numSteps;