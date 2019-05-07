function [Len,Obj,Info] = StepOutLengthFast(Len,y,Params,Opts)
  
  % function [Len,Obj,flag] = StepOutLength(Len,y,Params,Opts)
  %
  % Finds two length-scales either side of the optimum,
  % and evaluates the objective for a third points
  % between these two flankers. The FAST version
  % aggressively downsamples the data based on the
  % shortest length-scale (Len(1)) and the minimum
  % allowed length-scale (Opts.MinLength)

  flag = 1;  % zero when criteria is reached
    
  fprintf('\n\nFAST STEP OUT ALGORITHM')
  
  LenHist = []; % history of length scales
  ObjHist = []; % history of objective evaluations

  numSteps = 0; % number of steps
  
  while flag~=0 % loop while criteria not satisfied
    numSteps = numSteps+1; % update number of steps
    
    % Down sample to make things fast
    DS = max([floor(Len/Opts.MinLength),1]);
    LenDS = Len/DS;
    yDS = y(1:DS:end);
    
    % Get the objective for the current three points
    [LenDS,Obj,InfoLap] = GetLaplaceObjGPPADMulti(LenDS,yDS, ...
                                                  Params,Opts);
    Len = LenDS*DS;
    
    % Update the histories
    LenHist = [LenHist;Len];
    ObjHist = [ObjHist;Obj];
    
    % Display useful information
    fprintf(['\nFast step out number: ',num2str(numSteps)])   
    fprintf(['\nFast step out time-scales: ',num2str(Len)])   
    fprintf(['\nFast step out objectives: ',num2str(Obj)])
    fprintf(['\nFast step out down-sample rate: ',num2str(DS)])

    % Update the fraction by which we alter the
    % length-scales: This stops us moving backwards to
    % exactly where we've just come from. Another choice
    % would have been a random frac e.g. 
    % frac = rand*0.5+1.5;

    frac = 0.5*numSteps^-0.9+1.5;
    
    % Figure out whether the optimum is still at the
    % flankers or whether it is in the middle
    if InfoLap.MaxLoc ==0 % optimum in the middle - we're done
      flag = 0;
    elseif InfoLap.MaxLoc==-1 % optimum at the left flanker
      Len = [Len(1)/frac,Len(1),frac*Len(1)];
    elseif InfoLap.MaxLoc==1 % optimum at the right flanker
      Len = [Len(3)/frac,Len(3),frac*Len(3)];
    else % error if we are here
      disp('Error in StepOutLength')
      return;
    end
  end
    
  % return various bits of information
  Info.TCh = InfoLap.TCh;
  Info.TxCh = InfoLap.TxCh;
  Info.LenHist = LenHist;
  Info.ObjHist = ObjHist;
  Info.DS = DS;
  