function [Len,Obj,Info] = StepOutLength(Len,y,Params,Opts)
  
  % function [Len,Obj,flag] = StepOutLength(Len,y,Params,Opts)
  %
  % Finds two length-scales either side of the optimum,
  % and evaluates the objective for a third points
  % between these two flankers.

  flag = 1;  % zero when criteria is reached
    
  fprintf('\n\nSTEP OUT ALGORITHM')
  
  LenHist = []; % history of length scales
  ObjHist = []; % history of objective evaluations

  numSteps = 0; % number of steps
  
  while flag~=0 % loop while criteria not satisfied
    numSteps = numSteps+1; % update number of steps
    
    % Get the objective for the current three points
    [Len,Obj,InfoLap] = GetLaplaceObjGPPADMulti(Len,y, ...
                                             Params,Opts);
    % Update the histories
    LenHist = [LenHist;Len];
    ObjHist = [ObjHist;Obj];
    
    % Display useful information
    fprintf(['\nStep out number: ',num2str(numSteps)])   
    fprintf(['\nStep out time-scales: ',num2str(Len)])   
    fprintf(['\nStep out objectives: ',num2str(Obj)])

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

  