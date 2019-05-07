function [Len,Obj,Info] = BisectionLaplace(Len,Obj,y,Params,Opts)
  
  % function [len,Info] = BisectionLaplace(Len,Obj,y,Params,Opts)
  %
  % Finds a best estimate for the time-scales starting with
  % three time-scales and associated objective evaluations
  % (two of which flank the optima, and a third somewhere in
  % between), and then shrinking the flankers down by
  % bisecting either side.

  ParamsCur = Params;

  LenHist = [];
  ObjHist = [];

  fprintf('\n\nBISECTION ALGORITHM')
  
  if Obj(2)<Obj(1) || Obj(2)<Obj(3)
    fprintf('\aERROR - central objective not the largest')
    return
  end
  
  
  for it=1:Opts.NumItsBisect 

    LenHist = [LenHist;Len];
    ObjHist = [ObjHist;Obj];
    
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
    
    % compute the objective at the new bisection point
        
    if flag ==-1 % left bisection
      lenNew = (Len(1)+Len(2))/2;
      ParamsCur.len = lenNew;
      ObjNew = GetLaplaceObjGPPAD(y,ParamsCur,Opts);
      
      % Update the lengths and the objectives
      if ObjNew>Obj(2)
        Len(3) = Len(2);
        Len(2) = lenNew;
        Obj(3) = Obj(2);
        Obj(2) = ObjNew;
      else
        Len(1) = lenNew;
        Obj(1) = ObjNew;
      end
      
    else % right bisection
      lenNew = (Len(2)+Len(3))/2;
      ParamsCur.len = lenNew;
      ObjNew = GetLaplaceObjGPPAD(y,ParamsCur,Opts);
      
      % Update the lengths and the objectives
      if ObjNew>Obj(2)
        Len(1) = Len(2);
        Len(2) = lenNew;
        Obj(1) = Obj(2);
        Obj(2) = ObjNew;
      else
        Len(3) = lenNew;
        Obj(3) = ObjNew;
      end
    
    end
      DispBisectionInfo(Len,Obj,it,Opts);    
  end
    
  Info.LenHist = LenHist;
  Info.ObjHist = ObjHist;

  