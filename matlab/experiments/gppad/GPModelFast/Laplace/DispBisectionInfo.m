function  DispBisectionInfo(Len,Obj,it,Opts)

  % function DispBisectionInfo(Len,Obj)
  %
  % Displays information from the bisection function
    
  LenStr = '';
  ObjStr = '';
  
  for l=1:3
    LenStr = [LenStr,num2str(Len(l)),' '];
    ObjStr = [LenStr,num2str(round(Obj(l))),' '];
  end
  
  itstr = num2str(it);
  ITstr = num2str(Opts.NumItsBisect);
  
  fprintf(['\nBisection Progress: ',itstr,'/',ITstr,'\n'])  
  fprintf(['Time scales: ',LenStr])  