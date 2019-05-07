function [W,SX,u,v] = utp_ws(p,n)

switch p
  case 3
    [W,SX] = ut3_ws(n); u=[]; v=[];  
  case 5
    [W,SX,u] = ut5_ws(n); v=[];  
  case 7
    [W,SX,u,v] = ut7_ws(n);    
  case 9
    [W,SX,u,v] = ut9_ws(n);
  otherwise
    error('Not implemented')  
end
    
    
    