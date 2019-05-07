function [varc,varx,mux,Info] = GridMargParams(y,Opts);

% function [varc,varx,mux,Info] = GridMargParams(y,Opts);
%
% Performs a simple grid search to estimate parameters
% from the marginal data
  
varcRng = Opts.varcRng;  
varxRng = Opts.varxRng;  
muxRng = Opts.muxRng;  
N = Opts.NumGridPoints;  

varcs = logspace(log10(varcRng(1)),log10(varcRng(2)),N);
varxs = logspace(log10(varxRng(1)),log10(varxRng(2)),N);
muxs = linspace(muxRng(1),muxRng(2),N);

Obj = -inf;

for l=1:N
 %fprintf(['Grid Search Progress ',num2str(l),'/',num2str(N),'\r']) % WW suppressed
  for m=1:N
    for n=1:N
      
      ObjCur = GetObjNumericalInt(y,varxs(l), varcs(m),muxs(n),Opts);
    
      if ObjCur>Obj
        varx = varxs(l);
        varc = varcs(m);
        mux = muxs(n);
        Obj = ObjCur;
        Info.Obj = ObjCur;
        
      end
      
    end
  end
end

