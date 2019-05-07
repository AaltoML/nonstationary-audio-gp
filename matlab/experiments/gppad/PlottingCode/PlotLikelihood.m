function PlotLikelihood(Lens,Objs)
  
% Finding local optima
K = length(Objs);
  
if Objs(1)>Objs(2)
  LenMax = Lens(1);
  ObjMax = Objs(1);
else
  LenMax = [];
  ObjMax = [];
end

  for k=2:K-1
    if Objs(k)>Objs(k-1)&Objs(k)>Objs(k+1)
      LenMax = [LenMax;Lens(k)];
      ObjMax = [ObjMax;Objs(k)];
    end
  end

  if Objs(K)>Objs(K-1)
        LenMax = [LenMax;Lens(K)];
        ObjMax = [ObjMax;Objs(K)];
  end

figure
hold on
plot(Lens,Objs,'-k','linewidth',2)
temp = plot(LenMax,ObjMax,'o','markerfacecolor',[1,0,0],'markeredgecolor',[1,0,0],'markersize',5);
%set(gca,'xscale','log','xlim',[min(Lens),max(Lens)])
set(gca,'xlim',[min(Lens),max(Lens)])

