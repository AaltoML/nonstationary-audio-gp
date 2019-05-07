function PlotSingleChannelErrorBarsCompare(y1,y2,x1,xVar1,x2,xVar2,Params1,Params2,FSamp)

tol = 3;

FontName = 'Times'; FSlg = 12;
time = [1:length(y1)]/FSamp;

a1 = sqrt(Params1.varc)*log(1+exp(x1));
aUpper1 = sqrt(Params1.varc)*log(1+exp(x1+tol*sqrt(xVar1)));
aLower1 = sqrt(Params1.varc)*log(1+exp(x1-tol*sqrt(xVar1)));

a2 = sqrt(Params2.varc)*log(1+exp(x2));
aUpper2 = sqrt(Params2.varc)*log(1+exp(x2+tol*sqrt(xVar2)));
aLower2 = sqrt(Params2.varc)*log(1+exp(x2-tol*sqrt(xVar2)));

car1 = y1./a1;
car2 = y2./a2;

scale = sqrt(var(car1));
car1 = car1/scale;
car2 = car2/scale;

a1 = a1*scale;
a2 = a2*scale;

aLower1 = aLower1*scale;
aUpper1 = aUpper1*scale;

aLower2 = aLower2*scale;
aUpper2 = aUpper2*scale;

figure
hold on

temp1 = [time,time(end:-1:1)]';
temp2 = [aUpper1;aLower1(end:-1:1)];
temp3 = [aUpper2;aLower2(end:-1:1)];


plot(time,y2,'-k')
plot(time,y1,'-','color',[1,1,1]*0.5)

patch(temp1,temp2,[1,1/2,1/2],'facealpha',0.9,'linestyle','none')
patch(temp1,temp3,[1/2,1/2,1],'facealpha',0.9,'linestyle','none')

plot(time,a1,'-r','linewidth',2)
plot(time,a2,'-b','linewidth',2)
%legend('data','amplitude')
set(gca,'ylim',max(abs([aUpper1;y1;aUpper2;y2]))*[-1.025,1.025])
set(gca,'xlim',[min(time),max(time)])

drawnow;