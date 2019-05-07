function PlotSingleChannelErrorBars(y,x,xVar,Params,FSamp)

tol = 3;

FontName = 'Times'; FSlg = 12;
time = [1:length(y)]/FSamp;

a = sqrt(Params.varc)*log(1+exp(x));
aUpper = sqrt(Params.varc)*log(1+exp(x+tol*sqrt(xVar)));
aLower = sqrt(Params.varc)*log(1+exp(x-tol*sqrt(xVar)));

car = y./a;

scale = sqrt(var(car));
car = car/scale;
a = a*scale;
aLower = aLower*scale;
aUpper = aUpper*scale;

figure
subplot(2,1,1)
hold on

temp1 = [time,time(end:-1:1)]';
temp2 = [aUpper;aLower(end:-1:1)];

patch(temp1,temp2,[1,1/2,1/2],'facealpha',1/2,'linestyle','none')

plot(time,y,'-k')
plot(time,a,'-r','linewidth',2)
legend('data','amplitude')
set(gca,'ylim',max(abs([aUpper;y]))*[-1.025,1.025])


subplot(2,1,2)
hold on
plot(time,car,'-b')
set(gca,'ylim',max(abs([car]))*[-1.025,1.025])

drawnow;