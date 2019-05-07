function PlotSingleChannelHil(y,a,FSamp,Rng)

FontName = 'Times'; FSlg = 12;
time = [1:length(y)]/FSamp;

car = y./a;

scale = max(y)/max(a);

a = a*scale;



aHil = abs(hilbert(y));
cHil = y./aHil;

figure
%subplot(2,1,1)
hold on
plot(time,y,'-k')
plot(time,aHil,'-b','linewidth',2)
plot(time,a,'-r','linewidth',2)
legend('data','hilbert','pad')
set(gca,'xlim',Rng)

%subplot(2,1,2)
%hold on
%plot(time,car,'-b')
%drawnow;