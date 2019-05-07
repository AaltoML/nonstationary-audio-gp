function PlotSingleChannel(y,a,FSamp)

FontName = 'Times'; FSlg = 12;
time = [1:length(y)]/FSamp;

car = y./a;

scale = sqrt(var(car));
car = car/scale;
a = a*scale;

figure
subplot(2,1,1)
hold on
plot(time,y,'-k')
plot(time,a,'-r','linewidth',2)
legend('data','amplitude')

subplot(2,1,2)
hold on
plot(time,car,'-b')
drawnow;