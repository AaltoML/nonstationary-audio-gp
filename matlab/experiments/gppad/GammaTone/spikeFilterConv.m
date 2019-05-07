spike = [zeros(1,3) 0.1*blackman(64)' zeros(1,64)];
%spike = [zeros(1,3) 0.1 zeros(1,37)];
filtSpike1 = filter([0.15],[1 -0.95],spike);
filtSpike2 = filter([0.19],[1 -0.5 -0.45],spike);
filtSpike3 = filter([0.5 0.4],[1 -0.3 -0.2],spike);
filtSpike4 = filter([0.1 0.1 0.1 0.1 0.3],[1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1],spike);
filtSpike5 = filter([0.01 0 0 0 0 0 0 0 0 0 0 0.15],[1 -0.93],spike);
filtSpike6 = filter([0.07 0 0 0 0 0 0 0 0 0 0 0.09],[1 -0.45 -0.5],spike);
figure(99);clf
plot(spike,'k--')
hold on
plot(filtSpike1,'LineWidth',2)
plot(filtSpike2,'LineWidth',2)
plot(filtSpike3,'LineWidth',2)
plot(filtSpike4,'LineWidth',2)
plot(filtSpike5,'LineWidth',2)
plot(filtSpike6,'LineWidth',2)
leg=legend('input','output 1 - 1st order','output 2 - 2nd order','output 3 - 2nd order','output 4 - 5th order','output 5 - distributed lag','output 6 - distributed lag');
set(leg,'FontSize',14,'LineWidth',2)