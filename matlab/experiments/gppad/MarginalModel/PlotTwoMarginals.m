function PlotTwoMarginals(varx1,varc1,mux1,varx2,varc2,mux2,Opts) 

% function PlotTwoMarginals(varx1,varc1,mux1,varx2,varc2,mux2,Opts) 
% 
% Plots the log of the marginal probability of the data
% for two different distributions
  
tol = 1.5;  

yLim1 = tol*sqrt(varc1)*log(1+exp(mux1+tol*sqrt(varc1)));
yLim2 = tol*sqrt(varc2)*log(1+exp(mux2+tol*sqrt(varc2)));

yLim = yLim1*[-1,1];
%yLim = max([yLim1,yLim2])*[-1,1];

y = linspace(yLim(1),yLim(2),1000)';

logPY1 = GetLogPY(y,varx1,varc1,mux1,Opts);
logPY2 = GetLogPY(y,varx2,varc2,mux2,Opts);
  
figure
subplot(2,1,1)
hold on
plot(y,logPY1);
plot(y,logPY2,'-r');
legend('1','2')

subplot(2,1,2)
hold on
plot(y,exp(logPY1)/sum(exp(logPY1)));
plot(y,exp(logPY2)/sum(exp(logPY2)),'-r');
drawnow
