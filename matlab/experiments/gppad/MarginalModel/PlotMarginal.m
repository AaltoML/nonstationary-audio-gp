function PlotMarginal(varx,varc,mux,Opts) 

%  function PlotMarginal(varx,varc,mux,Opts) 
% 
% Plots the log of the marginal probability of the data
  
tol = 2;  
yLim = tol*sqrt(varc)*log(1+exp(mux+tol*sqrt(varc)))*[-1,1];
  
y = linspace(yLim(1),yLim(2),1000)';

logPY = GetLogPY(y,varx,varc,mux,Opts);
  
figure
subplot(2,1,1)
plot(y,logPY);

subplot(2,1,2)
plot(y,exp(logPY)/sum(exp(logPY)));
