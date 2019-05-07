function [y,a,x] = FMStudentAR(Dims,Params)
  
  % Draws from PAD with an AR student-t prior over
  % envelopes
  %
  % INPUTS  
  % Params = Structure of parameters including
  %    vary = observation noise (optional)
  %    varx = scale parameter of stundent-t
  %    lamx = AR dynamical parameter of stundent-t [tau,1]
  %    mux = mean parameter of stundent-t
  %    alphax = shape parameter of stundent-t
  %    varc = variance of the carriers  
  % Dims = structure of dimensions
  %    tau = AR process order
  %    T = length of time series to be generated
   
  [T,tau] = UnpackDims(Dims);  
  [varx,lamx,mux,alphax,varc,vary] = UnpackParams(Params); 
  
  Offset = 50*tau; % extra samples to make sure we sample
                   % from the stationary distribution of
                   % the process
  
  x = zeros(T+Offset,1);
    
  for t=tau+1:T+Offset
    x(t) = lamx'*x([t-1:-1:t-tau],1)+sqrt(varx)*trnd(alphax);  
  end  

  x = x(Offset+1:T+Offset,1)+mux;
  
  a = GetAmp(x);
  
  y = sqrt(varc)*randn(T,1).*a;