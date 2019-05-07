function [T,tau] = UnpackDims(Dims); 

  % function [T,tau] = UnpackDims(Dims); 
  %
  % Unpacks the dimensions structure into individual
  % variables
  %  
  % INPUTS  
  % Dims = Structure of dimensions including
  %    T = length of time series
  %    tau = order of autoregressive process
  % Outputs - as above, but as the component variables
    
  T = Dims.T;
  tau = Dims.tau;
