function [T,Tx] = UnpackDimsGPFast(Dims); 

  % function [T,Tx] = UnpackDimsGPFast(Dims); 
  %
  % Unpacks the dimensions structure into individual
  % variables
  %  
  % INPUTS  
  % Dims = Structure of dimensions including
  %    T = length of observations
  %    Tx = length of latents
  % Outputs - as above, but as the component variables
    
  T = Dims.T;
  Tx = Dims.Tx;
 
