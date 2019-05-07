function Dims = PackDimsGPFast(T,Tx); 

  % function Dims = PackDimsGPFast(T,Tx); 
  %
  % Unpacks the dimensions structure into individual
  % variables
  %  
  % INPUTS  
  % Dims = Structure of dimensions including
  %    T = length of observations
  %    Tx = length of latents
  % Outputs - as above, but as the component variables
    
  Dims.T = T;
  Dims.Tx = Tx;
 
