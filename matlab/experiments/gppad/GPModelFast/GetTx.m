function Tx = GetTx(T,tau);

  
  % function Tx = GetTx(T,tau);
  %
  % Gets the length of the latent variables by adding on
  % an extra bit for edge effects (tau) and rounding up
  % to the nearest power of 2
    
  Tx = 2^ceil(log2(T+tau));  
    