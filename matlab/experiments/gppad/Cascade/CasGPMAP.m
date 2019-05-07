function [A,X,c,Params,Info] = CasGPMAP(y,len,Opts)
  
  % function [A,X,c,Params,Info] = CasGPMAP(y,lenx,Opts)
  %
  % Decomposes the signal into a cascade of modulators
  % and a carrier
  %
  % INPUTS
  % y = signal [T,1]
  % len = time-scales to demodulate at measured in
  %       samples [M,1]
    
  % Opts = a structure containing various options required for
  %        inference. See CASGPMAPOpts.m for more details
  %
  % OUTPUTS
  % A = recovered modulators [T,M]
  % X = recovered transformed modulators [T,M]
  % c = carrier [T,1]
  % Params = structure of parameters including
  %     varc = variance of the carriers
  %     varx = variance of the transformed modulator [M,1]
  %     len = time-scales of the modulators (same as the
  %           input len) [M,1]
  %     mux = mean of the transformed modulators [M,1]  
  % Info = structure containing various information about
  %        inference and learning.  
  
  % rescale the signal to limit artifacts  
  sigy = sqrt(var(y));
  y = y/sigy;
  
  % Option for inference sub-routine to return the hidden
  % and observed signal
  Opts.Long = 1;
  
  % Sizes of the model
  M = length(len);
  T = length(y);
    
  % Figure out the missing region to add
  Tx = GetTx(T,max(len)*Opts.tol);
  Dims.T = T;  Dims.Tx = Tx;   Dims.M = M; 

  % Preallocate the modulators
  A = zeros(Tx,M);
  
 
  % Demodulate the signal at the shortest time-scale
  m = 1;
  fprintf(['\n','Level ',num2str(m),' \n'])

  Opts.Tx = Tx;
  [A(:,1),X,c,Params,Info_1] = FBGPMAP(y,len(1),Opts);

  Params.len = len;

  % Do recursive demodulation of the modulator
  for m=2:M
      fprintf(['\n','Level ',num2str(m),' \n'])

      % Estimate the parameters and the modulator
      [A(:,m),x,c,Params_m,Info_m] = FBGPMAP(A(1:Tx-len(m)* ...
                                               Opts.tol, ...
                                               m-1),len(m),Opts);

      % Update the previous modulator by dividing out the
      % new modulator
      A(:,m-1) = A(:,m-1)./A(:,m);
      
      % Update the parameters
      Params.varc = Params.varc*Params_m.varc;
      Params.varx = [Params.varx,Params_m.varx];
      Params.mux = [Params.mux,Params_m.mux];

  end

  % Fine tune by joint optimisation of all levels of the cascade
  fprintf(['\n','Fine tuning entire cascade. \n'])
  
  X = log(exp(A)-1);
  
  [X,Obj] = MAPGPCasFast(X,y,Params,Dims,Opts.NumItsCas);
  
  % Put the output into the correct format
  A = log(1+exp(X(1:T,:)));
  c = y./prod(A,2);
  X = X(1:T,:);
  
  Info.ObjFT = Obj;
  
  % rescale things
  Params.varc = Params.varc*sigy^2;
  