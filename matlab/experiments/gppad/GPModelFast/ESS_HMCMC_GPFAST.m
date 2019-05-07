function [X,Params,Info] = ESS_HMCMC_GPFAST(x,y,Params,Opts)
  
  % [x,Params,Info] = ESS_HMCMC_FAST(x,y,Params,Opts) 
  %
  % Carries out elliptical slice sampling and HMCMC for
  % the FAST GPPAD model
        
    ITS = Opts.ITS;
    NumIts = size(ITS,1);
    
    OptsHMCMC = Opts.OptsHMCMC;
    OptsESS = [];
    
    Info.HMCMCAcc = zeros(NumIts,1);
    Info.ESSRej = zeros(NumIts,1);
    
    % Duration of the sampling proceedures
    Info.timeESS =  zeros(NumIts,1);
    Info.timeHMCMC =  zeros(NumIts,1);
    Info.timeTOTAL =  zeros(NumIts,1);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X = zeros(length(x),NumIts);
    
    for it = 1:NumIts,
      
      % Display information
      DisESSHMCMC(Info,OptsHMCMC,OptsESS,it,ITS);
      
      % HMCMC Latent Update
      tic;
      OptsHMCMC.N = ITS(it,1);     
      
      % For randomising the step size epsilon:
      %[x,InfoHMCMC] = HMCMCGPFastV2(x,y,Params,OptsHMCMC);
      %OptsHMCMC.mean_step_size = InfoHMCMC.mean_step_size;
      
      [x,InfoHMCMC] = HMCMCGPFast(x,y,Params,OptsHMCMC);
      OptsHMCMC.epsilon = InfoHMCMC.epsilon;
      
      Info.HMCMCAcc(it) = InfoHMCMC.acc;
      Info.timeHMCMC(it)=toc;
      
      % ESS Latent Update
      tic;
      OptsESS.N = ITS(it,2);
      [x,InfoESS] = ESSMCMCGPFast(x,y,Params,OptsESS);
      Info.ESSRej(it) = InfoESS.NumRejectsPerSample;
      Info.timeESS(it)=toc;
      % Update samples
      X(:,it) = x;
      
      Info.timeTOTAL(it)=Info.timeESS(it)+Info.timeHMCMC(it);

end
    