function DisESSHMCMC(Info,OptsHMCMC,OptsESS,it,ITS);

% function DisInfESSML(Params,Info,OptsHMCMC,OptsESS,OptsML,it,ITS);
%
% Displays information about the function ESS_HMCMC_GPFAST
 
fprintf(['Sampling Progress ',num2str(it),'/',num2str(size(ITS, ...
                                                  1)),'\n']);
