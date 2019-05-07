% Plots draws from a student-t AR process with a SqE kernel
% for a variety of sparsenesses and length-scales

clear;

seed = 4; % random seed
addpath ..
LoadLocalPaths;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
T = 1000;

alphaxs = [1/2,1.5,8];
lenxs = [5,10,30];
tau = 7*max(lenxs); vary = 0;

NumAlphas = length(alphaxs);
NumLens= length(lenxs); 

mVarx = 1;
mux = -1; 
varc = 1; 

X = zeros(T,NumAlphas,NumLens);

for aa = 1:NumAlphas
  disp(['Progress ',num2str(aa),'/',num2str(NumAlphas)])
  for ll = 1:NumLens

    lenx = lenxs(ll);
    alphax = alphaxs(aa);
    
    [lamx,varx] = Cov2LamVar(mVarx,lenx,tau);
    lamx = lamx';

    Params = PackParams(varx,lamx,mux,alphax,varc,vary);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Dimensions
    
    tau = length(lamx);
    Dims = PackDims(T,tau);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Draw sample from the forward model
    randn('seed',seed)
    rand('seed',seed)  
    [y,a,X(:,aa,ll)] = FMStudentAR(Dims,Params);

  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hspace = 0.05;
vspace = 0.05;
top = 0.05;
bottom = 0.05;
left = 0.05;
right = 0.05;

width =( 1-left-right-hspace*(NumLens-1))/NumLens;
height =( 1-top-bottom-vspace*(NumAlphas-1))/NumAlphas;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

for aa=1:NumAlphas
  for ll = 1:NumLens
    posCur = [left+(ll-1)*(width+hspace),bottom+(aa-1)*(height+vspace), ...
              width,height];
    axCur = axes('position',posCur);
    hold on
    plot(X(:,aa,ll),'-k')
  
  end
end
