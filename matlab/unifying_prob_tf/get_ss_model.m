function [F,L,Qc,H,Pinf,dF,dQc,dPinf] = get_ss_model(w,kernel)
  

  % TODO - calculate these
  dF=[];
  dQc=[];
  dPinf=[];
  
  D = length(w)/3;
  omega = w(1:D);
  lengthScale = w(D+1:2*D);
  magnSigma2 = w(2*D+1:3*D);
  

  se_approx_order = 4;

  % Define the hyperparameters
  dt = 1;  % step size is 1 sample, regardless of sample rate
  w = omega;  % om
  
%   if strcmp(kernel,'exp')
%     lengthScale = 1 ./ lamx;
%     dl_dlam = @(lam) - lam^(-2);
%   elseif strcmp(kernel,'matern32')
%     lengthScale = sqrt(3) ./ lamx;
%     dl_dlam = @(lam) -sqrt(3) * lam^(-2);
%   elseif strcmp(kernel,'matern52')
%     lengthScale = sqrt(5) ./ lamx;
%     dl_dlam = @(lam) -sqrt(5) * lam^(-2);
%   end
  
  
  
  Qc1 = zeros(D,1);
  F1=[];L1=[];H1=[];Pinf1=[];dF1_dm=[];dF1_dl=[];dPinf1_dm=[];dPinf1_dl=[];
  dQc1_dm=[];dQc1_dl=[];
  for d=1:D
    % Construct the continuous-time SDE (from Solin+Sarkka AISTATS 2014)
    cf_to_ss = str2func(strcat('cf_',kernel,'_to_ss'));
    [F1d,L1d,Qc1d,H1d,Pinf1d,dF1d,dQc1d,dPinf1d] = cf_to_ss(magnSigma2(d), lengthScale(d), se_approx_order);
    F1 = blkdiag(F1, F1d);
    dF1_dm = blkdiag(dF1_dm, dF1d(:,:,1));
    dF1_dl = blkdiag(dF1_dl, dF1d(:,:,2));
    L1 = vertcat(L1,L1d);
    Qc1(d) = Qc1d;
    dQc1_dm = blkdiag(dQc1_dm, dQc1d(:,:,1));
    dQc1_dl = blkdiag(dQc1_dl, dQc1d(:,:,2));
    H1 = horzcat(H1,H1d);
    Pinf1 = blkdiag(Pinf1,Pinf1d);
    dPinf1_dm = blkdiag(dPinf1_dm,dPinf1d(:,:,1));
    dPinf1_dl = blkdiag(dPinf1_dl,dPinf1d(:,:,2));
  end
  tau1 = length(L1d); % tau = model order (1 for Exponential, 2 for Matern 3/2, etc.)
  dF1(:,:,:) = zeros([size(dF1_dm),2]);
  dF1(:,:,1) = dF1_dm;
  dF1(:,:,2) = dF1_dl;
  dQc1(:,:,:) = zeros([size(dQc1_dm),2]);
  dQc1(:,:,1) = dQc1_dm;
  dQc1(:,:,2) = dQc1_dl;
  dPinf1(:,:,:) = zeros([size(dPinf1_dm),2]);
  dPinf1(:,:,1) = dPinf1_dm;
  dPinf1(:,:,2) = dPinf1_dl;

    
  % Construct the continuous-time SDE: periodic (just a sinusoid) (from Solin+Sarkka AISTATS 2014)
  tau2=2; % real + imaginary
  F2=[];L2=[];Qc2=[];H2=[];dF2_dw=[];
  for d=1:D
    F2 = blkdiag(F2,[0 -w(d); w(d) 0]);
    dF2_dw = blkdiag(dF2_dw,[0 -1; 1 0]);
    L2 = blkdiag(L2,eye(tau2));
    Qc2 = blkdiag(Qc2,zeros(tau2));
    H2 = horzcat(H2,[1 0]);
  end
%   dF2(:,:,:) = zeros([size(dF2_dw),3]);
%   dF2(:,:,3) = dF2_dw;


%%
  % The product of the two (from Solin+Sarkka AISTATS 2014)
%   F = kron(F1,eye(tau2)) + kron(eye(tau1),F2);  % <- first order approach
%   L = kron(L1,eye(tau2));
%   Qc = kron(Qc1,eye(tau2));
  F2_kron=[];L=[];Qc=[];dF2_kron=[];
  for d=1:D  % for higher-order models we must iterate to stack the kronecker products along the diagonal
    idx1 = tau1*(d-1)+1:tau1*d;
    idx2 = tau2*(d-1)+1:tau2*d;
    F2d = F2(idx2,idx2);
    F2d_kron = kron(eye(tau1),F2d);
    F2_kron = blkdiag(F2_kron,F2d_kron);
    dF2d_dw = dF2_dw(idx2,idx2);
    dF2d_kron = kron(eye(tau1),dF2d_dw);
    dF2_kron = blkdiag(dF2_kron,dF2d_kron);
    L = blkdiag(L, kron(L1(idx1),L2(idx2,idx2)));
    Qc = blkdiag(Qc, kron(Qc1(d),L2(idx2,idx2)));
  end
  F = kron(F1,eye(tau2)) + F2_kron;
%   dF2_kron = zeros([size(F),3]);
%   for i=1:2
%       dF1(:,:,i)
%       size(kron(dF1(:,:,i),eye(tau2)))
%       size(dF2_kron(:,:,i))
%     dF(:,:,i) = kron(dF1(:,:,i),eye(tau2));
%   end
  dF(:,:,1) = kron(dF1(:,:,1),eye(tau2));
  dF(:,:,2) = kron(dF1(:,:,2),eye(tau2));
  dF(:,:,3) = dF2_kron;
  dQc(:,:,1) = kron(dQc1(:,:,1),eye(tau2));
  dQc(:,:,2) = kron(dQc1(:,:,2),eye(tau2));
  dQc(:,:,3) = zeros(size(dQc(:,:,1)));
  dPinf(:,:,1) = kron(dPinf1(:,:,1),eye(tau2));
  dPinf(:,:,2) = kron(dPinf1(:,:,2),eye(tau2));
  dPinf(:,:,3) = zeros(size(dPinf(:,:,1)));
  H = kron(H1,[1 0]);
  Pinf = kron(Pinf1,eye(tau2));
  
  
end