function [F,L,Qc,H,Pinf,dF,dQc,dPinf] = ss_modulators_nmf(w_subband,w_modulator,kernel1,kernel2)
  
  se_approx_order = 6;
  
  D = length(w_subband)/3;
  N = length(w_modulator)/2;
  sig1 = w_subband(1:D);
  len1 = w_subband(D+1:2*D);
  omega = w_subband(2*D+1:3*D);
  sig2 = w_modulator(1:N);
  len2 = w_modulator(N+1:2*N);
  
  cf_to_ss1 = str2func(strcat('cf_',kernel1,'_to_ss'));
  F1_temp = cf_to_ss1(1, 1, se_approx_order);
  tau1 = size(F1_temp,1); % tau = model order (1 for Exponential, 2 for Matern 3/2, etc.)
  
  tau2=2; % cosine kernel: real + imaginary
  
  cf_to_ss2 = str2func(strcat('cf_',kernel2,'_to_ss'));
  F2_temp = cf_to_ss2(1, 1, se_approx_order);
  tau3 = size(F2_temp,1);

%% periodic subband
  
  F1=[];L1=[];Qc1=[];H1=[];Pinf1=[];
  dF1_1 = zeros(tau1*D,tau1*D,2*D);
  dF1_2 = zeros(tau1*D,tau1*D,2*D);
  dQc1_1 = zeros(D,D,2*D);
  dQc1_2 = zeros(D,D,2*D);
  dPinf1_1 = zeros(tau1*D,tau1*D,2*D);
  dPinf1_2 = zeros(tau1*D,tau1*D,2*D);
  for d=1:D
    [F1d,L1d,Qc1d,H1d,Pinf1d,dF1d,dQc1d,dPinf1d] = cf_to_ss1(sig1(d), len1(d), se_approx_order);
    F1 = blkdiag(F1,F1d);
    L1 = vertcat(L1,L1d);
    Qc1 = blkdiag(Qc1,Qc1d);
    H1 = blkdiag(H1,H1d);
    Pinf1 = blkdiag(Pinf1,Pinf1d);
    dF1_1(1+(d-1)*tau1:d*tau1,1+(d-1)*tau1:d*tau1,d) = dF1d(:,:,1); % d_sig
    dF1_2(1+(d-1)*tau1:d*tau1,1+(d-1)*tau1:d*tau1,D+d) = dF1d(:,:,2); % d_len
    dQc1_1(d,d,d) = dQc1d(:,:,1); % d_sig
    dQc1_2(d,d,D+d) = dQc1d(:,:,2); % d_len
    dPinf1_1(1+(d-1)*tau1:d*tau1,1+(d-1)*tau1:d*tau1,d) = dPinf1d(:,:,1); % d_sig
    dPinf1_2(1+(d-1)*tau1:d*tau1,1+(d-1)*tau1:d*tau1,D+d) = dPinf1d(:,:,2); % d_len
  end
  dF1 = dF1_1 + dF1_2;
  dQc1 = dQc1_1 + dQc1_2;
  dPinf1 = dPinf1_1 + dPinf1_2;
  
  % Construct the continuous-time SDE: periodic (just a sinusoid) (from Solin+Sarkka AISTATS 2014)
  F_cos=[];L_cos=[];Qc_cos=[];H_cos=[];dF_cos=[];
  dF_cos_kron = zeros(tau1*tau2*D+N*tau3,tau1*tau2*D+N*tau3,D); % in NMF model change D*tau3 to N*tau3
  for d=1:D
    F_cos = blkdiag(F_cos,[0 -omega(d); omega(d) 0]);
    L_cos = blkdiag(L_cos,eye(tau2));
    Qc_cos = blkdiag(Qc_cos,zeros(tau2));
    H_cos = horzcat(H_cos,[1 0]);
    dF_cos = blkdiag(dF_cos,[0 -1; 1 0]);
  end
  
  % product of the two
  F_cos_kron=[];L_sm=[];Qc_sm=[];
  dF_cos_kron = zeros(tau1*tau2*D+N*tau3,tau1*tau2*D+N*tau3,D);
  for d=1:D  % for higher-order models we must iterate to stack the kronecker products along the diagonal
    idx1 = tau1*(d-1)+1:tau1*d;
    idx2 = tau2*(d-1)+1:tau2*d;
    F_cos_d = F_cos(idx2,idx2);
    F_cos_d_kron = kron(eye(tau1),F_cos_d);
    F_cos_kron = blkdiag(F_cos_kron,F_cos_d_kron);
    L_sm = blkdiag(L_sm, kron(L1(idx1),L_cos(idx2,idx2)));
    Qc_sm = blkdiag(Qc_sm, kron(Qc1(d,d),L_cos(idx2,idx2)));
    dF_cos_d = dF_cos(idx2,idx2);
    dF_cos_d_kron = kron(eye(tau1),dF_cos_d);
    dF_cos_kron(1+(d-1)*tau1*tau2:d*tau1*tau2,1+(d-1)*tau1*tau2:d*tau1*tau2,d) = dF_cos_d_kron;
  end
  F_sm = kron(F1,eye(tau2)) + F_cos_kron;
  H_sm = kron(H1,[1 0]);
  Pinf_sm = kron(Pinf1,eye(tau2));
  
  for d=1:D
    dF_sm_1(:,:,d) = blkdiag(kron(dF1(:,:,d),eye(tau2)),zeros(N*tau3));
    dF_sm_2(:,:,d) = blkdiag(kron(dF1(:,:,D+d),eye(tau2)),zeros(N*tau3));
    dQc_sm_1(:,:,d) = blkdiag(kron(dQc1(:,:,d),eye(tau2)),zeros(N));
    dQc_sm_2(:,:,d) = blkdiag(kron(dQc1(:,:,D+d),eye(tau2)),zeros(N));
    dPinf_sm_1(:,:,d) = blkdiag(kron(dPinf1(:,:,d),eye(tau2)),zeros(N*tau3));
    dPinf_sm_2(:,:,d) = blkdiag(kron(dPinf1(:,:,D+d),eye(tau2)),zeros(N*tau3));
  end
  dF_sm = cat(3,dF_sm_1,dF_sm_2,dF_cos_kron);
  dQc_sm = cat(3,dQc_sm_1,dQc_sm_2,zeros(size(dQc_sm_1)));
  dPinf_sm = cat(3,dPinf_sm_1,dPinf_sm_2,zeros(size(dPinf_sm_1)));
  
  
%% slow varying modulator 
  
  F2=[];L2=[];Qc2=[];H2=[];Pinf2=[];
  dF2_1 = zeros(tau3*N,tau3*N,2*N);
  dF2_2 = zeros(tau3*N,tau3*N,2*N);
  dQc2_1 = zeros(N,N,2*N);
  dQc2_2 = zeros(N,N,2*N);
  dPinf2_1 = zeros(tau3*N,tau3*N,2*N);
  dPinf2_2 = zeros(tau3*N,tau3*N,2*N);
  for d=1:N
    [F2d,L2d,Qc2d,H2d,Pinf2d,dF2d,dQc2d,dPinf2d] = cf_to_ss2(sig2(d), len2(d), se_approx_order);
    F2 = blkdiag(F2,F2d);
    L2 = blkdiag(L2,L2d);
    Qc2 = blkdiag(Qc2,Qc2d);
    H2 = blkdiag(H2,H2d);
    Pinf2 = blkdiag(Pinf2,Pinf2d);
    dF2_1(1+(d-1)*tau3:d*tau3,1+(d-1)*tau3:d*tau3,d) = dF2d(:,:,1); % d_sig
    dF2_2(1+(d-1)*tau3:d*tau3,1+(d-1)*tau3:d*tau3,N+d) = dF2d(:,:,2); % d_len
    dQc2_1(d,d,d) = dQc2d(:,:,1); % d_sig
    dQc2_2(d,d,N+d) = dQc2d(:,:,2); % d_len
    dPinf2_1(1+(d-1)*tau3:d*tau3,1+(d-1)*tau3:d*tau3,d) = dPinf2d(:,:,1); % d_sig
    dPinf2_2(1+(d-1)*tau3:d*tau3,1+(d-1)*tau3:d*tau3,N+d) = dPinf2d(:,:,2); % d_len
  end
  dF2 = dF2_1 + dF2_2;
  dQc2 = dQc2_1 + dQc2_2;
  dPinf2 = dPinf2_1 + dPinf2_2;
  
  for d=1:size(dF2,3)
    dF2_blk(:,:,d) = blkdiag(zeros(tau1*tau2*D),dF2(:,:,d));
    dQc2_blk(:,:,d) = blkdiag(zeros(tau2*D),dQc2(:,:,d));
    dPinf2_blk(:,:,d) = blkdiag(zeros(tau1*tau2*D),dPinf2(:,:,d));
  end
  
%% combine

  F = blkdiag(F_sm,F2);
  L = blkdiag(L_sm,L2);
  Qc = blkdiag(Qc_sm,Qc2);
  H = blkdiag(H_sm,H2);
  Pinf = blkdiag(Pinf_sm,Pinf2);
  
  dF = cat(3,dF_sm,dF2_blk);
  dQc = cat(3,dQc_sm,dQc2_blk);
  dPinf = cat(3,dPinf_sm,dPinf2_blk);
  
end