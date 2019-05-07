function [varargout] = gf_ep_mods_nmf_mixture(w,x,y,ss,mom,xt,kernel1,kernel2,J,ep_fraction,ep_damping,ep_itts)
% GF_EP_MODULATOR_NMF_CONSTRAINTS - Solve TF-NMF GP model by Power EP with
% constraints on optimisation
%
% Syntax:
%   [...] = gf_ep_modulator_nmf(w,x,y,k,xt)
%
% In:
%   w     - Log-parameters (nu, sigma2, theta)
%   x     - Training inputs
%   y     - Training outputs
%   ss    - State space model function handle, [F,L,Qc,...] = @(x,theta) 
%   mom   - Moment calculations for ADF
%   xt    - Test inputs (default: empty)
%
% Out (if xt is empty or not supplied):
%
%   e     - Negative log marginal likelihood
%   eg    - ... and its gradient
%
% Out (if xt is not empty):
%
%   Eft   - Predicted mean
%   Varft - Predicted marginal variance
%   Covft - Predicted joint covariance matrix (not used here)
%   lb    - 95% confidence lower bound
%   ub    - 95% confidence upper bound
%
% Description:
%   Consider the following GP model:
%
%       f ~ GP(0,k(x,x')),
%     y_i ~ p(y_i | f(x_i)),  i=1,2,...,n,
%
%   where k(x,x') is the prior covariance function. The state space model 
%   giving you linear time complexity in handling the latent function
%   is specified by the function handle 'ss' such that it returns the 
%   state space model matrices
%
%     [F,L,Qc,H,Pinf,dF,dQc,dPinf] = ss(x,theta),
%
%   where theta holds the hyperparameters. See the paper [1] for details.
%     The non-Gaussian likelihood is dealt with using single-sweep EP, also
%   known as assumed density filtering (ADF). See paper [2] for details.
%
%   NOTE: This code is proof-of-concept, not optimized for speed.
%
% References:
%
%   [1] Simo Sarkka, Arno Solin, Jouni Hartikainen (2013).
%       Spatiotemporal learning via infinite-dimensional Bayesian
%       filtering and smoothing. IEEE Signal Processing Magazine,
%       30(4):51-61.
%   [2] Hannes Nickisch, Arno Solin, and Alexander Grigorievskiy (2018). 
%       State space Gaussian processes with non-Gaussian likelihood. 
%       International Conference on Machine Learning (ICML). 
%
% Copyright:
%   2014-2018   Arno Solin
%
%  This software is distributed under the GNU General Public
%  License (version 3 or later); please refer to the file
%  License.txt, included with the software, for details.

%% Check defaults

  % Is there test data
  if nargin < 10, ep_fraction = 0.5; end
  if nargin < 11, ep_damping = 0.1; end
  if nargin < 12, ep_itts = 30; end
   
  
%% Figure out the correct way of dealing with the data

  % Combine observations and test points
  xall = [x(:); xt(:)];
  yall = [y(:); nan(numel(xt),1)];
    
  % Make sure the points are unique and in ascending order
  [xall,sort_ind,return_ind] = unique(xall,'first');
  yall = yall(sort_ind);
  
  % Only return test indices
  return_ind = return_ind(end-numel(xt)+1:end);
  
  
%% Set up model
  
  lik_param = w{1};
  F_z=[];L_z=[];Qc_z=[];H_z=[];Pinf_z=[];
  F_g=[];L_g=[];Qc_g=[];H_g=[];Pinf_g=[];
  Wnmf = [];
  D = 0; N = 0;
  for j=1:J
    param1 = w{2}{j};
    D_ = length(param1)/3;
    D = D + D_;
    param2 = w{3}{j};
    N_ = length(param2) / 2;
    N = N + N_;
    Wnmf = blkdiag(Wnmf,w{4}{j});
    
    cf_to_ss1 = str2func(strcat('cf_',kernel1{j},'_to_ss'));
    F1_temp = cf_to_ss1(1, 1, 6);
    tau1 = size(F1_temp,1); % tau = model order (1 for Exponential, 2 for Matern 3/2, etc.)
    tau2=2; % cosine kernel: real + imaginary
%     cf_to_ss2 = str2func(strcat('cf_',kernel2{j},'_to_ss'));
%     F2_temp = cf_to_ss2(1, 1, 6);
%     g_tau = size(F2_temp,1);
    z_tau = tau1*tau2;
    
    [F_j,L_j,Qc_j,H_j,Pinf_j,~,~,~] = ss(x,param1,param2,kernel1{j},kernel2{j});
    F_z = blkdiag(F_z,F_j(1:D_*z_tau,1:D_*z_tau));
    L_z = blkdiag(L_z,L_j(1:D_*z_tau,1:D_*tau2));
    Qc_z = blkdiag(Qc_z,Qc_j(1:D_*tau2,1:D_*tau2));
    H_z = blkdiag(H_z,H_j(1:D_,1:D_*z_tau));
    Pinf_z = blkdiag(Pinf_z,Pinf_j(1:D_*z_tau,1:D_*z_tau));
    F_g = blkdiag(F_g,F_j(D_*z_tau+1:end,D_*z_tau+1:end));
    L_g = blkdiag(L_g,L_j(D_*z_tau+1:end,D_*tau2+1:end));
    Qc_g = blkdiag(Qc_g,Qc_j(D_*tau2+1:end,D_*tau2+1:end));
    H_g = blkdiag(H_g,H_j(D_+1:end,D_*z_tau+1:end));
    Pinf_g = blkdiag(Pinf_g,Pinf_j(D_*z_tau+1:end,D_*z_tau+1:end));
  end
  F = blkdiag(F_z,F_g);
  L = blkdiag(L_z,L_g);
  Qc = blkdiag(Qc_z,Qc_g);
  H = blkdiag(H_z,H_g);
  Pinf = blkdiag(Pinf_z,Pinf_g);
%   keyboard
  
%% Prediction of test inputs (filtering and smoothing)

  % Check that we are predicting
  if ~isempty(xt)
    
    % Allocate space for results
    MS = zeros(size(F,1),size(yall,1));
    PS = zeros(size(F,1),size(F,1),size(yall,1));
    ttau = zeros(D+N,size(yall,1));
    tnu = zeros(D+N,size(yall,1));
    
    lZ = zeros(1,size(yall,1));
    R = zeros(D+N,size(yall,1));
    
    % Time step, dt
    dt = 1;
    
    % Discrete-time model
    [A,Q] = lti_disc(F,L,Qc,dt);
    
    % Iterate EP
    for itt = 1:ep_itts
        
      % Set initial state
      m = zeros(size(F,1),1);
      P = Pinf;
      
      % Track convergence (not very optimal)
      maxDiffP = 0;
      maxDiffM = 0;
      PSP = PS;
      MSP = MS;
      
      % ### Forward filter
      for k=1:numel(yall)
          
          % Prediction step
          if k > 1
            m = A * m;
            P = A * P * A' + Q;
          end
          
          % Update step
          if ~isnan(yall(k))
              
              % Latent marginal, cavity distribution
              fmu = H*m; W = P*H'; HPH = diag(H*P*H');
              
              % do normal update on first pass
              if itt == 1 || k == numel(yall)
                  
                  % Propagate moments through likelihood
                  [lZ(k),dlZ,d2lZ] = mom(lik_param,fmu,HPH,Wnmf,yall,k);
                  
                  % Perform moment matching
                  ttau(:,k) = (1-ep_damping)*ttau(:, k) + ep_damping/ep_fraction*(-d2lZ'./(1+d2lZ'.*HPH));
                  tnu(:,k) = (1-ep_damping)*tnu(:, k) + ep_damping/ep_fraction*((dlZ'-fmu.*d2lZ')./(1+d2lZ'.*HPH));
                  
                  % This is the equivalent measurement noise
                  R(:,k) = 1 ./ ttau(:,k); %-(1+d2lZ'.*HPH)./d2lZ';
                  
              end
              
              % Enforce positivity->lower bound ttau by zero
              ttau(:,k) = max(ttau(:,k),0);
              
              if min(ttau(:,k))==0
                  z = ttau(:,k).*HPH+1;
                  K = bsxfun(@times,W,(ttau(:,k)./z)');
                  v = ttau(:,k).*fmu - tnu(:,k);
                  m = m - W*(v./z);
                  P = P - K*W';
              else
                  K = bsxfun(@rdivide,W,(HPH+1./ttau(:,k))');
                  v = tnu(:,k)./ttau(:,k) - fmu;
                  m = m + K*v;
                  P = P - K*H*P;
              end
              
          end
          
          % Store filter estimate
          MS(:,k)   = m;
          PS(:,:,k) = P;
          
      end
      
      % Output debugging info
      if nargout>5
          out.tnu = tnu;
          out.ttau = ttau;
          out.lZ = lZ;
          out.R = R;
          out.MF = MS;
          out.PF = PS;
      end
      
      % ### Backward smoother and EP step
      
      % Rauch-Tung-Striebel smoother
      for k=size(MS,2)-1:-1:1
          
          % Smoothing step (using Cholesky for stability)
          PSk = PS(:,:,k);
          
          % Pseudo-prediction
          PSkp = A*PSk*A'+Q;
          
          % Solve the Cholesky factorization
          [L,notposdef] = chol(PSkp,'lower');
          
          % Numerical problems in Cholesky, retry with jitter
          if notposdef>0
              jitterSigma2 = 1e-4;
              jitter = sqrt(jitterSigma2)*diag(rand(size(A,1),1));
              L = chol(PSkp+jitter,'lower');
          end
          
          % Continue smoothing step
          G = PSk*A'/L'/L;
          
          % Do update
          m = MS(:,k) + G*(m-A*MS(:,k));
          P = PSk + G*(P-PSkp)*G';
          
          % Store estimate
          MS(:,k)   = m;
          PS(:,:,k) = P;
          
          if itt < ep_itts
              if ~isnan(yall(k))
          
                  % ### EP ###

                  % Get marginal
                  m_marginal = H*m;
                  v_marginal = diag(H*P*H');

                  % Compute cavity
                  v_cav = 1./(1./v_marginal - ep_fraction*ttau(:, k));
                  m_cav = v_cav.*(m_marginal./v_marginal - ep_fraction*tnu(:, k));
                  
                  % check which sites to update
                  update_idx = v_cav > 0;

                  % Compute gradients of the normalizer of the tilted dst
                  [~,dlZ_,d2lZ_] = mom(lik_param,m_cav,v_cav,Wnmf,yall,k);

                  % Moment matching
                  ttau(update_idx,k) = (1-ep_damping)*ttau(update_idx, k) + ep_damping/ep_fraction*(-d2lZ_'./(1+d2lZ_'.*v_cav));
                  tnu(update_idx,k) = (1-ep_damping)*tnu(update_idx, k) + ep_damping/ep_fraction*((dlZ_'-m_cav.*d2lZ_')./(1+d2lZ_'.*v_cav));

                  % This is the equivalent measurement noise
                  R(:,k) = 1 ./ ttau(:,k); %-(1+d2lZ_'.*v_cav)./d2lZ_';

                  % Max diff in m and P
                  maxDiffM = max(maxDiffM,max(max(abs(H*MSP(:,k)-H*m))));
                  maxDiffP = max(maxDiffP,max(max(abs(H*PSP(:,:,k)*H'-H*P*H'))));
                  
              end
          end
          
      end % end smoother iteration
    
      if itt < ep_itts
        fprintf('%.02i - max diff in m: %.6g - max diff in P: %.6g - nll: %.6g\n', ...
            itt,maxDiffM,maxDiffP,-sum(lZ))
      end
      
    end % end EP iteration
    
    % Output debugging info
    if nargout>5
      out.tnu = tnu;
      out.ttau = ttau;
      out.lZ = lZ;
      out.R = R;
    end
    
    % Estimate the joint covariance matrix if requested
    if nargout > 2
      
      % Allocate space for results
      %Covft = zeros(size(PS,3));
      Covft = [];
          
      % Lower triangular
      %{
      for k = 1:size(PS,3)-1
        GSS = GS(:,:,k);
        for j=1:size(PS,3)-k
          Covft(k+j,k) = H*(GSS*PS(:,:,k+j))*H';
          GSS = GSS*GS(:,:,k+j);
        end
      end
      %}
    
    end
    
    out.MS = MS;
    out.PS = PS;
  
    % These indices shall remain to be returned
    MS = MS(:,return_ind);
    PS = PS(:,:,return_ind);
    
    % Return mean
    Eft = H*MS;
    
    % Return variance
    if nargout > 1
        Varft = zeros(D+N,size(MS,2));
        for k=1:size(MS,2)
            Varft(:,k)  = diag(H*PS(:,:,k)*H');
        end
    end
    
    % Return values
    varargout = {Eft,Varft};
 
    % Also return joint covariance and upper/lower 95% bounds
    if nargout > 3
        
        % Join upper triangular and variances
        %if ~filteronly
        %    Covft = Covft(return_ind,return_ind);
        %    Covft = Covft+Covft'+diag(Varft(:));
        %else
        %    Covft = [];
        %end
        
        % The bounds
        lb = Eft - 1.96*sqrt(Varft);
        ub = Eft + 1.96*sqrt(Varft);
        varargout = {Eft,Varft,Covft,lb,ub,out};
        
    end

  end
  
  
%% Evaluate negative log marginal likelihood and its gradient

  if isempty(xt)
  
    error('this mixture script is not for training')
    
  end
  
  
  
  
  
  