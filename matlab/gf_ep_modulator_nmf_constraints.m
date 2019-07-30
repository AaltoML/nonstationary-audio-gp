function [varargout] = gf_ep_modulator_nmf_constraints(w,x,y,ss,mom,xt,kernel1,kernel2,num_lik_params,D,N,...
                                                       ep_fraction,ep_damping,ep_itts,constraints,w_fixed,tune_hypers)
% GF_EP_MODULATOR_NMF_CONSTRAINTS - Solve TF-NMF GP model by Power EP with
% constraints on optimisation
%
% Syntax:
%   [...] = gf_ep_modulator_nmf_constraints(w,x,y,ss,mom,xt,kernel1,kernel2,num_lik_params,D,N,ep_fraction,ep_damping,ep_itts)
%
% In:
%   w              - Log-parameters
%   x              - Training inputs
%   y              - Training outputs
%   ss             - State space model function handle, [F,L,Qc,...] = @(x,theta) 
%   mom            - Moment calculations for EP
%   xt             - Test inputs (default: empty)
%   kernel1        - kernel for subbands
%   kernel2        - kernel for modulators
%   num_lik_params - number of hypers. in likelihood model
%   D              - number of subbands
%   N              - number of modulators
%   ep_fraction    - the power in Power EP
%   ep_damping     - how much to damp the EP updates
%   ep_itts        - number of EP iterations
%   constraints    - a set of constraints
%   w_fixed        - the fixed model parameters
%   tune_hypers    - whether or not we are optimising the parameters
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
%   TODO
%
% References:
%
%   [1] William Wilkinson, Michael Riis Andersen, Josh Reiss, Dan Stowell, Arno Solin (2019) 
%       End-to-End Probabilistic Inference for Nonstationary Audio Analysis. 
%       International Conference on Machine Learning (ICML). 
%
%  This software is distributed under the Apache License, Version 2.0

%% Check defaults

  % Is there test data
  if nargin < 6, xt = []; end
  if nargin < 12, ep_fraction = 0.5; end
  if nargin < 13, ep_itts = 30; end
   
  
%% Figure out the correct way of dealing with the data

  % Combine observations and test points
  xall = [x(:); xt(:)];
  yall = [y(:); nan(numel(xt),1)];
    
  % Make sure the points are unique and in ascending order
  [~,sort_ind,return_ind] = unique(xall,'first');
  yall = yall(sort_ind);
  
  % Only return test indices
  return_ind = return_ind(end-numel(xt)+1:end);
  
  
%% Set up model
  % unwrap the optimisable and fixed parameters with constraints
  w_ind = 1; wf_ind = 1;  
  if tune_hypers(1)
    lik_param = w(1:num_lik_params);
    w_ind = w_ind + num_lik_params;
  else
    lik_param = w_fixed(1:num_lik_params);
    wf_ind = wf_ind + num_lik_params;
  end 
  param1 = [];
  param2 = [];
  for i=2:6
    if tune_hypers(i)
      if i<=4
        param1 = [param1; sigmoid(w(w_ind:w_ind+D-1),constraints(i-1,:))];
        w_ind = w_ind + D;
      else
        param2 = [param2; sigmoid(w(w_ind:w_ind+N-1),constraints(i-1,:))];
        w_ind = w_ind + N;
      end
    else
      if i<=4
        param1 = [param1; sigmoid(w_fixed(wf_ind:wf_ind+D-1),constraints(i-1,:))];
        wf_ind = wf_ind + D;
      else
        param2 = [param2; sigmoid(w_fixed(wf_ind:wf_ind+N-1),constraints(i-1,:))];
        wf_ind = wf_ind + N;
      end
    end
  end
  
  if tune_hypers(7)
    Wnmf = reshape(sigmoid(w(w_ind(1):end),constraints(6,:)),[D,N]);
  else
    Wnmf = reshape(sigmoid(w_fixed(wf_ind(1):end),constraints(6,:)),[D,N]);
  end
  
  % Form the state space model
  [F,L,Qc,H,Pinf,~,~,~] = ss(x,param1,param2,kernel1,kernel2);
  
  if true            % balance state space model for improved numerical stability
    [T,F] = balance(F); L = T\L; H = H*T;                     % balance F,L,Qc,H
    LL = T\chol(Pinf,'lower'); Pinf = LL*LL';                     % balance Pinf
%     for j=1:size(dF,3)                                    % balance dF and dPinf
%       dF(:,:,j) = T\dF(:,:,j)*T; dPinf(:,:,j) = T\dPinf(:,:,j)/T;
%     end
  end
  
  
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
    
    ep_damp = ep_damping(1);
    
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
                  [lZ(k),dlZ,d2lZ] = mom(lik_param,fmu,HPH,Wnmf,1,yall,k);
                  
                  % Perform moment matching
                  ttau(:,k) = (1-ep_damp)*ttau(:, k) + ep_damp*(-d2lZ'./(1+d2lZ'.*HPH));
                  tnu(:,k) = (1-ep_damp)*tnu(:, k) + ep_damp*((dlZ'-fmu.*d2lZ')./(1+d2lZ'.*HPH));
                  
                  % Enforce positivity->lower bound ttau by zero
                  ttau(:,k) = max(ttau(:,k),0);
                  
                  % This is the equivalent measurement noise
                  R(:,k) = 1 ./ ttau(:,k); %-(1+d2lZ'.*HPH)./d2lZ';
                  
              end
              
              % Commented out the old 1D update specific implementation
              %{
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
              %}
              
              % Pick which to be updated differently
              ii = (ttau(:,k)==0);
              
              % If any tau==0
              if any(ii)
                  %fprintf('Model %i hit bound.\n',find(ii))
                  z = ttau(ii,k).*HPH(ii)+1;
                  K = bsxfun(@times,W(:,ii),(ttau(ii,k)./z)');
                  v = ttau(ii,k).*fmu(ii) - tnu(ii,k);
                  m = m - W(:,ii)*(v./z);
                  P = P - K*W(:,ii)';
              end
              % If any normal cases
              if any(~ii)
                  K = bsxfun(@rdivide,W(:,~ii),(HPH(~ii)+1./ttau(~ii,k))');
                  v = tnu(~ii,k)./ttau(~ii,k) - fmu(~ii);
                  m = m + K*v;
                  P = P - K*H(~ii,:)*P;
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
      
      if itt < ep_itts
        ep_damp = ep_damping(itt+1);
      end
      
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
                  [lZ(k),dlZ,d2lZ] = mom(lik_param,m_cav,v_cav,Wnmf,ep_fraction,yall,k);

                  % Moment matching
                  ttau(update_idx,k) = (1-ep_damp)*ttau(update_idx, k) + ...
                                       ep_damp*ep_fraction*(-d2lZ(update_idx)'./(1+d2lZ(update_idx)'.*v_cav(update_idx)));
                  tnu(update_idx,k) = (1-ep_damp)*tnu(update_idx, k) + ...
                                      ep_damp*ep_fraction*((dlZ(update_idx)'-m_cav(update_idx).*d2lZ(update_idx)')./(1+d2lZ(update_idx)'.*v_cav(update_idx)));
                  
                  % Enforce positivity->lower bound ttau by zero
                  ttau(:,k) = max(ttau(:,k),0);
                  
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
  
    % Size of inputs
    d = size(F,1);
            
    % Allocate space for results
    gdata = zeros(1,length(w));
    
    % Set up
    dt = 1;
    
    ttau = zeros(D+N,size(yall,1));
    tnu = zeros(D+N,size(yall,1));
    lZ = zeros(1,size(yall,1));
    MS = zeros(d,size(yall,1));
    if ep_itts > 1
        PS = zeros(d,d,size(yall,1));
    end
    
    % Allocate space for expm results
%     A = expm(F*dt);
%     Q  = Pinf - A*Pinf*A';
    [A,Q] = lti_disc(F,L,Qc,dt);
    
    ep_damp = ep_damping(1);
    
    % Iterate EP
    for itt = 1:ep_itts
        
      % Track convergence (not very optimal)
%       maxDiffP = 0;
%       maxDiffM = 0;
%       PSP = PS;
%       MSP = MS;
        
      % Set initial state
      m = zeros(d,1);
      P = Pinf;
      
      if itt == 1 || itt < ep_itts
        % Loop over all observations
        for k=1:numel(yall)
            if k>1
                m = A*m;
                P = A*P*A' + Q;
            end
            
            if ~isnan(yall(k))
                % Latent marginal, cavity distribution
                fmu = H*m; W = P*H'; fs2 = diag(H*P*H');

                if min(fs2) <= 0
                    keyboard
                end

                % do normal update on first pass
                if itt == 1 || k == numel(yall)

                    % Propagate moments through likelihood
                    [lZ(k),dlZ,d2lZ] = mom(lik_param,fmu,fs2,Wnmf,1,yall,k);

                    % Perform moment matching
                    ttau(:,k) = (1-ep_damp)*ttau(:, k) + ep_damp*(-d2lZ'./(1+d2lZ'.*fs2));
                    tnu(:,k) = (1-ep_damp)*tnu(:, k) + ep_damp*((dlZ'-fmu.*d2lZ')./(1+d2lZ'.*fs2));

                end

                % Enforce positivity->lower bound ttau by zero
                ttau(:,k) = max(ttau(:,k),0);

                % Gauss: H*m~N(r,W^-1) ttau = W, tnu = W*r, r=y-m(x)
                if min(ttau(:,k))==0
                    z = ttau(:,k).*fs2+1;
                    K = bsxfun(@times,W,(ttau(:,k)./z)');
                    v = ttau(:,k).*fmu - tnu(:,k);
                    m = m - W*(v./z);
                    P = P - K*W';
                else
                    K = bsxfun(@rdivide,W,(fs2+1./ttau(:,k))');
                    v = tnu(:,k)./ttau(:,k) - fmu;
                    m = m + K*v;
                    P = P - K*H*P;
                end
            end
            
            if itt < ep_itts
                % Store filter estimate
                MS(:,k)   = m;
                PS(:,:,k) = P;
            end

        end
      end
      
      if itt < ep_itts
          % ### Backward smoother and EP step
          
            ep_damp = ep_damping(itt+1);

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
                    [lZ(k),dlZ,d2lZ] = mom(lik_param,m_cav,v_cav,Wnmf,ep_fraction,yall,k);

                    % Moment matching
                    ttau(update_idx,k) = (1-ep_damp)*ttau(update_idx, k) + ...
                                         ep_damp*ep_fraction*(-d2lZ(update_idx)'./(1+d2lZ(update_idx)'.*v_cav(update_idx)));
                    tnu(update_idx,k) = (1-ep_damp)*tnu(update_idx, k) + ...
                                        ep_damp*ep_fraction*((dlZ(update_idx)'-m_cav(update_idx).*d2lZ(update_idx)')./(1+d2lZ(update_idx)'.*v_cav(update_idx)));

                end

                % Max diff in m and P
    %             maxDiffM = max(maxDiffM,max(max(abs(H*MSP(:,k)-H*m))));
    %             maxDiffP = max(maxDiffP,max(max(abs(H*PSP(:,:,k)*H'-H*P*H'))));

            end % end smoother iteration

    %         fprintf('%.02i - max diff in m: %.6g - max diff in P: %.6g - nll: %.6g\n', ...
    %                 itt,maxDiffM,maxDiffP,-sum(lZ))
      end
      
    end % end EP iteration
    
    % Sum up the lik
    edata = -sum(lZ);
    
    % Account for log-scale
%     gdata = gdata.*exp(w);
    
    % Return negative log marginal likelihood and gradient
    varargout = {edata,gdata};
    
  end
  
  
  
  
  
  