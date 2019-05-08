function [varargout] = gf_giekf_modulator_nmf_constraints(w,x,y,ss,mom,xt,kernel1,kernel2,num_lik_params,D,N,...
                                                          g_iter,l_iter,constraints,w_fixed,tune_hypers,GradObj)
% gf_giekf_modulator_nmf - Solve NMF model by extended Kalman filtering
%
% Syntax:
%   [...] = gf_giekf_modulator_nmf(w,x,y,k,xt,...)
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
  if nargin < 6, xt = []; end
   
  
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
  [F,L,Qc,H,Pinf,dF,dQc,dPinf] = ss(x,param1,param2,kernel1,kernel2);

  % Measurement noise variance
  sigma2 = exp(lik_param);
  R = sigma2;
  
  % Concatenate derivatives
  dF    = cat(3,zeros(size(F)),dF);
  dQc   = cat(3,zeros(size(Qc)),dQc);
  dPinf = cat(3,zeros(size(Pinf)),dPinf);
  dR    = zeros(1,1,size(dF,3)); dR(1) = 1;
  
  % Measurement model
  linkf = @(x) log(1+exp(x));
  dlinkf = @(x) exp(x)./(exp(x)+1);
  d2linkf = @(x) dlinkf(x) .* (1 - dlinkf(x));
  handle = @(x,p) funh(x,H,linkf,D,N,Wnmf);
  dhandle = @(x,p) funhd(x,H,linkf,dlinkf,D,N,Wnmf);
  d2handle = @(x,p) funhd2(x,H,linkf,dlinkf,d2linkf,D,N,Wnmf);
  
  
%% Prediction of test inputs (filtering and smoothing)

  % Check that we are predicting
  if ~isempty(xt)
    
    % Allocate space for results
    MS = zeros(size(F,1),size(yall,1));
    PS = zeros(size(F,1),size(F,1),size(yall,1));
%     ttau = zeros(D+N,size(yall,1));
%     tnu = zeros(D+N,size(yall,1));
    
%     lZ = zeros(1,size(yall,1));
    R = zeros(D+N,size(yall,1));
    
    % Time step, dt
    dt = 1;
    
    % Discrete-time model
    [A,Q] = lti_disc(F,L,Qc,dt);
    
    % Iterate EKF
    for itt = 1:g_iter
        
      % Set initial state
      if itt==1         
        m = zeros(size(F,1),1);
      end
      P = Pinf;
      
 
      % Track convergence (not very optimal)
      maxDiffP = 0;
      PSP = PS;
      
      % ### Forward filter
      for k=1:numel(yall)
          
          % Prediction step
          if k>1
            m = A * m;
            P = A * P * A' + Q;
          end
          
          % Update step
          if ~isnan(yall(k))
              
              % Latent marginal, cavity distribution
              fmu = H*m; W = P*H'; HPH = diag(H*P*H');
              
              % Check derivative
              %[D0,D1] = der_check(@(x) handle(x,[]),@(x) dhandle(x,[]),1,10*randn(size(m)))
              %D0(:)-D1(:)
              
              
              % Do EKF update
%               [m,P] = ekf_update1(m,P,yall(k),dhandle,sigma2,handle,[],[]);
              [m,P] = iekf_update1(m,P,yall(k),dhandle,sigma2,handle,[],[],l_iter);

              
          end
          
          % Store filter estimate
          MS(:,k)   = m;
          PS(:,:,k) = P;
          
      end
      
      % Output debugging info
      if nargout>5
%           out.tnu = tnu;
%           out.ttau = ttau;
%           out.lZ = lZ;
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
          
          % Max diff in P
          maxDiffP = max(maxDiffP,max(max(abs(H*PSP(:,:,k)*H'-H*P*H'))));
          
      end % end smoother iteration
    
      maxDiffP
      %ll=-sum(lZ)
      
    end % end EKF iteration
    
    % Output debugging info
    if nargout>5
%       out.tnu = tnu;
%       out.ttau = ttau;
%       out.lZ = lZ;
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
    nparam = size(dF,3);
%     nparam = length(w);
    steps = numel(yall);
    
    % Allocate space for results
    edata = 0;
    gdata = zeros(1,length(w));
    
    % Set up
    Z  = zeros(d);
    m  = zeros(d,1);
    P  = Pinf;
    dm = zeros(d,nparam);
    dP = dPinf;

    % Equal time steps
    dt = 1;
    
    % Discrete-time model
%     [A,Q] = lti_disc(F,L,Qc,dt);
  
    if strcmp(GradObj,'off')
        nparam = 0;
    end
  
    % Allocate space for expm results
    AA = zeros(2*d,2*d,max(nparam,1));
    
    % Loop through all parameters (Kalman filter prediction step)
    for j=1:nparam
        
        % The first matrix for the matrix factor decomposition
        FF = [ F        Z;
              dF(:,:,j) F];
        
        % Solve the matrix exponential
        AA(:,:,j) = expm(FF*dt);
        
    end
    
%     A = AA(1:d,1:d,1);
%     Q  = Pinf - A*Pinf*A';
    A = expm(F*dt);
    Q = Pinf - A*Pinf*A';
    
    % Loop over all observations
    for k=1:steps
        
        % Loop through all parameters (Kalman filter prediction step)
        for j=1:nparam
           
            % Solve the differential equation
            dm(:,j) = AA(d+(1:d),:,j)*[m; dm(:,j)];
            
            % The derivatives of A and Q
            dA = AA(d+1:end,1:d,j);
            dAPinfAt = dA*Pinf*A';
            dQ = dPinf(:,:,j) - dAPinfAt - A*dPinf(:,:,j)*A' - dAPinfAt';
            
            % The derivatives of P
            dAPAt = dA*P*A';
            dP(:,:,j) = dAPAt + A*dP(:,:,j)*A' + dAPAt' + dQ;
            
        end
        
        % Set predicted m and P
        m = A*m;
        P = A*P*A' + Q;
        
        % Evaluate measurement model
        mu = handle(m,[]);
        JH = dhandle(m,[]);
        dJH = d2handle(m,[]);
        
        % Start the Kalman filter update step and precalculate variables
        
        S = JH*P*JH' + R;
        
        [LS,notposdef] = chol(S,'lower');
        
        % If matrix is not positive definite, add jitter
        if notposdef>0
            jitterSigma2 = 1e-4;
            jitter = jitterSigma2*diag(rand(size(S,1),1));
            [LS,notposdef] = chol(S+jitter,'lower');
            
            % Return nan to the optimizer
            if notposdef>0
                varargout = {nan*edata,nan*gdata};
                return;
            end
        end
        
        % Continue update
        HtiS = JH'/LS/LS';
        iS   = eye(size(S))/LS/LS';
        K    = P*HtiS;
        v    = yall(k) - mu;
        vtiS = v'/LS/LS';
        
        % Loop through all parameters (Kalman filter update step derivative)
        for j=1:nparam
            
            % Innovation covariance derivative
            if j <= nparam - D*N
                dmdJH = dm(:,j)'*dJH;
            else
                W_ = zeros(size(Wnmf));
                W_(j - nparam + D*N) = 1;
                dmdJH = funhd(m,H,linkf,dlinkf,D,N,W_);
%                 dmdJH'
%                 keyboard
            end
            
            dS = dmdJH*P*JH' + JH*dP(:,:,j)*JH' + JH*P*dmdJH' + dR(:,:,j);
            
            % Evaluate the energy derivative for j (optimized from above)
            gdata(j) = gdata(j) ...
                + .5*sum(iS(:).*dS(:)) ...
                - .5*(JH*dm(:,j))*vtiS' ...
                - .5*vtiS*dS*vtiS'     ...
                - .5*vtiS*(JH*dm(:,j));
            
            % Kalman filter update step derivatives
            dK        = dP(:,:,j)*HtiS + P*dmdJH'/LS/LS' - P*HtiS*dS/LS/LS';
            dm(:,j)   = dm(:,j) + dK*v - K*JH*dm(:,j);
            dKSKt     = dK*S*K';
            dP(:,:,j) = dP(:,:,j) - dKSKt - K*dS*K' - dKSKt';
            
        end
        
        % Evaluate the energy
        edata = edata + .5*size(S,1)*log(2*pi) + sum(log(diag(LS))) + .5*vtiS*v;
        
        % Finish Kalman filter update step
        m = m + K*v;
        P = P - K*S*K';        
        
%         keyboard
    end
    
    % Account for log-scale
%     ww = w(1:end-D*N);
%     ww = w(1:end);
%     gdata = gdata.*exp(ww(:)');
    % Return negative log marginal likelihood and gradient
    varargout = {edata,gdata};

  end
  

end

% For assumed density filtering methods
function [y] = funh(x,H,linkf,D,N,W)
  
  y = (H(1:D,:) * x)' * W * linkf(H(D+1:D+N,:) * x);
  
end

% For EKF
function [dy] = funhd(x,H,linkf,dlinkf,D,N,W)

  partials = [W * linkf(H(D+1:D+N,:) * x); ...
         diag((H(1:D,:) * x)' * W) * dlinkf(H(D+1:D+N,:) * x)];
  dy = partials'*H;
  
end
  
function [d2y] = funhd2(x,H,linkf,dlinkf,d2linkf,D,N,W)

  z = H(1:D,:) * x;
  g = H(D+1:D+N,:) * x;
  partials = [zeros(size(z,1)),            W.*(ones(D,1)*dlinkf(g)'); ...
          (W.*(ones(D,1)*dlinkf(g)'))', diag((z'*W)' .* d2linkf(g))];
  d2y = H'*partials*H;
  
end
