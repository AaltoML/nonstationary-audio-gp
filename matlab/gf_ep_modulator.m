function [varargout] = gf_ep_modulator(w,x,y,ss,mom,xt,kernel1,kernel2,num_lik_params)
% GF_EP_MODULATOR - Solve modulator GP model by Expectation Propagation
%
% Syntax:
%   [...] = gf_adf(w,x,y,k,xt)
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

  % Log transformed parameters
  lik_param = w(1:num_lik_params);
  param = exp(w(num_lik_params+1:end));
  
  % Form the state space model
  [F,L,Qc,H,Pinf,dF,dQc,dPinf] = ss(x,param,kernel1,kernel2);
  D = sum(H,[1,2]);
  
  
%% Prediction of test inputs (filtering and smoothing)

  % Check that we are predicting
  if ~isempty(xt)

    % Set initial state
    m = zeros(size(F,1),1);
    P = Pinf;
    
    % Allocate space for results
    MS = zeros(size(m,1),size(yall,1));
    PS = zeros(size(m,1),size(m,1),size(yall,1));
    MP = MS;
    PP = PS;
    ttau = zeros(D,size(yall,1));
    tnu = zeros(D,size(yall,1));
    lZ = zeros(1,size(yall,1));
    R = zeros(D,size(yall,1));
    fs2 = zeros(D,size(yall,1));
    fmus = zeros(D,size(yall,1));

    
    % Initial dt
%     dt = inf;
    dt = 1;
    
    % ### Forward filter
    [A,Q] = lti_disc(F,L,Qc,dt);
    
    % The filter recursion
    for k=1:numel(yall)
%         keyboard
        % Prediction step
        m = A * m;
        P = A * P * A' + Q;

        % Store predicted mean
        MP(:,k) = m;
        PP(:,:,k) = P;
        
        % Update step
        if ~isnan(yall(k))
            
            % Latent marginal, cavity distribution
            fmu = H*m; W = P*H'; HPH = diag(H*P*H');
            
            % Propagate moments through likelihood 
            [lZ(k),dlZ,d2lZ] = mom(lik_param,fmu,HPH,yall,k);
            
            % Perform moment matching (in natural parameter space)
            ttau(:,k) = -d2lZ'./(1+d2lZ'.*HPH); % tilted precision
            tnu(:,k) = (dlZ'-fmu.*d2lZ')./(1+d2lZ'.*HPH); % tilted precision-adjusted mean
            
            % This is the equivalent measurement noise
            R(:,k) = -(1+d2lZ'.*HPH)./d2lZ';
            
            fs2(:,k) = HPH;
            fmus(:,k) = fmu;
            
            % Enforce positivity->lower bound ttau by zero
            ttau(:,k) = max(ttau(:,k),0);
            
            if min(ttau(:,k))==0
%               warning('Moment matching hit bound.')
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
        
        % Store estimate
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
      out.fs2 = fs2;
      out.fmu = fmus;
    end
    
    % ### Backward smoother
      
    % Allocate space for storing the smoother gain matrix
    GS = zeros(size(F,1),size(F,2),size(yall,1));

    % Rauch-Tung-Striebel smoother
    for k=size(MS,2)-1:-1:1
        
      MSk = MS(:,k);
      MSkp = A*MSk;
      
      % Smoothing step (using Cholesky for stability)
      PSk = PS(:,:,k);
      
      % Pseudo-prediction
%       PSkp = A(:,:,k+1)*PSk*A(:,:,k+1)'+Q(:,:,k+1);
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
%       G = PSk*A(:,:,k+1)'/L'/L;
      G = PSk*A'/L'/L;
      
      % Do update
%       m = MS(:,k) + G*(m-A(:,:,k+1)*MS(:,k));
      m = MSk + G*(m-MSkp);
      P = PSk + G*(P-PSkp)*G';
      
      
      %%%%%% ATTEMPT AT BACKWARDS EP UPDATE %%%%%%
      %%%%%%  FOLLOWING S3.3 OF QI & MINKA  %%%%%%
      % convert to natural parameters
      post_tau = 1 ./ diag(H*P*H');
      post_nu = H*m .* post_tau;
      % remove site to obtain cavity
      cav_tau = post_tau - ttau(:,k);
      cav_nu = post_nu - tnu(:,k);
      % convert to moment parameters
      cav_v = 1 ./ cav_tau;
      cav_m = cav_nu ./ cav_tau;
      
      % Propagate moments through likelihood 
      [lZ_k,dlZ_k,d2lZ_k] = mom(lik_param,cav_m,cav_v,yall,k);
      
      % Perform moment matching (in natural parameter space)
      ttau(:,k) = -d2lZ_k'./(1+d2lZ_k'.*cav_v); % tilted precision
      tnu(:,k) = (dlZ_k'-cav_m.*d2lZ_k')./(1+d2lZ_k'.*cav_v); % tilted precision-adjusted mean

      % Enforce positivity->lower bound ttau by zero
      ttau(:,k) = max(ttau(:,k),0);
      
      W = P*H'; % is this right?
      % Is it correct to just run the same update as in the filtering
      % stage? Qi & Minka seems to suggest so.
      if min(ttau(:,k))==0
    %               warning('Moment matching hit bound.')
        z = ttau(:,k).*cav_v+1;
        K = bsxfun(@times,W,(ttau(:,k)./z)');
        v = ttau(:,k).*cav_m - tnu(:,k);
        m = m - W*(v./z);
        P = P - K*W';
      else
        K = bsxfun(@rdivide,W,(cav_v+1./ttau(:,k))');
        v = tnu(:,k)./ttau(:,k) - cav_m;
        m = m + K*v;
        P = P - K*H*P;
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % Store estimate
      MS(:,k)   = m;
      PS(:,:,k) = P;
      GS(:,:,k) = G;
%       keyboard
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
        Varft = zeros(D,size(MS,2));
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
%     nparam = numel(param);
    nparam = 0;  % TEMPORARY
    steps = numel(yall);
            
    % Allocate space for results
%     edata = 0;
    gdata = zeros(1,length(w));
    
    % Set up
    dt = -inf;
    
    ttau = zeros(D,size(yall,1));
    tnu = zeros(D,size(yall,1));
    lZ = zeros(1,size(yall,1));
%     dz =  zeros(size(yall,1),nparam);
        
    % Allocate space for expm results
    AA = zeros(2*d,2*d,nparam);
    
    % Loop over all observations
    for k=1:steps
        
        % The time discretization step length
        if (k>1)
            dt = xall(k)-xall(k-1);
        else
            dt = 0;
        end
        
        % Set predicted m and P
        if k>1
            m = A*m;
            P = A*P*A' + Q;
        end
        
        % Latent marginal, cavity distribution
        fmu = H*m; W = P*H'; fs2 = diag(H*P*H');

        % Propagate moments through likelihood
        [lZ(k),dlZ,d2lZ] = mom(lik_param,fmu,fs2,yall,k);
                
        % Perform moment matching
        ttau(:,k) = -d2lZ'./(1+d2lZ'.*fs2);
        tnu(:,k) = (dlZ'-fmu.*d2lZ')./(1+d2lZ'.*fs2);
        
        % Enforce positivity->lower bound ttau by zero
        ttau(:,k) = max(ttau(:,k),0);
      
        % Gauss: H*m~N(r,W^-1) ttau = W, tnu = W*r, r=y-m(x)
        z = ttau(:,k).*fs2+1;
        K = bsxfun(@times,W,(ttau(:,k)./z)');
        v = ttau(:,k).*fmu - tnu(:,k);
        
        % Do update
        m = m - W*(v./z);
        P = P - K*W';
        
    end
    
    % Sum up the lik
    edata = -sum(lZ);
    
    % Account for log-scale
    gdata = gdata.*exp(w);
    
    % Return negative log marginal likelihood and gradient
    varargout = {edata,gdata};

  end
  
end
  
  