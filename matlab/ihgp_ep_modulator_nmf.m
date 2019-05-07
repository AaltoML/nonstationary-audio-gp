function [varargout] = ihgp_ep_modulator_nmf(w,x,y,ss,mom,xt,kernel1,kernel2,num_lik_params,D,N,ep_fraction,ep_damping,ep_itts)
% IHGP_ADF - Infinite-horizon GP with single-sweep EP (ADF)
%
% Syntax:
%   [...] = ihgp_adf(w,x,y,k,xt)
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
%   Covft - Predicted joint covariance matrix (disabled output)
%   lb    - 95% confidence lower bound
%   ub    - 95% confidence upper bound
%
% Description:
%   Consider the general GP modelling problem:
%
%     f(t) ~ GP(0,k(x,x')),
%      y_i ~ p(y_i | f(x_i)),  i=1,2,...,n,
%
%   where k(x,x') = k_theta(x,x'). This function performs assumed density
%   filtering (ADF) / single-sweep expectation propagtaion for dealing with
%   non-Gaussian likelihoods (observation models).
%
%   The state space model is specified by the function handle 'ss' such
%   that it returns the state space model matrices
%
%     [F,L,Qc,H,Pinf,dF,dQc,dPinf] = ss(x,theta),
%
%   where theta holds the hyperparameters. See the paper for details.
%   This code is assuming a stationary covariance function even though the
%   methodology per se does not require it.
%
%   NOTE: This code is proof-of-concept, not optimized for speed.
%
% References:
%
%   [1] Arno Solin, James Hensman, and Richard E. Turner (2018). 
%       Infinite-horizon Gaussian processes. Advances in Neural 
%       Information Processing Systems (NIPS). Montreal, Canada. 
%
% Copyright:
%   2018 Arno Solin
%
%  This software is distributed under the GNU General Public
%  License (version 3 or later); please refer to the file
%  License.txt, included with the software, for details.

%%% Check defaults

  % Is there test data
  if nargin < 6, xt = []; end
  if nargin < 12, ep_fraction = 0.5; end
  if nargin < 14, ep_itts = 30; end
   
  
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
  param1 = exp(w(num_lik_params+1:num_lik_params+3*D));
  param2 = exp(w(num_lik_params+3*D+1:num_lik_params+3*D+2*N));
  Wnmf = reshape(exp(w(num_lik_params+3*D+2*N+1:end)),[D,N]);
  
  % Form the state space model
  [F,L,Qc,H,Pinf] = ss(x,param1,param2,kernel1,kernel2);
  
  
%% Solve a bunch of DAREs

    % Set A and Q (this part only works for stationary covariance functions)
    dt = 1; %xall(2)-xall(1);
    %A = expm(F*dt);
    %Q = Pinf - A*Pinf*A';
    [A,Q] = lti_disc(F,L,Qc,dt);
    Q = (Q+Q')/2;
     
    % Cell array of stationary covs
    PPlist = {};
    PPlisto = cell(N+D,1);%{};
    
    % Find indices    
    ilist = [find(sum(H,1)) size(H,2)+1];   

    % Set up each separate model
    for n=1:(N+D)
    
    ii = ilist(n):ilist(n+1)-1;
    
    % Set up interpolation setup
    ro = logspace(-2,4,32)';
    PPlisto{n} = nan(numel(ro),numel(Q(ii,ii)));
    for j=1:numel(ro)
      try
        [PP,~,~,~] = dare(A(ii,ii)',H(n,ii)',Q(ii,ii),ro(j));
        PPlisto{n}(j,:) = PP(:)';
      catch
        warning('Failed forward DARE calculation %i.',j)
        if false
          varargout = {nan,nan*w};
          return
        else
          ro(j) = nan;
        end          
      end
    end
    PPlisto{n}(isnan(ro),:)=[]; ro(isnan(ro)) = [];
    
    % Interpolation for lookup (cubic)
    r = logspace(-2,4,200)';
    U = apxGrid('interp',{ro},r,3);
    PPlist{n} = U*PPlisto{n}; %#ok
    end % End setting up
    
    % Output
    if nargout>5
      out.r = r;
      out.ro = ro;
      out.PPlist = PPlist;
    end
    
%% Prediction of test inputs (filtering and smoothing)

  % Check that we are predicting
  if ~isempty(xt)
    
    % *** Prepare smoother (steady-state) ***
%     PGlist = {};
    PGlist = cell(N+D,1);
    
    % Set up each separate model
    for n=1:(N+D)
    
    ii = ilist(n):ilist(n+1)-1;
    
    % Set up interpolation setup
    for j=1:numel(ro)
        PP = reshape(PPlisto{n}(j,:),size(Q(ii,ii)));
        S  = H(n,ii)*PP*H(n,ii)'+ro(j);
        K = PP*H(n,ii)'/S;
        P = PP-K*ro(j)*K';
        
        % Calculate cholesky
        [L,notpositivedefinite] = chol(A(ii,ii)*P*A(ii,ii)'+Q(ii,ii),'lower');
        if notpositivedefinite>0
            [V,DD]=eig(A(ii,ii)*P*Aii,ii'+Qii,ii); 
            ind=diag(DD)>0; APAQ = V(:,ind)*DD(ind,ind)*V(:,ind)';
            L = cholcov(APAQ)';
        end
        
        % Solve the associated DARE (with some stabilization tricks)
        G = P*A(ii,ii)'/L'/L;
        QQ = P-G*(PP)*G'; QQ = (QQ+QQ')/2;
        [V,DD]=eig(QQ); ind=diag(DD)>0; QQ = V(:,ind)*DD(ind,ind)*V(:,ind)';
        try
            PS2 = dare(G',0*G,QQ);
        catch
            PS2 = P*0;
            ro(j) = nan;
            fprintf('Failed smoother DARE %i\n',j)
        end
        PGlist{n}(j,:) = [PS2(:)' G(:)']; %%#ok
    end
    PGlist{n}(isnan(ro),:)=[]; ro(isnan(ro)) = [];
    
    % Interpolate
    U = apxGrid('interp',{ro},r,3);
    PGlist{n} = U*PGlist{n};
    
    end
      
    %for itt=1:ep_itts
    
    % *** Set-up filter ***
    
    % Set initial state
    m = zeros(size(F,1),1);
    P = Pinf;
    
    % Allocate space for results
    MS = zeros(size(m,1),size(yall,1));
%     PS = zeros(size(m,1),size(m,1),size(yall,1));
    ttau = zeros(D+N,size(yall,1));
    tnu = zeros(D+N,size(yall,1));
        
%     lZ = zeros(1,size(yall,1));
%     R = zeros(D+N,size(yall,1));
    R = exp(lik_param) .* ones(D+N,size(yall,1));

    ys = nan(D+N,size(yall,1));
    
    % Time step, dt
%     dt = 1;
    
    % ### Forward filter

    % Iterate EP
    for itt = 1:ep_itts
        
        ep_damp = ep_damping(itt);
        
        lZ = 0;
        
        % Track convergence (not very optimal)
        maxDiffM = 0;
        PSP = P;
        MSP = MS;
    
        % The filter recursion
        for k=1:numel(yall)
          if mod(k,16000) == 0, fprintf('filtering step %i / %i \n',k,numel(yall)); end
            % Look-up
            if k>1
              PP = zeros(size(P));
              for n=1:(N+D)
                [~,ind] = min(abs(r-R(n,k-1)));
                ii = ilist(n):ilist(n+1)-1;
                PP(ii,ii) = reshape(PPlist{n}(ind,:),size(A(ii,ii)));
              end
              %[~,ind] = min(abs(r-R(k-1)));
              %PP = reshape(PPlist(ind,:),size(A));
            else
              PP = Pinf;
            end

            % Latent marginal, cavity distribution
            fmu = H*A*m; W = PP*H'; HPH = diag(H*W);
            
            % do normal update on first pass
            if itt == 1 || k == numel(yall)
                  
                % Propagate moments through likelihood
                [lZ_k,dlZ,d2lZ] = mom(lik_param,fmu,HPH,Wnmf,1,yall,k);
                lZ = lZ + lZ_k;               
      
                % Perform moment matching
                ttau(:,k) = (1-ep_damp)*ttau(:, k) + ep_damp*(-d2lZ'./(1+d2lZ'.*HPH));
                tnu(:,k) = (1-ep_damp)*tnu(:, k) + ep_damp*((dlZ'-fmu.*d2lZ')./(1+d2lZ'.*HPH));

                % This is the equivalent measurement noise
                R(:,k) = 1 ./ ttau(:,k); %-(1+d2lZ'.*HPH)./d2lZ';
                
            end

            % Enforce positivity->lower bound ttau by zero
            ttau(:,k) = max(ttau(:,k),0);

            % Equivalent data
            ys(:,k) = tnu(:,k)./ttau(:,k);

            % Each separate update (we can get rid of this loop!)
            for n=1:(N+D)

                ii = ilist(n):ilist(n+1)-1;

                % Deal with special cases
                if ttau(n,k)==0
                    %warning('Moment matching hit bound!')
                    R(n,k) = inf;
                    m(ii) = A(ii,ii)*m(ii);
                    P(ii,ii) = PP(ii,ii);
                else

                    % Gain
                    K = W(ii,n)/(HPH(n)+R(n,k));

                    % Precalculate
                    AKHA = A(ii,ii)-K*H(n,ii)*A(ii,ii);

                    % The stationary filter recursion
                    m(ii) = AKHA*m(ii) + K*ys(n,k);             % O(m^2)
                    P(ii,ii) = PP(ii,ii)-K*R(n,k)*K';

                end        

            end

            % Store estimate
            MS(:,k)   = m;
    %         PS(:,:,k) = P;

        end


        % Output debugging info
        if nargout>5
          out.tnu = tnu;
          out.ttau = ttau;
          %out.z = z;
          out.lZ = lZ;
          out.R = R;
          out.MF = MS;
    %       out.PF = PS;
        end

        % ### Backward smoother
        %{
        % Final state
        try

            % Final state covariance 
            P = PS(:,:,end);

            % The gain and covariance
            [L,notpositivedefinite] = chol(A*P*A'+Q,'lower');
            if notpositivedefinite>0
                [V,DD]=eig(A*P*A'+Q); ind=diag(DD)>0; APAQ = V(:,ind)*DD(ind,ind)*V(:,ind)';
                L = cholcov(APAQ)';
            end
            G = P*A'/L'/L;

            QQ = P-G*(A*P*A'+Q)*G'; QQ = (QQ+QQ')/2;
            [V,DD]=eig(QQ); ind=diag(DD)>0; QQ = V(:,ind)*DD(ind,ind)*V(:,ind)';
            PS2 = dare(G',0*G,QQ);
            PS(:,:,end) = PS2;

        catch

            % Deal with this separately

            % Look-up
            [~,ind]=min(abs(r'-R(:,k)),[],2);
            ind(isinf(R(:,k))) = numel(r);

            % Assign values
            for n=1:(N+D)
              ii = ilist(n):ilist(n+1)-1;
              PG = PGlist{n}(ind(n),:);
              PS(ii,ii,end) = reshape(PG(1:end/2),size(P(ii,ii)));
            end
        end
        %}
        P = zeros(size(A));
        G = zeros(size(A));

        % Rauch-Tung-Striebel smoother
        for k=size(MS,2)-1:-1:1
            if mod(size(MS,2)-k,16000) == 0
                fprintf('smoothing step %i / %i \n',size(MS,2)-k,numel(yall));
            end
            
            % Look-up
            [~,ind]=min(abs(r'-R(:,k)),[],2);
            ind(isinf(R(:,k))) = numel(r);

            % Assign values
            for n=1:(N+D)
              ii = ilist(n):ilist(n+1)-1;
              PG = PGlist{n}(ind(n),:);
              P(ii,ii) = reshape(PG(1:end/2),size(P(ii,ii)));
              G(ii,ii) = reshape(PG(end/2+1:end),size(P(ii,ii)));
            end

            % Backward iteration
            m = MS(:,k) + G*(m-A*MS(:,k));       % O(m^2)

            % Store estimate
            MS(:,k)   = m;
    %         PS(:,:,k) = P;

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
                  
                  
                  % TODO: Do something more clever. Right now we are
                  % computing moments for all sites, even though we might
                  % only update a subset.
                
                  % Compute gradients of the normalizer of the tilted dst
                  [lZ_k,dlZ_,d2lZ_] = mom(lik_param,m_cav,v_cav,Wnmf,ep_fraction,yall,k);
                  if itt > 1, lZ = lZ + lZ_k; end
                  if ~isreal(dlZ_)
                      warning('moment matching resulted in complex values')
                      dlZ_ = real(dlZ_);
                      d2lZ_ = real(d2lZ_);
                  end

                  % Moment matching
                  ttau(update_idx, k) = (1-ep_damp)*ttau(update_idx, k) + ep_damp/ep_fraction*(-d2lZ_(update_idx)'./(1+d2lZ_(update_idx)'.*v_cav(update_idx)));
                  tnu(update_idx,k) = (1-ep_damp)*tnu(update_idx, k) + ep_damp/ep_fraction*((dlZ_(update_idx)'-m_cav(update_idx).*d2lZ_(update_idx)')./(1+d2lZ_(update_idx)'.*v_cav(update_idx)));

                  % This is the equivalent measurement noise
                  R(update_idx,k) = 1 ./ ttau(update_idx,k); %-(1+d2lZ_'.*v_cav)./d2lZ_';
                  % Max diff in m and P
                  
                  maxDiffM = max(maxDiffM,max(max(abs(H*MSP(:,k)-H*m))));

                end
            end
        
        end  % end smoother
        
        maxDiffP = max(max(abs(H*PSP*H'-H*P*H')));
        
        if itt < ep_itts
          fprintf('%.02i - max diff in m: %.6g - max diff in P: %.6g - nll: %.6g\n', ...
                  itt,maxDiffM,maxDiffP,-lZ)
        end
    
        
    end  % end EP iteration
    
    
    
    % Debugging information
    if nargout>5
        out.MS = MS;
%         out.PS = PS;
    end
      
    % Estimate the joint covariance matrix if requested
    if nargout > 2
      
      % Allocate space for results
      Covft = []; %zeros(size(PS,3));
          
      % Lower triangular (this output is not used)
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
  
    % These indices shall remain to be returned
    MS = MS(:,return_ind);
%     PS = PS(:,:,return_ind);
    
    % Return mean
    Eft = H*MS;
    
    % Return variance
    if nargout > 1
        Varft = repmat(diag(H*P*H'),[1,size(MS,2)]);
        varargout = {Eft,Varft};
    else
        varargout = {Eft};
    end
    
    
 
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
        if nargout > 5
            varargout = {Eft,Varft,Covft,lb,ub,out};
        else
            varargout = {Eft,Varft,Covft,lb,ub};
        end
        
    end

  end
  
  
%% Evaluate negative log marginal likelihood and its gradient

  if isempty(xt)
      
    ep_damp = ep_damping(1);
  
    % Size of inputs
    d = size(F,1);
    nparam = numel(w);
            
    % Allocate space for results
    gdata = zeros(1,nparam);
    
    % Set up
    m  = zeros(d,1);    
    ttau = zeros(D+N,size(yall,1));
    tnu = zeros(D+N,size(yall,1));
    lZ = zeros(1,size(yall,1));
    
    % Loop over all observations
    for k=1:numel(yall)
        
        % Look-up
        if k>1
          PP = zeros(size(P));
          for n=1:(N+D)
            [~,ind] = min(abs(r-R(n,k-1)));
            ii = ilist(n):ilist(n+1)-1;
            PP(ii,ii) = reshape(PPlist{n}(ind,:),size(A(ii,ii)));
          end
          %[~,ind] = min(abs(r-R(k-1)));
          %PP = reshape(PPlist(ind,:),size(A));
        else
          PP = Pinf;
        end
        
        % Latent marginal, cavity distribution
        fmu = H*A*m; W = PP*H'; HPH = diag(H*W);
        
        % Propagate moments through likelihood
        [lZ(k),dlZ,d2lZ] = mom(lik_param,fmu,HPH,Wnmf,1,yall,k);
        
        % Perform moment matching
        ttau(:,k) = (1-ep_damp)*ttau(:,k) + ep_damp*(-d2lZ'./(1+d2lZ'.*HPH));
        tnu(:,k) = (1-ep_damp)*tnu(:,k) + ep_damp*((dlZ'-fmu.*d2lZ')./(1+d2lZ'.*HPH));
        
        % This is the equivalent measurement noise
        R(:,k) = 1 ./ ttau(:,k); %-(1+d2lZ'.*HPH)./d2lZ';

        % Enforce positivity->lower bound ttau by zero
        ttau(:,k) = max(ttau(:,k),0);
        
        % Equivalent data
        ys(:,k) = tnu(:,k)./ttau(:,k);
        
        % Each separate update (we can ghet rid of this loop!)
        for n=1:(N+D)
            
            ii = ilist(n):ilist(n+1)-1;
            
            % Deal with special cases
            if ttau(n,k)==0
                %warning('Moment matching hit bound!')
                R(n,k) = inf;
                m(ii) = A(ii,ii)*m(ii);
                P(ii,ii) = PP(ii,ii);
            else
                
                % Gain
                K = W(ii,n)/(HPH(n)+R(n,k));
                
                % Precalculate
                AKHA = A(ii,ii)-K*H(n,ii)*A(ii,ii);
                
                % The stationary filter recursion
                m(ii) = AKHA*m(ii) + K*ys(n,k);             % O(m^2)
                P(ii,ii) = PP(ii,ii)-K*R(n,k)*K';
                
            end        
        
        end        
        
    end

    % Sum up
    edata = -sum(lZ);
    
    % Account for log-scale
    gdata = gdata.*exp(w);
    
    % Return negative log marginal likelihood and gradient
    varargout = {edata,gdata};

  end
  
  
  
end
