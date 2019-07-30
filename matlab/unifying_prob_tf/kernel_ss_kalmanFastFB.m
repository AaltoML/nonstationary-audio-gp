function  [lik,Xfin,Pfin,varargout] = kernel_ss_kalmanFastFB(A,Q,C,P0,K,vary,y,varargin)
% 
%                      
%                     
%% Some of the tweaks from the original function
                      
T = length(y);
Y = reshape(y,[1,1,T]);
lik=0;

if nargin<=7
  verbose = 0 ;
else
  verbose = varargin{1};
end

if nargin<=8
  KF = 0 ;
else
  KF = varargin{2};
end

if verbose==1
  if KF==1
    disp('Infnite-horizon Kalman Filtering')
  else
    disp('Infnite-horizon Kalman Smoothing')
  end
end

%% ARNO's re-write starts here

    tic

    % T / 2*number of progress values displayed
    CntInt=T/5; 

    % Set initial state and assign model matrices
    m = zeros(size(A,1),1);
%     P = P0;
    H = C;
    R = vary;
    
    
    % Check if this agrees with the predictive Riccati solution
    try
        
        % The gain matrix (inc. prediction)
        [PP,~,~,report] = dare(A',H',Q,R);
        
        % Innovation Cholesky (inc. prediction)
        S  = H*PP*H'+R;
        
    catch
        error('Unstable DARE solution!')
    end
    
    % Stationary gain
    K = PP*H'/S;
    
    % Precalculate
    AKHA = A-K*H*A;
    
    % Allocate space for results
    MS = zeros(size(m,1),size(Y,3));
    PS = zeros(size(m,1),size(m,1),size(Y,3));
    
    % ### Forward filter
    if verbose==1
        disp('filtering stage:')
    end
    
    % Stationary filter covariance
    PF2 = PP - K*H*PP;
    HA = H*A;
    
    % Constant part of lik
    lik = .5*log(2*pi)*size(Y,3) + .5*log(S)*size(Y,3);

    % The filter recursion
    for k=1:size(Y,3)
        
        % Progress
        if verbose==1 && mod(k-1,CntInt)==0
            fprintf(['Progress ',num2str(floor(50*k/size(Y,3))),'%%','\r'])
        end
        
        if ~isnan(Y(:,:,k))
        
            % Innovation mean
            v = Y(:,:,k) - HA*m;
        
            % The stationary filter recursion
            m = AKHA*m + K*Y(:,:,k);             % O(m^2)
        
            % Evaluate the energy (neg. log lik): Check this
            lik = lik + .5*v^2/S;
        
        else
            
            % Do prediction for unseen points
            m = A*m;
        
        end
        
        % Store estimate
        MS(:,k)   = m;
        PS(:,:,k) = PF2; % Waste of memory to store this for every step!
        
    end
%     PS = PF2;
    if verbose==1
        toc
    end
       
    % ### Backward smoother

    % Should we run the smoother?
    if KF==1
        
        % Only return filter  
        
    else
        if verbose==1
            disp('smoothing stage:')
        end
      
        % Solve backward smoother gain
        G = PF2*A'/PP;
    
        % Solve smoother state covariance
        QQ = PF2-G*PP*G';
        QQ = (QQ+QQ')/2;
        P = dare(G',zeros(size(QQ)),QQ);
        
        % Steady-state Rauch-Tung-Striebel smoother
        for k=size(MS,2)-1:-1:1
            
            % Progress
            if verbose==1 && mod(k+1,CntInt)==0
                fprintf(['Progress ',num2str(50+floor(50*(size(MS,2)-k+1)/size(MS,2))),'%%','\r'])
            end
            
            % Do smoother update
            m = MS(:,k) + G*(m-A*MS(:,k));
            
            % Store estimate
            MS(:,k)   = m;
            PS(:,:,k) = P;
            
        end
%         PS = P;
        if verbose==1
            toc
        end
    end
    
    % Return variables of interest
    lik = -lik;
    Xfin = reshape(MS,[1 size(MS)]);
    Pfin = PS;

%%
    
   if verbose==1
     fprintf('                                        \r')
   end    
    
