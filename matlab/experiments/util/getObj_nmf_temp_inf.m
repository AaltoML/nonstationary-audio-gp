function [obj,dobj] = getObj_nmf_temp_inf(logH,A,W,vary,varargin)

% Performs inference in the NMF model, but not learning. Used for
% missing data or denoising experiments.
%
% Vanilla NMF without temporal priors:
% function [obj,dobj] = getObj_nmf_temp_inf(logH,A,W,vary)
%
% NMF with temporal priors
% function [obj,dobj] = getObj_nmf_temp_inf(logH,A,W,vary,muinf,varinf,lam)
%
% NMF prior with (optional) temporal constraints (K = number of
% basis functions can be under/over/complete).
%
% A \approx Ahat = H*W
%
% The weights are normalised to have sum 1
%
% The cost function used in the likelihood is:
% \sum_{d,t} A_{d,t}./Ahat_{t,d} + \sum_{d,t} \log(Ahat_{t,d}) 
%
% Uses Inverse Gamma Markov chain priors on the temporal basis functions:
% h_{k,1} ~ IG(alp1,bet1)
% h_{k,t}|h_{k,t-1} ~ IG(alp_{k},bet_{k,t})
% alp_k = 2+lam_k^2+muinf_k^2/varinf_k
% bet_{t,k} = (alp_{t,k}-1)(lam_k*(h_{k,t-1}-muinf_k)+muinf_k)
%
% The steady state mean of the IG chain is muinf(k), the variance
% varinf(k), and the correlation coefficient between successive
% variables is lam(k).
%
% The first variable in the chain is drawn from an inverse Gamma
% distribution with mean and variance equal to the steady state.
%
% INPUTS
% logH = vectorised temporal basis matrix H a[T*K,1]
% A = data [T,D]
% W = spectral basis [K,D]
% vary = observation noise [T,D]
% muinf = steady state mean of the IGMC priors [1,K]
% varinf = steady state variance of the IGMC priors [1,K]
% lam = steady state correlation-coefficient for successive
%       variables in the IGMC priors [1,K]


[T,D] = size(A);
[K,D] = size(W);
logH = reshape(logH,[T,K]);

H = exp(logH);

Ahat = H*W+vary;

%%%

obj1 = sum(A(:)./Ahat(:)+log(Ahat(:)));

dobj1dA = (Ahat-A)./Ahat.^2;

dobj1dH = dobj1dA*W';

dobj1dlogh = dobj1dH(:).*H(:);

%%%
if nargin>6

  muinf = varargin{1};
  varinf = varargin{2};
  lam = varargin{3};
  
alp1 = muinf.^2./varinf+2;
bet1 = (alp1-1).*muinf;

alp = (2 + lam(:).^2 + muinf(:).^2./varinf(:))';

%a = alp+1;
%b = [0,0]

a = (alp-1).*lam;
b = (alp-1).*(1-lam).*muinf;
ahtpb = (ones(T-1,1)*a).*H(1:T-1,:)+ones(T-1,1)*b;

obj2 = sum(bet1./H(1,:))+sum((alp1+1).*logH(1,:)) ...
       - alp*sum(log(ahtpb),1)' ...
       + (alp+1)*sum(logH(2:T,:),1)' ...
       + sum(sum(ahtpb./H(2:T,:),1),2);


dobj2dH = zeros(T,K);

dobj2dH(1,:) = a./H(2,:)-bet1./H(1,:).^2+(alp1+1)./H(1,:)-alp.*a./ahtpb(1,:);

dobj2dH(2:T-1,:) = (ones(T-2,1)*a)./H(3:T,:)-ahtpb(1:T-2,:)./H(2:T-1,:).^2 ...
                   + (ones(T-2,1)*(alp+1))./H(2:T-1,:) ...
                   -(ones(T-2,1)*(alp.*a))./ahtpb(2:T-1,:);

dobj2dH(T,:) = -ahtpb(T-1,:)./H(T,:).^2 + (1+alp)./H(T,:);

dobj2dlogh = dobj2dH(:).*H(:);

obj = (obj1+obj2)/T;
dobj = [dobj1dlogh+dobj2dlogh]/T;

elseif nargin>4
    
  varinf = varargin{1};
  lam = varargin{2};
  
  obj2 = ((1+lam.^2)./(2*varinf))*sum(H(2:T-1,:).^2)' ...
         + (1./(2*varinf))*(H(T,:).^2)' +  (lam.^2./(2*varinf))*(H(1,:).^2)' ...
	 -(lam./varinf)*sum(H(2:T,:).*H(1:T-1,:),1)';  
  
  dobj2dH = zeros(T,K);
  
  dobj2dH(1,:) = (lam.^2./varinf).*H(1,:) - (lam./varinf).*H(2,:);
  
  dobj2dH(2:T-1,:) = ones(T-2,1)*((1+lam.^2)./varinf).*H(2:T-1,:) ...
                     - ones(T-2,1)*(lam./varinf).*(H(3:T,:)+H(1:T-2,:));
  
  dobj2dH(T,:) = 1./varinf.*H(T,:) -  (lam./varinf).*H(T-1,:);

  dobj2dlogh = dobj2dH(:).*H(:);
  

  obj = (obj1+obj2)/T;
  dobj = [dobj1dlogh+dobj2dlogh]/T;


else
  obj = (obj1)/T;
  dobj = [dobj1dlogh]/T;
end

%obj = obj1/T;
%dobj = [dobj1dlogh;dobj1dlogw]/T;

% obj = obj2;
% dobj = [dobj2dlogh;zeros(K*D,1)];




