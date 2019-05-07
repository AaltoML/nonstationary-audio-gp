function [Obj,varargout] = get_Obj_pSTFT_all(theta,vary,specTar,minVar,limOm,limLam,bet,kernel);

% Will Wilkinson 20181026
% calculates spectral density and its gradients for the spectral mixture GP
% with any kernel function. Based on Richard Turner's original
% code.


% function Obj =
% get_Obj_pSTFT_spec(theta,vary,specTar,minVar,limOm,limLam,bet);
%
% Optional:
%
% function [Obj,dObj] =
% get_Obj_pSTFT_spec(theta,vary,specTar,minVar,limOm,limLam,bet);
%
% For fitting a probabilistic spectrogram model to a signal.  
%
% x_{1,t,d} = lam_d x_{1,t-1,d} +\eta_{1,t,d} \varx_d^{1/2}
% x_{2,t,d} = lam_d x_{2,t-1,d} +\eta_{2,t,d} \varx_d^{1/2}
% \eta_{i,t,d} ~ \Norm(0,1)
% y_t = real(\sum_{d} exp(i om_d)*(x_{1,t,d}+i x_{2,t,d}))

% Fits the parameters (lam_d, varx_d,om_d) such that the spectra are
% matched to the target spectrum specified in specTar.
%
% The cost function is: \sum_t [ log(spec_t) + specTar_t/spec_t ]
% Which can be derived from the fact that the model is a Gaussian
% with a spectrum parameterised using the AR(1)/cosine parameters.
%
% The function parameterises the AR processes using the marginal
% variance of the processes rather than the conditional (marginal
% variance = conditional variance/(1-lam^2)) since the objective is
% more simply behaved in this space and the constraints on the
% parameters are simpler to enforce.
%
% INPUTS
% theta = parameters of the probabilistic spectrogram to be optimised:
%         first D components are the log marginal variances, the next
%         D components are the transformed sinusiod frequencies, the
%         final D parameters the transformed AR lambda
%         parameters. Size is therefore [3*D,1]
% vary = white noise variance (if this is set to zero it can result
%        in numerical problems arising from the division
%        spectAR2./specTar) a sensible setting is: vary=max(specTar)*1e-4
% specTar = target spectrum of the data to be fit, size [N,1]
% minVar = minimum of the marginal variances size [D,1]
% limOm = min/max centre-frequencies e.g. [0,1/2], size [D,2]
% limLam = min/max bandwidths e.g. [0,0.999], size [D,2]
% bet = strength of the gamma shrinkage prior on the marginal
%       variance parameter (set to zero for no pruning)
% 
% OUTPUTS
% Obj = objective
% OPTIONAL OUTPUT:
% dObj = derivative of the objective wrt the parameters, size [3*D,1]
% 
% See Chapter 5 of my thesis (Statistical Models for Natural Sounds
% by R.E.Turner) for more details about AR(2) Filter banks.

D = length(theta)/3;

dVar= exp(theta(1:D));
mVar = minVar+dVar;

om = limOm(:,1)+(limOm(:,2)-limOm(:,1))./(1+exp(-theta(D+1:2*D)));
lam = limLam(:,1)+(limLam(:,2)-limLam(:,1))./(1+exp(-theta(2*D+1:3*D)));

% cts model parameters:


N = length(specTar);
omegas = linspace(0,pi,ceil(N/2));
omegas = [omegas,-omegas([floor(N/2):-1:1])];

% conditional variance
% cVar = mVar .* (1 - lam.^2);


ss_func = str2func(strcat('cf_',kernel,'_to_ss'));  
if strcmp(kernel,'exp')
    len = 1./lam;
    dl_dlam = -1.*lam.^-2;
elseif strcmp(kernel,'matern32')
    len = sqrt(3)./lam;
    dl_dlam = -sqrt(3).*lam.^-2;
% elseif strcmp(kernel,'matern52')
else
%     if ~strcmp(kernel,'matern52')
%         disp('Warning - hyperparameter mapping may be incorrect')
%     end
    len = sqrt(5)./lam;
    dl_dlam = -sqrt(5).*lam.^-2;
end
num_hypers = 3;
   


% Get the component spectra
spec = ones(1,N)*vary; spec_om = zeros(1,N); I2 = eye(2);
  
for d=1:D
    % Construct the continuous-time SDE: Matern (from Solin+Sarkka AISTATS 2014)
    [F1,L1,Qc1,H1] = ss_func(mVar(d), len(d), 4);
    
    tau1 = length(L1);
    Itau1 = eye(tau1); Itau2 = eye(2*tau1);

    % Construct the continuous-time SDE: periodic (just a sinusoid)
    F2 = [0 -om(d); om(d) 0];    

    % The product of the two (from Solin+Sarkka AISTATS 2014)
    F = kron(F1,I2) + kron(Itau1,F2);
    L = kron(L1,I2);
    Qc = kron(Qc1,I2);
    H = kron(H1,[1 0]);

  % for the objective
%   alp1_c = lam(d).^2 + (omegas-om(d)).^2;
%   alp2_c = lam(d).^2 + (omegas+om(d)).^2;
%   spec = spec + cVar(d) .* lam(d) .* (alp1_c.^-1 + alp2_c.^-1);

    for i=1:N
        % can get small imaginery component due to numerical error
        spec_om(i) = real( H*((F - 1i*omegas(i)*Itau2)\L)*Qc*(H*((F - 1i*omegas(i)*Itau2)\L))' );
    end
    spec = spec + (1 - lam(d).^2) .* spec_om;
end




if nargout>1
  % Derivative of the components wrt parameters

  dObjdtransVar = zeros(D,1);
  dObjdtransOm = zeros(D,1);
  dObjdtransLam = zeros(D,1);
  
  dObjdspec = 1./spec-specTar'./(spec.^2);
  
  for d=1:D
      
    % Construct the continuous-time SDE: Matern (from Solin+Sarkka AISTATS 2014)
    [F1,L1,Qc1,H1,~,dF1_,dQc1_,~] = ss_func(mVar(d), len(d), 4);
    
    tau1 = length(L1);
    Itau1 = eye(tau1); Itau2 = eye(2*tau1);

    % Construct the continuous-time SDE: periodic (just a sinusoid)
    F2 = [0 -om(d); om(d) 0];    

    % The product of the two (from Solin+Sarkka AISTATS 2014)
    F = kron(F1,I2) + kron(Itau1,F2);
    L = kron(L1,I2);
    Qc = kron(Qc1,I2);
    H = kron(H1,[1 0]);

    
      % Derivative of the components wrt parameters
        dF1 = zeros([size(F1), num_hypers]);
        dF1(:, :, 1:2) = dF1_;
        dF1(:, :, 3) = zeros(size(F1));
        dQc1 = zeros([size(Qc1), num_hypers]);
        dQc1(:, :, 1:2) = dQc1_;
        dQc1(:, :, 3) = zeros(size(Qc1));


        dF2 = zeros([size(F2), num_hypers]);
        dF2(:, :, 1:2) = zeros([size(F2), 2]);
        dF2(:, :, 3) = [0 -1; 1 0];


        % gradients
        dF = zeros([size(F), num_hypers]);
        dQc = zeros([size(Qc), num_hypers]);
        for h=1:num_hypers
            dF(:,:,h) = kron(dF1(:,:,h),I2) + kron(Itau1,dF2(:,:,h));
            dQc(:,:,h) = kron(dQc1(:,:,h),I2);
        end

        % Gradient of spectral density
        dF_dvar = dF(:,:,1);
        dF_dl = dF(:,:,2);
        dF_dom = dF(:,:,3);

        dQc_dvar = dQc(:,:,1);
        dQc_dl = dQc(:,:,2);
        dQc_dom = dQc(:,:,3);
        
        B = L*L';
        
        dS_dvar = zeros(1,N); dS_dom = zeros(1,N); dS_dlam = zeros(1,N);
        for i=1:N
            % Evaluate analytical gradient
            G = F - 1i*omegas(i)*Itau2;
            
            J = H/G;
            K_dvar = dF_dvar*(G\B);
            K_dl = dF_dl*(G\B);
            K_dom = dF_dom*(G\B);
            % numerically more robust implementation
            JL = J*L;  

            dS_dvar(i) = dQc_dvar(1)*(JL*JL') - Qc(1)*trace((J'*J)*(K_dvar + K_dvar'));
            dS_dl = dQc_dl(1)*(JL*JL') - Qc(1)*trace((J'*J)*(K_dl + K_dl'));
            dS_dlam(i) = dS_dl*dl_dlam(d);
            dS_dom(i) = dQc_dom(1)*(JL*JL') - Qc(1)*trace((J'*J)*(K_dom + K_dom'));
        end
 
    % Derivative wrt transformed marginal variance
%     alp1_c = lam(d).^2 + (omegas-om(d)).^2;
%     alp2_c = lam(d).^2 + (omegas+om(d)).^2; 
%     dspecdtransVar = (mVar(d)-minVar(d)) .* (1 - lam(d).^2) .* lam(d) .* (alp1_c.^-1 + alp2_c.^-1);
    dspecdtransVar = (mVar(d)-minVar(d)) .* (1 - lam(d).^2) .* dS_dvar;
    dObjdtransVar(d) = sum(dObjdspec.*dspecdtransVar); 
    
    
    % Derivative wrt transformed centre frequency
%     dspecdom = 2 .* mVar(d) .* (1 - lam(d).^2) .* lam(d) * (alp1_c.^-2 .* (omegas-om(d)) - alp2_c.^-2 .* (omegas+om(d)));
    dspecdom = (1 - lam(d).^2) .* dS_dom;
    domdtransOm = (limOm(d,2)-limOm(d,1))*(1/4)./cosh(theta(D+d)/2).^2;
    dObjdtransOm(d) = sum(dObjdspec.*dspecdom)*domdtransOm; 

    
    % Derivative wrt transformed lambda
%     dspecdlam = mVar(d) .* (alp1_c.^-1 + alp2_c.^-1 - 2 .* lam(d).^2 * ...
%                           (alp1_c.^-2 + alp2_c.^-2));
%     dspecdlam = mVar(d) .* ((1 - 3.*lam(d).^2) .* (alp1_c.^-1 + alp2_c.^-1) ...
%                              - 2 .* lam(d).^2 .* (1 - lam(d).^2) * (alp1_c.^-2 + alp2_c.^-2) );                    
    dspecdlam = (1 - lam(d).^2) .* dS_dlam - 2 .* lam(d) .* spec_om;
    dlamdtransLam = (limLam(d,2)-limLam(d,1))*(1/4)./cosh(theta(2*D+d)/2).^2;
    dObjdtransLam(d) = sum(dObjdspec.*dspecdlam)*dlamdtransLam;  
  end
  
  
end


% Likelihood cost function
Obj1 = sum(log(spec))+sum(specTar'./spec);

% Prior cost function
Obj2 = bet*sum(mVar);

% Rescale to get nice units
Obj = (Obj1+Obj2)/N;

if nargout>1
    dObj1 = [dObjdtransVar;dObjdtransOm;dObjdtransLam;];

    dObj2 = [bet*dVar;zeros(2*D,1)];

    dObj = [dObj1+dObj2]/N;
    varargout{1}= dObj;
end

