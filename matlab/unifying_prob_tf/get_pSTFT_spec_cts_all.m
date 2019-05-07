function spec = get_pSTFT_spec_cts_all(freqs,lamx,varx,om,kernel)

% Will Wilkinson 20181026

% function spec = get_pSTFT_spec(freqs,lamx,varx,om)
%
% Computes the spectrum of a probabilistic STFT model
% x_{1,t} = lam x_{1,t-1} +\eta_{1,t} \varx^{1/2}
% x_{2,t} = lam x_{2,t-1} +\eta_{2,t} \varx^{1/2}
% \eta_{i,t} ~ \Norm(0,1)
% y_t = Real(exp(i om)*(x_{1,t}+i x_{2,t}))
%
% see getCompSpecPFB.m for the complex spectrum
%
% INPUTS
% freqs = frequencies at which to evaluate the spectrum, [1,N]
% lamx = dynamical AR parameters [D,1]
% varx = dynamical noise parameters [D,1]
% om = mean frequencies of the sinusoids [D,1]
%
% OUTPUTS
% spec = spectrum of the probabilistic STFT model
% 

D = length(lamx);
N = length(freqs);

spec = zeros(D,N);
omegas = 2*pi*freqs(:)';

% for d=1:D

%   spec(d,:) = 1/2*varx(d)./(1+lamx(d).^2-2*lamx(d)*cos(omegas-om(d))) ...
%             + 1/2*varx(d)./(1+lamx(d).^2-2*lamx(d)*cos(omegas+om(d))); 
%   spec(d,:) = varx(d) .* lamx(d) .* (lamx(d).^2 + (omegas-om(d)).^2).^-1 ...
%             + varx(d) .* lamx(d) .* (lamx(d).^2 + (omegas+om(d)).^2).^-1;
        
% end

ss_func = str2func(strcat('cf_',kernel,'_to_ss'));
if strcmp(kernel,'exp')
    len = 1./lamx;
elseif strcmp(kernel,'matern32')
    len = sqrt(3)./lamx;
% elseif strcmp(kernel,'matern52')
else
    if ~strcmp(kernel,'matern52')
        disp('Warning - hyperparameter mapping may be incorrect')
    end
    len = sqrt(5)./lamx;
end


spec_om = zeros(1,N); I2 = eye(2);
for d=1:D
    % Construct the continuous-time SDE: Matern (from Solin+Sarkka AISTATS 2014)
    [F1,L1,Qc1,H1] = ss_func(varx(d), len(d), 4);
    
    tau1 = length(L1);
    Itau1 = eye(tau1); Itau2 = eye(2*tau1);

    % Construct the continuous-time SDE: periodic (just a sinusoid)
    F2 = [0 -om(d); om(d) 0];    

    % The product of the two (from Solin+Sarkka AISTATS 2014)
    F = kron(F1,I2) + kron(Itau1,F2);
    L = kron(L1,I2);
    Qc = kron(Qc1,I2);
    H = kron(H1,[1 0]);

    for i=1:N
        % can get small imaginery component due to numerical error
        spec_om(i) = real( H*((F - 1i*omegas(i)*Itau2)\L)*Qc*(H*((F - 1i*omegas(i)*Itau2)\L))' );
    end
    spec(d,:) = spec_om;
    
end


