function [varargout] = likModulatorNMFPowerSq(link, hyp, y, mu, s2, W, p, ep_fraction, inf, i)

% INPUTS:
% link: link function, e.g. softplus
% hyp: observation noise (log-transformed)
% y: observation
% mu: mean
% s2: variance
% W: NMF weights
% p: order of cubature for Gaussian integral - cubature exact for pth-order polynomials
% alpha: Power EP alpha
% inf: inference method - only EP implemented

if nargin<4, varargout = {'1'}; return; end   % report number of hyperparameters

sn2 = exp(hyp);

if nargin<9                              % prediction mode if inf is not present
    error('prediction mode not implemented for modulator likelihood')
else
  switch inf 
  case 'infLaplace'
      error('infLaplace not implemented for modulator likelihood')

  case 'infEP'   
    if nargin<10                                             % no derivative mode
	  
      jitter = 1e-10;
      [D, N] = size(W);
      mu_z = mu(1:D); mu_g = mu(D+1:end); % z = subband, g = modulator
      s2_z = s2(1:D); s2_g = s2(D+1:end);
      if ismember(p,[3,5,7,9])
        [wn,xn_unscaled,~] = utp_ws(p,N); % Symmetric, pth-order cubature
        xn = (mu_g + diag(s2_g.^0.5)*xn_unscaled)';
        wn = wn';
      else
%         [wn,xn] = gh_ws(mu_g,diag(s2_g),p); % Gauss-Hermite quadrature
%         xn = xn';
%         wn = wn';
      [xn, wn] = mvhermgauss(mu_g,s2_g,p); % Gauss-Hermite quadrature - Will's implementation
      end
      
      % pre-calculate these terms:
      link_xn_W = sqrt(link(xn)*W'); % the sqrt here means we are modelling the square of the subband amplitudes, i.e. the spectrogram
      sn2_link_xn2_s2_z = sn2/ep_fraction + link_xn_W.^2 * s2_z;
      link_xn_mu_z = link_xn_W * mu_z;
      xn_mu_g_s2_g = (xn - mu_g') ./ s2_g';
      pEP_const = (2*pi*sn2)^(0.5*(1-ep_fraction)) * ep_fraction^(-0.5);
      
      normy_xn = normpdf(y, ...
                           link_xn_mu_z, ...
                           sqrt(sn2_link_xn2_s2_z));
      
      Z = pEP_const .* max(sum(wn .* normy_xn), jitter);
      Zinv = 1/Z;
      lZ = log(Z);
      if nargout>1
        dZ_integrand1 = link_xn_W .* (y - link_xn_mu_z) ...
                             ./ sn2_link_xn2_s2_z .* normy_xn;
                       
        dZ1 = sum(wn .* dZ_integrand1);
        dlZ(1:D) = Zinv .* pEP_const .* dZ1;
        
        dZ_integrand2 = xn_mu_g_s2_g .* normy_xn;
        
        dZ2 = sum(wn .* dZ_integrand2);
        dlZ(D+1:D+N) = Zinv .* pEP_const .* dZ2;
        
        if nargout>2
            d2Z_integrand1 = link_xn_W.^2 .* ( ...
                                   ( (y - link_xn_mu_z) ./ sn2_link_xn2_s2_z ).^2 ...
                                   - (1 ./ sn2_link_xn2_s2_z) ) ...
                                  .* normy_xn;
            d2lZ(1:D) = -dlZ(1:D).^2 + Zinv .* pEP_const .* sum(wn .* d2Z_integrand1);
            
            d2Z_integrand2 = ( xn_mu_g_s2_g.^2 ...
                                    - (1 ./ s2_g') ) ...
                                  .* normy_xn;
            d2lZ(D+1:D+N) = -dlZ(D+1:D+N).^2 + Zinv .* pEP_const .* sum(wn .* d2Z_integrand2);
        else
            d2lZ = {};
        end
      else
          dlZ = {};
      end
      % TODO: Do something more clever, rather than simply casting the
      % moments as reals before returning
      varargout = {real(lZ),real(dlZ),real(d2lZ)};
    else                                                       % derivative mode
        error('derivative mode not implemented for modulator likelihood')
%       dlZhyp = ((y-mu).^2./(sn2+s2)-1) ./ (1+s2./sn2);   % deriv. w.r.t. hyp.lik
%       varargout = {dlZhyp};
    end

  case 'infVB'
      error('infVB not implemented for modulator likelihood')
  end

  
end