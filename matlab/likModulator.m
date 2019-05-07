function [varargout] = likModulator(link, hyp, y, mu, s2, p, inf, i)

% INPUTS:
% link: link function, e.g. softplus
% hyp: observation noise (log-transformed)
% mu: mean
% s2: variance
% p: order of cubature for Gaussian integral - cubature exact for pth-order polynomials
% inf: inference method - only EP implemented

if nargin<4, varargout = {'1'}; return; end   % report number of hyperparameters

sn2 = exp(hyp);

if nargin<7                              % prediction mode if inf is not present
    error('prediction mode not implemented for modulator likelihood')
else
  switch inf 
  case 'infLaplace'
      error('infLaplace not implemented for modulator likelihood')

  case 'infEP'   
    if nargin<8                                             % no derivative mode
      
      jitter = 1e-8;
      D = length(mu)/2;
      mu_z = mu(1:D); mu_g = mu(D+1:end); % z = subband, g = modulator
      s2_z = s2(1:D); s2_g = s2(D+1:end);
      if ismember(p,[3,5,7,9])
        [wn,xn_unscaled,~] = utp_ws(p,D); % Symmetric, pth-order cubature
        xn = (mu_g + diag(s2_g.^0.5)*xn_unscaled)';
        wn = wn';
      else
%         [wn,xn] = gh_ws(mu_g,diag(s2_g),p); % Gauss-Hermite quadrature
%         xn = xn';
%         wn = wn';
        [xn, wn] = mvhermgauss(mu_g,s2_g,p); % Gauss-Hermite quadrature - Will's implementation
      end
      
      normy = @(x) normpdf(y, ...
                           link(x)*mu_z, ...
                           sqrt(sn2+(link(x).^2)*s2_z));
      normy_xn = normy(xn);
      link_xn = link(xn);
      sn2_link_xn2_s2_z = sn2 + link_xn.^2 * s2_z;
      link_xn_mu_z = link_xn * mu_z;
      xn_mu_g_s2_g = (xn - mu_g') ./ s2_g';

%       Z_integrand = @(x) normy(x);
%       Z = sum(wn .* Z_integrand(xn));
      Z = max(sum(wn .* normy_xn), jitter);
      Zinv = 1./Z;
      lZ = log(Z);
      if ~real(lZ)
          keyboard
      end
      if nargout>1
        dZ_integrand1 = link_xn .* (y-link_xn_mu_z) ...
                             ./ sn2_link_xn2_s2_z .* normy_xn;
%         dZ_integrand1 = @(x) link(x) .* (y-link(x)*mu_z) ...
%                              ./ (sn2+(link(x).^2)*s2_z) .* normy(x);
                       
%         dZ1 = sum(wn .* dZ_integrand1);
%         dZ1 = sum(wn .* dZ_integrand1(xn));
%         dlZ(1:D) = Zinv .* dZ1;
        dlZ_z = Zinv .* sum(wn .* dZ_integrand1);
        
%         dZ_integrand2 = @(x) (x - mu_g') ./ s2_g' .* normy(x);
        dZ_integrand2 = xn_mu_g_s2_g .* normy_xn;
        
%         dZ2 = sum(wn .* dZ_integrand2);
%         dlZ(D+1:2*D) = Zinv .* dZ2;
        dlZ_g = Zinv .* sum(wn .* dZ_integrand2);
        
        if nargout>2
            d2Z_integrand1 = link_xn.^2 .* ( ...
                                   ( (y-link_xn_mu_z) ./ sn2_link_xn2_s2_z ).^2 ...
                                   - (1 ./ sn2_link_xn2_s2_z) ) ...
                                  .* normy_xn;
%             d2Z_integrand1 = @(x) link(x).^2 .* ( ...
%                                    ( (y-link(x)*mu_z) ./ (sn2+(link(x).^2)*s2_z) ).^2 ...
%                                    - (1 ./ (sn2+(link(x).^2)*s2_z)) ) ...
%                                   .* normy(x);
            d2lZ_z = -dlZ_z.^2 + Zinv .* sum(wn .* d2Z_integrand1);
%             d2lZ(1:D) = -dlZ(1:D).^2 + Zinv .* sum(wn .* d2Z_integrand1(xn));
            
            d2Z_integrand2 = ( xn_mu_g_s2_g.^2 - (1 ./ s2_g') ) .* normy_xn;
%             d2Z_integrand2 = @(x) ( ((x - mu_g') ./ s2_g').^2 ...
%                                     - (1 ./ s2_g') ) ...
%                                   .* normy(x);
            d2lZ_g = -dlZ_g.^2 + Zinv .* sum(wn .* d2Z_integrand2);
%             d2lZ(D+1:2*D) = -dlZ(D+1:2*D).^2 + Zinv .* sum(wn .* d2Z_integrand2(xn));
        else
            d2lZ = {};
        end
      else
          dlZ = {};
      end
      varargout = {lZ,[dlZ_z,dlZ_g],[d2lZ_z,d2lZ_g]};
%       keyboard
    else                                                       % derivative mode
        error('derivative mode not implemented for modulator likelihood')
%       dlZhyp = ((y-mu).^2./(sn2+s2)-1) ./ (1+s2./sn2);   % deriv. w.r.t. hyp.lik
%       varargout = {dlZhyp};
    end

  case 'infVB'
      error('infVB not implemented for modulator likelihood')
  end

  
end