function [xn, wn] = mvhermgauss(mu, s2, N)
  %   Return the evaluation locations, and weights for multivariate
  %   Hermite-Gauss quadrature.
  %   INPUTS
  %   mu: means (Dx1)
  %   s2: variances (Dx1)
  %   N: Number of Gauss-Hermite evaluation points.
  %   OUTPUTS
  %   xn: eval_locations (N**DxD)
  %   wn: weights (N**D)
  if size(mu,1) > 1
      mu = mu'; s2 = s2'; % make row vector
  end
      
    D = length(mu);
    sig = diag(sqrt(s2));
    [t,w] = gauher(N);
    one2N = num2cell(repmat(1:N,D,1),2);  % 2 is general here
    [grid_loc{1:numel(one2N)}] = ndgrid(one2N{:});
    x_loc = cell2mat(cellfun(@(x) t(x(:)),grid_loc,'UniformOutput',false));
    w_loc = cell2mat(cellfun(@(x) w(x(:)),grid_loc,'UniformOutput',false));
    wn = prod(w_loc,2);  % N**D
    xn = x_loc * sig + repmat(mu,N^D,1);  % N**DxD
      
end