function [ W, SX, u ] = ut3_ws(n,kappa)
% [ W, SX, u ] = ut3_ws(n,[kappa])
%
% Return weights and sigma-points for 3rd order
% UT for dimension n with parameter kappa (default 1-n).

    if nargin < 2
        kappa = 1-n;
        kappa = 0; % Arno changed
    end
        
    % Weights
    W = zeros(1,2*n+1);
    for j=1:2*n+1
        if j==1
            wm = kappa / (n + kappa);
        else
            wm = 1 / (2 * (n + kappa));
        end
        W(j) = wm;
    end

    % Sigma points
    SX = [zeros(n,1) eye(n) -eye(n)];
    SX = sqrt(n + kappa)*SX;
    u  = sqrt(n + kappa);
end

