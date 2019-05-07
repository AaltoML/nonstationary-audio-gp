function [ W, SX, u ] = ut5_ws(n)
% [ W, SX, u ] = ut5_ws(n)
%
% Return weights and sigma-points for 5th order
% UT for dimension n.

    % The weights and sigma-point from McNamee & Stenger
    I0 = 1;
    I2 = 1;
    I4 = 3;
    I22 = 1;
    
    u = sqrt(I4/I2);
    
    A0 = I0 - n*(I2/I4)^2 * (I4 - 0.5 * (n-1)*I22);
    A1 = 0.5 * (I2/I4)^2*(I4 - (n-1)*I22);
    A11 = 0.25*(I2/I4)^2*I22;
 
    U0 = sym_set(n,[]);
    U1 = sym_set(n,[u]);
    U2 = sym_set(n,[u u]);
    
    SX = [U0 U1 U2];
    W  = [A0*ones(1,size(U0,2)) A1*ones(1,size(U1,2)) A11*ones(1,size(U2,2))];
    
%     % Check the matrix
%     B = [1 2*n 2*nupk(n,2); 
%          0 2*u^2 4*(n-1)*u^2;
%          0 2*u^4 4*(n-1)*u^4;
%          0 0     4*u^4];
%     B*[A0; A1; A11]
%     [I0; I2; I4; I22]
%     
    
end

function N=nupk(n,k) % see also: rising factorial
  N=prod((n-k+1):n);
end

function N=ndownk(n,k)
  N=nupk(n,k)/factorial(k);
end