function [ W, SX, u, v ] = ut7_ws(n)
% [ W, SX, u, v ] = ut7_ws(n)
%
% Return weights and sigma-points for 7th order
% UT for dimension n.

    % The weights and sigma-point from McNamee & Stenger
    I222 = 1;
    I22 = 1;
    I24 = 3;
    I2 = 1;
    I6 = 15;
    I4 = 3;
    I0 = 1;
    
    tmp = roots([I2^2-I0*I4 0 -(I2*I4-I0*I6) 0 (I4^2-I2*I6)]);
    tmp = tmp(tmp > 0);
    u = tmp(1);
    v = tmp(2);
    
    u2 = u*u;
    u4 = u2*u2;
    u6 = u4*u2;
    v2 = v*v;
    v4 = v2*v2;
    v6 = v4*v2;
    
    A111 = I222 / 8 / u6;
    tmp = 0.25 * ([u4 v4; u6 v6] \ ([I22; I24] - 8*(n-2)*[u4; u6] * A111));
    A11 = tmp(1);
    A22 = tmp(2);
    
    tmp = -2*(n-1)*[A11; A22] + 0.5 * ([u2 v2; u4 v4] \ ([I2; I4] - 8*(n-1)*(n-2)/2*[u2; u4]*A111));
    A1 = tmp(1);
    A2 = tmp(2);
    
    A0 = I0 - 2*n*(A1+A2) - 4*n*(n-1)/2*(A11+A22)-8*n*(n-1)*(n-2)/6*A111;
    
    U0 = sym_set(n,[]);
    U1 = sym_set(n,[u]);
    V1 = sym_set(n,[v]);
    U2 = sym_set(n,[u u]);
    V2 = sym_set(n,[v v]);
    U3 = sym_set(n,[u u u]);
    
    SX = [U0 U1 V1 U2 V2 U3];
    W  = [A0*ones(1,size(U0,2))  A1*ones(1,size(U1,2)) A2*ones(1,size(V1,2)) ...
              A11*ones(1,size(U2,2)) A22*ones(1,size(V2,2)) A111*ones(1,size(U3,2))];
end

