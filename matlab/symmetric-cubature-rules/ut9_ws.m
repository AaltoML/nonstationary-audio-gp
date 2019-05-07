function [ W, SX, u, v ] = ut9_ws(n)
% [ W, SX, u, v ] = ut9_ws(n)
%
% Return weights and sigma-points for 9th order
% UT for dimension n.
%
% See this reference for details:
%
%     @article{McNamee+Stenger:1967,
%           title = {Construction of fully symmetric numerical integration formulas},
%          author = {McNamee, John and Stenger, Frank},
%         journal = {Numerische Mathematik},
%          volume = {10},
%          number = {4},
%           pages = {327--344},
%            year = {1967},
%       publisher = {Springer}
%     }
% 

    % The weights and sigma-point from McNamee & Stenger
    
    % Calculate I_{i,j,k} such that
    %   I_{i,j} = \int \int x^i y^j N(x,y) dx dy
    % et cetara for the rest. 
    I2222 = 1;
    I224 = 3;
    I222 = 1; 
    I44 = 9;
    I26 = 15;
    I24 = 3;
    I22 = 1;
    I8  = 105; % 8th central moment
    I6 = 15;
    I4 = 3;
    I2 = 1;
    I0 = 1;
    
    % The unknowns u and v are the zeros of the polynomial
    tmp = roots([I4^2-I2*I6 0 -(I4*I6-I2*I8) 0 (I6^2-I4*I8)]);
    tmp = tmp(tmp > 0);
    u = tmp(1);
    v = tmp(2);
    
    % The powers of u and v
    u2 = u*u;
    u4 = u2*u2;
    u6 = u4*u2;
    u8 = u4*u4;
    v2 = v*v;
    v4 = v2*v2;
    v6 = v4*v2;
    v8 = v4*v4;
    
    % This is valid, provided that the following holds
    %   I26-(u2+v2)*I24+u2*v2*I22 = 0
    
    % The linear unknowns are given by
    A1111 = I2222 / 16 / u8;

    tmp = 1/8*([u6 v6; u8 v8] \ ([I222; I224]-16*(n-3)*A1111*[u6; u8]));
    A111 = tmp(1);
    A222 = tmp(2);

    A12 = (I26 - I44) / (4*u2*v2*(u2-v2)^2);
    
    tmp = -2*(n-2)*[A111; A222] + 1/4*([u6 v6; u8 v8] \ ...
      ([I24; I26] - 4*[u4*v2+u2*v4; u6*v2+u2*v6]*A12 - ...
      16*ndownk(n-2,2)*[u6; u8]*A1111));
    A11 = tmp(1);
    A22 = tmp(2);
    
    tmp = -2*(n-1)*[A11+A12; A22+A12] - 4*ndownk(n-1,2)*[A111; A222] + ...
      0.5*([u2 v2; u4 v4] \ ([I2; I4] - 16*ndownk(n-1,3)*[u2; u4]*A1111));
    A1 = tmp(1);
    A2 = tmp(2);
        
    A0 = I0 - 2*n*(A1+A2) - 4*ndownk(n,2)*(A11 + 2*A12 + A22) - ...
      -8*ndownk(n,3)*(A111+A222) - 16*ndownk(n,4)*A1111;
  
    % We could get rid of the ndownk..

    % The set of evaluation points
    U0 = sym_set(n,[]);
    U1 = sym_set(n,[u]);
    V1 = sym_set(n,[v]);
    U2 = sym_set(n,[u u]);
    UV = sym_set(n,[u v]);
    V2 = sym_set(n,[v v]);
    U3 = sym_set(n,[u u u]);
    V3 = sym_set(n,[v v v]);
    U4 = sym_set(n,[u u u u]);
    
    % Points
    SX = [U0 U1 V1 U2 UV V2 U3 V3 U4];
    
    % Corresponding weights
    W  = [A0*ones(1,size(U0,2))  A1*ones(1,size(U1,2)) A2*ones(1,size(V1,2)) ...
          A11*ones(1,size(U2,2)) A12*ones(1,size(UV,2)) A22*ones(1,size(V2,2)) ...
          A111*ones(1,size(U3,2)) A222*ones(1,size(V3,2)) A1111*ones(1,size(U4,2))];
   
end

function N=nupk(n,k) % see also: rising factorial
  N=prod((n-k+1):n);
end

function N=ndownk(n,k)
  N=nupk(n,k)/factorial(k);
end

%% Debug  
%     % <lazy>    
%     
%     % Define problem
%     B = [1 2*n  2*n  2^2*ndownk(n,2) 2^3*ndownk(n,2)   2^2*ndownk(n,3) 2^3*ndownk(n,3)      2^3*ndownk(n,3)      2^4*ndownk(n,4);
%          0 2*u2 2*v2 2^2*(n-1)*u2    2^2*(n-1)*(u2+v2) 2^2*(n-1)*v2    2^3*ndownk(n-1,2)*u2 2^3*ndownk(n-1,2)*v2 2^4*ndownk(n-1,3)*u2; 
%          0 2*u4 2*v4 2^2*(n-1)*u4    2^2*(n-1)*(u4+v4) 2^2*(n-1)*v4    2^3*ndownk(n-1,2)*u4 2^3*ndownk(n-1,2)*v4 2^4*ndownk(n-1,3)*u4; 
%          0 2*u6 2*v6 2^2*(n-1)*u6    2^2*(n-1)*(u6+v6) 2^2*(n-1)*v6    2^3*ndownk(n-1,2)*u6 2^3*ndownk(n-1,2)*v6 2^4*ndownk(n-1,3)*u6; 
%          0 2*u8 2*v8 2^2*(n-1)*u8    2^2*(n-1)*(u8+v8) 2^2*(n-1)*v8    2^3*ndownk(n-1,2)*u8 2^3*ndownk(n-1,2)*v8 2^4*ndownk(n-1,3)*u8;
%          0 0    0    2^2*u4          2^2*(u2*v2+u2*v2) 2^2*v4          2^3*(n-2)*u4         2^3*(n-2)*v4         2^4*ndownk(n-1,2)*u4;
%          0 0    0    2^2*u6          2^2*(u4*v2+u2*v4) 2^2*v6          2^3*(n-2)*u6         2^3*(n-2)*v6         2^4*ndownk(n-1,2)*u6;
%          0 0    0    2^2*u8          2^2*(u6*v2+u2*v6) 2^2*v8          2^3*(n-2)*u8         2^3*(n-2)*v8         2^4*ndownk(n-1,2)*u8;
%          0 0    0    2^2*u8          2^2*(u4*v4+u4*v4) 2^2*v8          2^3*(n-2)*u8         2^3*(n-2)*v8         2^4*ndownk(n-1,2)*u8;
%          0 0    0    0               0                 0               2^3*u6               2^3*v6               2^4*(n-3)*u6;
%          0 0    0    0               0                 0               2^3*u8               2^3*v8               2^4*(n-3)*u8;
%          0 0    0    0               0                 0               0                    0                    2^4*u8               ];
%     
%     II = [I0; I2; I4; I6; I8; I22; I24; I26; I44; I222; I224; I2222];
%     
%     % Solve
%     A = B\II;
%     % B([1 2 3 6 7 9 10 11 12],:)\II([1 2 3 6 7 9 10 11 12])
%     
%     A1111 = I2222 / 16 / u8;
% 
%     tmp = 1/8*([u6 v6; u8 v8] \ ([I222; I224]-16*(n-3)*A1111*[u6; u8]));
%     A111 = tmp(1);
%     A222 = tmp(2);
% 
%    % Set variables
%    A0 = A(1); A1 = A(2); A2 = A(3); A11 = A(4); A12 = A(5); A22 = A(6);
%    A111 = A(7); A222 = A(8); A1111 = A(9);
%    
%    % </lazy>
    


