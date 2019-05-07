function [xv,lmb] = LanczosGPPAD(d2pxhalf,d2py,vstart, ...
                                 epert,nev,jmax,tresh,tolconv)

% function [xv,lmb] = ...
%   LanczosGPPAD(d2pxhalf,d2py,vstart,epert,nev,jmax,tresh,tolconv)
  
% Runs Lanczos on a Hermitian matrix, for the GPPAD model
% selective orthogonalization as in Parlett et al Fortran cod
% Incoming parameters:
% n matrix order
% nev number of eigenvalues sought
% tolconv requested accuracy of eigenvalues, default at 100*eps
% jmax maximum size of basis, max number of matrix vector multiplies
% vstart starting vector default randn(n,1)
%
cpustart=cputime;
%tresh=sqrt(eps);
%global A R % COMMENTED BY RICH
reps=10*sqrt(eps);
eps1=eps;
wjm1=[0 0];
wj=[0 eps1 eps1];
wjp1=[];
wjps=[];
flagreort=0;
nconv=0;
r=vstart;
v=[];
alpha=[];
beta=[];
beta(1)=norm(r);
for j=1:jmax,
  % Basic recursion
  v(:,j)=r/beta(j);

 % r=AOP(v(:,j)); % Editted by Rich
 r = SigDSigTimesV(v(:,j),d2py,d2pxhalf);

  if j>1, r=r-v(:,j-1)*beta(j); end;
  alpha(j)=r'*v(:,j);
  r=r-v(:,j)*alpha(j);
  beta(j+1)=norm(r);
  % Estimate |A|_2 by |T|_1
  if j==1, anorm=abs(alpha(1)+beta(2));
     else
        anorm=max(anorm,beta(j)+abs(alpha(j))+beta(j+1));
     end;
  epert=anorm*eps;

  % Update recurrence to determine selective reorth
  wjp1(1)=0;
  jv=2:j;
  wjp1(jv)=beta(jv).*wj(jv+1)+(alpha(jv-1)-alpha(j)).*wj(jv)+beta(jv-1).*wj(jv-1)-beta(j)*wjm1(jv);
  wjp1(jv)=(wjp1(jv)+sign(wjp1(jv))*epert)/beta(j+1);
  wjp1(j+1)=epert/beta(j+1);
  wjp1(j+2)=eps1;
  wjps(j)=max(abs(wjp1));
  % Test if it is time to reorthogonalize
  if max(abs(wjp1)) > tresh,
    flagreort=1;
    vjj1=[v(:,j) r];
    h1=v(:,1:j-1)'*vjj1;
    vjj1=vjj1-v(:,1:j-1)*h1;
    if norm(h1)>reps,
      vjj1=vjj1-v(:,1:j-1)*(v(:,1:j-1)'*vjj1);
    end;
    r=vjj1(:,2)-vjj1(:,1)*(vjj1(:,1)'*vjj1(:,2));
    v(:,j)=vjj1(:,1);
    wjp1(2:j+1)=eps1;
    wj(2:j)=eps1;
  end;
  wjm1=wj;
  wj=wjp1;
  % Is it time to test for convergence
  if flagreort | j==jmax | beta(j+1)<epert,
    
    temp = diag(alpha)+diag(beta(2:j),1)+diag(beta(2:j),-1);
    if sum(isnan(temp(:)))>0||sum(isinf(temp(:)))>0
      lmb = NaN;
      xv = NaN;
      disp('NaN Lanczos Error')
      keyboard
      return;
    end
    [s d]=eig(diag(alpha)+diag(beta(2:j),1)+diag(beta(2:j),-1));
    bndv=abs(beta(j+1)*s(j,:));
    convd=bndv<tolconv*anorm;
    nconv=sum(convd);
%    fprintf('%4.0f %3.0f %8.2e %8.2e %15.8e %15.8e \n',j,nconv,cputime-cpustart,wjps(j),alpha(j),beta(j+1));
    flagreort=0;
  end;
  if beta(j+1)<epert | nconv>=nev , break; end;
end;
iconv=find(convd);
[lmb ilmb]=sort(-diag(d(iconv,iconv)));
lmb=-lmb;
xv=v*s(:,iconv(ilmb));
