function Varx = GetErrorBarGPPAD(x,y,Params,Dims,Opts)
  
  % function Varx = GetErrorBarGPPAD(x,y,Params,Dims,Opts)
  %
  % Approximate error-bars using Laplace's approximation for
  % the GPPAD model
  
[d,h] = GetHessian(x,y,Params,Dims);
[xv,lmb] = GetEigSigDSig(h,d,Opts);

lmb(lmb<0)=0;

Varx = ones(Dims.Tx,1)*Params.varx;

K = size(xv,2);

for k= 1:K
  b = ifft(h.^(1/2).*fft(xv(:,k)));
  Varx = Varx - abs(b).^2*lmb(k)/(lmb(k)+1);
end

if sum(Varx<0)>0

  clf
  subplot(4,1,1)
  plot(Varx)

  subplot(4,1,2)
  plot(lmb)
  
  subplot(4,1,3)
  a = log(1+exp(x));
  hold on
  plot(y/sqrt(var(y)),'-k')
  plot(a/sqrt(var(a)),'-r')
  
  subplot(4,1,4)
  plot(x)

  figure
  subplot(2,1,1)
  plot(d)
  
  subplot(2,1,2)
  plot(h)
 
  figure
  hold on
  plot(y(2000:4000)/sqrt(var(y)),'-k')
  plot(a(2000:4000)/sqrt(var(a)),'-r')

  keyboard
end