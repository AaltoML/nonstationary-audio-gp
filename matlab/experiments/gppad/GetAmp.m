function a = GetAmp(x)
  
  % function a = GetAmp(x)
  %
  % converts transformed envelopes into envelopes
  % avoiding numerical problems
   
  ind1 = x>500;
  ind2 = x<-30;
  
  a = log(1+exp(x));
  
  a(ind1) = x(ind1);
  
  a(ind2) = exp(x(ind2));