function x = GetX(a)

  % function x = GetX(a)
  %
  % converts amplitudes into transformed envelopes
    
  x = log(exp(a)-1);
  
  ind1 = a>500;
  
  x(ind1) = a(ind1);
  
  ind2 = a<exp(-20);
  x(ind2) = log(a(ind2));