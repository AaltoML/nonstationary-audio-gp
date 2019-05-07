function PlotFilterBank(CF,ERB,FSamp)
  
  % Plots the spectra of a filterbank of a gamma-tone filterbank
  
  % Inputs  
  % CF = centre frequencies [D,1]
  % ERB = equivalent rectangular bandwidths [D,1]
  % FSamp = sample rate  
    
  N = 4; % order of filters
  B = 1.019*24.7*(4.37e-3*CF+1);
  
  D = length(CF);
  NumPoints = 4000;
  Freqs = linspace(0,FSamp/2,NumPoints); 
  
  Spec = zeros(D,NumPoints);
  
  for d=1:D
    temp = factorial(N-1)*(B(d)./(2*pi*i*(CF(d)-Freqs)+2*pi*B(d))).^N;
    Spec(d,:) = temp.*conj(temp);
    Spec(d,:) = Spec(d,:)/sum(Spec(d,:));
  end
    
  fig1= figure;
  hold on;
  
  for d=1:D,
    plot(Freqs,Spec(d,:),'-k')
  end
  