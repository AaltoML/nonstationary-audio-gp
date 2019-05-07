function [Spec,Freq,ERB] = GetGammaToneSpec(CF,FSamp,Flim,NumPoints);

  % function [Spec,Freqs] = GetGammaToneSpec(CF,ERB,FSamp,Flim,NumPoints);

    Freq = linspace(Flim(1),Flim(2),NumPoints);
    N = 4; % order of filters
    ERB = 24.7*(4.37e-3*CF+1);
    B = 1.019*ERB; 

    temp = factorial(N-1)*(B./(2*pi*i*(CF-Freq*FSamp)+2*pi*B)).^N ...
           +factorial(N-1)*(B./(2*pi*i*(CF+Freq*FSamp)+2*pi*B)).^N;
    Spec = temp.*conj(temp);
    Spec = Spec/sum(Spec);
