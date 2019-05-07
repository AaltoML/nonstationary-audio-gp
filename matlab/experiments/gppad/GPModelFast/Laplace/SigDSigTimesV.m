function b = SigDSigTimesV(v,d2py,d2pxhalf);

TT = size(d2pxhalf,1);
T = size(d2py,1);

a = fft(v);
a = d2pxhalf.*a;
a = ifft(a);

b = zeros(TT,1);
b(1:T,1) = d2py.*a(1:T,1);
b = fft(b);
b = d2pxhalf.*b;
b = ifft(b);
