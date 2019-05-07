function [z,a,c] = InitPADSample(y,Params,WinSz);
  
  % Initialises PAD
    
T = length(y);    

if WinSz>=T-1
  WinSz = T-2;
end


% A = Params.A;
% len = Params.len;
% sh = Params.sh;
C = 1./Params.varc;
%vary = Params.vary;

% % Slow method
% yAbsExt = abs([y(1)*ones(WinSz,1);y;y(T)*ones(WinSz,1)]);
%
% a2 = zeros(T,1);
% for t=1:T
%   Points = WinSz + [t-WinSz:t+WinSz];
%   a2(t) = mean(yAbsExt(Points))*sqrt(C);
% end

% Optimised method

Win = [ones(WinSz,1);zeros(2*T-2*WinSz-2,1);ones(WinSz,1)]/(2*WinSz);
fftWin = fft(Win);
fftAbsy = fft([abs(y);zeros(T-2,1)]);

a = real(ifft(fftWin.*fftAbsy))*sqrt(C);
a = a(1:T,1);
tiny = 1e-8;
a(a<tiny) = tiny;
z = log(exp(a)-1);
c = y./a;

% figure
% hold on
% plot(a,'-r','linewidth',2)
% plot(a2(1:T),'-m')
% keyboard