% Tests the gammatone function and the way in which we
% want to use it.

% Tests the placing of the gammatone centre frequencies

A = 24.7*4.37e-3*1.019;
B = 24.7*1.019;
N = 10;

f = zeros(N,1);

Gamma = 2.5;

f(1) = 100;


for n=2:N
  f(n,1) = f(n-1,1)*(2+A*Gamma)/(2-A*Gamma)+2*B*Gamma/(2-A*Gamma);
end

df = A*f+B;

% figure
% hold on
% plot(f,ones(N,1),'.k')
% plot(f+df/2,ones(N,1),'.b')
% plot(f-df/2,ones(N,1),'.b')
  
% Test the gammatone function itself

SoundFiles = '/home/rich/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/';
CurrentFile = '03 - Birds & Squirrels.wav';
%CurrentFile = '39 - Mountain Stream.wav';

[y,FSamp,Bits] = wavread([SoundFiles,CurrentFile],'native');

Range = 100000+[1:10000];
yGT = gammatone_c(y(Range,1), FSamp, 1000);

%Range = [1000:4:10000];

figure
hold on
plot(y(Range))
plot(yGT,'-k')