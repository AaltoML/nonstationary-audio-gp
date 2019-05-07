%NOW modified so cf spacing emulates DSAM - outputs now look very similar
%also 23/08/2007

function [forward,feedback,cf,ERB,B]=GammaToneMake(fs,numChannels,loFreq,hiFreq,method)
% [forward,feedback,cf,ERB,B]=GammaToneMake(fs,numChannels,loFreq,hiFreq,method)

% computes the filter coefficients for a bank of Gammatone filters. These
% filters were defined by Patterson and Holdworth for simulating
% the cochlea. The results are returned as arrays of filter
% coefficients. Each row of the filter arrays (forward and feedback)
% can be passed to the MatLab "filter" function.

if loFreq < 75; loFreq = 75; end %line to stop errors when a very low fequency is used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I added these parameters to play with (Nick C 2007) from Apple TR #35, to
% make it easier to switch spacing method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5; method = 'Moore';warning('No method specified: using Moore/Glasberg'); end
switch lower(method)
    case{'lyon','stanley'}
        EarQ = 8;       % Lyon + Stanley Parameters (1988)
        minBW = 125;
        order = 2;
    case{'greenwood'}
        EarQ = 7.23824; % Greenwood Parameters (1990) as (nearly) in DSAM
        minBW = 22.8509;
        order = 1;
    case{'moore','glasberg'}
        EarQ = 9.26449; % Glasberg and Moore Parameters (1990)
        minBW = 24.7;
        order = 1;
    case{'wierddsam'}
        EarQ = 9.26; % Glasberg and Moore Parameters (1990)
        minBW = 15.719; %WORKS IF YOU SWCREW WITH THIS PARAMETER AS SO. . . 
        order = 1;
    otherwise
        warning('Method specified does not exist: using Moore/Glasberg');
        EarQ = 9.26449; % Glasberg and Moore Parameters (1990)
        minBW = 24.7;
        order = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All of the following expressions are derived in Apple TR #35, "An
% Efficient Implementation of the Patterson-Holdsworth Cochlear
% Filter Bank.", with exeption of the modification of the cf calculation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=1/fs;
%cf = -(EarQ*minBW)+exp((1:numChannels)'*( -log(hiFreq + EarQ*minBW) + log(loFreq + EarQ*minBW) )/numChannels)*(hiFreq + EarQ*minBW);

ERBlo = ((loFreq/EarQ)^order + minBW^order) ^ (1/order);
ERBhi = ((hiFreq/EarQ)^order + minBW^order) ^ (1/order);
overlap = (ERBhi/ERBlo)^(1/(numChannels-1));
ERB = ERBlo * (overlap.^(0:numChannels-1));
cf = EarQ*(((ERB.^order) - (minBW.^order)).^(1/order));

B=1.019*2*pi*ERB; %in rad here - note to self: some models require B in Hz (NC)
gain = abs((-2*exp(4*i*cf*pi*T)*T + 2*exp(-(B*T) + 2*i*cf*pi*T).*T.*(cos(2*cf*pi*T) - sqrt(3 - 2^(3/2))*sin(2*cf*pi*T))) .*(-2*exp(4*i*cf*pi*T)*T +2*exp(-(B*T) + 2*i*cf*pi*T).*T.*(cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) *sin(2*cf*pi*T))).*(-2*exp(4*i*cf*pi*T)*T +2*exp(-(B*T) + 2*i*cf*pi*T).*T.*(cos(2*cf*pi*T) -sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) .*(-2*exp(4*i*cf*pi*T)*T+2*exp(-(B*T) + 2*i*cf*pi*T).*T.*(cos(2*cf*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) ./(-2 ./ exp(2*B*T) - 2*exp(4*i*cf*pi*T) +2*(1 + exp(4*i*cf*pi*T))./exp(B*T)).^4);
feedback=zeros(length(cf),9);
forward=zeros(length(cf),5);
forward(:,1) = T^4 ./ gain;
forward(:,2) = -4*T^4*cos(2*cf*pi*T)./exp(B*T)./gain;
forward(:,3) = 6*T^4*cos(4*cf*pi*T)./exp(2*B*T)./gain;
forward(:,4) = -4*T^4*cos(6*cf*pi*T)./exp(3*B*T)./gain;
forward(:,5) = T^4*cos(8*cf*pi*T)./exp(4*B*T)./gain;
feedback(:,1) = ones(length(cf),1);
feedback(:,2) = -8*cos(2*cf*pi*T)./exp(B*T);
feedback(:,3) = 4*(4 + 3*cos(4*cf*pi*T))./exp(2*B*T);
feedback(:,4) = -8*(6*cos(2*cf*pi*T) + cos(6*cf*pi*T))./exp(3*B*T);
feedback(:,5) = 2*(18 + 16*cos(4*cf*pi*T) + cos(8*cf*pi*T))./exp(4*B*T);
feedback(:,6) = -8*(6*cos(2*cf*pi*T) + cos(6*cf*pi*T))./exp(5*B*T);
feedback(:,7) = 4*(4 + 3*cos(4*cf*pi*T))./exp(6*B*T);
feedback(:,8) = -8*cos(2*cf*pi*T)./exp(7*B*T);
feedback(:,9) = exp(-8*B*T);