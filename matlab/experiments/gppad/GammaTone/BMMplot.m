function BMMplot(stimIn, fc, sr, fcLabel) 
%function BMMplot(stimIn, fc, sr, fcLabel) 

% BASILAR MEMBRANE MOTION PLOTTER v0.1
% function by Nick Clark ~ nick@ihr.mrc.ac.uk

% This is a simple plotting routine which I put together to mimic the BMM
% plot types seen frequently in “Journal Of The Acoustical Society Of 
% America” (JASA) articles among others, as well as software such as AIM 
% and AMS.  This allows data from each channel to be viewed as stacked 
% line graphs.

stimOut = zeros(size(stimIn));
for n = 1:size(stimIn,1)
    stimOut(n, :) = n;
end
stimOut = stimIn./max(max(abs(stimIn))) + stimOut;
timeAx = [1/sr:1/sr:size(stimIn,2)*(1/sr)]*1e3;

plot(timeAx,stimOut','k');
set(gca, 'YTick',        fcLabel, ...
         'YTickLabel',   fc(fcLabel),...
         'TickDir',      'out',...
         'Xlim',         [timeAx(1) timeAx(end)],...
         'Ylim',         [0 size(stimIn,1)+1],...
         'Box',          'off');
axis on
xlabel('Time[ms]');
ylabel('Centre Frequency[Hz]');