function PlotSpectrogram(A,FSamp,CF,PlotOpts)

  
FontName = 'Times'; FSlg = 12;
time = [1:size(A,1)]/FSamp;
FLim = [min(CF),max(CF)];
T = size(A,1);

if isfield(PlotOpts,'PP')
    PP = PlotOpts.PP;
    DS = floor(T/PP);
    time = time(1:DS:T);
    A = A(1:DS:T,:); 
end

figure
hold on
title('Log amplitudes')
surf(time,CF,log(A'),'linestyle','none')
view(0,90)
xlabel('time /s','FontName',FontName,'FontSize',FSlg)
ylabel('frequency /Hz','FontName',FontName,'FontSize', FSlg)
set(gca,'yscale','log','ylim',FLim)
drawnow;


if isfield(PlotOpts,'Rng')
    Rng = PlotOpts.Rng;
    set(gca,'xlim',Rng);
else
  set(gca,'xlim',[min(time),max(time)]);
end
