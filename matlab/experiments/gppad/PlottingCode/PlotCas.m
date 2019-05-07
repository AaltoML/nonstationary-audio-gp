function PlotCas(A,y,FSamp,PlotOpts);

%LoadFigureSettings;  

% Figure settings 
FontName = 'Times'; FSlg = 12; FSsm = 10; 

LWthick = 1.5; LWthin = 1;

AmpCols = [1,0,0;
           1,0,1;
           1,1/2,0]; 


[T,M] = size(A);

% Size of plot areas
top = 0.05;
bottom = 0.1;
left = 0.1;
right = 0.03;
vspace = 0.01;

height = (1-top-bottom-M*vspace)/(M+1);
width = 1-left-right;

PosAx = [left,1-top-height,width,height];
Down = [0,height+vspace,0,0];

% Total envelope/times
a = prod(A,2);
ts = ([1:length(y)]-1)/FSamp;

visible = 'on';

% Plot each of the panels
fig1 = figure(1);

for m=1:M
  axm = axes('position',PosAx);
  hold on

  aCur = A(:,M-m+1);
  plot(ts,aCur,'-','linewidth',LWthick,'Color',AmpCols(M-m+1,:))
  set(gca,'ylim',[0 1]*max(aCur)*1.5)
  set(gca,'xlim',[0 max(ts)])
  set(gca,'fontname',FontName,'FontSize',FSsm)
  set(gca,'xticklabel','')
  mstr = num2str(M-m+1);

  PosAx = PosAx - Down;
end

ax = axes('position',PosAx);
hold on
plot(ts,y,'-k','linewidth',LWthin)
plot(ts,a,'-r','linewidth',LWthin)
plot(ts,-a,'-r','linewidth',LWthin)
set(gca,'ylim',[-1 1]*max(abs(y))*1.1)
set(gca,'xlim',[0 max(ts)])
set(gca,'fontname',FontName,'FontSize',FSsm)
xlabel('time /s','fontname',FontName,'FontSize',FSlg)

%set(gcf,'paperpositionmode','manual','paperposition',PP,'papersize',PS);
%print(gcf,'-depsc','-loose',[Location,Name,Num]);
%clf;
