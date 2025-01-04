%% Fig S5 Supplemental Figure Mean Field Theory

% set the color 
Bu = brewermap(11,'Blues');    % blues

dFolder = '../data';
saveFolder = '../figures';

% ===============================================================
% define the layout of the graphics, we have 3x3 subplots
% ===============================================================
errorMarkerSize  = 10;
LineWidth = 2;
labelFontSize = 28;
axisFontSize = 24;
ticketWidth = 1.5;

colorInx = 1;   % the color index in a 11 order of colors

figureSize = [0 0 13 4];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
[ha, pos] = tight_subplot(1,3,[.07 .07],[.25 .05],[.1 .03]);


fName = 'MFTsigmac_summ.mat';
load(fullfile(dFolder,fName));

% sp of W vs sigmac
axes(ha(1))
plot(sigall,sparsity,'Color',Bu(10,:))
% lg = legend('N = 100,n = 2','Location','southeast');
% legend boxoff
% xlim([1,5])
xlabel('$\sigma_{c}$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'XTick',1:1:5,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

% mean W vs sigmac
axes(ha(2))
plot(sigall,meanvalue,'Color',Bu(10,:))
% xlim([1,5])
xlabel('$\sigma_{c}$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'XTick',1:1:5,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

%std w vs sigmac
axes(ha(3))
plot(sigall,stdvalue,'Color',Bu(10,:))
% xlim([1,5])
xlabel('$\sigma_{c}$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'XTick',1:1:5,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

%save the figure
% prefix = ['MFT_figure5_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])
