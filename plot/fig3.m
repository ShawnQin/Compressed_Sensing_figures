%% Figure 3 summarize how sparisty of W change with sigma_c, n, M, and N

close all
clear

%% define parameters used in graphics
defaultGraphicsSetttings
%define some colors using brewermap
RdBu = brewermap(11,'RdBu');   % red and blue
Bu = brewermap(9,'Blues');    % blues
Gr = brewermap(9,'Greys');    % greys

myBu = [3, 110, 184]/256;
%myOr = [224, 135, 51]/256;
myOr = [255, 185, 85]/256;

myRd = [202, 39, 44]/256;
lBu = [96,166,223]/255; %light blue
dpBu = [63,114,183]/255; % deep blue
dkBu = [50,78,147]/255;   %dark blue
Or = [220,150,71]/255;  % orange
brickRd = [201,69,89]/255;  %brick red
green = [107,169,81]/255;  %green
purple = [113,49,119]/255;  % purple


%% load the summary data
% dFolder = '../data/lognorm_SigDp';
dFolder = '../data';
saveFolder = '../figures';

%% define the layout of the graphics, we have 3x3 subplots

% ===============================================================
% ===============================================================
errorBarSize = 1.5;
errorMarkerSize  = 10;
LineWidth = 1.5;
labelFontSize = 24;
axisFontSize = 20;
ticketWidth = 1;

% basic parameters
N = 50;
M = 13;
sp = 3;
sig = 2;

colorInx = 1;   % the color index in a 11 order of colors

% sigma-dependent
fName = 'gcmi_sigdp_N50M13Sp3_19-Oct-2018.mat';
load(fullfile(dFolder,fName));
summSigdp = dataSumm;

% N-dependent
% NFile = 'gcmi_Ndp_M13Sp3sig2_19-Oct-2018.mat';
NFile = 'gcmi_Ndp_M13Sp3sig2_24-Oct-2018.mat';
load(fullfile(dFolder,NFile))
summNdp = dataSumm;

% sp-dependent
nFile = 'gcmi_spWdp_N50M13Sp9_19-Oct-2018.mat';
load(fullfile(dFolder,nFile))
summSpdp = dataSumm;

% M-dependent
Mfile = 'gcmi_Mdp_N50Sp3sig2_15-Dec-2018.mat';
load(fullfile(dFolder,Mfile))
summMdp = dataSumm;

%% first row plot how rho_w change with sigma_c, n, N, M

% ========================================================
% first row plot how rho_w change with sigma_c, n, N, M
% ========================================================
colInx1 = 9;
colInx2 = 9;
figure
figureSize = [0 0 16 2.8];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
[ha, pos] = tight_subplot(1,4,[.01 .08],[.23 .04],[.08 .01]);


inx1 = 4:1:17;  %the index corresponding to N = 100

axes(ha(1))
errorbar(summSigdp.allSig(inx1)',summSigdp.meanSpW(inx1),summSigdp.stdSpW(inx1),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',Gr(colInx1,:),'Color',Gr(colInx2,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend(['n = ',num2str(sp)],'Location','northwest');
legend boxoff
pbaspect([1.28 1 1])
xlim([1.8, 6.2])
ylim([0.64, 0.81])

xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_w$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

axes(ha(2))
inx2 = 1:1:5;
errorbar(summSpdp.sp(inx2)', summSpdp.meanSpW(inx2),summSpdp.stdSpW(inx2),'^-','MarkerSize',...
      errorMarkerSize,'MarkerFaceColor',Gr(colInx1,:),'Color',Gr(colInx2,:),'CapSize',0,'LineWidth',LineWidth)
lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
legend boxoff
pbaspect([1.28 1 1])
ylim([0.4,0.85])
xlim([1.8, 6.2])
xlabel('$n$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')


% allN = [40,50,70,90,100,120];
odorInx = [1,2,4:1:8];   %odor selected
axes(ha(3))
errorbar(summNdp.allN(odorInx)',summNdp.meanSpW(odorInx),summNdp.stdSpW(odorInx),'square-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',Gr(colInx1,:),'Color',Gr(colInx2,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('n=3,\sigma_c = 2','Location','northeast');
legend boxoff
pbaspect([1.28 1 1])
ylim([0.57,0.72])
xlim([37, 155])
xlabel('$N$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

MInx = 2:1:8;
axes(ha(4))
errorbar(summMdp.allM(MInx)',summMdp.meanSpW(MInx),summMdp.stdSpW(MInx),'diamond-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',Gr(colInx1,:),'Color',Gr(colInx2,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('n=3,\sigma_c = 2','Location','northeast');
legend boxoff
pbaspect([1.28 1 1])
ylim([0.5,0.62])
xlim([18, 51])
set(ha(4),'XTick',20:10:50)
xlabel('$M$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\rho_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

% prefix = ['fig3_gcmi_summ_rhow_all',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])


%% second row plot how sigma_w change with sigma_c, n, N, M
% ========================================================
% second row plot how sigma_w change with sigma_c, n, N, M
% ========================================================
figure
figureSize = [0 0 14 3];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
[ha, pos] = tight_subplot(1,4,[.01 .07],[.23 .04],[.05 .01]);


inx1 = 4:1:17;  %the index corresponding to N = 100

axes(ha(1))
errorbar(summSigdp.allSig(inx1)',summSigdp.meanSigW(inx1),summSigdp.stdSigW(inx1),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend(['n = ',num2str(sp)],'Location','northwest');
legend boxoff
xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_w$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

axes(ha(2))
inx2 = 1:1:5;
errorbar(summSpdp.sp(inx2)', summSpdp.meanSigW(inx2),summSpdp.stdSigW(inx2),'o-','MarkerSize',...
      errorMarkerSize,'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),'CapSize',0,'LineWidth',LineWidth)
lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
legend boxoff
ylim([1,1.5])
xlabel('$n$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

% allN = [40,50,70,90,100,120];
odorInx = [1,2,4:1:8];   %odor selected
axes(ha(3))
errorbar(summNdp.allN(odorInx)',summNdp.meanSigW(odorInx),summNdp.stdSigW(odorInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('n=3,\sigma_c = 2','Location','northeast');
legend boxoff
ylim([0.9,1.55])
xlabel('$N$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

MInx = 2:1:8;
axes(ha(4))
errorbar(summMdp.allM(MInx)',summMdp.meanSigW(MInx),summMdp.stdSigW(MInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',purple,'Color',purple,'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('N = 50,n=3,\sigma_c = 2','Location','northeast');
legend boxoff
ylim([0.9,1.4])
xlabel('$M$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

% prefix = ['fig3_gcmi_summ_sigW_all',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%% A row show how mu_w change
figure
figureSize = [0 0 14 3];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
[ha, pos] = tight_subplot(1,4,[.01 .07],[.23 .04],[.05 .01]);


inx1 = 4:1:17;  %the index corresponding to N = 100

axes(ha(1))
errorbar(summSigdp.allSig(inx1)',summSigdp.meanW(inx1),summSigdp.stdW(inx1),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend(['n = ',num2str(sp)],'Location','northwest');
legend boxoff
xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\mu_w$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

axes(ha(2))
inx2 = 1:1:5;
errorbar(summSpdp.sp(inx2)', summSpdp.meanW(inx2),summSpdp.stdW(inx2),'o-','MarkerSize',...
      errorMarkerSize,'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),'CapSize',0,'LineWidth',LineWidth)
lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
legend boxoff
% ylim([1,1.5])
xlabel('$n$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

% allN = [40,50,70,90,100,120];
odorInx = [1,2,4:1:8];   %odor selected
axes(ha(3))
errorbar(summNdp.allN(odorInx)',summNdp.meanW(odorInx),summNdp.stdW(odorInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('n=3,\sigma_c = 2','Location','northeast');
legend boxoff
% ylim([0.9,1.55])
xlabel('$N$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

MInx = 2:1:8;
axes(ha(4))
errorbar(summMdp.allM(MInx)',summMdp.meanW(MInx),summMdp.stdW(MInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',purple,'Color',purple,'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend('N = 50,n=3,\sigma_c = 2','Location','northeast');
legend boxoff
% ylim([0.9,1.4])
xlabel('$M$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

% prefix = ['fig3_gcmi_summ_muW_all',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])


%%
% ================================================================
% mu_w and rho_w, different sigma_c
% ================================================================
figure
figureSize = [0 0 3.8 3.2];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
[ha, pos] = tight_subplot(2,1,[.0 .2],[0.2,0.03],[0.2,0.07]);

axes(ha(1))
errorbar(summSigdp.allSig(inx1)',summSigdp.meanW(inx1),summSigdp.stdW(inx1),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
legend boxoff
% xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

axes(ha(2))
errorbar(summSigdp.allSig(inx1)',summSigdp.meanSigW(inx1),summSigdp.stdSigW(inx1),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
% lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
% legend boxoff
xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')
% prefix = ['fig3_gcmi_muRhoW_sigc_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%%
% ================================================================
% mu_w and rho_w, different n
% ================================================================
figure
% figureSize = [0 0 3.8 3];
figureSize = [0 0 3.8 3.2];


set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
% [ha, pos] = tight_subplot(2,1,[.0 .2],[0.2,0.03],[0.2,0.07]);
[ha, pos] = tight_subplot(2,1,[.0 .25],[0.22,0.03],[0.22,0.07]);


axes(ha(1))
errorbar(summSpdp.sp(inx2)',summSpdp.meanW(inx2),summSpdp.stdW(inx2),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
% lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
% legend boxoff
% xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

axes(ha(2))
errorbar(summSpdp.sp(inx2)',summSpdp.meanSigW(inx2),summSpdp.stdSigW(inx2),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myOr(colorInx,:),'Color',myOr(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
% lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
% legend boxoff
ylim([1,1.5])
xlabel('$n$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')
% prefix = ['fig3_gcmi_muRhoW_n_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%%
% ================================================================
% mu_w and rho_w, different N
% ================================================================
figure
figureSize = [0 0 3.8 3.2];
[ha, pos] = tight_subplot(2,1,[.0 .25],[0.22,0.03],[0.22,0.07]);

set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
[ha, pos] = tight_subplot(2,1,[.0 .2],[0.2,0.03],[0.2,0.07]);

axes(ha(1))
errorbar(summNdp.allN(odorInx)',summNdp.meanW(odorInx),summNdp.stdW(odorInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
% lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
% legend boxoff
ylim([-1.6,-1])
% xlabel('$\sigma_c$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

axes(ha(2))
errorbar(summNdp.allN(odorInx)',summNdp.meanSigW(odorInx),summNdp.stdSigW(odorInx),'o-','MarkerSize',errorMarkerSize,...
        'MarkerFaceColor',myRd(colorInx,:),'Color',myRd(colorInx,:),'LineWidth',LineWidth,...
        'CapSize',0)
% lg = legend(['\sigma_c = ',num2str(sig)],'Location','northeast');
% legend boxoff
% ylim([1,1.5])
xlabel('$N$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')
% prefix = ['fig3_gcmi_muRhoW_N_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%%
% ================================================================
% mu_w and rho_w, different M
% ================================================================
figure
figureSize = [0 0 3.8 3.2];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the margins
[ha, pos] = tight_subplot(2,1,[.0 .25],[0.22,0.03],[0.22,0.07]);

MInx = 2:1:8;
axes(ha(1))
errorbar(summMdp.allM(MInx)',summMdp.meanW(MInx),summMdp.stdW(MInx),'o-',...
    'MarkerSize',errorMarkerSize,'MarkerFaceColor',purple,'Color',purple,...
    'LineWidth',LineWidth,'CapSize',0)

ylim([-1.4,-1])
ylabel('$\mu_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear','XTick',[])

axes(ha(2))
errorbar(summMdp.allM(MInx)',summMdp.meanSigW(MInx),summMdp.stdSigW(MInx),...
    'o-','MarkerSize',errorMarkerSize,'MarkerFaceColor',purple,...
    'Color',purple,'LineWidth',LineWidth,'CapSize',0)
ylim([0.9,1.3])
xlabel('$M$','Interpreter','latex','FontSize',labelFontSize)
ylabel('$\sigma_{w}$','Interpreter','latex','FontSize',labelFontSize)
set(gca,'FontSize',axisFontSize,'LineWidth',ticketWidth,'XScale','linear')

% prefix = ['fig3_gcmi_muRhoW_M_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])