% This program plot figure 5 and supporting figure 6


close all
clear

%% define parameters used in graphics
defaultGraphicsSetttings
%define some colors using brewermap
RdBu = brewermap(11,'RdBu');   % red and blue
Bu = brewermap(11,'Blues');    % blues
Gr = brewermap(11,'Greys');    % greys

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

%% Data and figure path
% dFolder = '/Users/shan/Documents/GoogleDrive/olfactoryCoding/code/decoding/data/reconsFig6';
% outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
% use relative path
dFolder = '../data/reconsFig5';
outFolder = '../figures';

%% Plot figure 5b
% load the data
dName = 'recons_tensor_N100M20sp2sig2ns0.05_loss_14-Oct-2018New.mat';
load(fullfile(dFolder,dName))

% load plastic reconstruction, first layer is "learned"
temp = load(fullfile(dFolder,'loss_plastic_100M20noise0.05_sp0.55_L2.mat'));
plsLoss = [mean(temp.testLost),std(temp.testLost)];

% load reconstruction for hardwired sensitivity matrix
entrData = 'entropyDistr_N100M20_sp2_sig2_ns0.05_03-Oct-2018.mat';
load(fullfile(dFolder,entrData))
% parameters
N = 100;
M = 20;
spar = 2;
sig = 2;
noise = 0.05;

% allSp = 1:-0.05:0;
inx  = 1:1:18;

figure
hold on
% errorbar(allSp(inx)',nanmean(log10(cmbData(inx,:)),2),nanstd(log10(cmbData(inx,:)),0,2),...
%         'o-','MarkerSize',10,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:),'LineWidth',1.5,'CapSize',0)
errorbar(1-allSp(inx)',nanmean(log10(allTestLoss(inx,:)),2),nanstd(log10(allTestLoss(inx,:)),0,2),...
        'o-','MarkerSize',10,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:),'LineWidth',1.5,'CapSize',0)
plot([0;1],log10([plsLoss(1);plsLoss(1)]),'k--','LineWidth',1.5)
% shadow the suboptimal entropy region
ah = gca;
% entRange = range(summH(:,1));
thd = 0.95;  % range of the entropy
TOL = 1e-3;
cs = spaps((0.05:0.05:1)',summH(:,1),TOL,1,3);
d1 = fnder(cs);
minSp = fnzeros(d1,[0.1,0.99]);
maxEntr = fnval(cs, minSp(1,1));
minEntr = summH(1,1);
refEntr = minEntr + (maxEntr - minEntr)*thd;
subOptRange = fnzeros(fncmb(cs,'-',refEntr),[0.1,0.99]);

% add two vertical lines  show the optimal entropy
yR = ah.YLim;
plot(ah,[subOptRange(1,1);subOptRange(1,1)],yR','--','LineWidth',1.5)
plot(ah,[subOptRange(1,2);subOptRange(1,2)],yR','--','LineWidth',1.5)
box on
set(gca,'YLim',[-0.65,0.15])
xlabel('$\rho_w$','interpreter','latex')
ylabel('$\log_{10}$(error)','interpreter','latex')

% add the second y axis to show the differential entropy
yyaxis right
ylabel('differential entropy','Color',Gr(8,:))
errorbar(allSp',summH(:,1),summH(:,2),'d-','MarkerSize',8,'Color',Gr(8,:),...
    'MarkerFaceColor',Gr(8,:),'LineWidth',1,'CapSize',0)


% add extral three lines to show sp = 0.1, 0.55,
plot(ah,[0.1;0.1],yR','--','LineWidth',1.5)
plot(ah,[0.6;0.6],yR','--','LineWidth',1.5)
plot(ah,[0.95;0.95],yR','--','LineWidth',1.5)
hold off


% figNamePref = ['fig5_recon_tf',num2str(N),'M',num2str(M),...
%     'sp',num2str(spar),'_nSig',num2str(noise),'_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


%% plot figure 5c
% compare reconstruction performance for three different rho_w      
dNames = {'loss_hardWire_100M20noise0.05_sp0.9_L2.mat',...
            'loss_hardWire_100M20noise0.05_sp0.4_L2.mat',...
            'loss_hardWire_100M20noise0.05_sp0.05_L2.mat'};

% parameters
N = 100; M = 20; sp = 2; sig = 2; H = 100; L=2;
noise = 0.05;
allSp = [0.1, 0.6, 0.95];
        

% =====================================================
% comparison of reconstructed odor vector with input,
% plot a scatter in a row
% =====================================================
figure
figureSize = [0 0 12 4.3];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% tight plot
[ha, pos] = tight_subplot(1,3,[.05 .03],[.25 .08],[.1 .03]);

for i0 = 1:3
    load(fullfile(dFolder,dNames{i0}))
    
    axes(ha(i0));
    temp1 = original;
    temp2 = prediction;
    nonZeorX0 = zeros(size(temp1,1),sp);
    nonZeorXh = zeros(size(temp2,1),sp);
    for j0 = 1:size(temp1,1)
        nonZeorX0(j0,:) = sort(temp1(j0,temp1(j0,:) ~=0));
        nonZeorXh(j0,:) = sort(temp2(j0,temp2(j0,:) ~=0));
    end
    scatter(nonZeorX0(:,1),nonZeorXh(:,1),10,Bu(7,:),'filled')
    hold on
    scatter(nonZeorX0(:,2),nonZeorXh(:,2),10,Bu(7,:),'filled')
    lg = legend(['Error = ',num2str(round(testLost,2))],'Location','northwest');
    set(lg,'FontSize',16)
    legend boxoff
    xlim([-5.5,7])
    ylim([-5.5,7])
    xl = ha(i0).XLim;
    plot(xl',xl','--','LineWidth',1.5,'Color',Gr(8,:))
    hold off
    box on
    xlabel('$\ln(c)$','interpreter','latex')
    if i0==1
        ylabel('$\ln(\hat{c})$','interpreter','latex')
    else
        set(gca,'YTick',[])
    end
end
% figNamePref = ['fig5_recons_tf_N',num2str(N),'M',num2str(M),...
%     'sp',num2str(sp),'_nSig',num2str(noise),'_scatter_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


%% plot figure S6
% =============================================
% only compare reconstruction heatmap
% =============================================
% set graph size
figure
figureSize = [0 0 13 6];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% tight plot
[ha, pos] = tight_subplot(2,3,[.05 .08],[.15 .02],[.1 .08]);
inx = 501:600;  %select partial data

for i0 = 1:3
    load(fullfile(dFolder,dNames{i0}))
    
    axes(ha(i0))
    temp1 = original(inx,:); temp1(temp1 ==0) = -100;
    imagesc(temp1,[-6,6])
    set(gca,'XTick',[])
    ylabel('trial')

    axes(ha(i0+3))
    temp2 = prediction(inx,:);temp2(temp2 == 0) = -100;
    imagesc(temp2,[-6,6])
    ylabel('trial')
    xlabel('odorant index')

end
% figNamePref = ['figS6_recons_tf_N',num2str(N),'M',num2str(M),...
%     'sp',num2str(sp),'_nSig',num2str(noise),'_heatmap_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])