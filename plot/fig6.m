% This program plot figure 6, the classification performance

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
dFolder = '../decoding/data/classFig6';
outFolder = '../figures';


classData = 'LDA_N100M10sp3_ns0.05_nType2_P100_H500_dS0.1_27-Dec-2018.mat';
entrData = 'entropyDistr_N100M10_sp3_sig2_ns0.05_26-Sep-2018.mat';
load(fullfile(dFolder,classData))
% load('/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/decoding/LDA_N100M10_diffCluster/LDA_N100M10sp3_ns0.05_nType2_P100_H500_dS0.1_27-Dec-2018.mat')
load(fullfile(dFolder,entrData))

% parameters
N = 100;
M = 10;
sp = 3;
sig = 2;
noiseSig = 0.05;

% allSp = 0.05:0.05:1;
allSp = 0.1:0.05:1;

% selected data
ix = 4;   % corresponds to 100 patterns
% rawData = squeeze(summData.errorRate(:,:,ix));
rawData = allTestError;

fh = figure;
hold on
Bu = brewermap(11,'Blues');    % blues
errorbar(allSp',mean(rawData,2),std(rawData,0,2),...
        'o-','MarkerSize',10,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:),'LineWidth',1.5,'CapSize',0)
% lg = legend('N = 100,M=10,\sigma_n = 0.05, P = 100');
lg = legend('classification');
legend boxoff
set(lg, 'FontSize',16)

xlabel('$\rho_w$','interpreter','latex')
ylabel('classification error')

% add the second y axis to show the differential entropy
yyaxis right
ylabel('differential entropy','Color',Gr(8,:))
errorbar(allSp(1:19)',summH(:,1),summH(:,2),'d-','MarkerSize',8,'Color',Gr(8,:),...
    'MarkerFaceColor',Gr(8,:),'LineWidth',1,'CapSize',0)

ah = gca;
y2h = ah.YAxis(2);
y2h.Color = 'k';
box on


% shadow the suboptimal entropy region
entRange = range(summH(:,1));
thd = 0.95;  % range of the entropy
TOL = 1e-3;
cs = spaps((0.05:0.05:0.95)',summH(:,1),TOL,1,3);
% figure
% hold on
% fnplt(cs)
% scatter(allSp(1:19)',summH(:,1),'o')
% hold off
% first order derivative
d1 = fnder(cs);
minSp = fnzeros(d1,[0.1,0.99]);
maxEntr = fnval(cs, minSp(1,1));
minEntr = summH(1,1);
refEntr = minEntr + (maxEntr - minEntr)*thd;
subOptRange = fnzeros(fncmb(cs,'-',refEntr),[0.1,0.99]);

% add two vertical lines 
yR = ah.YLim;
plot(ah,[subOptRange(1,1);subOptRange(1,1)],yR','--','LineWidth',1.5)
plot(ah,[subOptRange(1,2);subOptRange(1,2)],yR','--','LineWidth',1.5)
hold off

% figNamePref = ['fig6_class_LDA_N',num2str(N),'M',num2str(M),...
%     'sp',num2str(sp),'_nSig',num2str(noiseSig),'_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])