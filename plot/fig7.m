% This program plot figure 7, optimal sensitivity matrix when both
% excitation and inhibtion are considered

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

%%
% load the summary data
dFolder = '../data/inhi_fig7';
outFolder = '../figures';

%% load data
dFile = 'gcmi_inhi_N50M10sp2_Acti_diffBasal__13-Oct-2018.mat';
load(fullfile(dFolder,dFile))

N = 50;  %in the main figure only plot one N
% inx = find(nRecp ==N);

% we also need some simulation on when r0 = 0;
load(fullfile(dFolder,'int_N50_R10_S2_sig2_2018-07-17.mat'))
f0 = mean(-allfmin);
stdf0= std(-allfmin);

%% figure 7e inhibitory fraction
% ==========================================
% plot how inhibitory fraction change with r0
% ==========================================
figure
hold on
errorbar(R0',meanRatio,stdRatio,'o-','MarkerSize',12,'MarkerFaceColor',Bu(9,:),...
    'MarkerEdgeColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
plot([0;1],[0;1],'k--','LineWidth',2)
hold off
box on
xlabel('$r_0$','Interpreter','latex','FontSize',28)
ylabel('inhibitory fraction','FontSize',28)
set(gca,'XTick',0:0.2:1,'FontSize',24,'LineWidth',1.5,'XScale','linear')
% figNamePref = ['fig7_gcmi_bothExciInhi_ratio_basal_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% figure 7d 
% ========================================================================
% plot how differential entorpy change with r_0
% ========================================================================
figure
hold on
errorbar(R0',meanFmin,stdFmin,'o-','MarkerSize',12,'MarkerFaceColor',Bu(9,:),...
    'MarkerEdgeColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
hold off
box on
xlabel('$r_0$','Interpreter','latex','FontSize',28)
ylabel('differential entropy','FontSize',28)
set(gca,'XTick',0:0.2:1,'FontSize',24,'LineWidth',1.5,'XScale','linear')
% figNamePref = ['gcmi_bothExciInhi_fmin_basal_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% Merge entropy and inhibitory fraction

% set the figure size
figure
figureSize = [0 0 8.8 9];
set(gcf,'Units','centimeters','Position',figureSize,...
'PaperPositionMode','auto');
[ha, pos] = tight_subplot(2,1,[0.02 0],[.15 .01],[.2 .03]);


axes(ha(1))
errorbar(R0',meanRatio,stdRatio,'o-','MarkerSize',8,'MarkerFaceColor',Bu(9,:),...
    'MarkerEdgeColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',1.5,'CapSize',0)
hold on
plot([0;1],[0;1],'k--','LineWidth',1)
xlim([0,1])
ylim([0,1])
hold off
box on
set(gca,'YTick',0:0.5:1,'XTick',[],'LineWidth',0.5,'FontSize',16);
% xlabel('$r_0$','Interpreter','latex','FontSize',28)
ylabel('$\rho_i$','Interpreter','latex','FontSize',20)


axes(ha(2))
errorbar(R0',meanFmin,stdFmin,'o-','MarkerSize',8,'MarkerFaceColor',Bu(9,:),...
    'MarkerEdgeColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',1.5,'CapSize',0)
box on
ylim([-0.85,0])
xlabel('$r_0$','Interpreter','latex','FontSize',20)
ylabel('$I$','Interpreter','latex','FontSize',20)
set(gca,'XTick',0:0.2:1,'FontSize',16,'LineWidth',0.5,'XScale','linear')

% figNamePref = ['fig7_merge_Entropy_InhiFrac',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


%% figure 7c histogram
% ========================================================================
% plot the example histogram
% ========================================================================
% exampFile = 'int_N50_R10_S2_sig2_alp0.18_2018-06-22.mat';
% load([dFolder,filesep,'gcmi_inhi',filesep,exampFile])

selFolder = [dFolder,filesep,'N50M10S2sig2_1013'];
exampFile = fullfile(selFolder,'gcmiInhi_N50_R10_S2_sig2_alp0.18_frac_2018-10-13.mat');
load(exampFile)

% essential parameters
nOdor = 50;
nRecp = 10;
sig = 2;
sp = 2;
r0 = 0.18;

inx = 1;
w = allMat(:,1);
we = log(w(allSign(:,1)>0));
wi = log(w(allSign(:,1)<0));

% histogram
figure
figureSize = [0 0 9 9];
set(gcf,'Units','centimeters','Position',figureSize,...
'PaperPositionMode','auto');
[ha, pos] = tight_subplot(2,1,[0 0],[.18 .01],[.23 .03]);

axes(ha(1))
h1 = histfit(we,15);
xlim([-4,3.5])
ylim([0,100])
set(gca,'XTick',[],'YTick',0:50:100,'FontSize',16);


axes(ha(2))
h2 =  histfit(wi,15);
xlim([-4,3.5])
ylim([0,120])
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('count','FontSize',16)
set(gca,'XTick',-4:2:2,'YTick',0:50:100,'FontSize',16);

% set the color
h1(1).FaceColor = Bu(5,:);%h1(1).EdgeColor = Gr(5,:);
h2(1).FaceColor = RdBu(4,:); %h2(1).EdgeColor = Gr(5,:);
h1(2).Color = Bu(10,:);
h2(2).Color =RdBu(2,:);
h1(2).LineWidth = 3;h2(2).LineWidth = 3;

% figNamePref = ['fig7_ExciInhi_example_histo_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

%% figure 7b  heatmap
% ========================================================================
% plot the  heatmap
% ========================================================================
% allBlue = brewermap(128,'Blues');  % here we make the color following "common rule"
% allRd = brewermap(128,'Reds');
allBlue = brewermap(128,'Reds');  % here we make the color following "common rule"
allRd = brewermap(128,'Blues');

% w = reshape(w,[nRecp,nOdor]);
sign = reshape(allSign(:,1),[nRecp,nOdor]);

[inx_Y,inx_X] = find(sign >0);
[inx_Yinhi,inx_Xinhi] = find(sign <0);
% symbolExci = zeros(nRecp,nOdor);

%map value to color
LB = -2;
UB = 3;
we(we < LB) = LB;
we(we > UB) = UB;

LB_inhi = -3;
UB_inhi = 2;
wi(wi < LB_inhi) = LB_inhi;
wi(wi > UB_inhi) = UB_inhi;
% colorExci = zeros(length(we),3);
colorExci = allBlue(round((we - LB)/(UB-LB)*127)+1,:);
colorInhi = allRd(round((wi - LB_inhi)/(UB_inhi- LB_inhi)*127)+1,:);

maxMarker = 20;
allSizeExci = (we - LB)/(UB-LB)*(maxMarker - 1)+1;
allSizeInhi = (wi - LB_inhi)/(UB_inhi-LB_inhi)*(maxMarker - 1)+1;

figure
figureSize = [0 0 15 4];
set(gcf,'Units','centimeters','Position',figureSize,...
'PaperPositionMode','auto');
hold on
for i0 = 1:length(we)
    plot(inx_X(i0)-0.5,inx_Y(i0)-0.5,'o','MarkerSize',allSizeExci(i0),'MarkerFaceColor',...
        colorExci(i0,:),'MarkerEdgeColor',colorExci(i0,:),'LineWidth',0.5)
end

for i0 = 1:length(wi)
     plot(inx_Xinhi(i0)-0.5,inx_Yinhi(i0)-0.5,'o','MarkerSize',allSizeInhi(i0),'MarkerFaceColor',...
        colorInhi(i0,:),'MarkerEdgeColor',colorInhi(i0,:),'LineWidth',0.5)
end
hold off
set(gca,'XTick',10:10:50,'YTick',[1,5,10])
% grid on
% set(gca,'XTick',[])
% daspect([1 1 1])
% figNamePref = ['fig7_ExciInhi_examp_heatmap_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


%% figure 7f distribution of excitatory and inhibitory elements
% prepare the data from Hallem and Carlson's 2006 cell paper
file1 = '../data/Carlson2006TableS1.xlsx';
[NUM1,~,~]=xlsread(file1);
file2 = '../data/CarslonORNdigit.xlsx';
[NUM2,~,~]=xlsread(file2);  %digital matrix defines strong excitation and inhibition

r0 = max(abs(min(NUM1,[],1)),[],1);  %adding 0.1 for stability

% here, spontaneous firing rate is defined as the maximum inhibitory
% ineractions
% spRate = abs(min(NUM1,[],1));
spMat = ones(110,1)*r0;
adjM = NUM1(1:110,:) + r0;  % adjust after considering the basal activity

% default hillCoef
h = 1;  % we can tune this to see how the reuslts change

% we assumed the same maxium spiking rate
Rmax = max(adjM(:)) + 1;  %to avoid numerical unstability, so we add 1

% alpha parameter
alp = Rmax./r0 - 1;  % each receptor has one alpha
alpMat = ones(110,1)*alp;  %this matrix is used for future use

%strong excitation based on the digital matrix
% allM = adjM(1:end-1,:);
allM = adjM;
strongExi = allM(NUM2 > 0);
strongW = (alpMat(NUM2 > 0)./(Rmax./allM(NUM2 > 0) - 1) - 1).^(1/h);

% consider all excitation
allInx = allM > spMat;
allExci = allM(allInx);
exciW = (alpMat(allInx)./(Rmax./allM(allInx) - 1) - 1).^(1/h);

%strong inhibition
strongInhi = max(allM(NUM2 < 0),1);  % for stability

strongInhiW = ((Rmax./strongInhi-1)./alpMat(NUM2 < 0) - 1).^(1/h);

% consider all excitation
allInx = allM < spMat;
allInhi = max(allM(allInx),1);   % for stability

inhiW = ((Rmax./allInhi-1)./alpMat(allInx) - 1).^(1/h);

% ============================================================
% merged tight histogram of all excitation and inhibition and lognormal fit
% ============================================================
% histogram
figureSize = [0 0 9 9];
set(gcf,'Units','centimeters','Position',figureSize,...
'PaperPositionMode','auto');
% [ha, pos] = tight_subplot(2,1,[0 0],[.2 .01],[.2 .03]);
[ha, pos] = tight_subplot(2,1,[0 0],[.18 .01],[.23 .03]);

axes(ha(1))
h1 = histfit(log(exciW),20);
xlim([-5,8])
ylim([0,350])
ylabel('count','FontSize',16)
set(gca,'XTick',[],'YTick',0:150:350,'FontSize',16);

axes(ha(2))
h2 =  histfit(log(inhiW(inhiW>0)),20);
xlim([-5,8])
ylim([0,350])
xlabel('$\ln(w)$','Interpreter','latex','FontSize',16)
ylabel('count','FontSize',16)
set(gca,'XTick',-5:5:5,'YTick',0:150:300,'FontSize',16);

% set the color
h1(1).FaceColor = Bu(5,:);h2(1).FaceColor = RdBu(4,:);
h1(2).Color = Bu(10,:);h2(2).Color =RdBu(2,:);
h1(2).LineWidth = 3;h2(2).LineWidth = 3;

% figNamePref = ['fig7_ExciInhi_expermient_histo_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


%% figure add:  histogram of mosquito data
% prepare the data from Hallem and Carlson's 2006 cell paper
file1 = '../data/CareyNature2010.xlsx';
[NUM1,~,~] = xlsread(file1);
file2 = '../data/CareyNature2010Digit.xlsx';
[NUM2,~,~] = xlsread(file2);  %digital matrix defines strong excitation and inhibition

r0 = max([NUM1(end,:);abs(min(NUM1,[],1))],[],1);  %adding 0.1 for stability
spMat = ones(110,1)*r0;
adjM = NUM1 + r0;  % adjust after considering the basal activity

% default hillCoef
h = 1;  % we can tune this to see how the reuslts change

% we assumed the same maxium spiking rate
Rmax = max(adjM(:)) + 0.1;  %to avoid numerical unstability, so we add 1

% alpha parameter
alp = Rmax./r0 - 1;  % each receptor has one alpha
alpMat = ones(110,1)*alp;  %this matrix is used for future use

%strong excitation based on the digital matrix
% allM = adjM(1:end-1,:);
allM = adjM(1:110,:);
strongExi = allM(NUM2 > 0);
strongW = (alpMat(NUM2 > 0)./(Rmax./allM(NUM2 > 0) - 1) - 1).^(1/h);

% consider all excitation
allInx = allM > spMat;
allExci = allM(allInx);
exciW = (alpMat(allInx)./(Rmax./allM(allInx) - 1) - 1).^(1/h);

%strong inhibition
strongInhi = max(allM(NUM2 < 0),0.1);  % for stability

strongInhiW = ((Rmax./strongInhi-1)./alpMat(NUM2 < 0) - 1).^(1/h);

% consider all excitation
allInx = allM < spMat;
allInhi = max(allM(allInx),1);   % for stability

inhiW = ((Rmax./allInhi-1)./alpMat(allInx) - 1).^(1/h);

% ============================================================
% merged tight histogram of all excitation and inhibition and lognormal fit
% ============================================================
% histogram all the data are used
figureSize = [0 0 9 9];
set(gcf,'Units','centimeters','Position',figureSize,...
'PaperPositionMode','auto');
% [ha, pos] = tight_subplot(2,1,[0 0],[.2 .01],[.2 .03]);
[ha, pos] = tight_subplot(2,1,[0 0],[.18 .01],[.23 .03]);

axes(ha(1))
h1 = histfit(log(exciW),20);
xlim([-7,7])
ylim([0,500])
ylabel('count','FontSize',16)
set(gca,'XTick',[],'YTick',0:150:500,'FontSize',16);

axes(ha(2))
h2 =  histfit(log(inhiW),20);
xlim([-7,7])
ylim([0,500])
xlabel('$\ln(w)$','Interpreter','latex','FontSize',16)
ylabel('count','FontSize',16)
set(gca,'XTick',-5:5:7,'YTick',0:150:500,'FontSize',16);

% set the color
h1(1).FaceColor = Bu(5,:);h2(1).FaceColor = RdBu(4,:);
h1(2).Color = Bu(10,:);h2(2).Color =RdBu(2,:);
h1(2).LineWidth = 3;h2(2).LineWidth = 3;

% figNamePref = ['fig7_ExciInhi_mosquito_histo_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])



% ----------------------------------------------------------------
% histogram, only use the strong excitatory and inhibitory interactions
% --------------------------------------------------------
figureSize = [0 0 11 9];
set(gcf,'Units','centimeters','Position',figureSize,...
'PaperPositionMode','auto');
% [ha, pos] = tight_subplot(2,1,[0 0],[.2 .01],[.2 .03]);
[ha, pos] = tight_subplot(2,1,[0 0],[.18 .01],[.23 .03]);

axes(ha(1))
h1 = histfit(log(strongW),20);
xlim([-5,7])
ylim([0,150])
ylabel('count','FontSize',16)
set(gca,'XTick',[],'YTick',0:50:150,'FontSize',16);

axes(ha(2))
h2 =  histfit(log(strongInhiW(strongInhiW>0)),20);
xlim([-5,7])
ylim([0,150])
xlabel('$\ln(w)$','Interpreter','latex','FontSize',16)
ylabel('count','FontSize',16)
set(gca,'XTick',-5:5:5,'YTick',0:50:150,'FontSize',16);

% set the color
h1(1).FaceColor = Bu(5,:);h2(1).FaceColor = RdBu(4,:);
h1(2).Color = Bu(10,:);h2(2).Color =RdBu(2,:);
h1(2).LineWidth = 3;h2(2).LineWidth = 3;

% figNamePref = ['fig7_ExciInhi_mosquito_Strong_histo_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


%% figure added P. Siccifolium data from Missbach2014 eLife
% notice that the rows are ORs and columns are different odorants
dataFile = "../data/ Psiccifolium_ORNresp_Missbach2014.xlsx";
[NUM3,~,~] = xlsread(dataFile);
[ROW,COL] = size(NUM3(:,1:(end-1)));

% the last column is the baseline firing rate
sp0 = NUM3(:,end);
estSp0 = abs(min(NUM3,[],2)); % estimated spontaneous firing rate
adjSp0 = max([sp0,estSp0],[],2); % adjusted spontaneous firing rate
rawFR = NUM3(:,1:(end-1));   % raw firing rate
spMat = adjSp0*ones(1,COL);  % matrix of spontaneous firing rate
adjM = NUM3(:,1:(end-1)) + spMat; % adding 0.1 to avoid infinity latter


% default hillCoef
h = 1;  % we can tune this to see how the reuslts change

% we assumed the same maxium spiking rate
Rmax = max(adjM(:)) + 1;  %to avoid numerical unstability, so we add 1

% alpha parameter
alp = Rmax./adjSp0 - 1;  % each receptor has one alpha
alpMat = alp*ones(1,COL);  %this matrix is used for future use


% consider all excitation
allInx = rawFR > 0;
allExci = adjM(allInx);
exciW = (alpMat(allInx)./(Rmax./adjM(allInx) - 1) - 1).^(1/h);


% consider all inhibition
allInx = rawFR < 0;
allInhi = max(adjM(allInx),1);   % for stability

inhiW = ((Rmax./allInhi-1)./alpMat(allInx) - 1).^(1/h);

% ============================================================
% merged tight histogram of all excitation and inhibition and lognormal fit
% ============================================================
% histogram all the data are used
figureSize = [0 0 9 9];
set(gcf,'Units','centimeters','Position',figureSize,...
'PaperPositionMode','auto');
% [ha, pos] = tight_subplot(2,1,[0 0],[.2 .01],[.2 .03]);
[ha, pos] = tight_subplot(2,1,[0 0],[.18 .01],[.23 .03]);

axes(ha(1))
h1 = histfit(log(exciW),20);
xlim([-6,8])
ylim([0,200])
ylabel('count','FontSize',16)
set(gca,'XTick',[],'YTick',0:100:200,'FontSize',16);

axes(ha(2))
h2 =  histfit(log(inhiW),20);
xlim([-6,8])
ylim([0,200])
xlabel('$\ln(w)$','Interpreter','latex','FontSize',16)
ylabel('count','FontSize',16)
set(gca,'XTick',-5:5:5,'YTick',0:100:200,'FontSize',16);

% set the color
h1(1).FaceColor = Bu(5,:);h2(1).FaceColor = RdBu(4,:);
h1(2).Color = Bu(10,:);h2(2).Color =RdBu(2,:);
h1(2).LineWidth = 3;h2(2).LineWidth = 3;

% figNamePref = ['fig7_ExciInhi_siccifolium_histo_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-dpdf',[outFolder,filesep,figNamePref,'.pdf'])
