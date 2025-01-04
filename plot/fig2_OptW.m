% This program generate figure 2 and related supplemental figures

close all
clear


% define the colors
Bu = [143,195,237]/256;
Or = [243,152,0]/256;

RdBu = brewermap(11,'RdBu');   % red and blue
% Bu = brewermap(11,'Blues');    % blues
Gr = brewermap(11,'Greys');    % greys
seqColor = brewermap(21,'YlGnBu');

lBu = [96,166,223]/255; %light blue
dpBu = [63,114,183]/255; % deep blue
dkBu = [50,78,147]/255;   %dark blue
Or = [220,150,71]/255;  % orange
brickRd = [201,69,89]/255;  %brick red
green = [107,169,81]/255;  %green
purple = [113,49,119]/255;  % purple

% data folder and folder to save the figures
dFolder = '../data/GcmiN100R30-10-03';
saveFolder = '../figures';

%% define some parameters

% the parameters of the system
N = 100;
M = 30;
sp = 2;
sig =2;

repeats = 100;  %number of permuation


%%  load and prepare the data
% load the data used
fName = 'N100_R30_S2_sig2_2018-10-03.mat';
load(fullfile(dFolder,fName));

% select one matrix and show its sparsity
ix0 = 1;   % slect one matrix, default the first
w = reshape(allMat(:,ix0),[M,N]);  %the example optimal sensitivity matrix

% truncates w to be in the range [-4*sigma, 4*sigma], for better look of
% the heatmap
trucate = true; % default
if trucate
   trc = 'truc'; %trucation or not
   w(w < -4*sig) = -4*sig;
   w(w > 4*sig) = 4*sig;
else
    trc = 'NOtruc'; %trucation or not
end

%%  Fig2 a
% ====================================================
% plot the heatmap of an optimal w
% ====================================================
figure
figureSize = [0 0 11 4];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

set(gcf,'renderer','Painters')
imagesc(w,[-4,4])
set(gca,'TickLength',[0,0])

colormap(jet);
c = colorbar; 
% c.TickLabels{1} = 'NaN'; 
% c.Label.String = '-log10(EC50)';
% c.FontSize = 16;
xlabel('Odorant')
ylabel('Receptor')
% prefix = ['fig2_gcmi_eplW_N',num2str(N),'M',num2str(M),'sig',num2str(sig),....
%     'sp',num2str(sp),'_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%% Fig2 b
% ====================================================
% plot the histogram of optimal w
% ====================================================
% need to reload the matrix, here we need all the elements
w = reshape(allMat(:,ix0),[M,N]);

Wid = 8;
Hei =  Wid*0.85;
figure
set(gcf,'Units','centimeters','Position',[0,0,Wid,Hei])
% figure
hold on
set(gcf,'renderer','Painters')
h1 = histogram(w(:),60,'Normalization','probability');
h1.FaceColor = lBu; h1.FaceAlpha = 0.4; h1.EdgeColor = 'none';
stairs([h1.BinEdges(1),h1.BinEdges,h1.BinEdges(end)],...
    [0,h1.Values,h1.Values(end),0],'Color',dpBu,'LineWidth',1)

box on
hold off
set(gca,'Layer','top')
xlim([-120,10])
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('probability')
% prefix = ['fig2_examp_AllW_N',num2str(N),'M',num2str(M),'sig',num2str(sig),....
%     'sp',num2str(sp),'_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% enlarged, only plot the active elements
figure
set(gcf,'renderer','Painters')
set(gcf,'Units','centimeters','Position',[0,0,Wid,Hei])

hold on
h2 = histogram(w(abs(w)<5),20,'Normalization','pdf');
h2.FaceColor = lBu; h2.FaceAlpha = 0.4; h2.EdgeColor = 'none';
stairs([h2.BinEdges(1),h2.BinEdges,h2.BinEdges(end)],...
    [0,h2.Values,h2.Values(end),0],'Color',dpBu,'LineWidth',1)

% get the normal fit parameter
pd = fitdist(w(abs(w)<6),'normal');
X = -6:0.05:6;
Y = normpdf(X,pd.mean,pd.sigma);
plot(X,Y,'Color',Or,'LineWidth',3)
hold off
box on
legend('active w','Gassian fit')
legend boxoff
set(gca,'XLim',[-6,6])
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
% prefix = ['fig2_examp_ActiveW_N',num2str(N),'M',num2str(M),'sig',num2str(sig),....
%     'sp',num2str(sp),'_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%% Fig 2 c
% Comparison of the correlation
% at this part we need to trucate W before calculating the rank correlation
w = reshape(allMat(:,ix0),[M,N]);
trucate = true;
if trucate
   trc = 'truc'; %trucation or not
   w(w < -4*sig) = -4*sig;
   w(w > 4*sig) = 4*sig;
else
    trc = 'NOtruc'; %trucation or not
end

%==================================================================
%correlation coefficients after shuffling, row-wise
%=================================================================
corrType = 'Kendall';  % or Spearman

% first, compare original and shuffuled correlation
L = size(w,1)*size(w,2);

t1 = corr(w','type',corrType);
y1 = triu(t1,1);
rowMean = mean(y1(y1~=0));
rowStd = std(y1(y1~=0));

newW = reshape(w(randperm(L)),size(w,1),size(w,2));
t2 = corr(newW','type',corrType);
y2 = triu(t2,1);

[bc1, edg1] = histcounts(y1(y1~=0),15);
[bc2,edg2] = histcounts(y2(y2~=0),15);

figure
hold on
Y1 = [bc1;bc1]/sum(bc1)/(edg1(2)- edg1(1)); %normalized, probability
ar1 = area(sort(edg1([1:end-1 2:end])),Y1(1:end));
ar1.FaceAlpha = 0.3;
ar1.FaceColor = lBu;
ar1.EdgeColor = dpBu;
ar1.LineWidth = 2;

Y2 =  [bc2;bc2]/sum(bc2)/(edg2(2)- edg2(1));
ar2 = area(sort(edg2([1:end-1 2:end])),Y2(1:end));
ar2.FaceAlpha = 0.4;
ar2.FaceColor = Gr(4,:);
ar2.EdgeColor = Gr(8,:);
ar2.LineWidth = 2;
set(gca,'Layer','top')
hold off
xlim([-1,1])
lg = legend('original','shuffled');
set(lg,'FontSize',16)
legend boxoff
box on
xlabel('$\tau$','Interpreter','latex')
ylabel('pdf')
% xlim([-1,1])
set(gca,'XTick',-1:0.5:1)

% prefix = ['fig2_exampCorrComp_RowShuff_',corrType,'_',trc,'_N',...
%     num2str(N),'M',num2str(M),'sig',num2str(sig),....
%     'sp',num2str(sp),'_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])


% statistics after shuffling
repeats = 100;
% L = size(w,1)*size(w,2);
allMeanCorr = zeros(repeats,1);
allStdCorr = zeros(repeats,1);

for i0 = 1:repeats
    newW = reshape(w(randperm(L)),size(w,1),size(w,2));
    temp = corr(newW','type',corrType);
    y3 = triu(temp,1);
    allMeanCorr(i0) = mean(y3(y3~=0));
    allStdCorr(i0) = std(y3(y3~=0));
end
figure
hold on
scatter(allMeanCorr,allStdCorr,20,[0.5,0.5,0.5],'filled')
xlim([-0.05,0.05])
ylim([0.04,0.11])
set(gca,'XTick',-0.05:0.05:0.05)
scatter(rowMean,rowStd,80,dkBu,'filled')
hold off
box on
lg = legend('shuffled','original');
xlabel('$\langle \tau \rangle$','Interpreter','latex')
ylabel('$\sigma_{\tau}$','Interpreter','latex')
% prefix = ['fig2_exampCorrComp_RowShuffScatter_',corrType,'_',trc,'_N',...
%     num2str(N),'M',num2str(M),'sig',num2str(sig),....
%     'sp',num2str(sp),'_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-painters','-dpdf',[saveFolder,filesep,prefix,'.pdf'])

%% Fig 2d and Fig SI, Guangwei Si's data
% plot the fly larvae and mouse OR sensitivity data

fileLavae = '../data/flyLarva_Guangwei_LogEC50_New.xlsx';
[NUM3,TEXT,~]=xlsread(fileLavae);
allW = NUM3(abs(NUM3) > 0);
newM = -NUM3;

receptorLabel = TEXT(1,2:end);
odorLabel = TEXT(2:end,1);

sr = '[\w-\s,]{3,30}';
for i0 = 1:length(receptorLabel)
    receptorLabel{i0} = char(regexp(receptorLabel{i0},sr,'match'));
end

for i0 = 1:length(odorLabel)
    odorLabel{i0} = char(regexp(odorLabel{i0},sr,'match'));
end

% ===============================================================
% Fit the data with a skewed gaussian
% ===============================================================
gaussian = @(x) (1/sqrt((2*pi))*exp(-x.^2/2));
skewStand = @(x,alpha) 2*gaussian(x).*normcdf(alpha*x);
skewPdf =  @(x,xi,omi,alp) 2/omi*gaussian((x-xi)/omi).*normcdf(alp*(x-xi)/omi);

dVec = newM(~isnan(newM));
% ecdf
[Y,X] = ecdf(dVec);

fun = @(x,xdata) normcdf((xdata-x(1))/x(2)) - 2*myOwenT((xdata-x(1))/x(2),x(3));
x0 = [mean(dVec),std(dVec),3];
lb = [0,1,0];
ub = [10,10,10];
% least square fit
optParam = lsqcurvefit(fun,x0,X(2:end),Y(2:end),lb,ub);

xi = optParam(1); omi = optParam(2); alp = optParam(3);

% ================================================
% fit histogram and fit with Gaussian distribution
% ================================================
% 
figure
figureSize = [0 0 5 4.5];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');
hold on
set(gcf,'renderer','Painters')
hh = histogram(abs(-allW),15,'Normalization','pdf',...
    'EdgeColor','none','FaceColor',lBu,'FaceAlpha',0.4);

pd = fitdist(abs(allW),'normal');
X = 0:0.05:10;
Y = normpdf(X,pd.mean,pd.sigma);
plot(X,Y,'LineWidth',3,'Color',Or)
Y2 = skewPdf(X,xi,omi,alp);
plot(X,Y2,'LineWidth',3,'Color',RdBu(2,:))
box on
ylim([0,0.5])
hold off
lg = legend('experiment','Gaussian fit','Skewed Gaussian fit');
set(lg,'FontSize',16)
% legend boxoff
xlabel('$\log_{10}(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
set(gca,'Layer', 'top')
% prefix = ['fig2d_GuangweiStrongW_fit_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% ================================
% simple histogram without fitting
% ================================
figure
set(gcf,'Units','centimeters','Position',[0,0,7,7*0.9])
hold on
h1 = histogram(abs(-allW),15,'Normalization','pdf','FaceColor',lBu,'FaceAlpha',...
    0.4,'EdgeColor','none');
stairs([h1.BinEdges(1),h1.BinEdges,h1.BinEdges(end)],...
    [0,h1.Values,h1.Values(end),0],'Color',dpBu,'LineWidth',0.75)
box on
ylim([0,0.4])
xlim([0,10])
hold off
xlabel('$\log_{10}(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
set(gca,'Layer', 'top')
% prefix = ['fig2d_GuangweiStrongW_hist_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%% Fig 2d  and Fig SI, Saito 2009 data
% in this data set contain a panel of 63 odorants, 53 mouse Or and 10 Human
% Ors
dataFile = "../data/Saito2009s3.xls";
[NUM4,~,~] = xlsread(dataFile);

% get rido of humman receptor index
humanInx = [1,5,14,17,18,20,23,34,35,40];
allInx = 1:1:62;
selInx = setdiff(allInx,humanInx);
mouseMat = NUM4(selInx,:);  % mouse-only data
mouseVec = mouseMat(mouseMat ~= 0);
allVec = NUM4(NUM4 ~=0);

% ================================
% statistical test, test if data is Gaussian
% ================================
[H, pValue, SWstatistic] = swtest(allVec, 0.01);

gaussian = @(x) (1/sqrt((2*pi))*exp(-x.^2/2));
skewStand = @(x,alpha) 2*gaussian(x).*normcdf(alpha*x);
skewPdf =  @(x,xi,omi,alp) 2/omi*gaussian((x-xi)/omi).*normcdf(alp*(x-xi)/omi);

dVec = abs(allVec);
% ecdf
[Y,X] = ecdf(dVec);

fun = @(x,xdata) normcdf((xdata-x(1))/x(2)) - 2*myOwenT((xdata-x(1))/x(2),x(3));
x0 = [mean(dVec),std(dVec),3];
lb = [0,1,0];
ub = [10,10,10];
% least square fit
optParam = lsqcurvefit(fun,x0,X(2:end),Y(2:end),lb,ub);

xi = optParam(1); omi = optParam(2); alp = optParam(3);

% ==================================================================
% plot the histogram of mouse only data and add a Gaussian fit
% ========================================================================
figure
set(gcf,'Units','centimeters','Position',[0,0,7,7*0.9])

% set(gcf,'Units','inches','Position',figureSize,...
% 'PaperPositionMode','auto');
hold on
set(gcf,'renderer','Painters')
h1 = histogram(-mouseVec,25,'Normalization','pdf','FaceColor',lBu,'FaceAlpha',...
    0.4,'EdgeColor','none');
stairs([h1.BinEdges(1),h1.BinEdges,h1.BinEdges(end)],...
    [0,h1.Values,h1.Values(end),0],'Color',dpBu,'LineWidth',1)
% xlim([-2,2])
% fit a lognormal distribution
pd = fitdist(-mouseVec,'normal');
X = 2:0.05:8;
Y = normpdf(X,pd.mean,pd.sigma);
plot(X,Y,'LineWidth',4,'Color',Or)

Y2 = skewPdf(X,xi,omi,alp);
plot(X,Y2,'LineWidth',4,'Color',RdBu(2,:))
hold off
lg = legend('Experiment','Normal Fit','Skewed Fit','Location','northwest');
set(lg,'FontSize',16)
legend boxoff
box on

% legend('w','Gassian fit')
% legend boxoff
xlabel('$\log_{10}(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
% prefix = ['fig2d_Saito_mouse_fit_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%% Merge fly larvae and mouse data together


figureSize = [0 0 13 5.5];
set(gcf,'Units','centimeters','Position',figureSize,...
'PaperPositionMode','auto');
set(gcf, 'Renderer', 'painters');

% [ha, pos] = tight_subplot(2,1,[0 0],[.2 .01],[.2 .03]);
[ha, pos] = tight_subplot(1,2,[0 0],[.25 .01],[.13 .03]);

axes(ha(1))
h1 = histogram(abs(-allW),15,'Normalization','pdf','FaceColor',lBu,'FaceAlpha',...
    0.4,'EdgeColor','none');
hold on
stairs([h1.BinEdges(1),h1.BinEdges,h1.BinEdges(end)],...
    [0,h1.Values,h1.Values(end),0],'Color',dpBu,'LineWidth',0.75)
box on
ylim([0,0.65])
xlim([0,10.5])
hold off
xlabel('$\log_{10}(w)$','Interpreter','latex','FontSize',16)
ylabel('$p_s(w)$','Interpreter','latex','FontSize',16)
set(gca,'Layer', 'top')


axes(ha(2))
h1 = histogram(-mouseVec,25,'Normalization','pdf','FaceColor',lBu,'FaceAlpha',...
    0.4,'EdgeColor','none');
hold on
stairs([h1.BinEdges(1),h1.BinEdges,h1.BinEdges(end)],...
    [0,h1.Values,h1.Values(end),0],'Color',purple,'LineWidth',1)
hold off
xlim([1,8])
ylim([0,0.65])
box on
xlabel('$\log_{10}(w)$','Interpreter','latex','FontSize',16)
% ylabel('pdf','Interpreter','latex','FontSize',16)
set(gca,'YTick',[]);

% prefix = ['fig2_Saito_mouse_histo_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%% Missbach2014 P. Siccifolium data
% notice that the rows are ORs and columns are different odorants
dataFile = "../data/Psiccifolium_ORNresp_Missbach2014.xlsx";
[NUM5,~,~] = xlsread(dataFile);
[ROW,COL] = size(NUM5(:,1:(end-1)));

% the last column is the baseline firing rate
sp0 = NUM5(:,end);
estSp0 = abs(min(NUM5,[],2)); % estimated spontaneous firing rate
adjSp0 = max([sp0,estSp0],[],2); % adjusted spontaneous firing rate
rawFR = NUM5(:,1:(end-1));   % raw firing rate
FRmat = NUM5(:,1:(end-1)) + adjSp0*ones(1,COL) + 0.1; % adding 0.1 to avoid infinity latter

exciEle = FRmat(FRmat>2*adjSp0*ones(1,COL));
rho_e = sum(exciEle)/COL/ROW;   % sparsity of excitatory interaction (strong criterion)

% estimate the sensitivity for strong excitation

%% Fig S1
%==================================================================
%correlation coefficients after shuffling, cloumn-wise
%=================================================================
% first, compare original and shuffuled correlation
t1 = corr(w,'type',corrType);
y1 = triu(t1,1);
rowMean = mean(y1(y1~=0));
rowStd = std(y1(y1~=0));

newW = reshape(w(randperm(L)),size(w,1),size(w,2));
t2 = corr(newW,'type',corrType);
y2 = triu(t2,1);

[bc1, edg1] = histcounts(y1(y1~=0),25);
[bc2,edg2] = histcounts(y2(y2~=0),25);

figure
hold on

Y1 = [bc1;bc1]/sum(bc1)/(edg1(2)- edg1(1)); %normalized, probability
ar1 = area(sort(edg1([1:end-1 2:end])),Y1(1:end));
ar1.FaceAlpha = 0.3;
ar1.FaceColor = lBu;
ar1.EdgeColor = dpBu;
ar1.LineWidth = 2;

% Y2 = [bc2;bc2]*2/size(w,1)/(size(w,1)-1); %normalized, probability
Y2 =  [bc2;bc2]/sum(bc2)/(edg2(2)- edg2(1));
ar2 = area(sort(edg2([1:end-1 2:end])),Y2(1:end));
ar2.FaceAlpha = 0.4;
ar2.FaceColor = Gr(4,:);
ar2.EdgeColor = Gr(8,:);
ar2.LineWidth = 2;
set(gca,'Layer','top')
hold off
xlim([-1,1])
lg = legend('original','shuffled');
set(lg,'FontSize',16)
legend boxoff
box on
xlabel('$\tau$','Interpreter','latex')
ylabel('pdf')
% xlim([-1,1])
set(gca,'XTick',-1:0.5:1)
hold off

% prefix = ['fig2_CorrComp_columnShuff_',corrType,'_',trc,'_N',num2str(N),'M',...
%     num2str(M),'sig',num2str(sig),....
%     'sp',num2str(sp),'_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-dpdf',[saveFolder,filesep,prefix,'.pdf'])

% statistics after shuffling
repeats = 100;
L = size(w,1)*size(w,2);
allMeanCorr = zeros(repeats,1);
allStdCorr = zeros(repeats,1);

for i0 = 1:repeats
    newW = reshape(w(randperm(L)),size(w,1),size(w,2));
    temp = corr(newW,'type',corrType);
    y3 = triu(temp,1);
    allMeanCorr(i0) = mean(y3(y3~=0));
    allStdCorr(i0) = std(y3(y3~=0));
end
figure
hold on
box on
scatter(allMeanCorr,allStdCorr,20,[0.5,0.5,0.5],'filled')
xlim([-0.01,0.01])
ylim([0.12,0.18])
set(gca,'XTick',-0.01:0.01:0.01)
scatter(rowMean,rowStd,80,dkBu,'filled')
hold off
lg = legend('shuffled','original');
xlabel('$\langle \tau \rangle$','Interpreter','latex')
ylabel('$\sigma_{\tau}$','Interpreter','latex')
% prefix = ['fig2_exampCorrComp_ColumnShuffScatter_',corrType,'_',trc,'_N',...
%     num2str(N),'M',num2str(M),'sig',num2str(sig),....
%     'sp',num2str(sp),'_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%%  Fig S2 a
% ========================================================
% Compare the eigenvalues of optimal and permutated matrix
% Row-wise
% ========================================================
repeats = 500;  %number of permuation

% select one matrix and show its sparsity
w = reshape(allMat(:,ix0),[M,N]);

% truncates w to be in the range [-4*sigma, 4*sigma]
trucate = true;   %modified on 12/06/2018
if trucate
   trc = 'truc'; %trucation or not
   w(w < -4*sig) = -4*sig;
   w(w > 4*sig) = 4*sig;
else
    trc = 'NOtruc'; %trucation or not
end

corrType = 'Kendall';  %or Spearman

C1 = corr(w','type',corrType);
eigenOpt = eig(C1);

% random permutation
eigenPerm = zeros(M, repeats);
aveOpt = [mean(eigenOpt),std(eigenOpt)];

for i0 = 1:repeats
    newW = reshape(w(randperm(M*N)),size(w,1),size(w,2));
    temp = corr(newW','type',corrType);
    eigenPerm(:,i0) = eig(temp);
end
avePerm = [mean(eigenPerm,1);std(eigenPerm,0,1)];

figure
set(gcf,'Renderer','painters')
hold on
histogram(avePerm(2,:),'FaceColor',Bu,'Normalization','probability','LineWidth',0.5);
ah = gca;
ylim = ah.YLim;
plot([aveOpt(2); aveOpt(2)],ylim,'r--','LineWidth',2)
box on
set(gca,'XLim',[0.25, 0.5])
legend('shuffle','original')
xlabel('$\sigma_{\lambda}$','Interpreter','latex')
ylabel('probability')
% prefix = ['fig2_CorrComp_Row_eigen_',corrType,'_',trc,'_N',num2str(N),'M',...
%     num2str(M),'sig',num2str(sig),....
%     'sp',num2str(sp),'_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% CDF compare
figure
hold on
for i0 = 1:repeats
    [f,x] = ecdf(eigenPerm(:,i0));
    plot(x,f,'LineWidth',1,'Color',Gr(6,:))
end

[f,x] = ecdf(eigenOpt);
plot(x,f,'LineWidth',2,'Color','red')
hold off
box on
lg = legend('shuffle','original');
legend boxoff
xlabel('$\lambda $','Interpreter','latex')
ylabel('ECDF')
% prefix = ['fig2_CorrComp_Row_CDF_',corrType,'_',trc,'_N',num2str(N),'M',...
%     num2str(M),'sig',num2str(sig),...
%     'sp',num2str(sp),'_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

%% figure S2 b
% ==================================================
% column-wise correlation matrix, eigenvalues comparison
% ==================================================
thd = 1e-5;  % threshold used to determine if a eigenvalue is zero
C1 = corr(w,'type',corrType);
eigenOpt = eig(C1);

% random permutation
eigenPerm = zeros(N, repeats);
aveOpt = [mean(eigenOpt(72:end)),std(eigenOpt(72:end))];

for i0 = 1:repeats
    newW = reshape(w(randperm(M*N)),size(w,1),size(w,2));
    temp = corr(newW,'type',corrType);
    eigenPerm(:,i0) = eig(temp);
end
avePerm = [mean(eigenPerm,1);std(eigenPerm,0,1)];

figure
set(gcf,'Renderer','painters')
hold on
histogram(avePerm(2,72:end),'FaceColor',Bu,'Normalization','probability','LineWidth',0.5);
ah = gca;
ylim = ah.YLim;
plot([aveOpt(2); aveOpt(2)],ylim,'r--','LineWidth',2)
box on
set(gca,'XLim',[0.9, 1.6])
legend('shuffle','original')
xlabel('$\sigma_{\lambda}$','Interpreter','latex')
ylabel('probability')
% prefix = ['fig2_CorrComp_Column_eigen_',corrType,'_',trc,'_N',num2str(N),...
%     'M',num2str(M),'sig',num2str(sig),....
%     'sp',num2str(sp),'_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])


figure
set(gcf,'Renderer','painters')
hold on
for i0 = 1:repeats
    [f,x] = ecdf(eigenPerm(72:end,i0));
    plot(x,f,'LineWidth',1,'Color',Gr(6,:))
end

[f,x] = ecdf(eigenOpt(72:end));
plot(x,f,'LineWidth',2,'Color','red')
hold off
lg = legend('shuffle','original');
legend boxoff
xlabel('$\lambda $','Interpreter','latex')
ylabel('ECDF')
% prefix = ['fig2_CorrComp_Column_CDF_',corrType,'_',trc,'_N',num2str(N),'M',...
%     num2str(M),'sig',num2str(sig),....
%     'sp',num2str(sp),'_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])


%% figS2 d
% parameters
nOdor = 100;
nRecp = 30;
spar = 2;
sig =2;
% looad data
dFolder = '../data';
load(fullfile(dFolder,'gcmi_diffEnt_OptW_randomShuffle_N100M30sig2sp2_03-Oct-2018_2.mat'))
% set graphics 
% defaultGraphicsSetttings
%define some colors using brewermap
Bu = brewermap(11,'Blues');    % blues
Gr = brewermap(11,'Greys');    % greys

figure
hold on
histogram(-mean(allH,2),20,'Normalization','probability','FaceColor',Bu(9,:))
Yrange = get(gca,'Ylim');
plot(-[mean(allH0);mean(allH0)],Yrange','--')
legend('shuffle','original')
legend boxoff
box on
set(gca,'XLim',[-15,-12.5])
xlabel('differential entropy')
ylabel('probability')
% prefix = ['figS2_gcmi_entrComp_OptW_Shuff','_N',num2str(nOdor),'M',num2str(nRecp),'sig',num2str(sig),....
%     'sp',num2str(spar),'_',date,'_2'];
% saveas(gcf,[outFolder,filesep,prefix,'.fig'])
% print('-painters','-depsc',[saveFolder,filesep,prefix,'.eps'])

%%
% ==============================================
% map the elements of original matrix into a stardard 
% Gaussian distribution and compare the eigenvalues
% ==============================================
w = reshape(allMat(:,ix0),[M,N]);

[f,x] = ecdf(w(:));
f(end) = 1-1e-4;
[~, inx] = sort(w(:));   %index
newW = norminv(f(2:end));
w0 = zeros(size(w,1),size(w,2));
w0(inx) = newW;

C0 = corr(w0','type',corrType);
eigenOpt = eig(C0);
eigenPerm = zeros(M, repeats);
aveOpt = [mean(eigenOpt),std(eigenOpt)];

for i0 = 1:repeats
    newW = reshape(w0(randperm(M*N)),size(w0,1),size(w0,2));
    temp = corr(newW','type',corrType);
    eigenPerm(:,i0) = eig(temp);
end
avePerm = [mean(eigenPerm,1);std(eigenPerm,0,1)];
largest = max(eigenPerm);  %largest eigen values

% compare the standard deviation of eigenvalues
figure
hold on
histogram(avePerm(2,:),'FaceColor',Bu(9,:),'Normalization','probability');
ah = gca;
ylim = ah.YLim;
plot([aveOpt(2); aveOpt(2)],ylim,'r--','LineWidth',2)
box on
legend('shuffle','original')
xlabel('$\sigma_{\lambda}$','Interpreter','latex')
ylabel('probability')

% compare the maxium eigenvalue
figure
hold on
histogram(largest,'FaceColor',Bu(9,:),'Normalization','probability');
ah = gca;
ylim = ah.YLim;
plot([max(eigenOpt); max(eigenOpt)],ylim,'r--','LineWidth',2)
box on
legend('shuffle','original')
xlabel('$\lambda_{max}$','Interpreter','latex')
ylabel('probability')

% ecdf of data and Gaussian
figure
plot(x,f)
xlabel('w')
ylabel('ecdf')

% standard Gaussian
pd = makedist('Normal');
X = -3:.1:3;
cdf_normal = cdf(pd,X);
figure
plot(X,cdf_normal)
xlabel('x')
ylabel('cdf')
