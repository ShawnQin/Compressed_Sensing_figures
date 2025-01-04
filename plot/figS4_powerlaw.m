% This program plot figures  related to power-law odor concentration
% distribution

% last revised 1/1/2019

%% graphical settings
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
saveFolder = '../figures';
% ===============================================================
% heatmap and histogram of an example 
% ===============================================================
selFolder = '../data';
selFile = 'Gcmi_power_N100_R20_S5_sig1.5_2018-10-15.mat';
load(fullfile(selFolder,selFile))

N = 100;
M = 20;
sp = 5;
alp = 1.5;

ix = 9;  % select on of the matrix
w = reshape(allMat(:,ix),[M,N]);

% heatmap
figure
set(gcf,'renderer','Painters')
imagesc(w,[-6,3])
set(gca,'TickLength',[0,0])
colormap(jet);
c = colorbar; 
xlabel('Odorant')
ylabel('Receptor')
% prefix = ['figS4_powerlaw_exampW_N',num2str(N),'M',num2str(M),'alp',num2str(alp),....
%     'sp',num2str(sp),'_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% histogram
figure
hold on
set(gcf,'renderer','Painters')
h1 = histogram(w(:),60,'Normalization','probability');
h1.FaceColor = lBu; h1.FaceAlpha = 0.4; h1.EdgeColor = 'none';
stairs([h1.BinEdges(1),h1.BinEdges,h1.BinEdges(end)],...
    [0,h1.Values,h1.Values(end),0],'Color',dpBu,'LineWidth',2)

box on
hold off
set(gca,'Layer','top')
xlim([-80,5])
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('probability')
% prefix = ['figS4_power_exampW_N',num2str(N),'M',num2str(M),'alp',num2str(alp),....
%     'sp',num2str(sp),'_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% enlarged
figure
set(gcf,'renderer','Painters')
hold on
h2 = histogram(w(abs(w)<8),20,'Normalization','pdf');
h2.FaceColor = lBu; h2.FaceAlpha = 0.4; h2.EdgeColor = 'none';
stairs([h2.BinEdges(1),h2.BinEdges,h2.BinEdges(end)],...
    [0,h2.Values,h2.Values(end),0],'Color',dpBu,'LineWidth',2)

% get the normal fit parameter
% pd = fitdist(w(abs(w)<6),'normal');
% X = -6:0.05:6;
% Y = normpdf(X,pd.mean,pd.sigma);
% plot(X,Y,'Color',Or,'LineWidth',4)
hold off
box on

% set(gca,'XLim',[-6,2])
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex')
% prefix = ['figS4_power_exampActiveW_N',num2str(N),'M',num2str(M),'alp',num2str(alp),....
%     'sp',num2str(sp),'_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

