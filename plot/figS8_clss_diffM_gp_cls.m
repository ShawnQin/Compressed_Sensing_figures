% This script plot how calssification error change with different
% parameters, such as the number of receptors, number of groups (labels)
% and the number of clusters

% last revised on 12/27/2018

%% graphics settings
close all
clear

defaultGraphicsSetttings

%define some colors using brewermap
RdBu = brewermap(11,'RdBu');   % red and blue
Bu = brewermap(11,'Blues');    % blues
Gr = brewermap(11,'Greys');    % greys
seqColor = brewermap(21,'YlGnBu');

% define my own colors
CN = 255;
myRed = [202,39,44]/CN;
myBrickRed = [242,65,95]/CN;
myOrang = [224,135,51]/CN;
myYellow = [255,241,0]/CN;
myGreen = [10,188,102]/CN;
myCyan = [84,195,241]/CN;
myLightPurple = [232,82,15]/CN;
myDarkPurple = [164,11,93]/CN;
myBlue = [3,110,184]/CN;
myBlueGreen = [76,216,196]/CN;

lBu = [96,166,223]/255; %light blue
dpBu = [63,114,183]/255; % deep blue
dkBu = [50,78,147]/255;   %dark blue
Or = [220,150,71]/255;  % orange
brickRd = [201,69,89]/255;  %brick red
green = [107,169,81]/255;  %green
purple = [113,49,119]/255;  % purple


errorBarSize = 1.5;
errorMarkerSize = 10;
LineWidth = 1.5;
labelFontSize = 24;
axisFontSize = 20;
ticketWidth = 1.5;

%%
outFolder = '../figures';

% ================================================================
% Classification with different groups
% ================================================================
dFolder = '../decoding/data/LDA_N100M10DiffG';

allFile = dir(fullfile(dFolder,filesep,'*.mat'));
filesRaw = {allFile.name}';

N = 100;
M = 10;
sp = 3;
sig = 2;
H = 500; %hidden layer size
nSig = 0.1;  %noise std

str_marker = '27-Dec-2018';  % only select data in the same batch
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(filesRaw,str,'once'));
files = filesRaw(FIND(str_marker));

allGroup = 2:1:6;     % all different groups
allSp = 0.1:0.05:1;   % the sparisty
meanError = zeros(length(allSp),length(allGroup));
stdError = zeros(length(allSp),length(allGroup));
minError = zeros(length(allGroup),2); % position and value
bestPerf = zeros(length(allGroup),2);  %store information of the optimal sparsity

optParam = zeros(length(files),1);
bestPerform = zeros(length(files),1);
TOL = 2e-5;
for i0 = 1:length(files)
    s1 = '(?<= *nType)[\d]+(?=_)';
    group = str2num(char(regexp(files{i0},s1,'match')));
    inx = find(allGroup == group);
    
    load(char(fullfile(dFolder,filesep,files{i0})));
    
    meanError(:,inx) = mean(allTestError,2);
    stdError(:,inx) = std(allTestError,0,2);
    
    [minX,minY] = sort(meanError(:,inx));
    minError(i0,2) = minX(1);  %minimum error
    minError(i0,1) = minY(1);  %minimum position
    
    % fit and find the minium error
%     for j0 = 1:size(allTestError,2)
%         figure
%         hold on
%         plot(allSp',allTestError(:,j0),'o')
%         cs = spaps(allSp, allTestError(:,j0),TOL,1,3);
%         fnplt(cs)
%         hold off
%         d1 = fnder(cs);
%         minSp = fnzeros(d1,[0.1,0.99]);
%         optParam(inx,j0) = minSp(1,1);
%         bestPerform(inx,j0) = fnval(cs,optParam(i0,1));
%         
%     end
%     close all
    figure
    hold on
    plot(allSp',meanError(:,inx),'o')
    cs = spaps(allSp', meanError(:,inx),TOL,1,3);
    fnplt(cs)
    hold off
    d1 = fnder(cs);
    minSp = fnzeros(d1,[0.1,0.99]);
    optParam(inx) = minSp(1,1);
    bestPerform(inx) = fnval(cs,optParam(inx,1));
end

% plot the figure
figure
hold on
for i0 = 1:length(allGroup)
    plot(allSp',meanError(:,i0),'Color',Bu(2*i0 + 1,:))
    [ix, iy] = sort(meanError(:,i0));
    bestPerf(i0,1) = allSp(iy(1));
    bestPerf(i0,2) = ix(1);
end
hold off
box on
lg = legend('group = 2','group = 3','group = 4','group = 5','group = 6');
set(lg,'FontSize',16)
xlabel('$\rho_w$','interpreter','latex')
ylabel('Error')
% figNamePref = ['figS8_class_Error_N',num2str(N),'M',num2str(M),...
%     'sp',num2str(sp),'_nSig',num2str(sig),'_groupDpd_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% plot how the best performance change with group numbers
figure
% errorbar(allGroup',mean(optParam,2),std(optParam,0,2),'o-','MarkerSize',errorMarkerSize,...
%         'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',LineWidth,...
%         'CapSize',0)
plot(allGroup',bestPerform,'o-','MarkerSize',12,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:))
ylim([0.0,0.05])
xlabel('Number of groups')
ylabel('minimum error')
% figNamePref = ['figS8_class_MinError_N',num2str(N),'M',num2str(M),...
%     'sp',num2str(sp),'_nSig',num2str(sig),'_groupDpd_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% plot the optimal sparsity to achieve best performance
figure
plot(allGroup',optParam,'o-','MarkerSize',12,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:))
ylim([0,1])
set(gca,'XTick',2:1:6)
xlabel('groups')
ylabel('$\rho_w^*$','interpreter','latex')
% figNamePref = ['figS8_class_MinErrorPosi_N',num2str(N),'M',num2str(M),...
%     'sp',num2str(sp),'_nSig',num2str(sig),'_groupDpd_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


%%
% ================================================================
% Classification performance with different receptors
% ================================================================
% % dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/decoding/LDA_N100DiffM';
dFolder = '../decoding/data/LDA_N100DiffM';
allFile = dir(fullfile(dFolder,filesep,'*.mat'));
filesRaw = {allFile.name}';

allM = [5,8,10,15,20,25,30];
N = 100;
sp = 3;
sig = 2;
H = 500;
P = 200;


str_marker = '27-Dec-2018';  % only select data in the same batch
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(filesRaw,str,'once'));
files = filesRaw(FIND(str_marker));

allSp = 0.1:0.05:1;   % the sparisty
meanError = zeros(length(allSp),length(allM));
stdError = zeros(length(allSp),length(allM));
minError = zeros(length(allM),2); % position and value
% bestPerf = zeros(length(allM),2);  %store information of the optimal sparsity

optParam = zeros(length(files),1);
bestPerform = zeros(length(files),1);
TOL = 2e-6;
for i0 = 1:length(files)
    s1 = '(?<= *100M)[\d]+(?=sp)';
    M = str2num(char(regexp(files{i0},s1,'match')));
    inx = find(M == allM);
    
    load(char(fullfile(dFolder,filesep,files{i0})));
    
    meanError(:,inx) = mean(allTestError,2);
    stdError(:,inx) = std(allTestError,0,2);
    
    [minX,minY] = sort(meanError(:,inx));
    minError(i0,2) = minX(1);  %minimum error
    minError(i0,1) = minY(1);  %minimum position

    figure
    hold on
    plot(allSp',meanError(:,inx),'o')
    cs = spaps(allSp', meanError(:,inx),TOL,1,3);
    fnplt(cs)
    hold off
    d1 = fnder(cs);
    minSp = fnzeros(d1,[0.1,0.99]);
    optParam(inx) = minSp(1,1);
    bestPerform(inx) = fnval(cs,optParam(inx,1));
end


% plot how classificaiton error change with sparsity, different M
figure
hold on
for i0 = 1:size(meanError,2)
    plot(allSp,meanError(:,i0),'Color',seqColor(4 + 2*i0,:))
    [ix, iy] = sort(meanError(:,i0));
end
lg = legend('M=5','M=8','M=10','M=15','M=20','M=25','M=30','Location','eastoutside');
set(lg,'FontSize',18)
legend boxoff

hold off
box on
xlabel('$\rho_w$','interpreter','latex')
ylabel('classification error')
% figNamePref = ['figS8_class_Error_N',num2str(N),...
%     'sp',num2str(sp),'_nSig',num2str(sig),'_diffM_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% best peformance, minimum error
figure
plot(allM',bestPerform,'o-','MarkerSize',12,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:))
box on
xlabel('$M$','interpreter','latex')
ylabel('minimum error')
% figNamePref = ['figS8_class_MinmumError_N',num2str(N),...
%     'sp',num2str(sp),'_nSig',num2str(sig),'_MDpd_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% sparsity that achieve best performance
figure
plot(allM',optParam(:,1),'o-','MarkerSize',12,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:))
box on
ylim([0.4,0.6])
xlabel('$M$','interpreter','latex')
ylabel('$\rho_w^*$','interpreter','latex')
% figNamePref = ['figS8_class_MinErrorPosi_N',num2str(N),...
%     'sp',num2str(sp),'_nSig',num2str(sig),'_MDpd_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


%%
% ================================================================
% Classification performance with different groups of labels
% ================================================================
% dFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData/decoding/LDA_N100M10_diffCluster';
dFolder = '../decoding/data/LDA_N100M10_diffCluster';
% allFile = dir(fullfile(dFolder,filesep,'*.mat'));
allFile = dir(fullfile(dFolder,filesep,'*.mat'));

filesRaw = {allFile.name}';

allP = [20,50,100,200,500,1000];
N = 100;
sp = 3;
sig = 2;
H = 500;
% P = 200;


str_marker = '27-Dec-2018';  % only select data in the same batch
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(filesRaw,str,'once'));
files = filesRaw(FIND(str_marker));

allSp = 0.1:0.05:1;   % the sparisty
meanError = zeros(length(allSp),length(allP));
stdError = zeros(length(allSp),length(allP));
minError = zeros(length(allP),2); % position and value
% bestPerf = zeros(length(allM),2);  %store information of the optimal sparsity

optParam = zeros(length(files),1);
bestPerform = zeros(length(files),1);
TOL = 2e-6;
for i0 = 1:length(files)
    s1 = '(?<= *P)[\d]+(?=_)';
    P = str2num(char(regexp(files{i0},s1,'match')));
    inx = find(P == allP);
    
    load(char(fullfile(dFolder,filesep,files{i0})));
    
    meanError(:,inx) = mean(allTestError,2);
    stdError(:,inx) = std(allTestError,0,2);
    
    [minX,minY] = sort(meanError(:,inx));
    minError(inx,2) = minX(1);  %minimum error
    minError(inx,1) = minY(1);  %minimum position

    figure
    hold on
    plot(allSp',meanError(:,inx),'o')
    cs = spaps(allSp(3:end)', meanError(3:end,inx),TOL,1,3);
    fnplt(cs)
    hold off
    d1 = fnder(cs);
    minSp = fnzeros(d1,[0.2,0.99]);
    optParam(inx) = minSp(1,1);
    bestPerform(inx) = fnval(cs,optParam(inx));
end


% plot how classificaiton error change with sparsity, different M
figure
hold on
for i0 = 1:size(meanError,2)
    plot(allSp,meanError(:,i0),'Color',seqColor(4 + 2*i0,:))
    [ix, iy] = sort(meanError(:,i0));
end
lg = legend('20','50','100','200','500','1000','Location','eastoutside');
set(lg,'FontSize',18)
legend boxoff

hold off
box on
% set(gca,'YScale','log')
xlabel('$\rho_w$','interpreter','latex')
ylabel('classification error')
% figNamePref = ['figS8_class_Error_N',num2str(N),...
%     'sp',num2str(sp),'_nSig',num2str(sig),'_diffCluster_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


% best peformance, minimum error
figure
plot(allP,bestPerform,'o-','MarkerSize',12,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:))
box on
set(gca,'XScale','log')
xlabel('number of cluster')
ylabel('minimum error')
% figNamePref = ['figS8_class_MinError_N',num2str(N),...
%     'sp',num2str(sp),'_nSig',num2str(sig),'_ClusterDpd_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

% sparsity that achieve best performance
figure
plot(allP,optParam,'o-','MarkerSize',12,'Color',Bu(9,:),'MarkerFaceColor',Bu(9,:))
box on
ylim([0.4,0.7])
xlabel('number of cluster')
ylabel('$\rho_w^*$','interpreter','latex')
% figNamePref = ['figS8_class_MinErrorPosi_N',num2str(N),...
%     'sp',num2str(sp),'_nSig',num2str(sig),'_ClusterDpd_',date];
% saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
% print('-depsc',[outFolder,filesep,figNamePref,'.eps'])

