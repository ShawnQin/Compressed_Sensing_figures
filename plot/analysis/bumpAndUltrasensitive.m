% this program extrat the fraction of ultra-sensitive elements in the
% sensitivity matrix
% In our previous simulation we observe that a little "bump" appear in the
% right sie of the W
% This might be a coding strategy

close all
clear

% default folder
dFolder = '/Users/shan/Documents/GoogleDrive/olfactoryCoding/code/data';
figFolder = '../figures';  % folder to store the figures
selpath = uigetdir(dFolder);

allFile = dir(fullfile(selpath,filesep,'*.mat'));
files = {allFile.name}';

s1 = '(?<= *N)[\d.]+(?=_)';
s2 = '(?<= *_R)[\d]+(?=_)';
s3 = '(?<= *_S)[\d]+(?=_)';
s5 = '(?<= *sig)[\d.]+(?=_)';
s6 = '(?<= *alp)[-,\d]+(?=_)';

allParaName  = {'N','M','n','$\sigma_c$','$\hat{\alpha}$'};
numOr = zeros(length(files),1);
numRecp = zeros(length(files),1);
spInput = zeros(length(files),1);
allSig = zeros(length(files),1);
allAlp = zeros(length(files),1);


% store the fraction of ultrasensitive elements
ultrSensFrac = zeros(length(files),2);
allThr = zeros(length(files),2);
allSensi = cell(length(files),1);

lcut = -15;  % only show the elements that are larger than thr
hcut = 10;   % only show elements that are samller than this value

for i0 = 1:length(files)
    load(fullfile(selpath,files{i0}))   
    numOr(i0) = str2num(char(regexp(files{i0},s1,'match')));

    numRecp(i0) = str2num(char(regexp(files{i0},s1,'match')));
    numRecp(i0) = str2num(char(regexp(files{i0},s2,'match')));
    spInput(i0) = str2num(char(regexp(files{i0},s3,'match')));
    allSig(i0) = str2num(char(regexp(files{i0},s5,'match')));
    allAlp(i0) = str2num(char(regexp(files{i0},s6,'match')));
    
    allActW = nan(size(allMat,2),1);
   
    
    % plot and set the upper and lower limit
    figure
    figureSize = [4 4 16 5];
    set(gcf,'Units','inches','Position',figureSize,'PaperPositionMode','auto');
    histogram(allMat(allMat > lcut & allMat < hcut))
    % select the lower bound
    disp('select the lower threshould by clicking!')
    [LB,~] = getpts;
    % select the upper bound
    disp('select the upper threshould by clicking!')
    [UB,~] = getpts;
    allThr(i0) = UB;
    allUltr = sum(allMat > UB,1)/size(allMat,1);
    ultrSensFrac(i0,:) = [mean(allUltr),std(allUltr)];
    allSensi{i0} = allMat(allMat > LB);
    
%     for j0 = 1:size(allMat,2)
%         temp = allMat(:,j0);
%         t = temp(temp > LB);
%         actData = t(t<=UB);
%         pd = fitdist(actData,'normal');
%         allActW(j0) = pd.mean;
%     end
       
end

% check which parameter is the variable
allParam = [numOr,numRecp,spInput,allSig,allAlp];
uniqPara = [length(unique(numOr)),length(unique(numRecp)),...
    length(unique(spInput)),length(unique(allSig)),length(unique(allAlp))];
ix = find(uniqPara > 1);

%% plot the figure
defaultGraphicsSetttings
myBu = brewermap(11,'Blues');
colorInx = 8;

% ================================================
% plot the fraction of ultrasensitive elements 
% ================================================

% sort the index of parameter
[~,newInx] = sort(allParam(:,ix));
figure
errorbar(allParam(newInx,ix),ultrSensFrac(newInx,1),ultrSensFrac(newInx,2),'o-','MarkerSize',10,...
        'MarkerFaceColor',myBu(colorInx,:),'Color',myBu(colorInx,:),'LineWidth',2,...
        'CapSize',0)
ylim([0.0,0.05])
xlabel(allParaName(ix),'interpreter','latex')
ylabel('$\rho_w^u$','interpreter','latex')
prefix = ['ultrSensi_SkewGauss_alp_N100M20_',date];
saveas(gcf,[figFolder,filesep,prefix,'.fig'])
print('-depsc',[figFolder,filesep,prefix,'.eps'])


% ================================================
% plot a tight histogram in 2 by 4 array
% ================================================
figure
figureSize = [0 0 12 10];
set(gcf,'Units','inches','Position',figureSize,...
'PaperPositionMode','auto');

% here we used "tightplot" to set the marginsre
row = 3;
col = 3;
% [ha, pos] = tight_subplot(row,col,[.02 .01],[.02 .05],[.08 .05]);

% only show the sensitive part
for i0 = 1:length(allSensi)
    subplot(row,col,i0)
    temp = allSensi{i0};
    histogram(temp(temp <10),'Normalization','pdf')
    if mod(i0,col) ==1
        ylabel('pdf')
    end
    
    if i0 > (row-1)*col
        xlabel('$\ln w$','interpreter','latex')
    end
end
