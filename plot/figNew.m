% this program show the distribution of sensitive W depends on input
% distribution
% we compare (1) lognoraml and skewed gaussian;(2) lognormal and stretched
% exponential
% To plot the figures, we need three data set

close all
clear

%% color setting and default folder to load and store the figures

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
lRd = [255, 141,139]/255; % light red
brickRd = [201,69,89]/255;  %brick red
green = [107,169,81]/255;  %green
purple = [113,49,119]/255;  % purple

% load the summary data
dFolder = '../data';
saveFolder = '../figures';


%% load data, lognormal distribution
% basic parameters
N = 50;
M = 13;
sp = 3;
sig = 2;
num = 40;    % number of repeats
%sigma_c dependent
sigFolder = '../data/lognorm_SigDp/GcmiDifcons10-15';
% fName = 'gcmi_sigdp_N50M13Sp3_19-Oct-2018.mat';
% load(fullfile(dFolder,fName));
% summSigdp = dataSumm;

allFile = dir(fullfile(sigFolder,filesep,'*.mat'));
files1 = {allFile.name}';
s5 = '(?<= *sig)[\d.]+(?=_)';
allSig = zeros(length(files1),1);
skGauss_allSk = nan(length(files1),num);
allCDF_ln = cell(length(files1),1);   % store all the ecdf

LB = -10;       %lower bound to register sensitive 
UB = 2.5;       % upper bound cut off

% define the skewed gaussian used to fit the senstive elements
for i0 = 1: length(files1)
    load(fullfile(sigFolder,files1{i0}))
    allSig(i0) = str2num(char(regexp(files1{i0},s5,'match')));
    
    % calculate the skewness by fitting a skewed Gaussian or direct
    % calculae
    [f,X] = ecdf(allMat(allMat > LB));
    allCDF_ln{i0} = [f,X];
    for j0 = 1:size(allMat,2)
        temp = allMat(:,j0);
        sensiW = temp(temp > LB & temp < UB);
        skGauss_allSk(i0,j0) = skewness(sensiW);          
    end  
end

[~,ix] = sort(allSig);
plot(allSig(ix),nanmean(skGauss_allSk(ix),2))


%% streched power law
% newFolder = '/Users/shan/Documents/GoogleDrive/olfactoryCoding/code/data/GcmiExp20190201-1';
newFolder = '../data/ExpSmallBeta0216';
allFile = dir(fullfile(newFolder,filesep,'*.mat'));
files = {allFile.name}';
s5 = '(?<= *sig)[\d.]+(?=_)';  % here sigma is just beta that coresponding to same variance of sigma

% calculate the skewness by fitting a skewed Gaussian distribution
% store the summary results
num = 40;   % repeat of simulation
fmin = nan(length(files),2);
allSk = nan(length(files),num);
fitSk = nan(length(files),num);
allSpW = nan(length(files),num);
sig = zeros(length(files),1);   % store the beta that correspond to same sigma_c in lognormal
allCDF_pl = cell(length(files),1);   % store all the ecdf


LB = -10;     % lower bound to register sensitivity
UB = 5;       % upper bound cut off

% define the skewed gaussian used to fit the senstive elements
for i0 = 1: length(files)
    load(fullfile(newFolder,files{i0}))
    sig(i0) = str2num(char(regexp(files{i0},s5,'match')));
    
    LB = -3*sqrt(2)/sig(i0);  % threshold 
    % calculate the skewness by fitting a skewed Gaussian or direct
    % calculae
    fmin(i0,:) = [mean(-allfmin),std(-allfmin)];
    [f,X] = ecdf(allMat(allMat > LB));
    allCDF_pl{i0} = [f,X];

    for j0 = 1:size(allMat,2)
        temp = allMat(:,j0);
        sensiW = temp(temp > LB & temp < UB);
        allSk(i0,j0) = skewness(sensiW);
        allSpW(i0,j0) = sum(temp>LB)/length(temp);  % sparsity
            
        % fit for individual matrix
%         [fitPara, fitSk(i0,j0)] = fitSkewedGauss(sensiW);
    end
    
    % fit the sensitive elements with a skewed Gaussian, use all the data
%     allSen = allMat(allMat > LB & allMat < UB);
%     [fitPara, fitSk(i0)] = fitSkewedGauss(allSen,'plot');
end

% ====================================
% plot how skewness change with beta
% ====================================
figure
hold on
[~,inx] = sort(round(sig*100));
errorbar(sig(inx(3:end)),nanmean(allSk(inx(3:end),:),2),nanstd(allSk(inx(3:end),:),0,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Gr(11,:),'Color',Gr(11,:),'LineWidth',2,'CapSize',0)

ylim([-3,1])
% set(gca,'FontSize',28)
xlabel('$\beta$','interpreter','latex','FontSize',20)
ylabel('Skewness','FontSize',20)
box on
set(gca,'FontSize',20,'LineWidth',1)
% prefix = ['Skewness_powerlaw_alp_N50M13_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% ====================================================
% plot how sparsity of optimal W change with \beta together 
% with the skewness using double y-axis
% ====================================================

figure
hold on
[~,inx] = sort(round(sig*100));
errorbar(sig(inx(3:end)),nanmean(allSpW(inx(3:end),:),2),nanstd(allSpW(inx(3:end),:),0,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Gr(11,:),'Color',Gr(11,:),'LineWidth',2,'CapSize',0)

ylim([0.6,0.75])
set(gca,'XTick',0.2:0.2:1)
xlabel('$\beta$','interpreter','latex')
ylabel('$\rho_w$','interpreter','latex')
box on
% prefix = ['sparsityW_powerlaw_alp_N50M13_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% ========================================================
% plot how skewness and sparisity of W change with beta
% using double y-axis
% ========================================================
figure
hold on
[~,inx] = sort(round(sig*100));
selInx = inx(3:end);

xlim([0.25,1.05])
% xlabel('$\sigma_c$','interpreter','latex')

% add the second y axis to show the differential entropy
errorbar(sig(selInx),nanmean(allSpW(selInx,:),2),nanstd(allSpW(selInx,:),0,2),'diamond-','MarkerSize',12,...
    'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
box on
ylim([0.6,0.75])
xlabel('$\beta$','interpreter','latex')
ylabel('$\rho_w$','interpreter','latex')


yyaxis right
errorbar(sig(selInx),nanmean(allSk(selInx,:),2),nanstd(allSk(selInx,:),0,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Gr(11,:),'Color',Gr(11,:),'LineWidth',2,'CapSize',0)
ylim([-3,1])
ylabel('Skewness')

ax = gca;
ax.YAxis(2).Color = 'k';

% prefix = ['Skewness_spW_powerlaw_alp_N50M13_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% ================================================================================
% plot distribution of W with two different beta
% ================================================================================
bt1 = 0.3;
bt2 = 0.7;
ix1 = find(round(sig*100) == bt1*100);
ix2 = find(round(sig*100) == bt2*100);

% load the data 
load(fullfile(newFolder,files{ix1}));
data1 = allMat(:);
sens1 = data1(data1 > -3*sqrt(2)/bt1 & data1 < UB);

load(fullfile(newFolder,files{ix2}));
data2 = allMat(:);
sens2 = data2(data2 > -3*sqrt(2)/bt2 & data2 < UB);


% plot the pdf of the two
% compare their empirical pdf
% [f1,X1] = ksdensity(sens1*bt1/sqrt(2));
% [f2,X2] = ksdensity(sens2*bt2/sqrt(2));
[f1,X1] = ksdensity(sens1);
[f2,X2] = ksdensity(sens2);
figure
hold on
plot(X1,f1)
plot(X2,f2)
box on
lg = legend('\beta = 0.3','\beta = 0.7','Location','northeast');
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('emperical pdf')
% xlim([-10,2])
set(gca,'YTick',10.^(-6:2:0))
set(gca,'YScale','log')

prefix = ['figNew_ecdf_examp_symmPower_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])
% ylabel('$\ln(pdf)$','Interpreter','latex') 

%% show the differenence between the sensitive elements when input distribution changes

% first, compare the histogram difference between Gaussian and Skewed
% select a sigma to plot
sigSel = 3.2;
ix1 = find(round(allSig*100) == sigSel*100);
% ix2 = find(round(sig*100) == sigSel*100);
ix2 = find(round(sig*100) == 0.7*100);

% load first example data
load(fullfile(sigFolder,files1{ix1}));
dataSkGauss = allMat;
sens1 = dataSkGauss(dataSkGauss > LB & dataSkGauss < UB);

load(fullfile(newFolder,files{ix2}));
datapowerlaw = allMat;
sens2 = datapowerlaw(datapowerlaw > LB & datapowerlaw < UB);

figure
hold on
h1 = histogram(sens1,50,'Normalization','pdf','FaceColor',lBu,'FaceAlpha',...
0.4,'EdgeColor','none');
stairs([h1.BinEdges(1),h1.BinEdges,h1.BinEdges(end)],...
[0,h1.Values,h1.Values(end),0],'Color',dpBu,'LineWidth',2)

h2 = histogram(sens2,30,'Normalization','pdf','FaceColor',lRd,'FaceAlpha',...
0.4,'EdgeColor','none');
stairs([h2.BinEdges(1),h2.BinEdges,h2.BinEdges(end)],...
[0,h2.Values,h2.Values(end),0],'Color',myRd,'LineWidth',2)


hold off
lg = legend('log-normal','exp(-\beta|ln(c)|)','Location','northwest');
set(lg,'FontSize',16)
legend boxoff
box on
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('pdf','Interpreter','latex') 


% compare their empirical pdf
[f1,X1] = ksdensity(sens1);
[f2,X2] = ksdensity(sens2);
figure
hold on
plot(X1,f1)
plot(X2,f2)
box on
lg = legend('log-normal','exp(-\beta|ln(c)|)','Location','northwest');
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('emperical pdf')
xlim([-10,2])
prefix = ['figNew_ecdf_examp_log',num2str(sigSel),'_',date];
saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
print('-depsc',[saveFolder,filesep,prefix,'.eps'])
% ylabel('$\ln(pdf)$','Interpreter','latex') 
set(gca,'YScale','log')


% ========================================================
% calculate their K-L divergence
% ========================================================
shift = true;   % true or false, shift or not
mult = 1;
% co = DKL_kNN_k_initialization(mult);
co = DKL_vME_initialization(mult);
sameSig = intersect(allSig,sig);  % same sigma
[sortedSig,~] = sort(sameSig);
% KLdiv = zeros(length(sortedSig),1);
KLdiv = zeros(length(sortedSig),num);


UB = 3;
X = LB:0.01:5;
for i0 = 1:length(sortedSig)
    ix1 = find(allSig == sortedSig(i0));
    ix2 = find(sig == sortedSig(i0));
     
    load(fullfile(sigFolder,files1{ix1}));
    allMat1 = allMat;
%     dataSkGauss = allMat;
%     sens1 = dataSkGauss(dataSkGauss > LB & dataSkGauss < UB);
%     f1 = ksdensity(sens1,X);  % emperical pdf
    
    load(fullfile(newFolder,files{ix2}));
    allMat2 = allMat;
%     datapowerlaw = allMat;
%     sens2 = datapowerlaw(datapowerlaw > LB & datapowerlaw < UB);
%     f2 = ksdensity(sens2,X);  % emperical pdf
%     
%     KLdiv(i0) = DKL_vME_estimation(f1,f2,co);
    
    % calculate each matrix
    
    COL = min(size(allMat1,2),size(allMat2,2));  % matrix from two data set may have different columns
    for j0 = 1:COL
        temp1 = allMat1(:,j0);
        sens1 = temp1(temp1 > LB & temp1 < UB);
       
        temp2 = allMat2(:,j0);
        sens2 = temp2(temp2 > LB & temp2 < UB);
        if shift
            % first find the middle point
            [f,X] = ecdf(temp1(temp1 > LB & temp1 < UB));
            ecdf_ln = [f,X];
            mIdx1 = find(ecdf_ln >= 0.5,1,'first');
            mediamW1 = ecdf_ln(mIdx1,2);
            
            
            [f,X] = ecdf(temp2(temp2 > LB & temp2 < UB));
            ecdf_pl = [f,X];
            mIdx2 = find(ecdf_pl >= 0.5,1,'first');
            mediamW2 = ecdf_pl(mIdx2,2);
            
            f1 = ksdensity(sens1-mediamW1,X-mediamW1);  % emperical pdf 
            f2 = ksdensity(sens2-mediamW2,X-mediamW2);  % emperical pdf 
            
        else
            f1 = ksdensity(sens1,X);  % emperical pdf 
            f2 = ksdensity(sens2,X);  % emperical pdf
 
        end
             
        
        KLdiv(i0,j0) = DKL_vME_estimation(f1,f2,co);
    end
end

% sort the results
[selSig,sortInx] = sort(sortedSig);

%plot how KL-devergence change with sigma_c
figure
KLplot = KLdiv;
KLplot(KLplot==0) = nan;
meanKL = nanmean(KLplot,2);
stdKL = nanstd(KLplot,0,2);
errorbar(selSig,meanKL(sortInx),stdKL(sortInx),'o-','MarkerSize',12,...
        'MarkerFaceColor',myBu,'Color',myBu,'LineWidth',2,'CapSize',0)
% plot(selSig,KLdiv(sortInx),'o-','MarkerSize',12,'MarkerFaceColor',lBu,...
%     'MarkerEdgeColor',dpBu,'LineWidth',2)
xlabel('$\sigma_c$','Interpreter','latex')
ylabel('K-L divergence')
ylim([0,2])

if shift
    prefix = ['figNew_KL_lognorm_diffSig_shift',date];
else
    prefix = ['figNew_KL_lognorm_diffSig',date];
end
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% ===============================
% compare the empirical cdf
% ===============================

figure
hold on
% for i0 = 1:length(alp)
%     plot(allCDF{inx(i0)}(:,2),allCDF{inx(i0)}(:,1))
% end
mIdx1 = find(allCDF_ln{ix1}(:,1) >= 0.5,1,'first');
mediamW1 = allCDF_ln{ix1}(mIdx1,2);
mIdx2 = find(allCDF_pl{ix2}(:,1) >= 0.5,1,'first');
mediamW2 = allCDF_pl{ix2}(mIdx2,2);
plot(allCDF_ln{ix1}(:,2)-mediamW1,allCDF_ln{ix1}(:,1))
plot(allCDF_pl{ix2}(:,2)-mediamW2,allCDF_pl{ix2}(:,1))

box on
hold off
lg = legend('log-normal','exp(-\beta|ln(c)|)','Location','northwest');
set(lg,'FontSize',16)
legend boxoff
box on
xlim([-5,UB])
xlabel('$\Delta\ln(w)$','Interpreter','latex')
ylabel('empirical cdf')

% prefix = ['figNew_ecdf_lognorm_power_shift',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])



%% skewed Gaussian
skewFolder = '../data/gcmi_skew/N100R20sig3_noReg';
allFile = dir(fullfile(skewFolder,filesep,'*.mat'));
fileNameCell = {allFile.name}';

num = 40;   % repeat of simulation
str_marker = 'N100_R20';   %folder with this string contains the data we need
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fileNameCell,str,'once'));
files = fileNameCell(FIND(str_marker));

N = 100;
M = 20;
sig = 3;

% threshold to register the active elements
LB = -6;
UB = 4;

s5 = '(?<= *sig)[\d.]+(?=_)';
s6 = '(?<= *alp)[-\d.]+(?=_)';


allSk = nan(length(files),num);
allKt = nan(length(files),num);
spW = nan(length(files),num);     % store the sparsity information
KLdiv = nan(length(files),1);     % store the K-L divergence
rawSensi = cell(length(files),1); % store all the sensitive elements
allCDF = cell(length(files),1);   % store all the ecdf
fmin = nan(length(files),2); 


% select two example matrices for future use
alp_slt = [-2,2];  %select data from these two alp
W_ix = [2,2];
seltW = cell(2,1);

alp = zeros(length(files),1);  % store the skewness parameters
for i0 = 1:length(files)
    load(fullfile(skewFolder,files{i0}))

    alp(i0) = str2num(char(regexp(files{i0},s6,'match')));
    
    fmin(i0,:) = [mean(-allfmin),std(-allfmin)];
    [f,X] = ecdf(allMat(allMat > LB));
    allCDF{i0} = [f,X];
    rawSensi{i0} = allMat(allMat > LB & allMat < UB);
    
    for j0 = 1:size(allMat,2)
        temp = allMat(:,j0);
        sensiW = temp(temp > LB & temp <= UB);  % here the upper limit should be checked,2.5
        spW(i0,j0) = length(temp(temp > LB))/length(temp);
        allSk(i0,j0) = skewness(sensiW);
        allKt(i0,j0) = 	kurtosis(sensiW)-3;  
        
        if ismember(alp(i0),alp_slt)
            ix = find(alp_slt==alp(i0));
            if j0 == W_ix(ix)
                seltW{ix} = sensiW;
            end
        end
    end
end

%plot the summary of the results
simSk = [nanmean(allSk,2),nanstd(allSk,0,2)];
simKt = [nanmean(allKt,2),nanstd(allKt,0,2)];

% ================================================
% plot how skewness change with beta
% ================================================
figure
hold on
[~,inx] = sort(round(alp*100));
errorbar(alp(inx),simSk(inx,1),simSk(inx,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
% errorbar(alp(inx),gaussSk(inx,1),gaussSk(inx,2),'o-','MarkerSize',12,...
%     'MarkerFaceColor',Gr(7,:),'Color',Gr(7,:),'LineWidth',2,'CapSize',0)
% lg = legend('simulation','Gaussian');
% set(lg,'FontSize',20)
% ylim([-0.7,0.7])
xlabel('$\hat{\alpha}$','interpreter','latex')
ylabel('Skewness')
box on
% prefix = ['Skewness_SkewGauss_alp_N',num2str(N),'M',num2str(M),'sig',num2str(sig),'_',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% ===============================
% compare the empirical cdf
% ===============================
% [orderAlp,inx] = sort(round(alp*100));
selAlp =  [-4,0,4];  % only plot these three
figure
hold on
% for i0 = 1:length(alp)
%     plot(allCDF{inx(i0)}(:,2),allCDF{inx(i0)}(:,1))
% end
for i0 = 1:length(selAlp)
    ix = find(alp == selAlp(i0));
    % shift the position according to medium value
    mIdx = find(allCDF{ix}(:,1) >= 0.5,1,'first');
    mediamW = allCDF{ix}(mIdx,2);
%     plot(allCDF{ix}(:,2) - mediamW,allCDF{ix}(:,1))
    plot(allCDF{ix}(:,2),allCDF{ix}(:,1))
end
box on
xlim([-4,3])
lg = legend('\hat{\alpha} = -4', '\hat{\alpha} = 0', '\alpha = 4','Location','northwest');
set(lg,'FontSize',16)
legend boxoff

hold off
ylabel('emperical cdf')
xlabel('$\Delta\ln(w)$','interpreter','latex')

% prefix = ['figNew_cdf_skewGauss_diffSig_noShift',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% ===============================
% compare the empirical pdf
% ===============================
[orderAlp,inx] = sort(round(alp*100));
selAlp =  [-2,2];  % only plot these three

figure
hold on

for i0 = 1:length(selAlp)
    ix = find(alp == selAlp(i0));
    load(fullfile(skewFolder,files{ix}))
    sensi = allMat(allMat > LB & allMat < UB);
    [f,loc] = ksdensity(sensi);
    
    plot(loc,f)
end
box on
xlim([-5,3])
lg = legend(['\hat{\alpha} =',num2str(selAlp(1))], ['\hat{\alpha} =',num2str(selAlp(1))],'Location','northwest');
set(lg,'FontSize',16)
legend boxoff

hold off
ylabel('emperical cdf')
xlabel('$\ln(w)$','interpreter','latex')

% prefix = ['figNew_pdf_skewGauss_diffSig_noShift',date];
prefix = ['figNew_pdf_skewGauss_noShift',num2str(N),'M',num2str(M),'sig',num2str(sig),'_',date];

% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])


% ====================================
% K-L divergence distance
% ====================================
X = LB:0.01:UB;
absAlp = abs(alp);
alpPair = unique(absAlp(absAlp >0)); % all the alpha pair
for i0 = 1:length(alpPair)
    ix1 = find(alp == alpPair(i0));   % positive 
    ix2 = find(alp == -alpPair(i0));  % negative 
    
    % fit a ksdensity
    temp1 = rawSensi{ix1};
    sens1 = temp1(temp1 > LB & temp1 < UB);
    f1 = ksdensity(sens1,X);  % emperical pdf        
        
    temp2 = rawSensi{ix2};
    sens2 = temp2(temp2 > LB & temp2 < UB);
    f2 = ksdensity(sens2,X);  % emperical pdf   
    
    KLdiv(i0) = DKL_vME_estimation(f1,f2,co);
    
end

%plot the KL divergence with different alpha
figure
plot(absAlp,KLdiv,'o-','MarkerSize',12,'MarkerEdgeColor',myBu,...
    'MarkerFaceColor',myBu,'LineWidth',2)
ylim([0,0.7])
xlabel('$|\hat{\alpha}|$','Interpreter','latex')
ylabel('K-L divergence')

% prefix = ['figNew_KLdiv_skewGauss_diffAlp',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])


%% comparision of input distribution
sig = 3;
bt = sqrt(2)/sig;
gaussian = @(x) (1/sqrt((2*pi))/sig*exp(-x.^2/2/sig^2));
symmetricPower = @(x) bt/2.*exp(-bt*(abs(x)));

X = -10:0.05:10;

figure
hold on
plot(X,gaussian(X))
plot(X,symmetricPower(X))
hold off
xlabel('$\ln(w)$','Interpreter','latex')
ylabel('pdf')
box on
% legend('lognormal','symmetric powerlaw')
% prefix = ['fig4new_lognoram_symmetricPowerlaw',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])



