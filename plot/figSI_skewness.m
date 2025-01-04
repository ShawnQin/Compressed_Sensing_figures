%% Show the skewness of W change with input skewness
% We used a skewed Gaussian concentration distribution and observed that
% the distribution of optimal sensitivity matrix elements follow skewed
% gaussian 

% last revised 1/5/2019

clear
clc

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
dataFolder = '../data/gcmi_skew/N100R20sig3_noReg';
figFolder = '../figures';

allFile = dir(fullfile(dataFolder,filesep,'*.mat'));
fileNameCell = {allFile.name}';

num = 40;   % repeat of simulation
% this set of data happened to with the same alpha
% str_marker = '(?<= *_sig)[2]+(?=_)';   %folder with this string contains the data we need
str_marker = 'N100_R20';   %folder with this string contains the data we need
FIND = @(str) cellfun(@(c) ~isempty(c), regexp(fileNameCell,str,'once'));
files = fileNameCell(FIND(str_marker));

% threshold to register the active elements
LB = -6; % threshold of non-zero elements
UB = 4;  % cut off of the ultrasensitive elements

s1 = '(?<= *N)[\d.]+(?=_)';
s2 = '(?<= *_R)[\d]+(?=_)';
s3 = '(?<= *_S)[\d]+(?=_)';
s5 = '(?<= *sig)[\d.]+(?=_)';
s6 = '(?<= *alp)[-\d.]+(?=_)';


GaussSk = nan(length(files),num);
GaussKt = nan(length(files),num);
allSk = nan(length(files),num);
allKt = nan(length(files),num);
spW = nan(length(files),num); % store the sparsity information

% select two example matrices for future use
alp_slt = [-2,2];  %select data from these two alp
W_ix = [2,2];
seltW = cell(2,1);

alp = zeros(length(files),1);  % store the skewness parameters
for i0 = 1:length(files)
    load(fullfile(dataFolder,files{i0}))

    alp(i0) = str2num(char(regexp(files{i0},s6,'match')));    
    
    for j0 = 1:size(allMat,2)
        temp = allMat(:,j0);
        sensiW = temp(temp > LB & temp <=UB);  % here the upper limit should be checked,2.5
        spW(i0,j0) = length(temp(temp>LB))/length(temp);
        allSk(i0,j0) = skewness(sensiW);
        allKt(i0,j0) = 	kurtosis(sensiW)-3;  
    
    % generate Gaussian data
        gd = std(sensiW)*randn(length(sensiW),1) + mean(sensiW);
        GaussSk(i0,j0) = skewness(gd);
        GaussKt(i0,j0) = kurtosis(gd)-3;
        
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
gaussSk = [nanmean(GaussSk,2),nanstd(GaussSk,0,2)];
simKt = [nanmean(allKt,2),nanstd(allKt,0,2)];
GaussKt = [nanmean(GaussKt,2),nanstd(GaussKt,0,2)];

% plot how skewness change with beta
figure
hold on
[~,inx] = sort(round(alp*100));
errorbar(alp(inx),simSk(inx,1),simSk(inx,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
% errorbar(alp(inx),gaussSk(inx,1),gaussSk(inx,2),'o-','MarkerSize',12,...
%     'MarkerFaceColor',Gr(7,:),'Color',Gr(7,:),'LineWidth',2,'CapSize',0)
% lg = legend('simulation','Gaussian');
% set(lg,'FontSize',20)
ylim([-0.7,1.8])
xlabel('$\hat{\alpha}$','interpreter','latex')
ylabel('Skewness')
box on
prefix = ['Skewness_SkewGauss_alp_N50M13_',date];
%saveas(gcf,[figFolder,filesep,prefix,'.fig'])
%print('-depsc',[figFolder,filesep,prefix,'.eps'])


% ========================================================================
% Fit a skewed Gaussian distribution
% ==========================================================
gaussian = @(x) (1/sqrt((2*pi))*exp(-x.^2/2));
skewStand = @(x,alpha) 2*gaussian(x).*normcdf(alpha*x);
skewPdf =  @(x,xi,omi,alp) 2/omi*gaussian((x-xi)/omi).*normcdf(alp*(x-xi)/omi);
fun = @(x,xdata) normcdf((xdata-x(1))/x(2)) - 2*myOwenT((xdata-x(1))/x(2),x(3));


fitedParm = zeros(2,3);
lb = [-10,0.1,-10];
ub = [10,10,10];
alpSelt = [-2,2];

for i0 = 1:2
    [Y,X] = ecdf(seltW{i0}); %the ecdf
    x0 = [mean(seltW{i0}),std(seltW{i0}),0];
    optParam = lsqcurvefit(fun,x0,X(2:end),Y(2:end),lb,ub);
    fitedParm(i0,:) = optParam;
    
    
    % plot
    figure
    figureSize = [0 0 5 4.5];
    set(gcf,'Units','inches','Position',figureSize,'PaperPositionMode','auto');
    hold on
    
    h1 = histogram(seltW{i0},25,'Normalization','pdf','FaceColor',lBu,'FaceAlpha',...
    0.4,'EdgeColor','none');
    stairs([h1.BinEdges(1),h1.BinEdges,h1.BinEdges(end)],...
    [0,h1.Values,h1.Values(end),0],'Color',dpBu,'LineWidth',2)
    Y2 = skewPdf(X,optParam(1),optParam(2),optParam(3));
    plot(X,Y2,'LineWidth',4,'Color',RdBu(2,:))
    hold off
    lg = legend('Experiment','Skewed Fit','Location','northwest');
    set(lg,'FontSize',16)
    legend boxoff
    box on
    xlabel('$\ln(w)$','Interpreter','latex')
    ylabel('pdf','Interpreter','latex')
%     prefix = ['SI_skew_alp',num2str(alpSelt(i0)),'_',date];
%     saveas(gcf,[figFolder,filesep,prefix,'.fig'])
%     print('-painters','-dpdf',[figFolder,filesep,prefix,'.pdf'])
end


% ===================================================
%  plot the sparsity of the W
% ===================================================
figure
[~,inx] = sort(alp);
errorbar(alp(inx),nanmean(spW(inx),2),nanstd(spW(inx),0,2),'o-','MarkerSize',12,...
    'MarkerFaceColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
ylim([0.4,0.7])
xlabel('$\hat{\alpha}$','interpreter','latex')
ylabel('$\rho_w$','interpreter','latex')
box on
% prefix = ['sparW_SkewGauss_alp_N50M13_',date];
% saveas(gcf,[figFolder,filesep,prefix,'.fig'])
% print('-depsc',[figFolder,filesep,prefix,'.eps'])

