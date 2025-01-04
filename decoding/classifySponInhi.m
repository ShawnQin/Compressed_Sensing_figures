function summData = classifySponInhi(numTypes, NP, HiddenSize, noiseSig, dS, varargin)
% this program comparing the classification accuracy for different basal
% activity, modified from program "classifyRndPNtoKC.m"
% We use the optimized sensitivity matrix derived for different r0
% the input is just the OR response pattern and related tag (odor valence), 
% to incorporate the expansion from PN to KC, we limit the matrix from PN
% to KC be fixed, but only alow the learning from KC to MBON
% we also use global inhibition at the KC level to generate sparse
% representation. This is achieved through "normalization" at the KC level
% the last classifier is smiple LDA

% After reading Abbott's 2017 nature papper, add truncated Gaussian noise

% last revised on 12/29/2018
% ========================================================================

% numTypes      number of categories to be classified, default 2
% NP            number of odor clsuters, or centroids
% HiddenSize    number of KC cells, default 200
% noiseSig      standard deviation of noise at the ORN response level,0.05
% dS            the relative standard deviation within a cluster, 0.1
% varargin{1}   extral input, such as the random number seed

%parpool(4)

% parse input parameter, whether use random number seed
if nargin > 5
    rngSeed = varargin{1};
    rng(rngSeed)
else
    rng('shuffle')
end

% specify the input data type, random or cluster
inputType = 'cluster';   % generating clustered odor clouds
classMethod = 'LDA';     % classification methods, can be logistic,svm,LDA,or multiClass
KCmethod = 'relu';       % activation function at KC level, can be relu, binary, linear

% add the search path, this is only useful when runing on the cluster
addpath(genpath('/lustre1/tangc_pkuhpc/ssqin/olfactoryCoding/optiInterMatrix/myFunc'))

% the spontaneous activity that we have used in our previous search
spAll = [0.05, 0.08, 0.1:0.02:0.18,0.25,0.3:0.1:0.7,0.82,0.84,0.88,0.9,0.95];

% ====== parameters of the network======================

% we use optimized sensitivity matrix with both exciation and inhibition,
% here we need to loard previous data
dFolder = './data/inhi_N50M10S2sig2_1013';

hiddenLayerSize = HiddenSize;      %set the hidden layer size
nOdor = 50;
nRecep = 10;
spar = 2;
sig = 2;
% ========================================================

% the parameters of the input
sparsity = spar;
param.nRecep = nRecep;
param.nOdor = nOdor;
param.lSig = sig;
param.lMu = 0;
param.nSamp = 1e4;
param.sparsity = true;
param.spRatio = sparsity/param.nOdor;
param.dType = 'lognorm';
param.eig = [];
param.corr = false;
param.regularize = false;
param.spType = 'absolute';          %type of sparisty, could be absolute or average
param.noiseSig = noiseSig;              % introduce small noise on the response
        
param.nPattern = NP;        % number of pattern, centroids
param.withinPattSamp = 50;  % number of samples in each pattern
param.patternStd = dS;
param.sigw = 0.1;           %the distribution of PN-KC weight

spPoject = 6;             % number of PN  each KC receive
trainRatio = 0.8;         % fraction of data used to train
KCsp = 0.1;               % KC response sparsity

% initialization
errRate = nan(length(spAll),40);  % store the results
for i0 = 1:length(spAll)    
    r0 = spAll(i0);
    dFile = ['gcmiInhi_N50_R10_S2_sig2_alp',num2str(r0),'_frac_2018-10-13.mat'];
    load(fullfile(dFolder,dFile))
    
    L = length(allfmin);
    for j0 = 1:L

        % direct read the matrix and sign from input     
        W = allMat(:,j0);
        Sign = allSign(:,j0);
        
        if strcmp(inputType,'cluster')
            [trainData, targets] = genOdorClouds(param,spar,numTypes,1);
            NS = param.nPattern*param.withinPattSamp;
            param.nSamp = NS;
        else
            error('the input data type can be only  random or cluster!')
        end
        
        % generate noise on the input or not, using random truncated
        % Gaussian noise
        if ~isempty(param.noiseSig)
            pd = makedist('Normal');
            t = truncate(pd,0,inf);
            noiseAdd = param.noiseSig*random(t,param.nRecep,NS);
        else
            noiseAdd = 0;
        end
        
        % set the value of inhibitory and excitatory elements
        w_in = zeros([param.nRecep,param.nOdor]);
        w_ac = w_in;
        w_in(Sign==-1) = W(Sign==-1);
        w_ac(Sign==1) = W(Sign==1);
        
        alpha = (1-r0)/r0;  % determine the basal activity
        resp = (1+alpha*(1+w_in*trainData)./(1+w_ac*trainData)).^(-1) + noiseAdd;
        
        % define the matrix from PN to KC, a random projection matrix
        % response at the KC level,introduce global inhibition     
        W1 = PNtoKC(param.nRecep,hiddenLayerSize,spPoject,'gaussian');

        rhat = mean(resp,2)/norm(mean(resp,2));      % average vector
        normResp = W1*(resp - rhat*(rhat'*resp));    % normalization, from Abbott 2010 PNAS
        
        % generate sparse representation at KC level
        KCout = genKCresp(normResp,KCmethod,1-KCsp);
        
        % remove useless KCs (those have all 0s)
        KCsel = KCout(sum(KCout > 0,2) > 0, :);
        
        % divide the data into training and testing set
        allLabel = targets(1,:);  % binary label
        cvp = cvpartition(param.nSamp,'Holdout',1-trainRatio);
        idxTrn = training(cvp);   % Training set indices
        idxTest = test(cvp);      % Test set indices
       
        
        trainSet = KCsel(:,idxTrn)';
        trainLabel = allLabel(idxTrn)';
        testSet = KCsel(:,idxTest)';
        testLabel = allLabel(idxTest)';
        
        % test the performance on testing set
        if strcmp(classMethod,'svm')
            % training classifier
            SVMModel = fitcsvm(trainSet,trainLabel,'KernelFunction','linear','Standardize',true,...
            'ClassNames',[0,1]);
       
            % test performance
            [label,score] = predict(SVMModel,testSet);
            errRate(i0,j0) = sum(abs(testLabel - label))/length(testLabel);
        elseif strcmp(classMethod,'logistic')
            % training
            mdl = fitglm(trainSet,trainLabel,'Distribution','binomial','Link','logit');
            ypred = predict(mdl,testSet);
            errRate(i0,j0) = sum(abs((ypred > 0.5) - testLabel'))/length(testSet);

        elseif strcmp(classMethod,'multiClass')
%             [allLabel,~] = find(targets==1);
            [B,dev,stats] = mnrfit(trainSet,trainLabel);
            pihat = mnrval(B,testSet);
            
            
        elseif strcmp(classMethod,'LDA')
            MdlLinear = fitcdiscr(trainSet,trainLabel,'DiscrimType','pseudolinear');
            [label,score,cost] = predict(MdlLinear,testSet);
            confMat = confusionmat(testLabel,label);
            errRate(i0,j0) = 1- sum(diag(confMat))/length(label);
        end       

    end
end


% put all the relevant data into a struct
summData = struct('hiddenSize',hiddenLayerSize,'noiseSig',noiseSig,...
   'errorRate',errRate,'allSp',spAll,'group',numTypes);

% save the final trainging results
sName = ['LDA_inhi_r0_','_N',num2str(nOdor),'M',num2str(nRecep),'H',...
    num2str(hiddenLayerSize),'_sp',num2str(spar),'_G',num2str(numTypes),'_ns',...
    num2str(noiseSig),'_',date,'.mat'];
save(sName,'summData');
end

%% generate KC response 
function KCout = genKCresp(resp,method,varargin)
% return the sparse KC output
% resp should be the total input to KC, linear function of PN/OSN
% method is a string, specified in "wta",'relu','linear'


% winner takes all
if strcmp(method,'wta')  %winner takes all
    KCout = zeros(size(resp,1),size(resp,2));
    if nargin > 2      
        thd = varargin{1}; %threhold for the rectified linear output, 0 ~1
    else
        thd = 0.95;         % default, only select the 5% strongest response
    end    
    
    tempThd = thd;
    for i0 = 1:size(resp,1)
        [F,X] = ecdf(resp(i0,:));
        inx = find(F>tempThd,1,'first');
        thd = X(inx);
        KCout(i0,resp(i0,:)>thd) = 1; 
    end
    
elseif strcmp(method,'relu')
    if nargin > 2      
        thd = varargin{1};   %threhold for the rectified linear output, 0 ~1
    else
        thd = 0.95;          %default, only select the 5% strongest response
    end    
        % tune the threshold such that only 10% KC response
    [F,X] = ecdf(resp(:));
    inx = find(F > thd,1,'first');
    thd = X(inx);
    
    KCout = max(0,resp - thd);
        
elseif strcmp(method,'linear')
    KCout = resp;       %just return the original data
else
    error('method not support yet, has to be wta, relu, or linear!')
end
end


function plotFigure(summData)

% folder to save the figure
outFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';

% default graphic setting
defaultGraphicsSetttings

% color
Bu = brewermap(11,'Blues');    % blues

errorbar(summData.allSp',mean(summData.errorRate,2),std(summData.errorRate,0,2),...
    'o-','MarkerSize',12,'MarkerFaceColor',Bu(10,:),'Color',Bu(10,:),'LineWidth',2,...
    'CapSize',0)
lg = legend('N=100,M=10,sp=3,\sigma_c=2');
legend boxoff
set(lg,'FontSize',16)
xlabel('sparsity of W')
ylabel('classification error')

figNamePref = ['clasErr_PN_KCrand_N100M20_sp3_sig2_H500_diffSpW_logistic',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])
end

%% connection from PN to KC
function W = PNtoKC(nRecp,hiddenLayerSize,sp,weightType,varargin)
% nRecp                 number of receptors, PN
% hiddenLayerSize       number of KCs, an integer
% sp                    average number of PNs projected to KC
% weightType            could be gaussian, or uniform random, or binary
% varargin{1}           exact connection numbers each KC receive

% initialize the W
W = zeros(hiddenLayerSize,nRecp);
if nargin > 4
    nn = varargin{1};
    if nn > nRecp
        error('number of connection must be smaller than receptors!')
    end
    allConn = ones(hiddenLayerSize,1)*nn;
else
    allConn = min(nRecp,binornd(10,(sp-1)/11,[hiddenLayerSize,1]) + 1); % revised on 08/10/2018
end

% different connection strength
if strcmp(weightType,'binary')
    for i0 = 1:hiddenLayerSize
        W(i0,:) = ones(allConn(i0),1);
    end
elseif strcmp(weightType,'gaussian')
    % truncated Gaussian noise
    pd = makedist('Normal');
    t = truncate(pd,-2,2);
    for i0 = 1:hiddenLayerSize
        W(i0,randperm(nRecp,allConn(i0))) = (random(t,allConn(i0),1) + 2)/4;
    end
elseif strcmp(weightType,'uniform')
    for i0 = 1:hiddenLayerSize
        W(i0,:) = rand(allConn(i0),1);
    end
end
end

%% generate centroid odor mixture that are more uniform in the space
function [trainData, targets] = genOdorClouds(param,spar,nType,nLabel,varargin)
% param      a struct specify all the parameters needed
% spar       sparsity of odor
% nType      number of clusters
% nLabel     how to specify the label, scaler or a vector
% varargin{1}  spcecify the centroid concentration type, identical, or
% random

if nargin > 4
    cType = varargin{1};
else
    cType = 'random';       % default
end

NP = param.nPattern;        % number of pattern, centroids
NS = param.withinPattSamp;  % number of samples in each pattern
ds = param.patternStd;      % variation level within pattern

% set the centroid, which is chose uniformly cover the space
% step 1, get the index
inx = zeros(NP,spar);
for i0 = 1:NP
    inx(i0,:) = randperm(param.nOdor,spar);
end

% step 2, set the concentration
centroid = zeros(param.nOdor,NP);
if strcmp(cType,'identical')
    for i0 = 1:NP
        centroid(inx(i0,:),i0) = 1;
    end
elseif strcmp(cType,'random')
    for i0 = 1:NP
        centroid(inx(i0,:),i0) = exp(randn(1,spar)*param.lSig);
    end
    
end

% step 3 randomly assign a label to these centroid
if nLabel > 1
centroidLabel = zeros(nType,NP);
eachSamp = ceil(NP/nType);
for i0 = 1:nType-1
    centroidLabel((i0-1)*eachSamp+1 : i0*eachSamp,i0) = 1;
end

centroidLabel((nType-1)*eachSamp+1 : NP,nType) = 1;
else
    centroidLabel = randi(nType,1,NP);
end
    
% step 4, generate clouds around centroids
trainData = zeros(param.nOdor,NP*NS);


% there are different ways to store the labels depending on different
% decoder used in MATLAB, such as LDA, SVM
if nLabel > 1    % only need a scalar
    targets = zeros(nType,NP*NS);
    for i0 = 1:NP
        trainData(centroid(:,i0) > 0,(i0-1)*NS+1:1:i0*NS) = (1+ds*randn(spar,NS)).*centroid(centroid(:,i0) > 0,i0);
        targets(:,(i0-1)*NS+1:1:i0*NS) = repmat(centroidLabel(:,i0),1,NS);
    end
else
    trainData = zeros(param.nOdor,NP*NS);
    targets = zeros(1,NP*NS);
    
    for i0 = 1:NP
        trainData(centroid(:,i0) > 0,(i0-1)*NS+1:1:i0*NS) = exp((1+ds*randn(spar,NS))).*centroid(centroid(:,i0) > 0,i0);
        targets((i0-1)*NS+1:1:i0*NS) = repmat(centroidLabel(i0),1,NS);
    end
end

end 