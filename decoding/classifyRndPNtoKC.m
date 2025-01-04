% this program comparing the classification accuracy for random odor-Or
% interaction and optimized, maximum entropy encoding
% the input is just the Or response pattern and related tag (odor valence), 
% to incorporate the expansion from PN to KC, we limit the matrix from PN
% to KC fixed, but only alow the learning from KC to MBON
% we also use global inhibition at the KC level to generate sparse
% representation
% logistic regression model is used
% After reading Abbott's 2017 nature papper, add truncated Gaussian noise

% this script is modified from "classifier"
% last revised on 09/23/2018

%%
function summData = classifyRndPNtoKC(nOdor,nRecp,spar,sig,numTypes,NP,HiddenSize, noiseSig,dS, varargin)
% nOdor         number of odorants, an integer
% nRecp         number of receptors, an integer
% sig           odorant concentration, defaut 2
% numTypes      groups of lables, an interger larger than 1
% NP            number of patterns
% varargin{1}   extral input, such as the random number seed

%parpool(4)

% parse input parameter, whether use random number seed
if nargin > 9
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
spAll = 0.1:0.05:1;     %the range of sparsity of W    

% random number generator seeds

% ====== parameters of the network======================
hiddenLayerSize = HiddenSize;      %set the hidden layer size
L = 40;                     %repeats
% ======================================================


% =======================================================
% parameters on the measurement matrix, should be adjusted according to the
% number of odorant and ORN
% [muW,sigW,rhoW,fmin] = selectMeanSig(nOdor,spar);
muW = -1;
sigW = 2;
% ========================================================

% the parameters of the input
sparsity = spar;
param.nRecep = nRecp;
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
errRate = zeros(length(spAll),L);  % store the results
for j0 = 1:L
    for i0 = 1:length(spAll)     
        sp = spAll(i0);  % input sparsity
        
        w = normrnd(muW,sigW,[nRecp,nOdor]);
        w_zero = rand(nRecp,nOdor);
        w(w_zero>sp) = -inf;
        
        if strcmp(inputType,'random')
            eigVal = specifyEig(param.nOdor,param.eig);
            corrCoefMat = randCorrCoef('buildin',eigVal);
            trainData = genTrainData(param,corrCoefMat);
            NS = param.nSamp;
            targets = zeros(2,param.nSamp);  
            for k0 = 1:param.nSamp
                inx = find(trainData(:,k0) > 0);
                if all(inx <= param.nOdor/2)
                    targets(1,k0) = 1;
                elseif all(inx > param.nOdor/2)
                    targets(2,k0) = 1;
                else
                    if trainData(inx(1),k0) > trainData(inx(2),k0)
                        targets(1,k0) = 1;
                    else
                        targets(2,k0) = 1;
                    end
                end
              
            end
        elseif strcmp(inputType,'cluster')
            [trainData, targets] = genOdorClouds(param,spar,numTypes,1);
            NS = param.nPattern*param.withinPattSamp;
            param.nSamp = NS;
        else
            error('the input data type can be only  random or cluster!')
        end
        
        % add some noise on the input or not
        % random truncated Gaussian noise
        pd = makedist('Normal');
        t = truncate(pd,0,inf);
        if ~isempty(param.noiseSig)
            resp = reshape(exp(w),nRecp,nOdor)*trainData./(1+reshape(exp(w),nRecp,nOdor)*trainData)...
                + param.noiseSig*random(t,param.nRecep,NS);
        else
            resp = reshape(exp(w),nRecp,nOdor)*trainData./(1+reshape(exp(w),nRecp,nOdor)*trainData);
        end
        
        % define the matrix from PN to KC, a random projection matrix
        % response at the KC level,introduce global inhibition
        
%         W1 = zeros(hiddenLayerSize,nRecp);
        W1 = PNtoKC(nRecp,hiddenLayerSize,spPoject,'gaussian');
%         W1(randperm(hiddenLayerSize*nRecp,round(hiddenLayerSize*nRecp*spPoject))) ...
%             = normrnd(1,param.sigw,round(hiddenLayerSize*nRecp*spPoject),1);
%         W1(randperm(hiddenLayerSize*nRecp,round(hiddenLayerSize*nRecp*spPoject))) ...
%             = rand(round(hiddenLayerSize*nRecp*spPoject),1);
        rhat = mean(resp,2)/norm(mean(resp,2));
        normResp = W1*(resp - rhat*(rhat'*resp));  % from Abbott 2010 PNAS
%         normResp = (W1 - mean(W1(W1 > 0))*spPoject)*resp;  % Abbott 2017 Neuron
%         KCinput = W1*resp;
        
        % generate sparse representation at KC level
        KCout = genKCresp(normResp,KCmethod,1-KCsp);
        
        % remove useless KCs
        KCsel = KCout(sum(KCout > 0,2) > 0, :);
% %         KCout = KCout(sum(KCout > 0,2) > 0,:); % index of useless KC
%         KCsel = [];
%         for k0 = 1:size(KCout,1)
%             temp = KCout(k0,KCout(k0,:)>0);
%             if range(temp) > 0.01
%                 KCsel = [KCsel;KCout(k0,:)];
%             end
%         end
%         KCsel = KCsel(:,sum(KCsel > 0,2));
%         KCout = KCout(:,sum(KCout > 0,1) > 0); % index of useless KC
        
        % divide the data into training and testing set
        allLabel = targets(1,:);  % binary label
        cvp = cvpartition(param.nSamp,'Holdout',1-trainRatio);
        idxTrn = training(cvp); % Training set indices
        idxTest = test(cvp);    % Test set indices
        
        
        % get out rows with predictor that is
%         temp1 = KCout(:,idxTrn)';
%         ix1 = sum(temp1>0,1)>0;
%         temp1 = temp1(:,ix1);
%         
%         temp2 = KCout(:,idxTest)';
%         ix2 = sum(temp2>0,1)>0;
%         temp2 = temp2(:,ix2);
        
        trainSet = KCsel(:,idxTrn)';
        trainLabel = allLabel(idxTrn)';
%         tblTrn = array2table(temp1);
%         tblTrn.group = allLabel(idxTrn)';
        testSet = KCsel(:,idxTest)';
        testLabel = allLabel(idxTest)';
        
%         trainInx = randperm(param.nSamp,round(param.nSamp*trainRatio));
%         testInx = setdiff(1:1:param.nSamp,trainInx);
%         trainSet = KCsel(:,trainInx);
%         testSet = KCsel(:,testInx);
%         trainLabel = allLabel(trainInx);
%         testLabel = allLabel(testInx);
        
        % this part is used for debuging
%         myTrain = array2table(KCsel');
%         myTrain.Group = allLabel';
%         save('testClass.mat','myTrain');
       
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
%             [model, llh] = logitBin(trainSet, trainLabel, 0, 0.05);
%             [y, p] = logitBinPred(model, testSet);
%             sum(abs(y-testLabel))
%             B1 = logisticReg(trainSet',trainLabel');
%             ypred = glmval(B1,X,'logit');
%             errRate(i0,j0) = sum(abs((ypred > 0.5) - testLabel'))/length(testSet);
        elseif strcmp(classMethod,'multiClass')
%             [allLabel,~] = find(targets==1);
            [B,dev,stats] = mnrfit(trainSet,trainLabel);
            pihat = mnrval(B,testSet);
            
            
        elseif strcmp(classMethod,'LDA')
            MdlLinear = fitcdiscr(trainSet,trainLabel,'DiscrimType','pseudolinear');
            [label,score,cost] = predict(MdlLinear,testSet);
            confMat = confusionmat(testLabel,label);
            errRate(i0,j0) = 1- sum(diag(confMat))/length(label);
%             MdlLinear.ClassNames([2 3])
%             K = MdlLinear.Coeffs(2,3).Const;
%             L = MdlLinear.Coeffs(2,3).Linear;
        end       

    end
end


summData = struct('hiddenSize',hiddenLayerSize,'noiseSig',noiseSig,...
   'errorRate',errRate,'allSp',spAll,'group',numTypes);

% save the final trainging results
sName = ['classify_LDA','_N',num2str(nOdor),'M',num2str(nRecp),'H',...
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
if strcmp(method,'wta')
    KCout = zeros(size(resp,1),size(resp,2));
    if nargin > 2      
        thd = varargin{1}; %threhold for the rectified linear output, 0 ~1
    else
        thd = 0.95;   % default only select the 5% strongest response
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
        thd = varargin{1}; %threhold for the rectified linear output, 0 ~1
    else
        thd = 0.95;   % default only select the 5% strongest response
    end    
        % tune the threshold such that only 10% KC response
    [F,X] = ecdf(resp(:));
    inx = find(F>thd,1,'first');
    thd = X(inx);
    
%     MAXRESP = max(resp(:)) - thd;
%     KCout = max(0,resp - thd)/MAXRESP; %normalized to 0 ~ 1
    KCout = max(0,resp - thd);
        
elseif strcmp(method,'linear')
    KCout = resp;  %just return the original data
else
    error('method not support yet, has to be wta, relu, or linear!')
end
end
%% generate clustered input
function [trainData, targets] = genClusterInput(param,spar,numTypes,varargin)
% param      a struct specify all the parameters needed

        
NP = param.nPattern;       % number of pattern, centroids
NS = param.withinPattSamp;  % number of samples in each pattern
ds = param.patternStd;      % variation level within pattern


% set the centroid;
% trainCentroid = zeros(param.nOdor,NP);
eigVal = specifyEig(param.nOdor,param.eig);
corrCoefMat = randCorrCoef('buildin',eigVal);
trainCentroid = genTrainData(param,corrCoefMat);

% randomly assign label to these centroid
centroidLabel = zeros(numTypes,NP);
for i0 = 1:NP
    centroidLabel(randi(numTypes),i0) = 1;
end
% centroidLabel = binornd(1,1/numTypes,[1,NP]);
% centroidLabel = [centroidLabel;~centroidLabel];

% sample and assign label within each pattern
trainData = zeros(param.nOdor,NP*NS);
targets = zeros(numTypes,NP*NS);
for i0 = 1:NP
%     inx = find(trainCentroid(:,i0) > 0);
    trainData(trainCentroid(:,i0) > 0,(i0-1)*NS+1:1:i0*NS) = (1+ds*randn(spar,NS)).*trainCentroid(trainCentroid(:,i0) > 0,i0);
    targets(:,(i0-1)*NS+1:1:i0*NS) = repmat(centroidLabel(:,i0),1,NS);
end

% if multiple categories are specified
if nargin > 3
    classMethod = varargin{1};
    if strcmp(classMethod,'multiClass')
        [targets, ~] = find(targets == 1);
        targets = targets';  % transform it, to be consistent
    end
end



end

%% logistic regression with regularization
function B1 = logisticReg(X,Ybool)
[B,FitInfo] = lassoglm(X,Ybool,'binomial',...
    'NumLambda',25,'CV',10);
%%
lassoPlot(B,FitInfo,'PlotType','CV');
legend('show','Location','best') % show legend

indx = FitInfo.Index1SE;
B0 = B(:,indx);
nonzeros = sum(B0 ~= 0);
cnst = FitInfo.Intercept(indx);
B1 = [cnst;B0];

end


%% plot the sparsity dependent classification error
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

%% select the optimal parameters for sparse matrix
function [muW,sigW,rhoW,fmin] = selectMeanSig(N,sp)
% this function return the optimal parameter of parameterized W

% load the data
% dFoler = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';  % run local
dFoler = './';  % run on cluster
dName= fullfile(dFoler,'gcmi_distri_summData_16-Jul-2018.mat');
load(dName)
% allN = allN;
% allSp = allSp;
inx1 = find(allN ==N);
inx2 = find(allSp ==sp);

muW = meanW(inx1,inx2);
sigW = meanSigW(inx1,inx2);
rhoW = meanRho(inx1,inx2);
fmin = meanFmin(inx1,inx2);
end

function [trainData, targets] = genTwoSimple(param,spar,varargin)
% param      a struct specify all the parameters needed

        
NP = param.nPattern;       % number of pattern, centroids
NS = param.withinPattSamp;  % number of samples in each pattern
ds = param.patternStd;      % variation level within pattern


% set the centroid;
eigVal = specifyEig(param.nOdor,param.eig);
corrCoefMat = randCorrCoef('buildin',eigVal);
trainCentroid = genTrainData(param,corrCoefMat);

% assign the label based on the position of appeared odorants
centroidLabel = zeros(2,NP);
for i0 = 1:NP
    inx = find(trainCentroid(:,i0) > 0);
    if (all(inx <= round(param.nOdor/2))) 
        centroidLabel(1,i0) = 1;
    elseif (all(inx > round(param.nOdor/2)))
        centroidLabel(2,i0) = 1;
    else
        if trainCentroid(inx(1),i0) > trainCentroid(inx(2),i0)
            centroidLabel(1,i0) = 1;
        else
            centroidLabel(2,i0) = 1;
        end
    end
end
% centroidLabel = binornd(1,1/numTypes,[1,NP]);
% centroidLabel = [centroidLabel;~centroidLabel];

% sample and assign label within each pattern
trainData = zeros(param.nOdor,NP*NS);
targets = zeros(2,NP*NS);
for i0 = 1:NP
%     inx = find(trainCentroid(:,i0) > 0);
    trainData(trainCentroid(:,i0) > 0,(i0-1)*NS+1:1:i0*NS) = (1+ds*randn(spar,NS)).*trainCentroid(trainCentroid(:,i0) > 0,i0);
    targets(:,(i0-1)*NS+1:1:i0*NS) = repmat(centroidLabel(:,i0),1,NS);
end


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
    cType = 'random';   % default
end

        
NP = param.nPattern;        % number of pattern, centroids
NS = param.withinPattSamp;  % number of samples in each pattern
ds = param.patternStd;      % variation level within pattern


% set the centroid, which is chose uniformly cover the space
% step 1, get the index
% inx = ceil(lhsdesign(NP,spar)*param.nOdor);
inx = zeros(NP,spar);
for i0 = 1:NP
    inx(i0,:) = randperm(param.nOdor,spar);
end

% step 2, set the concentration, range -2sig ~ 2 sig
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

% step 3 randomly assign label to these centroid
if nLabel > 1
centroidLabel = zeros(nType,NP);
eachSamp = ceil(NP/nType);
for i0 = 1:nType-1
    centroidLabel((i0-1)*eachSamp+1 : i0*eachSamp,i0) = 1;
end
centroidLabel((nType-1)*eachSamp+1 : NP,nType) = 1;
%     for i0 = 1:NP-1
%         centroidLabel(randi(nType),i0) = 1;
%     end
else
    centroidLabel = randi(nType,1,NP);
end
    
% step 4, generate clouds around centroids
trainData = zeros(param.nOdor,NP*NS);

if nLabel > 1
    targets = zeros(nType,NP*NS);
    for i0 = 1:NP
%     inx = find(trainCentroid(:,i0) > 0);
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