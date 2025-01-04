% This program plot the heatmap of the sensitivity matrix for partial
% inhibition and how the inhibitory fraction reduce to the zero elements in
% the excitation-only situation

clear
close all

%% graphics settings
defaultGraphicsSetttings
%define some colors using brewermap
RdBu = brewermap(11,'RdBu');   % red and blue
Bu = brewermap(11,'Blues');    % blues

%% load the data
dFolder = '../data';
outFolder = '../figures';

file = 'gcmi_partInhi_summData_N50M20sp2sig2frac_0.5_25-Jun-2019.mat';
load(fullfile(dFolder,file));

%% basic parameters
N = 50;
M = 20;
sp =2;
frac = 0.5;

rmin = 0.03;  % lower limit of basal activity
rmax = 0.4;   % upper limit of basal activity
NI = ceil(M*frac); %number of ORNs that have basal activity
allr0 = rmin + (rmax-rmin)/NI*(0:(NI-1));

%% How the fraction of inhibitory interaction change with the basal activty
figure
hold on
errorbar(allr0',summData.meanInhiRatio,summData.stdInhiRatio,'o-','MarkerSize',12,'MarkerFaceColor',Bu(9,:),...
    'MarkerEdgeColor',Bu(9,:),'Color',Bu(9,:),'LineWidth',2,'CapSize',0)
plot([0;0.5],[0;0.5],'k--','LineWidth',2)
errorbar(0,summData.menaRho,summData.stdRho,'o-','MarkerSize',12,'MarkerFaceColor',RdBu(2,:),...
    'MarkerEdgeColor',RdBu(2,:),'Color',RdBu(2,:),'LineWidth',2,'CapSize',0)
hold off
box on
xlim([-0.02,0.5])
ylim([-0.02,0.5])
xlabel('$r_0$','Interpreter','latex','FontSize',28)
ylabel('inhibitory fraction','FontSize',28)
set(gca,'XTick',0:0.2:1,'FontSize',24,'LineWidth',1.5,'XScale','linear')
figNamePref = ['gcmi_partialInhi_ratio_basal_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])


%% heatmap of a typical sensitivity matrix
allBlue = brewermap(128,'Blues');  % here we make the color following "common rule"
allRd = brewermap(128,'Reds');

% select one matrix to plot
inx = 1;
inhiW = reshape(log(allInhiW(:,inx)),[],N);
exciW = reshape(log(allExciW(:,inx)),[],N);
sign = reshape(allSign(:,inx),[],N);  % sign of the inhibitory part

[exciInx_inhiW_Y,exciInx_inhiW_X] = find(sign >0);  % index of excitation in inhibitory part of the matrix
[inhiInx_inhiW_Y,inhiInx_inhiW_X] = find(sign <0);  % index of inhibitione
inhiW_exci = inhiW(sign >0);
inhiW_inhi = inhiW(sign <0);


[exciInx_exciW_Y,exciInx_exciW_X] = find(exciW >-8); 
exciW_exci = exciW(exciW>-8);


%map value to color, excitatory
LB = -8;
UB = 5;
exciW_exci(exciW_exci < LB) = LB;
exciW_exci(exciW_exci > UB) = UB;

% map value to color, inhibitory
LB_inhi = -4;
UB_inhi = 1;
inhiW_inhi(inhiW_inhi < LB_inhi) = LB_inhi;
inhiW_inhi(inhiW_inhi > UB_inhi) = UB_inhi;

inhiW_exci(inhiW_exci < LB) = LB;
inhiW_exci(inhiW_exci > UB) = UB;

% colorExci = zeros(length(we),3);
colorExci_exci = allBlue(round((exciW_exci - LB)/(UB-LB)*127)+1,:);
colorInhi_exci = allBlue(round((inhiW_exci - LB)/(UB-LB)*127)+1,:);

colorInhi = allRd(round((inhiW_inhi - LB_inhi)/(UB_inhi- LB_inhi)*127)+1,:);

maxMarker = 20;
allSizeExci = (exciW_exci - LB)/(UB-LB)*(maxMarker - 1)+1;
allSizeExci_inhW = (inhiW_exci - LB)/(UB-LB)*(maxMarker - 1)+1;
allSizeInhi = (inhiW_inhi - LB_inhi)/(UB_inhi-LB_inhi)*(maxMarker - 1)+1;

figure
figureSize = [0 0 15 7];
set(gcf,'Units','centimeters','Position',figureSize,...
'PaperPositionMode','auto');
hold on

% excitatory elements in inhibitory part
for i0 = 1:length(inhiW_exci)
    plot(exciInx_inhiW_X(i0)-0.5,exciInx_inhiW_Y(i0) + size(inhiW,1)-0.5,'o','MarkerSize',allSizeExci_inhW(i0),'MarkerFaceColor',...
        colorInhi_exci(i0,:),'MarkerEdgeColor',colorInhi_exci(i0,:),'LineWidth',0.5)
end

% excitatory elements excitaotry part
for i0 = 1:length(exciW_exci)
    plot(exciInx_exciW_X(i0) -0.5,exciInx_exciW_Y(i0)-0.5,...
        'o','MarkerSize',allSizeExci(i0),'MarkerFaceColor',...
        colorExci_exci(i0,:),'MarkerEdgeColor',colorExci_exci(i0,:),'LineWidth',0.5)
end

% inhibitory elements
for i0 = 1:length(inhiW_inhi)
     plot(inhiInx_inhiW_X(i0)-0.5,inhiInx_inhiW_Y(i0) + size(inhiW,1)-0.5,'o','MarkerSize',allSizeInhi(i0),'MarkerFaceColor',...
        colorInhi(i0,:),'MarkerEdgeColor',colorInhi(i0,:),'LineWidth',0.5)
end
hold off
set(gca,'XTick',10:10:50,'YTick',[5,10,15,20])
box on
%grid on
% set(gca,'XTick',[])
% daspect([1 1 1])
figNamePref = ['R2_ExciInhi_diffr0_heatmap_',date];
saveas(gcf,[outFolder,filesep,figNamePref,'.fig'])
print('-depsc',[outFolder,filesep,figNamePref,'.eps'])
