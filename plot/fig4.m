%% Figure 4 Mean Field Theory with N = 2
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
dFolder = '../data';
saveFolder = '../figures';

% dFolder = '/Users/shan/Dropbox/olfactionProject/data/tempfigure';
% saveFolder = '/Users/shan/Documents/MATLAB/theoNeurosci/olfaction/figureData';
dName = 'data_fit_MFT.mat'; %more direct and accuarate, see SI
% dName = 'MFT_fit.mat';  % Yuhai's ansatz
load(fullfile(dFolder, dName))

inx = [1 3 6];
allsig = 1.1:0.2:2.9; % all the sigma

[Y1,Y2] = sort(all,2,'descend');
largest = zeros(size(all,1),2);
for i0 = 1:size(largest,1)
    largest(i0,1) = q_all(Y2(i0,1));
end
largest(:,2) = Y1(:,1);

% ==================================================
% plot three I2(m) with different sigma_;c
% ==================================================
figure
hold on
for i0 = 1:length(inx)
    plot(q_all',all(inx(i0),:)','Color',Bu(2+i0*3,:),'LineWidth',3)
    ah = gca;
    yl = ah.YLim;
    plot([largest(inx(i0),1);largest(inx(i0),1)],yl','--','LineWidth',2,'Color',Bu(2+i0*3,:))
end
ah.YLim = [0 5];
lg = legend('\sigma_c = 1.1','\sigma_c = 1.5','\sigma_c = 2.1');
set(lg,'FontSize',16)
legend boxoff
box on

xlabel('$\rho_w$','interpreter','latex')
ylabel('$I_2$','interpreter','latex')
hold off
% prefix = ['fig4_MFT_N2_entr_rhoW',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])

% ==================================================
% plot peak position of I2 with differnt sigma_c
% ==================================================
figure
% plot(allsig',largest(:,1),'o-','MarkerSize',12,'MarkerFaceColor',myOr,...
%     'MarkerEdgeColor','k','LineWidth',1)
plot(allsig',largest(:,1),'o-','MarkerSize',15,'MarkerFaceColor',Bu(10,:),...
    'MarkerEdgeColor',Bu(10,:),'LineWidth',2)
ylim([0.6 0.75])
xlabel('$\sigma_c$','Interpreter','latex')
ylabel('$\rho_w*$','Interpreter','latex')
% prefix = ['fig4_MFT_N2_rhow_sigc',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])