%% Figure 1C nonlinear ORN
% saveFolder = '../decoding/data/reconsFig6';
saveFolder = '../figures';


C = 10.^(-2:0.01:2);
Km = 1;
n = 1;
f = @(Km,C) (C/Km).^n./(1+(C/Km).^n);

c1 = 1/9;  % start, 10%
c2 = 9;    % end, 90% 
Wid = 8;   % figure width
Hei = Wid*0.85; % height of figure

figure
set(gcf,'Units','centimeters','Position',[0,0,Wid,Hei])
hold on
plot(C,f(Km,C),'k','LineWidth',3)
plot([c1,c1],[0,1],'k--','LineWidth',1)
plot([c2,c2],[0,1],'k--','LineWidth',1)
hold off
box on
set(gca,'xscale','log','XTick',10.^(-2:1:2),'LineWidth',1)
xlabel('$\sum_{j=1}^NW_{ij}c_j$','Interpreter','latex')
ylabel('$(r_i)$','Interpreter','latex')
% 
% prefix = ['fig1_nonlinearORN_h1',date];
% saveas(gcf,[saveFolder,filesep,prefix,'.fig'])
% print('-depsc',[saveFolder,filesep,prefix,'.eps'])