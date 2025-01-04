function [optParam,fitSkew] = fitSkewedGauss(data,varargin)
% this function fit a skewed gaussian distribution of data, and return a
% parameter struct
% data      a vector
% 
% define the essential paramters
gaussian = @(x) (1/sqrt((2*pi))*exp(-x.^2/2));
skewStand = @(x,alpha) 2*gaussian(x).*normcdf(alpha*x);
skewPdf =  @(x,xi,omi,alp) 2/omi*gaussian((x-xi)/omi).*normcdf(alp*(x-xi)/omi);
fun = @(x,xdata) normcdf((xdata-x(1))/x(2)) - 2*myOwenT((xdata-x(1))/x(2),x(3));

fitedParm = zeros(2,3);
lb = [-10,0.1,-10];
ub = [10,10,10];

[Y,X] = ecdf(data); %the ecdf
x0 = [mean(data),std(data),0];
optParam = lsqcurvefit(fun,x0,X(2:end),Y(2:end),lb,ub);

% calculate the skewness of fitted skewed Gaussian distribution
delta = optParam(3)/sqrt((1+optParam(3)^2));
fitSkew = (4-pi)*delta^3/2*(2/(pi - 2*delta^2))^(3/2);

% plot the figure or not
if nargin > 1
    pltFlag = varargin{1}; % a string to specify whether plot or 
    if strcmp(pltFlag,'plot')
        
        % define some colors
        lBu = [96,166,223]/255; %light blue
        RdBu = brewermap(11,'RdBu');   % red and blue
        dpBu = [63,114,183]/255; % deep blue

        
        figure
        hold on
        h1 = histogram(data,25,'Normalization','pdf','FaceColor',lBu,'FaceAlpha',...
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
    end
end

end