function plotbalance(X1, YMatrix1, YMatrix2,time)
%CREATEFIGURE(X1,YMATRIX1,YMATRIX2)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data
%  YMATRIX2:  matrix of y data

%  Auto-generated by MATLAB on 12-Mar-2013 14:17:46

% Create figure
figure2 = figure;

% Create subplot
subplot1 = subplot(1,2,1,'Parent',figure2,'YScale','log','YMinorTick','on',...
    'FontSize',16,...
    'FontName','Hoefler Text');
box(subplot1,'on');
hold(subplot1,'all');


% Create multiple lines using matrix input to semilogy
semilogy1 = semilogy(X1,YMatrix1,'Parent',subplot1,'LineWidth',2);
set(semilogy1(1),'Color',[1 0 0],'DisplayName',sprintf('T = %d',time));
set(semilogy1(2),'Color',[0 1 0],'DisplayName','T = 0');
set(semilogy1(3),'Marker','o','LineStyle','none','Color',[0 0 1],...
    'DisplayName','Pack',...
    'LineWidth',0.5);

% Create xlabel
xlabel('Diameter','FontSize',16,'FontName','Hoefler Text');

% Create ylabel
ylabel('Fractional Area','FontSize',16,'FontName','Hoefler Text');

% Create title
title('A(d)','FontSize',16,...
    'FontName','Hoefler Text');

% Create legend
legend(subplot1,'show');

% Create subplot
subplot2 = subplot(1,2,2,'Parent',figure2,'XScale','log','XMinorTick','on',...
    'FontSize',16,...
    'FontName','Hoefler Text');
box(subplot2,'on');
hold(subplot2,'all');

% Create multiple lines using matrix input to semilogx
semilogx1 = semilogx(X1,YMatrix2,'Parent',subplot2,'LineWidth',2);
set(semilogx1(1),'Color',[0 0.498039215803146 0],'DisplayName','Melt');
set(semilogx1(2),'Color',[1 0 0],'DisplayName','Redist');
set(semilogx1(3),'Color',[0 1 0],'DisplayName','Advect');
set(semilogx1(4),'Color',[0 0.749019622802734 0.749019622802734],...
    'DisplayName','Swell');

% Create xlabel
xlabel('Diameter','FontSize',16,'FontName','Hoefler Text');

% Create ylabel
ylabel('\Delta A','FontSize',16,'FontName','Hoefler Text');

% Create title
title({'Steady State Balance'},'FontSize',16,...
    'FontName','Hoefler Text');

% Create legend
legend1 = legend(subplot2,'show');
set(legend1,...
    'Position',[0.754894586518775 0.712786523059929 0.0749436513899324 0.146322378716745]);
if saveplots == 1
hgexport(figure2,'Output/Balances.fig')
hgexport(figure2,'Output/Balances.eps')
end