function createfigure(X1, YMatrix1)
%CREATEFIGURE(X1,YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 21-Mar-2013 20:18:04

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
    'XScale','log',...
    'XMinorTick','on',...
    'FontSize',24,...
    'FontName','Hoefler Text');
box(axes1,'on');
hold(axes1,'all');
2
% Create multiple lines using matrix input to loglog
loglog1 = loglog(X1,YMatrix1,'LineWidth',2);
set(loglog1(1),'Color',[1 0 0],'DisplayName','Long-Term Distribution');
set(loglog1(2),'Color',[0 1 0],'DisplayName','Initial Distribution');

% Create xlabel
xlabel('Diameter','FontSize',16,'FontName','Hoefler Text');

% Create ylabel
ylabel({'Area Distribution',''},'FontSize',16,'FontName','Hoefler Text');
