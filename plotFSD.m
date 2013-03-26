function plotFSD(D,A,A0,A1,stats,Save,Timer)

FSDmat = [A;A0;A1];

Relmag = [0*Timer' + 1;stats(:,3)';stats(:,5)';stats(:,7)'];
Mag = [stats(:,1)';stats(:,2)';stats(:,4)';stats(:,6)'];

size(FSDmat)
size(Mag)
size(Relmag)

%% Create figure
figure1 = figure;

%% 
% Create subplot
subplot1 = subplot(4,1,1,'Parent',figure1,'YScale','log','YMinorTick','on');
box(subplot1,'on');
hold(subplot1,'all');

%% 
% Create multiple lines using matrix input to semilogy
semilogy1 = semilogy(D',FSDmat,'Parent',subplot1,'LineWidth',2);
set(semilogy1(1),'Color',[1 0 0],'DisplayName','Area dist at T =  1e+04');
set(semilogy1(2),'Color',[0 1 0],'DisplayName','Area dist at T =      0');
set(semilogy1(3),'Marker','o','LineStyle','none','Color',[0 0 1],...
    'DisplayName','Pack Area dist',...
    'LineWidth',0.5);

% Create xlabel
xlabel('Diameter');

% Create ylabel
ylabel('Fractional Area');

% Create title
title('Area dist at T =      0 and T =  1e+04','FontSize',16,...
    'FontName','Hoefler Text');

% Create legend
legend(subplot1,'show');

%% 
% Create subplot
subplot2 = subplot(4,1,2,'Parent',figure1,'YScale','log','YMinorTick','on');
box(subplot2,'on');
hold(subplot2,'all');

% Create multiple lines using matrix input to semilogy

semilogy2 = semilogy(Timer,Relmag,'Parent',subplot2,'LineWidth',2);

set(semilogy2(1),'Color',[0 1 0],'DisplayName','1');
set(semilogy2(2),'LineStyle','--','Color',[0 0 1],...
    'DisplayName','Melting/Advection');
set(semilogy2(3),'LineStyle','-.','Color',[1 0 0],...
    'DisplayName','Redistribution/Advection');
set(semilogy2(4),'LineStyle','-.','Color',[0 1 1],...
    'DisplayName','Swell/Advection');

% Create xlabel
xlabel('Time');

% Create ylabel
ylabel('Ratio of Area Change');

% Create title
title('Relative Area Change over time with respect to advection',...
    'FontSize',16,...
    'FontName','Hoefler Text');

% Create legend
legend(subplot2,'show');

%% 
% Create subplot
subplot3 = subplot(4,1,3,'Parent',figure1,'YScale','log','YMinorTick','on');
box(subplot3,'on');
hold(subplot3,'all');

% Create multiple lines using matrix input to semilogy
semilogy3 = semilogy(Timer,Mag,'Parent',subplot3,'LineWidth',2);
set(semilogy3(1),'Color',[0 1 0],'DisplayName','Advection');
set(semilogy3(2),'LineStyle','--','Color',[0 0 1],'DisplayName','Melting');
set(semilogy3(3),'LineStyle','-.','Color',[1 0 0],...
    'DisplayName','Redistribution');
set(semilogy3(4),'LineStyle','-.','Color',[0 1 1],'DisplayName','Swell');

% Create xlabel
xlabel('Time');

% Create ylabel
ylabel({'Area Change',''});

% Create title
title('Total Area change over time','FontSize',16,'FontName','Hoefler Text');

% Create legend
legend(subplot3,'show');

%%
% Create subplot
subplot4 = subplot(4,1,4,'Parent',figure1);
hold(subplot4,'all');

% Create multiple lines using matrix input to plot
plot1 = plot(D',Save,'Parent',subplot4,'LineWidth',2);
set(plot1(1),'LineStyle','--','Color',[0 1 0],...
    'DisplayName','Advection Area Changes');
set(plot1(2),'LineStyle','-.','Color',[0 0 1],'DisplayName','Melting');
set(plot1(3),'LineStyle','-.','Color',[1 0 0],...
    'DisplayName','Redistribution');
set(plot1(4),'LineStyle','--','Color',[0 1 1],'DisplayName','Swell');

% Create xlabel
xlabel({'Diameter',''});

% Create ylabel
ylabel({'Integrated Area Change',''});

% Create title
title('Integrated Changes in Area due to various processes, for each gridbox',...
    'FontSize',16,...
    'FontName','Hoefler Text');

% Create legend
legend(subplot4,'show');





end