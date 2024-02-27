colorMat=[...
    .000 .447 .741;
    .850 .325 .098;
    .929 .694 .125;
    .494 .184 .556;
    .466 .674 .188;
    .301 .745 .933;
    .635 .078 .184;
    .1235 .078 .484;
    .8235 .3278 .384;
    .5235 .7278 .184;
    .11 .66 .55;
    .51 .26 .85;
    .51 .96 .05;
    .51 .556 .15;
    .151 .356 .75;
    0 0 0    
    ];


lw=3;
lw2=2;
fs=15;
legendfs=12;

legendText={...
    'PrEP, B',...
    '$1\to2$, B',...
    '$2\to3$ (at diag), B',...
    '$2\to3$ (after diag), B',...
    '$3\to4$, B',...
    '$4\to3$, B',...
    '$4\to5$, B',...
    '$5\to4$, B',...
    'PrEP, H',...
    '$1\to2$, H',...
    '$2\to3$ (at diag), H',...
    '$2\to3$ (after diag), H',...
    '$3\to4$, H',...
    '$4\to3$, H',...
    '$4\to5$, H',...
    '$5\to4$, H',...
};

% load('./Indices_833Points/Indices_833Points/Tot_Sob_g1_All.mat');
% load('./Indices_833Points/Indices_833Points/Tot_Sob_g2_All.mat');
% load('./Indices_833Points/Indices_833Points/Tot_Sob_g3_All.mat');
% load('./Indices_833Points/Indices_833Points/Tot_Sob_g4_All.mat');
% load('./Indices_833Points/Indices_833Points/Tot_Tot_Sob_g1_All.mat');
% load('./Indices_833Points/Indices_833Points/Tot_Tot_Sob_g2_All.mat');
% load('./Indices_833Points/Indices_833Points/Tot_Tot_Sob_g3_All.mat');
% load('./Indices_833Points/Indices_833Points/Tot_Tot_Sob_g4_All.mat');


load('./Sob_g1_All.mat');
load('./Sob_g2_All.mat');
load('./Sob_g3_All.mat');
load('./Sob_g4_All.mat');
load('./Sob_g5_All.mat');
load('./Tot_Sob_g1_All.mat');
load('./Tot_Sob_g2_All.mat');
load('./Tot_Sob_g3_All.mat');
load('./Tot_Sob_g4_All.mat');
load('./Tot_Sob_g5_All.mat');

Sec_Sob_g1_All=Tot_Sob_g1_All-Sob_g1_All;
Sec_Sob_g2_All=Tot_Sob_g2_All-Sob_g2_All;
Sec_Sob_g3_All=Tot_Sob_g3_All-Sob_g3_All;
Sec_Sob_g4_All=Tot_Sob_g4_All-Sob_g4_All;
Sec_Sob_g5_All=Tot_Sob_g5_All-Sob_g5_All;

posMat=.5*[0 0 1600 700];
%yrs_Plot=(2017:1:2026)';
%yrs_Plot=(2017:1:(2017+num_Yrs-1))';
yrs_Plot=(2023:1:2034)';


%num_Yrs=4;
%generate some plots

tp=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
%subplot(2,2,1)
nexttile

plot(yrs_Plot, Sec_Sob_g1_All(1,:) ,':s','color',colorMat(1,:),'LineWidth',lw);
%plot(yrs_Plot, Vi1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
hold on
plot(yrs_Plot, Sec_Sob_g1_All(3,:) ,':s','color',colorMat(2,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g1_All(5,:) ,':s','color',colorMat(3,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g1_All(7,:) ,':s','color',colorMat(4,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g1_All(9,:) ,':s','color',colorMat(5,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g1_All(11,:) ,':s','color',colorMat(6,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g1_All(13,:) ,':s','color',colorMat(7,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g1_All(15,:) ,':s','color',colorMat(8,:),'LineWidth',lw);

plot(yrs_Plot, Sec_Sob_g1_All(2,:) ,'-o','color',colorMat(1,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g1_All(4,:) ,'-o','color',colorMat(2,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g1_All(6,:) ,'-o','color',colorMat(3,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g1_All(8,:) ,'-o','color',colorMat(4,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g1_All(10,:) ,'-o','color',colorMat(5,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g1_All(12,:) ,'-o','color',colorMat(6,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g1_All(14,:) ,'-o','color',colorMat(7,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g1_All(16,:) ,'-o','color',colorMat(8,:),'LineWidth',lw2);
xlim([yrs_Plot(1) yrs_Plot(end)])
%ylim([0 1])

%plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
%xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
%ylabel('Variance', 'Interpreter','latex','FontSize',50);
legend(legendText,'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',2);
%for k=2--size(survPrevRange,1)
%    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
%end
title('Intervention Effects on incidence (Black)', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;


nexttile
plot(yrs_Plot, Sec_Sob_g2_All(1,:) ,':s','color',colorMat(1,:),'LineWidth',lw);
%plot(yrs_Plot, Vi1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
hold on
plot(yrs_Plot, Sec_Sob_g2_All(3,:) ,':s','color',colorMat(2,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g2_All(5,:) ,':s','color',colorMat(3,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g2_All(7,:) ,':s','color',colorMat(4,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g2_All(9,:) ,':s','color',colorMat(5,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g2_All(11,:) ,':s','color',colorMat(6,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g2_All(13,:) ,':s','color',colorMat(7,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g2_All(15,:) ,':s','color',colorMat(8,:),'LineWidth',lw);

plot(yrs_Plot, Sec_Sob_g2_All(2,:) ,'-o','color',colorMat(1,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g2_All(4,:) ,'-o','color',colorMat(2,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g2_All(6,:) ,'-o','color',colorMat(3,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g2_All(8,:) ,'-o','color',colorMat(4,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g2_All(10,:) ,'-o','color',colorMat(5,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g2_All(12,:) ,'-o','color',colorMat(6,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g2_All(14,:) ,'-o','color',colorMat(7,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g2_All(16,:) ,'-o','color',colorMat(8,:),'LineWidth',lw2);

xlim([yrs_Plot(1) yrs_Plot(end)])
%ylim([0 1])

%plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
%ylabel('Variance', 'Interpreter','latex','FontSize',50);
legend(legendText,'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',2);
%for k=2--size(survPrevRange,1)
%    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
%end
title('Intervention Effects on incidence (Hisp.)', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;



nexttile
plot(yrs_Plot, Sec_Sob_g3_All(1,:) ,':s','color',colorMat(1,:),'LineWidth',lw);
%plot(yrs_Plot, Vi1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
hold on
plot(yrs_Plot, Sec_Sob_g3_All(3,:) ,':s','color',colorMat(2,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g3_All(5,:) ,':s','color',colorMat(3,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g3_All(7,:) ,':s','color',colorMat(4,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g3_All(9,:) ,':s','color',colorMat(5,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g3_All(11,:) ,':s','color',colorMat(6,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g3_All(13,:) ,':s','color',colorMat(7,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g3_All(15,:) ,':s','color',colorMat(8,:),'LineWidth',lw);

plot(yrs_Plot, Sec_Sob_g3_All(2,:) ,'-o','color',colorMat(1,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g3_All(4,:) ,'-o','color',colorMat(2,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g3_All(6,:) ,'-o','color',colorMat(3,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g3_All(8,:) ,'-o','color',colorMat(4,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g3_All(10,:) ,'-o','color',colorMat(5,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g3_All(12,:) ,'-o','color',colorMat(6,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g3_All(14,:) ,'-o','color',colorMat(7,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g3_All(16,:) ,'-o','color',colorMat(8,:),'LineWidth',lw2);
xlim([yrs_Plot(1) yrs_Plot(end)])
%ylim([0 1])

%plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
%ylabel('Variance', 'Interpreter','latex','FontSize',50);
legend(legendText,'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',2);
%for k=2--size(survPrevRange,1)
%    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
%end
title('Intervention Effects on IRR (Black)', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;



nexttile
plot(yrs_Plot, Sec_Sob_g4_All(1,:) ,':s','color',colorMat(1,:),'LineWidth',lw);
%plot(yrs_Plot, Vi1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
hold on
plot(yrs_Plot, Sec_Sob_g4_All(3,:) ,':s','color',colorMat(2,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g4_All(5,:) ,':s','color',colorMat(3,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g4_All(7,:) ,':s','color',colorMat(4,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g4_All(9,:) ,':s','color',colorMat(5,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g4_All(11,:) ,':s','color',colorMat(6,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g4_All(13,:) ,':s','color',colorMat(7,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g4_All(15,:) ,':s','color',colorMat(8,:),'LineWidth',lw);

plot(yrs_Plot, Sec_Sob_g4_All(2,:) ,'-o','color',colorMat(1,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g4_All(4,:) ,'-o','color',colorMat(2,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g4_All(6,:) ,'-o','color',colorMat(3,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g4_All(8,:) ,'-o','color',colorMat(4,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g4_All(10,:) ,'-o','color',colorMat(5,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g4_All(12,:) ,'-o','color',colorMat(6,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g4_All(14,:) ,'-o','color',colorMat(7,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g4_All(16,:) ,'-o','color',colorMat(8,:),'LineWidth',lw2);
xlim([yrs_Plot(1) yrs_Plot(end)])
%ylim([0 1])

%plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
%ylabel('Variance', 'Interpreter','latex','FontSize',50);
legend(legendText,'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',2);
%for k=2--size(survPrevRange,1)
%    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
%end
title('Intervention Effects on IRR (Hisp.)', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;


% 
nexttile
plot(yrs_Plot, Sec_Sob_g5_All(1,:) ,':s','color',colorMat(1,:),'LineWidth',lw);
%plot(yrs_Plot, Vi1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
hold on
plot(yrs_Plot, Sec_Sob_g5_All(3,:) ,':s','color',colorMat(2,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g5_All(5,:) ,':s','color',colorMat(3,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g5_All(7,:) ,':s','color',colorMat(4,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g5_All(9,:) ,':s','color',colorMat(5,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g5_All(11,:) ,':s','color',colorMat(6,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g5_All(13,:) ,':s','color',colorMat(7,:),'LineWidth',lw);
plot(yrs_Plot, Sec_Sob_g5_All(15,:) ,':s','color',colorMat(8,:),'LineWidth',lw);

plot(yrs_Plot, Sec_Sob_g5_All(2,:) ,'-o','color',colorMat(1,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g5_All(4,:) ,'-o','color',colorMat(2,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g5_All(6,:) ,'-o','color',colorMat(3,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g5_All(8,:) ,'-o','color',colorMat(4,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g5_All(10,:) ,'-o','color',colorMat(5,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g5_All(12,:) ,'-o','color',colorMat(6,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g5_All(14,:) ,'-o','color',colorMat(7,:),'LineWidth',lw2);
plot(yrs_Plot, Sec_Sob_g5_All(16,:) ,'-o','color',colorMat(8,:),'LineWidth',lw2);
xlim([yrs_Plot(1) yrs_Plot(end)])

%plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Principal sobol index', 'Interpreter','latex','FontSize',50);
%ylabel('Variance', 'Interpreter','latex','FontSize',50);
legend(legendText,'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',2);
%for k=2--size(survPrevRange,1)
%    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
%end
title('Intervention Effects on spending', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;
xlim([yrs_Plot(1) yrs_Plot(end)])
ylim([0 1])


set(gcf,'Position',posMat);
