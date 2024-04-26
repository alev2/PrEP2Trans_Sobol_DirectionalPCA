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

colorMat=[...
    .494 .184 .556;
    .929 .694 .125;
    .850 .325 .098;
    .000 .447 .741;
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
legendfs=14;


legendText={...
    %'Mixing level',...
    'PrEP, HETM',...
    'PrEP, HETF',...
    'PrEP, MSM',...
    'PrEP, PWID',...
};
outcomeText={...
    'Incidence, MSM',...
    'Incidence, HETF',...
    'Incidence, HETM',...
    'Incidence, PWID',...
    'Incidence, Total',...
    'Num. on PrEP, MSM',...
    'Num. on PrEP, HETF',...
    'Num. on PrEP, HETM',...
    'Num. on PrEP, PWID',...
    'Coverage, MSM',...
    'Coverage, HETF',...
    'Coverage, HETM',...
    'Coverage, PWID',...
};




load('./Sob_g1_All.mat');
load('./Sob_g2_All.mat');
load('./Sob_g3_All.mat');
load('./Sob_g4_All.mat');
load('./Sob_g5_All.mat');
%load('./Sob_g6_All.mat');
%load('./Sob_g7_All.mat');

load('./Tot_Sob_g1_All.mat');
load('./Tot_Sob_g2_All.mat');
load('./Tot_Sob_g3_All.mat');
load('./Tot_Sob_g4_All.mat');
load('./Tot_Sob_g5_All.mat');
%load('./Tot_Sob_g6_All.mat');
%load('./Tot_Sob_g7_All.mat');

posMat=[0 -200 1100 750];
yrs_Plot=yearRange(yearInds);


interventionCount=1;
tp=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile

%plot(yrs_Plot, Tot_Sob_g1_All(1,:) ,':s','color',colorMat(1,:),'LineWidth',lw);
plot(yrs_Plot, Tot_Sob_g1_All(2,:) ,':s','color',colorMat(2,:),'LineWidth',lw);
hold on
plot(yrs_Plot, Tot_Sob_g1_All(3,:) ,':s','color',colorMat(3,:),'LineWidth',lw);
plot(yrs_Plot, Tot_Sob_g1_All(4,:) ,':s','color',colorMat(4,:),'LineWidth',lw);
plot(yrs_Plot, Tot_Sob_g1_All(5,:) ,':s','color',colorMat(5,:),'LineWidth',lw);

xlim([yrs_Plot(1) yrs_Plot(end)])
ylim([0 1])

xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
%ylabel('Variance', 'Interpreter','latex','FontSize',50);
legend(legendText,'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',2);
title(strcat('Effect on\,', outcomeText{interventionCount}), 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;


nexttile
interventionCount=interventionCount+1;
%plot(yrs_Plot, Tot_Sob_g2_All(1,:) ,':s','color',colorMat(1,:),'LineWidth',lw);
plot(yrs_Plot, Tot_Sob_g2_All(2,:) ,':s','color',colorMat(2,:),'LineWidth',lw);
hold on
plot(yrs_Plot, Tot_Sob_g2_All(3,:) ,':s','color',colorMat(3,:),'LineWidth',lw);
plot(yrs_Plot, Tot_Sob_g2_All(4,:) ,':s','color',colorMat(4,:),'LineWidth',lw);
plot(yrs_Plot, Tot_Sob_g2_All(5,:) ,':s','color',colorMat(5,:),'LineWidth',lw);

xlim([yrs_Plot(1) yrs_Plot(end)])
ylim([0 1])

xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
%ylabel('Variance', 'Interpreter','latex','FontSize',50);
legend(legendText,'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',2);
title(strcat('Effect on\,', outcomeText{interventionCount}), 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;



nexttile
interventionCount=interventionCount+1;
%plot(yrs_Plot, Tot_Sob_g3_All(1,:) ,':s','color',colorMat(1,:),'LineWidth',lw);
plot(yrs_Plot, Tot_Sob_g3_All(2,:) ,':s','color',colorMat(2,:),'LineWidth',lw);
hold on
plot(yrs_Plot, Tot_Sob_g3_All(3,:) ,':s','color',colorMat(3,:),'LineWidth',lw);
plot(yrs_Plot, Tot_Sob_g3_All(4,:) ,':s','color',colorMat(4,:),'LineWidth',lw);
plot(yrs_Plot, Tot_Sob_g3_All(5,:) ,':s','color',colorMat(5,:),'LineWidth',lw);

xlim([yrs_Plot(1) yrs_Plot(end)])
ylim([0 1])

xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
%ylabel('Variance', 'Interpreter','latex','FontSize',50);
legend(legendText,'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',2);
title(strcat('Effect on\,', outcomeText{interventionCount}), 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;



nexttile
interventionCount=interventionCount+1;
%plot(yrs_Plot, Tot_Sob_g4_All(1,:) ,':s','color',colorMat(1,:),'LineWidth',lw);
plot(yrs_Plot, Tot_Sob_g4_All(2,:) ,':s','color',colorMat(2,:),'LineWidth',lw);
hold on
plot(yrs_Plot, Tot_Sob_g4_All(3,:) ,':s','color',colorMat(3,:),'LineWidth',lw);
plot(yrs_Plot, Tot_Sob_g4_All(4,:) ,':s','color',colorMat(4,:),'LineWidth',lw);
plot(yrs_Plot, Tot_Sob_g4_All(5,:) ,':s','color',colorMat(5,:),'LineWidth',lw);

xlim([yrs_Plot(1) yrs_Plot(end)])
ylim([0 1])

xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
%ylabel('Variance', 'Interpreter','latex','FontSize',50);
legend(legendText,'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',2);
title(strcat('Effect on\,', outcomeText{interventionCount}), 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;


% % % 
% % % 
% nexttile
% interventionCount=interventionCount+1;
% plot(yrs_Plot, Tot_Sob_g5_All(1,:) ,':s','color',colorMat(1,:),'LineWidth',lw);
% hold on
% plot(yrs_Plot, Tot_Sob_g5_All(2,:) ,':s','color',colorMat(2,:),'LineWidth',lw);
% plot(yrs_Plot, Tot_Sob_g5_All(3,:) ,':s','color',colorMat(3,:),'LineWidth',lw);
% plot(yrs_Plot, Tot_Sob_g5_All(4,:) ,':s','color',colorMat(4,:),'LineWidth',lw);
% plot(yrs_Plot, Tot_Sob_g5_All(5,:) ,':s','color',colorMat(5,:),'LineWidth',lw);
% 
% xlim([yrs_Plot(1) yrs_Plot(end)])
% ylim([0 1])
% 
% xlabel('Year', 'Interpreter','latex','FontSize',50);
% ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
% %ylabel('Variance', 'Interpreter','latex','FontSize',50);
% legend(legendText,'Interpreter','latex','FontSize',legendfs,'Location','best',...
%     'AutoUpdate','off','NumColumns',2);
% title(strcat('Effect on\,', outcomeText{interventionCount}), 'Interpreter','latex','FontSize',50);
% set(gca,'FontSize',fs);
% ax=gca;
% ax.YAxis.Exponent=0;

set(gcf,'Position',posMat);