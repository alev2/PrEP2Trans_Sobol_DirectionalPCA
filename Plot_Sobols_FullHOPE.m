colorMat=[...
    .000 .447 .741;
    .850 .325 .098;
    .929 .694 .125;
    .494 .184 .556;
    .466 .674 .188;
    .301 .745 .933;
    .635 .078 .184;
    0 0 0    
    ];
lw=3;
fs=15;
legendfs=12;

legendText={...
    'PrEP, B',...
    '$5\to4$, B',...
    '$4\to5$, B',...
    '$4\to3$, B',...
    'PrEP, H',...
    '$5\to4$, H',...
    '$4\to5$, H',...
    '$4\to3$, H'...
};

% load('./Indices_833Points/Indices_833Points/Sob_g1_All.mat');
% load('./Indices_833Points/Indices_833Points/Sob_g2_All.mat');
% load('./Indices_833Points/Indices_833Points/Sob_g3_All.mat');
% load('./Indices_833Points/Indices_833Points/Sob_g4_All.mat');
% load('./Indices_833Points/Indices_833Points/Tot_Sob_g1_All.mat');
% load('./Indices_833Points/Indices_833Points/Tot_Sob_g2_All.mat');
% load('./Indices_833Points/Indices_833Points/Tot_Sob_g3_All.mat');
% load('./Indices_833Points/Indices_833Points/Tot_Sob_g4_All.mat');


load('./HistDepInds_833Pts/Sob_g1_All.mat');
load('./HistDepInds_833Pts/Sob_g2_All.mat');
load('./HistDepInds_833Pts/Sob_g3_All.mat');
load('./HistDepInds_833Pts/Sob_g4_All.mat');
load('./HistDepInds_833Pts/Tot_Sob_g1_All.mat');
load('./HistDepInds_833Pts/Tot_Sob_g2_All.mat');
load('./HistDepInds_833Pts/Tot_Sob_g3_All.mat');
load('./HistDepInds_833Pts/Tot_Sob_g4_All.mat');

posMat=.5*[0 0 1600 700];
%yrs_Plot=(2017:1:2026)';
%yrs_Plot=(2017:1:(2017+num_Yrs-1))';
yrs_Plot=(2023:1:2031)';


%num_Yrs=4;
%generate some plots

subplot(2,2,1)
plot(yrs_Plot, Sob_g1_All(1,:) ,':o','color',colorMat(1,:),'LineWidth',lw);
%plot(yrs_Plot, Vi1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
hold on
plot(yrs_Plot, Sob_g1_All(3,:) ,':o','color',colorMat(3,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g1_All(5,:) ,':o','color',colorMat(5,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g1_All(7,:) ,':o','color',colorMat(7,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g1_All(2,:) ,':s','color',colorMat(2,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g1_All(4,:) ,':s','color',colorMat(4,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g1_All(6,:) ,':s','color',colorMat(6,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g1_All(8,:) ,':s','color',colorMat(8,:),'LineWidth',lw);
xlim([yrs_Plot(1) yrs_Plot(end)])

%plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
%xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Principal sobol index', 'Interpreter','latex','FontSize',50);
%ylabel('Variance', 'Interpreter','latex','FontSize',50);
legend(legendText,'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',2);
%for k=2:size(survPrevRange,1)
%    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
%end
title('Intervention Effects on incidence (Black)', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;


subplot(2,2,2)
plot(yrs_Plot, Sob_g2_All(1,:) ,':o','color',colorMat(1,:),'LineWidth',lw);
%plot(yrs_Plot, Vi1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
hold on
plot(yrs_Plot, Sob_g2_All(3,:) ,':o','color',colorMat(3,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g2_All(5,:) ,':o','color',colorMat(5,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g2_All(7,:) ,':o','color',colorMat(7,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g2_All(2,:) ,':s','color',colorMat(2,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g2_All(4,:) ,':s','color',colorMat(4,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g2_All(6,:) ,':s','color',colorMat(6,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g2_All(8,:) ,':s','color',colorMat(8,:),'LineWidth',lw);

xlim([yrs_Plot(1) yrs_Plot(end)])

%plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Principal sobol index', 'Interpreter','latex','FontSize',50);
%ylabel('Variance', 'Interpreter','latex','FontSize',50);
legend(legendText,'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',2);
%for k=2:size(survPrevRange,1)
%    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
%end
title('Intervention Effects on incidence (Hisp.)', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;



subplot(2,2,3)
plot(yrs_Plot, Sob_g3_All(1,:) ,':o','color',colorMat(1,:),'LineWidth',lw);
%plot(yrs_Plot, Vi1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
hold on
plot(yrs_Plot, Sob_g3_All(3,:) ,':o','color',colorMat(3,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g3_All(5,:) ,':o','color',colorMat(5,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g3_All(7,:) ,':o','color',colorMat(7,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g3_All(2,:) ,':s','color',colorMat(2,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g3_All(4,:) ,':s','color',colorMat(4,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g3_All(6,:) ,':s','color',colorMat(6,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g3_All(8,:) ,':s','color',colorMat(8,:),'LineWidth',lw);

xlim([yrs_Plot(1) yrs_Plot(end)])

%plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Principal sobol index', 'Interpreter','latex','FontSize',50);
%ylabel('Variance', 'Interpreter','latex','FontSize',50);
legend(legendText,'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',2);
%for k=2:size(survPrevRange,1)
%    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
%end
title('Intervention Effects on IRR (Black)', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;



subplot(2,2,4)
plot(yrs_Plot, Sob_g4_All(1,:) ,':o','color',colorMat(1,:),'LineWidth',lw);
%plot(yrs_Plot, Vi1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
hold on
plot(yrs_Plot, Sob_g4_All(3,:) ,':o','color',colorMat(3,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g4_All(5,:) ,':o','color',colorMat(5,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g4_All(7,:) ,':o','color',colorMat(7,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g4_All(2,:) ,':s','color',colorMat(2,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g4_All(4,:) ,':s','color',colorMat(4,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g4_All(6,:) ,':s','color',colorMat(6,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g4_All(8,:) ,':s','color',colorMat(8,:),'LineWidth',lw);
xlim([yrs_Plot(1) yrs_Plot(end)])

%plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Principal sobol index', 'Interpreter','latex','FontSize',50);
%ylabel('Variance', 'Interpreter','latex','FontSize',50);
legend(legendText,'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',2);
%for k=2:size(survPrevRange,1)
%    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
%end
title('Intervention Effects on IRR (Hisp.)', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;



set(gcf,'Position',posMat);
