colorMat=[...
    .000 .447 .741;
    .850 .325 .098;
    .929 .694 .125;
    .494 .184 .556;
    .466 .674 .188;
    .301 .745 .933;
    .635 .078 .184];
lw=3;
fs=15;
legendfs=12;

load('./Indices_833Points/Indices_833Points/Sob_g1_All.mat');
load('./Indices_833Points/Indices_833Points/Sob_g2_All.mat');
load('./Indices_833Points/Indices_833Points/Sob_g3_All.mat');
load('./Indices_833Points/Indices_833Points/Sob_g4_All.mat');
load('./Indices_833Points/Indices_833Points/Tot_Sob_g1_All.mat');
load('./Indices_833Points/Indices_833Points/Tot_Sob_g2_All.mat');
load('./Indices_833Points/Indices_833Points/Tot_Sob_g3_All.mat');
load('./Indices_833Points/Indices_833Points/Tot_Sob_g4_All.mat');

posMat=.5*[0 0 1600 700];
%yrs_Plot=(2017:1:2026)';
%yrs_Plot=(2017:1:(2017+num_Yrs-1))';
yrs_Plot=(2023:1:2031)';


%num_Yrs=4;
%generate some plots

subplot(1,2,1)
plot(yrs_Plot, Sob_g1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
%plot(yrs_Plot, Vi1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
hold on
plot(yrs_Plot, Sob_g1_All(2,:) ,':','color',colorMat(2,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g1_All(3,:) ,':','color',colorMat(3,:),'LineWidth',lw*2);
plot(yrs_Plot, Sob_g1_All(4,:) ,':','color',colorMat(4,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g1_All(5,:) ,':','color',colorMat(5,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g1_All(6,:) ,':','color',colorMat(6,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g1_All(7,:) ,':','color',colorMat(5,:),'LineWidth',lw);
plot(yrs_Plot, Sob_g1_All(8,:) ,':','color',colorMat(6,:),'LineWidth',lw);
%plot(yrs_Plot, Vi2_All(2,:) ,':','color',colorMat(2,:),'LineWidth',lw);
%plot(yrs_Plot, Vi2_All(3,:) ,':','color',colorMat(3,:),'LineWidth',lw*2);
xlim([yrs_Plot(1) yrs_Plot(end)])

%plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Principal sobol index', 'Interpreter','latex','FontSize',50);
legend(['Testing'],['PrEP'],['Maintain VLS'],['Initiate/reinitiate ART'],['Maintain ART (Not VLS)'],['Maintain ART (VLS)'],'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off');
%for k=2:size(survPrevRange,1)
%    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
%end
title('Relative effect on HIV incidence, interventions', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;


subplot(1,2,2)
plot(yrs_Plot, Tot_Sob_g1_All(1,:)-Sob_g1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
%plot(yrs_Plot, Vi1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
hold on
plot(yrs_Plot, Tot_Sob_g1_All(2,:)-Sob_g1_All(2,:) ,':','color',colorMat(2,:),'LineWidth',lw);
plot(yrs_Plot, Tot_Sob_g1_All(3,:)-Sob_g1_All(3,:) ,':','color',colorMat(3,:),'LineWidth',lw*2);
plot(yrs_Plot, Tot_Sob_g1_All(4,:)-Sob_g1_All(4,:) ,':','color',colorMat(4,:),'LineWidth',lw);
plot(yrs_Plot, Tot_Sob_g1_All(5,:)-Sob_g1_All(5,:) ,':','color',colorMat(5,:),'LineWidth',lw);
plot(yrs_Plot, Tot_Sob_g1_All(6,:)-Sob_g1_All(6,:) ,':','color',colorMat(6,:),'LineWidth',lw);
%plot(yrs_Plot, Vi2_All(2,:) ,':','color',colorMat(2,:),'LineWidth',lw);
%plot(yrs_Plot, Vi2_All(3,:) ,':','color',colorMat(3,:),'LineWidth',lw*2);
xlim([yrs_Plot(1) yrs_Plot(end)])

%plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Interaction sobol index', 'Interpreter','latex','FontSize',50);
legend(['Testing'],['PrEP'],['Maintain VLS'],['Initiate/reinitiate ART'],['Maintain ART (Not VLS)'],['Maintain ART (VLS)'],'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off');
%for k=2:size(survPrevRange,1)
%    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
%end
title('Relative effect on HIV incidence, interventions', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;


% 
% subplot(2,3,3)
% 
% plot(yrs_Plot, Sob_g3_All(1,:),':','color',colorMat(1,:),'LineWidth',lw);
% %plot(yrs_Plot, Vi3_All(1,:),':','color',colorMat(1,:),'LineWidth',lw);
% hold on
% plot(yrs_Plot, Sob_g3_All(2,:),':','color',colorMat(2,:),'LineWidth',lw);
% plot(yrs_Plot, Sob_g3_All(3,:) ,':','color',colorMat(3,:),'LineWidth',lw);
% %plot(yrs_Plot, Vi3_All(2,:),':','color',colorMat(2,:),'LineWidth',lw);
% %plot(yrs_Plot, Vi3_All(3,:) ,':','color',colorMat(3,:),'LineWidth',lw);
% xlim([yrs_Plot(1) yrs_Plot(end)])
% %plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
% xlabel('Year', 'Interpreter','latex','FontSize',50);
% ylabel('Principal sobol index', 'Interpreter','latex','FontSize',50);
% legend(['PrEP increase, MSM'],['PrEP increase, HetF'],['PrEP increase, HetM'],'Interpreter','latex','FontSize',legendfs,'Location','best',...
%     'AutoUpdate','off');
% %for k=2:size(survPrevRange,1)
% %    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
% %end
% title('Spillover of PrEP on HETM incidence', 'Interpreter','latex','FontSize',50);
% set(gca,'FontSize',fs);
% ax=gca;
% ax.YAxis.Exponent=0;
% 
% 
% 
% subplot(2,3,4)
% plot(yrs_Plot, Tot_Sob_g1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
% %plot(yrs_Plot, VTi1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
% hold on
% plot(yrs_Plot, Tot_Sob_g1_All(2,:) ,':','color',colorMat(2,:),'LineWidth',lw);
% plot(yrs_Plot, Tot_Sob_g1_All(3,:) ,':','color',colorMat(3,:),'LineWidth',lw*2);
% %plot(yrs_Plot, VTi1_All(2,:) ,':','color',colorMat(2,:),'LineWidth',lw);
% %plot(yrs_Plot, VTi1_All(3,:) ,':','color',colorMat(3,:),'LineWidth',lw*2);
% xlim([yrs_Plot(1) yrs_Plot(end)])
% 
% %plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
% xlabel('Year', 'Interpreter','latex','FontSize',50);
% ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
% legend(['PrEP increase, MSM'],['PrEP increase, HetF'],['PrEP increase, HetM'],'Interpreter','latex','FontSize',legendfs,'Location','best',...
%     'AutoUpdate','off');
% %for k=2:size(survPrevRange,1)
% %    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
% %end
% title('Spillover of PrEP on MSM incidence', 'Interpreter','latex','FontSize',50);
% set(gca,'FontSize',fs);
% ax=gca;
% ax.YAxis.Exponent=0;
% 
% 
% 
% subplot(2,3,5)
% plot(yrs_Plot, Tot_Sob_g2_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
% %plot(yrs_Plot, VTi2_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
% hold on
% plot(yrs_Plot, Tot_Sob_g2_All(2,:) ,':','color',colorMat(2,:),'LineWidth',lw);
% plot(yrs_Plot, Tot_Sob_g2_All(3,:) ,':','color',colorMat(3,:),'LineWidth',lw*2);
% %plot(yrs_Plot, VTi2_All(2,:) ,':','color',colorMat(2,:),'LineWidth',lw);
% %plot(yrs_Plot, VTi2_All(3,:) ,':','color',colorMat(3,:),'LineWidth',lw*2);
% xlim([yrs_Plot(1) yrs_Plot(end)])
% 
% %plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
% xlabel('Year', 'Interpreter','latex','FontSize',50);
% ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
% legend(['PrEP increase, MSM'],['PrEP increase, HetF'],['PrEP increase, HetM'],'Interpreter','latex','FontSize',legendfs,'Location','best',...
%     'AutoUpdate','off');
% %for k=2:size(survPrevRange,1)
% %    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
% %end
% title('Spillover of PrEP on HETF incidence', 'Interpreter','latex','FontSize',50);
% set(gca,'FontSize',fs);
% ax=gca;
% ax.YAxis.Exponent=0;
% 
% 
% subplot(2,3,6)
% 
% plot(yrs_Plot, Tot_Sob_g3_All(1,:),':','color',colorMat(1,:),'LineWidth',lw);
% %plot(yrs_Plot, VTi3_All(1,:),':','color',colorMat(1,:),'LineWidth',lw);
% hold on
% plot(yrs_Plot, Tot_Sob_g3_All(2,:),':','color',colorMat(2,:),'LineWidth',lw);
% plot(yrs_Plot, Tot_Sob_g3_All(3,:) ,':','color',colorMat(3,:),'LineWidth',lw);
% %plot(yrs_Plot, VTi3_All(2,:),':','color',colorMat(2,:),'LineWidth',lw);
% %plot(yrs_Plot, VTi3_All(3,:) ,':','color',colorMat(3,:),'LineWidth',lw);
% xlim([yrs_Plot(1) yrs_Plot(end)])
% %plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
% xlabel('Year', 'Interpreter','latex','FontSize',50);
% ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
% legend(['PrEP increase, MSM'],['PrEP increase, HetF'],['PrEP increase, HetM'],'Interpreter','latex','FontSize',legendfs,'Location','best',...
%     'AutoUpdate','off');
% %for k=2:size(survPrevRange,1)
% %    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
% %end
% title('Spillover of PrEP on HETM incidence', 'Interpreter','latex','FontSize',50);
% set(gca,'FontSize',fs);
% ax=gca;
% ax.YAxis.Exponent=0;


set(gcf,'Position',posMat);
