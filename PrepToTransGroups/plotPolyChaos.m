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


interventionChoices={...
    'Mixing',...
    'PrEP, HETM',...
    'PrEP, HETF',...
    'PrEP, MSM',...
    'PrEP, PWID',...
};

yr1String=num2str(yearRange(yearInds(yr1)));
yr2String=num2str(yearRange(yearInds(yr2)));
yrDist=(yr2-yr1)+1;
interventionString=interventionChoices{idx_dx};




posMat=[-1863,166,1071,703];
tp=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile
output_idx=5;
plot(relevantInterval, output_Of_Interest(:,output_idx),':','color',colorMat(5,:),'LineWidth',lw);
%ylim([5000 40000])
hold on
%plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
xlabel(interventionString, 'Interpreter','latex','FontSize',50);
ylabel(strcat(outcomeText{output_idx},' ',yr1String,'-',yr2String), 'Interpreter','latex','FontSize',50);
legend({'Total'},'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',1);
%plot(Sr.knots*100, output_Of_Interest_KnotPts(:,7),'o','color',colorMat(5,:),'LineWidth',lw);
%for k=2--size(survPrevRange,1)
%    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
%end
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;

nexttile
plot(relevantInterval, output_Of_Interest(:,1),':','color',colorMat(1,:),'LineWidth',lw);
hold on
plot(relevantInterval, output_Of_Interest(:,2),':','color',colorMat(2,:),'LineWidth',lw);
plot(relevantInterval, output_Of_Interest(:,3),':','color',colorMat(3,:),'LineWidth',lw);
plot(relevantInterval, output_Of_Interest(:,4),':','color',colorMat(4,:),'LineWidth',lw);
plot(relevantInterval, output_Of_Interest(:,5),':','color',colorMat(5,:),'LineWidth',lw);
%xlim([yrs_Plot(1) yrs_Plot(end)])
%ylim([1000 13000])
%plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
xlabel(interventionString, 'Interpreter','latex','FontSize',50);
ylabel(strcat('Incidence ', yr1String,'-',yr2String), 'Interpreter','latex','FontSize',50);
legend(outcomeText{1:5},'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',1);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;


nexttile
plot(relevantInterval, output_Of_Interest(:,6)/yrDist,':','color',colorMat(1,:),'LineWidth',lw);
hold on
plot(relevantInterval, output_Of_Interest(:,7)/yrDist,':','color',colorMat(2,:),'LineWidth',lw);
plot(relevantInterval, output_Of_Interest(:,8)/yrDist,':','color',colorMat(3,:),'LineWidth',lw);
plot(relevantInterval, output_Of_Interest(:,9)/yrDist,':','color',colorMat(4,:),'LineWidth',lw);
%xlim([yrs_Plot(1) yrs_Plot(end)])
%ylim([1000 13000])
%plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
xlabel(interventionString, 'Interpreter','latex','FontSize',50);
ylabel(strcat('Avg num. on PrEP ', yr1String,'-',yr2String), 'Interpreter','latex','FontSize',50);
legend(outcomeText{6:9},'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',1);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;


nexttile
plot(relevantInterval, output_Of_Interest(:,10)/yrDist,':','color',colorMat(1,:),'LineWidth',lw);
hold on
plot(relevantInterval, output_Of_Interest(:,11)/yrDist,':','color',colorMat(2,:),'LineWidth',lw);
plot(relevantInterval, output_Of_Interest(:,12)/yrDist,':','color',colorMat(3,:),'LineWidth',lw);
plot(relevantInterval, output_Of_Interest(:,13)/yrDist,':','color',colorMat(4,:),'LineWidth',lw);
%xlim([yrs_Plot(1) yrs_Plot(end)])
%ylim([1000 13000])
%plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
xlabel(interventionString, 'Interpreter','latex','FontSize',50);
ylabel(strcat('Avg \% elig. on PrEP ', yr1String,'-',yr2String), 'Interpreter','latex','FontSize',50);
legend(outcomeText{10:13},'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'AutoUpdate','off','NumColumns',1);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;

% 
% 
% nexttile
% 
% plot(relevantInterval*100, 100*output_Of_Interest(:,8),':','color',colorMat(3,:),'LineWidth',lw);
% hold on
% %plot(yrs_Plot, Vi1_All(1,:) ,':','color',colorMat(1,:),'LineWidth',lw);
% %xlim([yrs_Plot(1) yrs_Plot(end)])
% ylim([70 100])
% 
% %plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
% xlabel(interventionString, 'Interpreter','latex','FontSize',50);
% ylabel('2034 \% on ART [-]', 'Interpreter','latex','FontSize',50);
% %ylabel('Variance', 'Interpreter','latex','FontSize',50);
% legend({'\% on ART'},'Interpreter','latex','FontSize',legendfs,'Location','best',...
%     'AutoUpdate','off','NumColumns',1);
% 
% %plot(Sr.knots*100, 100*output_Of_Interest_KnotPts(:,8),'o','color',colorMat(3,:),'LineWidth',lw);
% 
% %for k=2--size(survPrevRange,1)
% %    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
% %end
% %title('Intervention effects on incidence', 'Interpreter','latex','FontSize',50);
% set(gca,'FontSize',fs);
% ax=gca;
% ax.YAxis.Exponent=0;
% 
% 
% nexttile
% 
% plot(relevantInterval*100, 100*output_Of_Interest(:,9),':','color',colorMat(5,:),'LineWidth',lw);
% hold on
% plot(relevantInterval*100, 100*output_Of_Interest(:,10),':','color',colorMat(1,:),'LineWidth',lw);
% plot(relevantInterval*100, 100*output_Of_Interest(:,11),':','color',colorMat(2,:),'LineWidth',lw);
% plot(relevantInterval*100, 100*output_Of_Interest(:,12),':','color',colorMat(4,:),'LineWidth',lw);
% %xlim([yrs_Plot(1) yrs_Plot(end)])
% ylim([50 100])
% 
% %plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
% xlabel(interventionString, 'Interpreter','latex','FontSize',50);
% ylabel('2034 \% VLS among  [-]', 'Interpreter','latex','FontSize',50);
% %ylabel('Variance', 'Interpreter','latex','FontSize',50);
% legend({'Total','Black/AA','Hisp/Latino','Other'},'Interpreter','latex','FontSize',legendfs,'Location','best',...
%     'AutoUpdate','off','NumColumns',1);
% %plot(Sr.knots*100, 100*output_Of_Interest_KnotPts(:,9),'o','color',colorMat(5,:),'LineWidth',lw);
% %plot(Sr.knots*100, 100*output_Of_Interest_KnotPts(:,10),'o','color',colorMat(1,:),'LineWidth',lw);
% %plot(Sr.knots*100, 100*output_Of_Interest_KnotPts(:,11),'o','color',colorMat(2,:),'LineWidth',lw);
% %plot(Sr.knots*100, 100*output_Of_Interest_KnotPts(:,12),'o','color',colorMat(4,:),'LineWidth',lw);
% 
% %for k=2--size(survPrevRange,1)
% %    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
% %end
% %title('Intervention effects on incidence', 'Interpreter','latex','FontSize',50);
% set(gca,'FontSize',fs);
% ax=gca;
% ax.YAxis.Exponent=0;
% 


set(gcf,'Position',posMat);
