function [] = plotInterventionDirection(baselineDirection, inputDirection,interventionNames,tileDim)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

    tp=tiledlayout(tileDim(1),tileDim(2),'TileSpacing','compact','Padding','compact');

    fs=12;

%     nexttile
%     bar(interventionNames,baselineDirection)
%     title('BaselineDirection','Interpreter','latex','FontSize',50);
%     set(gca,'FontSize',fs);
%     ax=gca;
%     ax.YAxis.Exponent=0;
%     ylim([-.1 1.2])


    for i=1:size(inputDirection,2)

        nexttile
        bar(interventionNames,.1*inputDirection(:,i)/max(abs(inputDirection(:,i))))
        title(strcat('{\boldmath$\tau_{',num2str(i),'}$}'),'Interpreter','latex','FontSize',70);
        set(gca,'FontSize',fs);
        ax=gca;
        ax.YAxis.Exponent=0;
       ylim([-.12 .12])

    end

   posMat=[-1651,167,1044,619];
%    posMat=[-1651,210,1045,576];
    set(gcf,'Position',posMat);
end