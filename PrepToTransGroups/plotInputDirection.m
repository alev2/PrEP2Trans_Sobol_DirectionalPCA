function [] = plotInputDirection(inputDirection,interventionNames,tileDim)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

    tp=tiledlayout(tileDim(1),tileDim(2),'TileSpacing','compact','Padding','compact');

    fs=12;

    for i=1:size(inputDirection,2)

        nexttile
        bar(interventionNames,.1*(inputDirection(:,i)/max(abs(inputDirection(:,i)))))
        title(strcat('InputPC:\,\,',num2str(i)),'Interpreter','latex','FontSize',50);
        set(gca,'FontSize',fs);
        ax=gca;
        ax.YAxis.Exponent=0;
        ylim([-.125 .125])

    end

   
    posMat=[-1845,251,1628,441];
    %set(gcf,'Position',posMat);
end