function [] = plotIdentified_InterventionDirection(baselineDirection, inputDirection,interventionNames,tileDim)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

    tp=tiledlayout(tileDim(1),tileDim(2),'TileSpacing','compact','Padding','compact');

    fs=21;
%     legString={...        
%     '\mbox{\boldmath$\theta^0$} + $\nabla$\mbox{\boldmath$\widetilde\Phi$} $($\mbox{\boldmath$\theta^0$}$)^{\dag}U_{1:4}^T$\mbox{\boldmath$y$}$_{\,Target\,1}$',...
%     '\mbox{\boldmath$\theta^0$} + $\nabla$\mbox{\boldmath$\widetilde\Phi$} $($\mbox{\boldmath$\theta^0$}$)^{\dag}U_{1:4}^T$\mbox{\boldmath$y$}$_{\,Target\,2}$',...
%     '\mbox{\boldmath$\theta^0$} + $\nabla$\mbox{\boldmath$\widetilde\Phi$} $($\mbox{\boldmath$\theta^0$}$)^{\dag}U_{1:4}^T$\mbox{\boldmath$y$}$_{\,Target\,3}$',...
%     '\mbox{\boldmath$\theta^0$} + $\nabla$\mbox{\boldmath$\widetilde\Phi$} $($\mbox{\boldmath$\theta^0$}$)^{\dag}U_{1:4}^T$\mbox{\boldmath$y$}$_{\,Target\,4}$',...
%     };
   legString={};
    for i=1:size(inputDirection,2)
        legString=[legString ...
            strcat(...
                     '\mbox{\boldmath$\theta^0$} + $\nabla$\mbox{\boldmath$\widetilde\Phi$} $($\mbox{\boldmath$\theta^0$}$)^{\dag}U_{1:4}^T$\mbox{\boldmath$y$}$_{\,Target\,',...
                      num2str(i),...
                     '}$'...
                   )...
                   ];
    end

    %        nexttile
        
        %inputDirection(5:7,:)=1-inputDirection(5:7,:);
        %inputDirection(5:7,2)=1-inputDirection(5:7,2);
        bar(interventionNames,inputDirection*100)
        title('Input-output ROM: Identified intervention comparison','Interpreter','latex','FontSize',50);
        
        legend(legString,'Interpreter','latex','FontSize',fs);
        ylabel('\% increase from baseline level [-]','FontSize',fs,'Interpreter','Latex')
        set(gca,'FontSize',fs);
        ax=gca;
        ax.YAxis.Exponent=0;
      ylim([-25 100])

    %end

   
    posMat=[-1651,167,1044,619];
    set(gcf,'Position',posMat);
end