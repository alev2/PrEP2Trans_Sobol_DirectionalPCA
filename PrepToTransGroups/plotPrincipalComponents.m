function [outputArg1,outputArg2] = plotPrincipalComponents(S)
colorMat=[...
    .000 .447 .741;
    .850 .325 .098;
    .929 .694 .125;
    .494 .184 .556;
    .466 .674 .188;
    .301 .745 .933;
    .635 .078 .184;
    .1235 .078 .484;
    .8235 .3278 .384];
posMat=[336,322,1328,492];
lw=2;
ms=10;
fs=18;
legendfs=18;
% A simple plotting tool that can show off the principal components.

    %if you give me a matrix rather than a vector
    if(size(S,2)>1)
    
        %if you give me mtx S from the SVD
        if(norm(diag(S))==norm(S,'fro'))
            S=diag(S);
        %otherwise compute its decomposition
        else
            [U,S,V]=svd(S);
            S=diag(S);
        end
        
    end

    %otherwise we were given the vector of singular values.

cumulativeSingularMass=cumsum(S)/sum(S);    

tp=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile
plot(S,'--d','LineWidth',lw,'MarkerSize',ms,'Color',colorMat(1,:));
xlabel('Principal component', 'Interpreter','latex','FontSize',50);
ylabel('PC Value', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;

nexttile
plot(S/sum(S),'--*','LineWidth',lw,'MarkerSize',ms,'Color',colorMat(2,:));
hold on
plot(cumulativeSingularMass,'--s','LineWidth',lw,'MarkerSize',ms,'Color',colorMat(3,:));
xlabel('Principal component', 'Interpreter','latex','FontSize',50);
ylabel('\% of total variance', 'Interpreter','latex','FontSize',50);
legend({'Individual','Cumulative'},'Interpreter','latex','FontSize',legendfs,'Location','best',...
    'NumColumns',1);
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;


set(gcf,'Position',posMat);

end