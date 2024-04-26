%build the table
%readAndBuildTable

fs=14;
fs2=14;
posMat=[-1845,251,1628,441];

interventionNames=categorical({...
    'Mixing level',...
    'PrEP, HETM',...
    'PrEP, HETF',...
    'PrEP, MSM',...
    'PrEP, PWID',...
});


interventionNames=reordercats(interventionNames,{...
    'Mixing level',...
    'PrEP, HETM',...
    'PrEP, HETF',...
    'PrEP, MSM',...
    'PrEP, PWID',...
} );


outputNames=categorical({...
    'NumPrEPMSM',...
    'NumPrEPHETF',...
    'NumPrEPHETM',...
    'NumPrEPPWID',...
 });

outputNames=reordercats(outputNames,{...
    'NumPrEPMSM',...
    'NumPrEPHETF',...
    'NumPrEPHETM',...
    'NumPrEPPWID',...
});


%low-rank approximation - number of "meta variables" we want to keep
apx_Rank=4;

%normalize the model outputs
rowMeans=mean(modelOutputs,2)*ones(1,size(modelOutputs,2));
stdMat=diag(std(modelOutputs,0,2));
modelOutputsNorm=stdMat\(modelOutputs-rowMeans);

%perform SVD on model outputs to get the principal components
[UFull,SFull,VFull]=svd(modelOutputsNorm);

%low rank approximation
UApx=UFull(:,1:apx_Rank); 
SApx=SFull(1:apx_Rank,1:apx_Rank);
VApx=VFull(:,1:apx_Rank);

modelOutputsNorm_Apx=UApx*SApx*VApx';
modelOutputs_Apx=stdMat*modelOutputsNorm_Apx + rowMeans;

%
reduced_Outputs=(UApx')*modelOutputsNorm_Apx;




%[Sob_r1,Tot_Sob_r1,m1,v1] = compute_sobol_indices_from_sparse_grid(S,Sr,reduced_Outputs(1,:),domain,'legendre'); 
%[Sob_r2,Tot_Sob_r2,m2,v2] = compute_sobol_indices_from_sparse_grid(S,Sr,reduced_Outputs(2,:),domain,'legendre'); 
%[Sob_r3,Tot_Sob_r3,m3,v3] = compute_sobol_indices_from_sparse_grid(S,Sr,reduced_Outputs(3,:),domain,'legendre'); 
[Sob_r1,Tot_Sob_r1,m1,v1] = compute_variance_decomposition_from_sparse_grid(S,Sr,reduced_Outputs(1,:),domain,'legendre'); 
[Sob_r2,Tot_Sob_r2,m2,v2] = compute_variance_decomposition_from_sparse_grid(S,Sr,reduced_Outputs(2,:),domain,'legendre'); 
[Sob_r3,Tot_Sob_r3,m3,v3] = compute_variance_decomposition_from_sparse_grid(S,Sr,reduced_Outputs(3,:),domain,'legendre'); 
[Sob_r4,Tot_Sob_r4,m4,v4] = compute_variance_decomposition_from_sparse_grid(S,Sr,reduced_Outputs(4,:),domain,'legendre'); 
%[Sob_r5,Tot_Sob_r5,m5,v5] = compute_sobol_indices_from_sparse_grid(S,Sr,reduced_Outputs(5,:),domain,'legendre'); 

% [Sob_r1,Tot_Sob_r1,m1,v1] =compute_variance_decomposition_from_sparse_grid(S,Sr,reduced_Outputs(1,:),domain,'legendre');
% [Sob_r2,Tot_Sob_r2,m2,v2] =compute_variance_decomposition_from_sparse_grid(S,Sr,reduced_Outputs(2,:),domain,'legendre');
% [Sob_r3,Tot_Sob_r3,m3,v3] =compute_variance_decomposition_from_sparse_grid(S,Sr,reduced_Outputs(3,:),domain,'legendre');
% [Sob_r4,Tot_Sob_r4,m4,v4] =compute_variance_decomposition_from_sparse_grid(S,Sr,reduced_Outputs(4,:),domain,'legendre');


colorMat=[...
    .000 .447 .741;
    .850 .325 .098;
    .929 .694 .125;
    .494 .184 .556;
    .466 .674 .188;
    .301 .745 .933;
    .635 .078 .184;
    ];




%tp=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
tp=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile
bar(outputNames,UApx(:,1)*SApx(1,1));
ylabel('PC Loading', 'Interpreter','latex','FontSize',50);
ylim([-25 25])
%xticklabels(interventionNames);
title('PC1', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;

nexttile
bar(outputNames,UApx(:,2)*SApx(2,2));
ylabel('PC Loading', 'Interpreter','latex','FontSize',50);
ylim([-25 25])
%xticklabels(interventionNames);
title('PC2', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;


nexttile
bar(outputNames,UApx(:,3)*SApx(3,3));
ylabel('PC Loading', 'Interpreter','latex','FontSize',50);
ylim([-25 25])
%xticklabels(interventionNames);
title('PC3', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;

nexttile
bar(outputNames,UApx(:,4)*SApx(4,4));
ylabel('PC Loading', 'Interpreter','latex','FontSize',50);
ylim([-25 25])
%xticklabels(interventionNames);
title('PC4', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs);
ax=gca;




set(gcf,'Position',posMat);
figure




%tp=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
tp=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile

bar(interventionNames, Tot_Sob_r1)
hold on
%bh=bar([13 15], Tot_Sob_r1([13 15]),'FaceColor',colorMat(2,:))

ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
ylim([0 4])



%xticklabels(interventionNames);
title('PC1', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs2);
ax=gca;

nexttile

bar(interventionNames,Tot_Sob_r2)
hold on
%bh=bar([2 14 16], Tot_Sob_r2([2 14 16]),'FaceColor',colorMat(2,:))

ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
ylim([0 1])


%xticklabels(interventionNames);
title('PC2', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs2);
ax=gca;

nexttile


bar(interventionNames,Tot_Sob_r3)
hold on
%bh=bar([1], Tot_Sob_r3([1]),'FaceColor',colorMat(2,:))


ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
ylim([0 1])


%xticklabels(interventionNames);
title('PC3', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs2);
ax=gca;



nexttile

bar(interventionNames,Tot_Sob_r4)
ylabel('Total sobol index', 'Interpreter','latex','FontSize',50);
ylim([0 1])


%xticklabels(interventionNames);
title('PC4', 'Interpreter','latex','FontSize',50);
set(gca,'FontSize',fs2);
ax=gca;



set(gcf,'Position',posMat);