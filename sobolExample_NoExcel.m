close all
clear

% initialize paramaters that we do not vary for HIV problem, as well as
% the sparse grid library.
sobolExamples_Setup

fprintf('setup done\n');

% Sobol indices computation

N = 1; % number of parameters...

%%PrEP_Params (multiplier, default is 1 for both):  tt_PrEPInitRateOrPctOnPrEP_Set5
%prepMult_Black=[0, 1]; %inds 1, 4, 7, 10
%prepMult_Hisp=[0, 1]; %inds 2, 5, 8, 11

%LoseVLS Params  tt_dropOutRate_VLSToANV_5
%move5to4_Black=[.02*.2502 .2502];
%move5to4_Hisp=[.02*.137 .137];

%Get VLS params: tt_ARTInitRateFromANV_5
%move4to5_Black=[.531, .98];
%move4to5_Hisp=[.473, .98];


%Drop off ART params: tt_dropOutRate_ANVToCare_5
%move4to3_Black=[.02*.246, .246];
%move4to3_Hisp=[.02*.252, .252];

%long ART params
longARTPct = [0 1];

num_Yrs=length(yearInds);

% CHANGE HERE! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int_p1 = longARTPct; % interval of variation for parameter 1 
%int_p2 = prepMult_Hisp; % interval of variation for parameter 2
%int_p3 = move5to4_Black; % interval of variation for parameter 3
%int_p4 = move5to4_Hisp; % interval of variation for parameter 3
%int_p5 = move4to5_Black; % interval of variation for parameter 3
%int_p6 = move4to5_Hisp; % interval of variation for parameter 3
%int_p7 = move4to3_Black; % interval of variation for parameter 3
%int_p8 = move4to3_Hisp; % interval of variation for parameter 3

domain = [int_p1',int_p2',int_p3', int_p4', int_p5', int_p6' int_p7' int_p8']; % define the N-dim. domain - we need this in the function that computes the Sobol indices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% type of knots to be used for the construction of the sparse grid
% (we use for symmetric Leja points for the moment, this should be ok!)
knots_1 = @(n) knots_leja(n,int_p1(1),int_p1(2),'sym_line');
%knots_2 = @(n) knots_leja(n,int_p2(1),int_p2(2),'sym_line');
%knots_3 = @(n) knots_leja(n,int_p3(1),int_p3(2),'sym_line');
%knots_4 = @(n) knots_leja(n,int_p4(1),int_p4(2),'sym_line');
%knots_5 = @(n) knots_leja(n,int_p5(1),int_p5(2),'sym_line');
%knots_6 = @(n) knots_leja(n,int_p6(1),int_p6(2),'sym_line');
%knots_7 = @(n) knots_leja(n,int_p7(1),int_p7(2),'sym_line');
%knots_8 = @(n) knots_leja(n,int_p8(1),int_p8(2),'sym_line');

% combine the knots for the three parameters in one - we need this as input to the function that creates the sparse grid
knots = {knots_1};%,knots_2,knots_3, knots_4, knots_5, knots_6, knots_7, knots_8}; 
fprintf('all sparsegrid prelim stuff done\n');

% level-to-knots function chosen in accordance with the type of knots above
lev2knots = @ lev2knots_2step;

% type of underlying polynomial approximation 
rule = @(i) sum(i-1);

% level (basically dictating the accuracy level of the sparse grid approximation)
w = 3; 

% create the sparse grid
fprintf('Creating sparse grid... ');

S = create_sparse_grid(N,w,knots,lev2knots,rule);
fprintf('done.\n');
% visualization: plot the sparse grid
%plot_sparse_grid(S,[1 2]) % first two parameters
%plot_sparse_grid(S,[2 3]) % second and third parameters

% create the reduced sparse grid to get a list of knots with no repetition
fprintf('Reducing sparse grid... ');
Sr = reduce_sparse_grid(S); 
fprintf('done.\n');


fprintf('Extracting points... ');
% extract the list of points 
fprintf('done.\n');
sg_pts = Sr.knots; % this is a matrix of dimension N x nb. of sparse grid points (each column of the matrix is one point) 
nb_pts = size(sg_pts,2); % variable to know how many points we have at hand


% DO SOMETHING HERE! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Create a .xls file for each of this points and collect the results of
% your model run in a table 
% - !! Order is important!! Store the model run according to the ordering of
% the points in sg_pts. 
% - Assuming 3 groups and results for 5 years, the table contains 15 data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outputRoot='./datasetInputTest_AllStuff/simulation_';
%outputRoot='./datasetInput/simulation_';
%outputRoot='./dataset20Years/simulation_';
Alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
fprintf(strcat('Num points: ',num2str(nb_pts),'.\n'));


min_MixMat=...
    [
        .99999; %loHETM -> HETF
        9.99e-6; %loHETM -> PWID
        .99999; %loHETF -> HETM
        9.99e-6; %loHETF -> MSM
         0.00; %loHETF -> PWID
        .975; %hiHETM -> HETF
        .025; %hiHETM -> PWID
        .9964764333; %hiHETF -> HETM
        .00087084249099998; %hiHETF -> MSM
        .002652724209;%hiHETF -> PWID
        .0914023621; %MSM -> HETF
        .906105747187; %MSM -> MSM
        .00249189071300004; %MSM -> PWID

    ];

max_MixMat=...
    [
        .997; %loHETM -> HETF
        .003; %loHETM -> PWID
        .997; %loHETF -> HETM
        .003; %loHETF -> MSM
         0.00; %loHETF -> PWID
        .935; %hiHETM -> HETF
        .065; %hiHETM -> PWID
        .9564764333; %hiHETF -> HETM
        .040871; %hiHETF -> MSM
        .002652724209;%hiHETF -> PWID
        .1714023621; %MSM -> HETF
        .826105747187; %MSM -> MSM
        .00249189071300004; %MSM -> PWID
    ];


mix_Delta=max_MixMat-min_MixMat;



%outcomePP=HIVEpiModel_Sobol('HOPE Model V10_05_LB20230224_2_20231005_Fresh.xlsm', sg_pts);

load('ExcelValues_AllParameters.mat');
load('ExcelValues_Populations.mat')
%Plot_Sobols

for i=1:nb_pts
    tic

    Parameters_Cur=ExcelValues_AllParameters;
    Parameters_Cur(indsDiscRate)=0.;
    mixing_Factor=sg_pts(1,i);
    
    mixing_Matrix=min_MixMat+mixing_Factor*mix_Delta;
    %prepMult_Hisp_SimVal=sg_pts(2,i);
    

    fprintf('Changing mixing matrix... ');
    
    Parameters_Cur(mixingInds)=mixing_Matrix;
    


    fprintf('Running HOPE... ');
    [outTest, Params]=HIVEpiModel_Sobol(Parameters_Cur, ExcelValues_Populations);
    fprintf('done.\n ');
    
    
    ann_NewInfections_MSM=sum(outTest.ann_NewInfectionsMSM,2);
    ann_NewInfections_PWID=sum(outTest.ann_NewInfectionsIDU,2);
    ann_NewInfections_HETM=sum(outTest.ann_NewInfectionsHETM,2);
    ann_NewInfections_HETF=sum(outTest.ann_NewInfectionsHETF,2);
    ann_NewInfections_Total=IncidenceTotal=sum(Outcome.ann_TotalNewInfections,2);
    
    

    spending=outTest.ann_ARTCarePrEPTransitionAndSEPCost_Disc;

    outputMatrix=[ann_NewInfections_MSM ann_NewInfections_HETF ann_NewInfections_HETM ann_NewInfections_PWID ann_NewInfections_Total];       
       
    incidenceOut=array2table([outputMatrix(yearInds,:)]);% annualDeaths(1:end,:) pctAware(1:end,:)]);
    
    fprintf('Writing table... ');    
    writetable(incidenceOut,strcat(outputRoot,num2str(i),'.xls'),'Range',...
       strcat(Alphabet(1),num2str(1),':',Alphabet(size(outputMatrix,2)),num2str(size(incidenceOut,1)+1)));%,'WriteVariableNames',false);
    fprintf('done.\n ');

    toc
end



% % CHANGE HERE! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%num_Yrs = 5; 
% read from your files and arrange the results in a convenient way
% we store the results for each group in one table 
values_g1 = 0*ones(num_Yrs,nb_pts); % incidence black
values_g2 = 0*ones(num_Yrs,nb_pts); % incidence hisp
values_g3 = 0*ones(num_Yrs,nb_pts); % irr black
values_g4 = 0*ones(num_Yrs,nb_pts); % irr hisp
values_g5 = 0*ones(num_Yrs,nb_pts); % spending

% read from file and copy the results in the tables
fprintf('Reading table... ');
%for k=1:nb_pts
for k=1:nb_pts
    
%    values_g1(:,k) = xlsread('xxx.xls', ...); % if you give another input to xlsread you can select columns in the file 
%    values_g2(:,k) = xlsread('xxx.xls', ...);
%    values_g3(:,k) = xlsread('xxx.xls', ...);
%    tableIn=table2array(readtable(strcat(outputRoot,num2str(k),'.xls')),...
%        'Range',strcat(Alphabet(1),num2str(1),':',Alphabet(1),num2str(size(incidenceOut,1)+1)));
    tableIn=table2array(readtable(strcat(outputRoot,num2str(k),'.xls')),...
        'Range',strcat(Alphabet(1),num2str(1),':',Alphabet(1),num2str(length(yearInds)+1)));

    values_g1(:,k) = tableIn(:,1); % if you give another input to xlsread you can select columns in the file 
    values_g2(:,k) = tableIn(:,2);
    values_g3(:,k) = tableIn(:,3);
    values_g4(:,k) = tableIn(:,4);
    values_g5(:,k) = tableIn(:,5);    
end
fprintf('done.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % compute Sobol indices at each year for each group
    
Sob_g1_All=[];
Sob_g2_All=[];
Sob_g3_All=[];
Sob_g4_All=[];
Sob_g5_All=[];

Tot_Sob_g1_All=[];
Tot_Sob_g2_All=[];
Tot_Sob_g3_All=[];
Tot_Sob_g4_All=[];
Tot_Sob_g5_All=[];

m1_All=[];
m2_All=[];
m3_All=[];
m4_All=[];
m5_All=[];


v1_All=[];
v2_All=[];
v3_All=[];
v4_All=[];
v5_All=[];

v1_new=[];
v2_new=[];
v3_new=[];
v4_new=[];
v5_new=[];

VTi1_All=[];
VTi2_All=[];
VTi3_All=[];
VTi4_All=[];
VTi5_All=[];

Vi1_All=[];
Vi2_All=[];
Vi3_All=[];
Vi4_All=[];
Vi5_All=[];

vti1=[];
vti2=[];
vti3=[];
vti4=[];
vti5=[];

vi1=[];vi2=[];vi3=[];vi4=[];vi5=[];

fprintf('Computing variance decomposition... ');
for j = 1:num_Yrs
      %quasi static
%   [Sob_g1,Tot_Sob_g1,m1,v1] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g1(j,:),domain,'legendre'); 
%   [Sob_g2,Tot_Sob_g2,m2,v2] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g2(j,:),domain,'legendre'); 
%   [Sob_g3,Tot_Sob_g3,m3,v3] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g3(j,:),domain,'legendre'); 
%   [Sob_g4,Tot_Sob_g4,m4,v4] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g4(j,:),domain,'legendre'); 
%   [Sob_g5,Tot_Sob_g5,m5,v5] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g5(j,:),domain,'legendre'); 


 %   [Sob_g1,Tot_Sob_g1,m1,v1,vi1,vti1] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi1_All, VTi1_All, v1_new, values_g1(j,:),domain,'legendre'); 
 %   [Sob_g2,Tot_Sob_g2,m2,v2,vi2,vti2] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi2_All, VTi2_All, v2_new, values_g2(j,:),domain,'legendre'); 
 %   [Sob_g3,Tot_Sob_g3,m3,v3,vi3,vti3] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi3_All, VTi3_All, v3_new, values_g3(j,:),domain,'legendre'); 
 %  [Sob_g4,Tot_Sob_g4,m4,v4,vi4,vti4] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi4_All, VTi4_All, v4_new, values_g4(j,:),domain,'legendre'); 
 %  [Sob_g5,Tot_Sob_g5,m5,v5,vi5,vti5] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi5_All, VTi5_All, v5_new, values_g5(j,:),domain,'legendre'); 

    [Sob_g1,Tot_Sob_g1,m1,v1] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g1(j,:),domain,'legendre'); 
    [Sob_g2,Tot_Sob_g2,m2,v2] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g2(j,:),domain,'legendre'); 
    [Sob_g3,Tot_Sob_g3,m3,v3] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g3(j,:),domain,'legendre'); 
    [Sob_g4,Tot_Sob_g4,m4,v4] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g4(j,:),domain,'legendre'); 
    [Sob_g5,Tot_Sob_g5,m5,v5] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g5(j,:),domain,'legendre'); 

    Sob_g1_All=[Sob_g1_All Sob_g1];
    Sob_g2_All=[Sob_g2_All Sob_g2];
    Sob_g3_All=[Sob_g3_All Sob_g3];
    Sob_g4_All=[Sob_g4_All Sob_g4];
    Sob_g5_All=[Sob_g5_All Sob_g5];
        
    Tot_Sob_g1_All=[Tot_Sob_g1_All Tot_Sob_g1];
    Tot_Sob_g2_All=[Tot_Sob_g2_All Tot_Sob_g2];
    Tot_Sob_g3_All=[Tot_Sob_g3_All Tot_Sob_g3];
    Tot_Sob_g4_All=[Tot_Sob_g4_All Tot_Sob_g4];    
    Tot_Sob_g5_All=[Tot_Sob_g5_All Tot_Sob_g5];   

    m1_All=[m1_All m1];
    m2_All=[m2_All m2];
    m3_All=[m3_All m3];
    m4_All=[m4_All m4];
    m5_All=[m5_All m5];

    v1_new=v1;
    v2_new=v2;
    v3_new=v3;
    v4_new=v4;
    v5_new=v5;

    v1_All=[v1_All v1];
    v2_All=[v2_All v2];
    v3_All=[v3_All v3];
    v4_All=[v4_All v4];
    v5_All=[v5_All v5];

    VTi1_All=vti1;%[VTi1_All vti1];
    VTi2_All=vti2;%[VTi2_All vti2];
    VTi3_All=vti3;%[VTi3_All vti3];
    VTi4_All=vti4;%[VTi3_All vti3];
    VTi5_All=vti5;%[VTi3_All vti3];

    Vi1_All=vi1;
    Vi2_All=vi2;%[Vi2_All vi2];
    Vi3_All=vi3;%[Vi3_All vi3];
    Vi4_All=vi4;%[Vi3_All vi3];
    Vi5_All=vi5;%[Vi3_All vi3];
    
    
end

save('Sob_g1_All','Sob_g1_All')
save('Sob_g2_All','Sob_g2_All')
save('Sob_g3_All','Sob_g3_All')
save('Sob_g4_All','Sob_g4_All')
save('Sob_g5_All','Sob_g5_All')

save('Tot_Sob_g1_All','Tot_Sob_g1_All')
save('Tot_Sob_g2_All','Tot_Sob_g2_All')
save('Tot_Sob_g3_All','Tot_Sob_g3_All')
save('Tot_Sob_g4_All','Tot_Sob_g4_All')
save('Tot_Sob_g5_All','Tot_Sob_g5_All')




fprintf('done.\n');


