close all
clear

% initialize paramaters that we do not vary for HIV problem, as well as
% the sparse grid library.
sobolExamples_Setup

fprintf('setup done\n');

% Sobol indices computation

N = 16; % number of parameters...

%%PrEP_Params (multiplier, default is 1 for both):  tt_PrEPInitRateOrPctOnPrEP_Set5
prepMult_Black=[0, 1]; %inds 1, 4, 7, 10
prepMult_Hisp=[0, 1]; %inds 2, 5, 8, 11


%testing (an RR)
move1to2_Black=[1 2];
move1to2_Hisp=[.864 2];


%LTC at diag
move2to3_atDiag_Black=[.79 .98];
move2to3_atDiag_Hisp=[.83 .98];

%LTC after diag
move2to3_afterDiag_Black=[.15 .3];
move2to3_afterDiag_Hisp=[.16 .3];

%initiate ART (given as RR)
move3to4_Black=[1.0 2.0];
move3to4_Hisp=[1.0 2.0];

%Drop off ART params: tt_dropOutRate_ANVToCare_5
move4to3_Black=[.02, .246];
move4to3_Hisp=[.02, .252];

%Get VLS params: tt_ARTInitRateFromANV_5
move4to5_Black=[.531, .98];
move4to5_Hisp=[.473, .98];

%LoseVLS Params  tt_dropOutRate_VLSToANV_5
move5to4_Black=[.02, .2502];
move5to4_Hisp=[.02, .137];

num_Yrs=length(yearInds);

% CHANGE HERE! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int_p1 = prepMult_Black; % interval of variation for parameter 1 
int_p2 = prepMult_Hisp; % interval of variation for parameter 2
int_p3 = move1to2_Black; % interval of variation for parameter 3
int_p4 = move1to2_Hisp; % interval of variation for parameter 4
int_p5 = move2to3_atDiag_Black; % interval of variation for parameter 5
int_p6 = move2to3_atDiag_Hisp; % interval of variation for parameter 6
int_p7 = move2to3_afterDiag_Black; % interval of variation for parameter 7
int_p8 = move2to3_afterDiag_Hisp; % interval of variation for parameter 8
int_p9 = move3to4_Black; % interval of variation for parameter 9
int_p10 = move3to4_Hisp; % interval of variation for parameter 10
int_p11 = move4to3_Black; % interval of variation for parameter 11
int_p12 = move4to3_Hisp; % interval of variation for parameter 12
int_p13 = move4to5_Black; % interval of variation for parameter 13
int_p14 = move4to5_Hisp; % interval of variation for parameter 14
int_p15 = move5to4_Black; % interval of variation for parameter 15
int_p16 = move5to4_Hisp; % interval of variation for parameter 16


domain = [int_p1',int_p2',int_p3', int_p4', int_p5', int_p6' int_p7' int_p8',...
          int_p9' int_p10' int_p11' int_p12' int_p13' int_p14' int_p15' int_p16' ]; % define the N-dim. domain - we need this in the function that computes the Sobol indices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% type of knots to be used for the construction of the sparse grid
% (we use for symmetric Leja points for the moment, this should be ok!)
knots_1 = @(n) knots_leja(n,int_p1(1),int_p1(2),'sym_line');
knots_2 = @(n) knots_leja(n,int_p2(1),int_p2(2),'sym_line');
knots_3 = @(n) knots_leja(n,int_p3(1),int_p3(2),'sym_line');
knots_4 = @(n) knots_leja(n,int_p4(1),int_p4(2),'sym_line');
knots_5 = @(n) knots_leja(n,int_p5(1),int_p5(2),'sym_line');
knots_6 = @(n) knots_leja(n,int_p6(1),int_p6(2),'sym_line');
knots_7 = @(n) knots_leja(n,int_p7(1),int_p7(2),'sym_line');
knots_8 = @(n) knots_leja(n,int_p8(1),int_p8(2),'sym_line');
knots_9 = @(n) knots_leja(n,int_p9(1),int_p9(2),'sym_line');
knots_10 = @(n) knots_leja(n,int_p10(1),int_p10(2),'sym_line');
knots_11 = @(n) knots_leja(n,int_p11(1),int_p11(2),'sym_line');
knots_12 = @(n) knots_leja(n,int_p12(1),int_p12(2),'sym_line');
knots_13 = @(n) knots_leja(n,int_p13(1),int_p13(2),'sym_line');
knots_14 = @(n) knots_leja(n,int_p14(1),int_p14(2),'sym_line');
knots_15 = @(n) knots_leja(n,int_p15(1),int_p15(2),'sym_line');
knots_16 = @(n) knots_leja(n,int_p16(1),int_p16(2),'sym_line');


% combine the knots for the three parameters in one - we need this as input to the function that creates the sparse grid
knots = {knots_1,knots_2,knots_3, knots_4, knots_5, knots_6, knots_7, knots_8...
         knots_9 knots_10 knots_11 knots_12 knots_13 knots_14 knots_15 knots_16}; 
fprintf('all sparsegrid prelim stuff done\n');

% level-to-knots function chosen in accordance with the type of knots above
lev2knots = @ lev2knots_2step;

% type of underlying polynomial approximation 
rule = @(i) sum(i-1);

% level (basically dictating the accuracy level of the sparse grid approximation)
w = 2; 

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

outputRoot='./datasetInput/simulation_';
%outputRoot='./datasetInput/simulation_';
%outputRoot='./dataset20Years/simulation_';
Alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
fprintf(strcat('Num points: ',num2str(nb_pts),'.\n'));


%outcomePP=HIVEpiModel_Sobol('HOPE Model V10_05_LB20230224_2_20231005_Fresh.xlsm', sg_pts);
prepIncreases_Black=...
    [
        .0628;
        .061;
        .4817;
        .0608;
    ];

prepIncreases_Hisp=...
    [
        .0152;
        .0151;
        .6809;
        .0293;
    ];

load('ExcelValues_AllParameters.mat');
load('ExcelValues_Populations.mat')
%Plot_Sobols



for i=1:nb_pts
    tic

    Parameters_Cur=ExcelValues_AllParameters;

    prepMult_Black_SimVal=sg_pts(1,i);
    prepMult_Hisp_SimVal=sg_pts(2,i);
    

    %change prep
    fprintf('Changing PrEP parameters... ');
    Parameters_Cur(indsPrEPBlack)=...
            ExcelValues_AllParameters(indsPrEPBlack)+prepIncreases_Black*prepMult_Black_SimVal;
    Parameters_Cur(indsPrEPHisp)=...
            ExcelValues_AllParameters(indsPrEPHisp)+prepIncreases_Hisp*prepMult_Hisp_SimVal;
    fprintf('done.\n');
    

    %change testing
    fprintf('Changing 1 to 2 parameters (testing)... ');    
    move1to2_Black_SimVal=sg_pts(3,i);
    move1to2_Hisp_SimVal=sg_pts(4,i);

    Parameters_Cur(inds1to2(1))=move1to2_Black_SimVal;
    Parameters_Cur(inds1to2(2))=move1to2_Hisp_SimVal;
    fprintf('done.\n');


    %change 2 to 3 at diag
    fprintf('Changing 2 to 3 at diag parameters... ');    
    move2to3_atDiag_Black_SimVal=sg_pts(5,i);
    move2to3_atDiag_Hisp_SimVal=sg_pts(6,i);

    Parameters_Cur(inds2to3_atDiag(1))=move2to3_atDiag_Black_SimVal;
    Parameters_Cur(inds2to3_atDiag(2))=move2to3_atDiag_Hisp_SimVal;
    fprintf('done.\n');


    %change 2 to 3 after diag
    fprintf('Changing 2 to 3 after diag parameters... ');    
    move2to3_afterDiag_Black_SimVal=sg_pts(7,i);
    move2to3_afterDiag_Hisp_SimVal=sg_pts(8,i);

    Parameters_Cur(inds2to3_afterDiag(1))=move2to3_afterDiag_Black_SimVal;
    Parameters_Cur(inds2to3_afterDiag(2))=move2to3_afterDiag_Hisp_SimVal;
    fprintf('done.\n');


    %change 3 to 4
    fprintf('Changing 3 to 4 parameters... ');    
    move3to4_Black_SimVal=sg_pts(9,i);
    move3to4_Hisp_SimVal=sg_pts(10,i);

    Parameters_Cur(inds3to4(1))=move3to4_Black_SimVal;
    Parameters_Cur(inds3to4(2))=move3to4_Hisp_SimVal;
    fprintf('done.\n');


    %change 4 to 3
    fprintf('Changing 4 to 3 parameters... ');    
    move4to3_Black_SimVal=sg_pts(11,i);
    move4to3_Hisp_SimVal=sg_pts(12,i);

    Parameters_Cur(inds4to3(1))=move4to3_Black_SimVal;
    Parameters_Cur(inds4to3(2))=move4to3_Hisp_SimVal;
    fprintf('done.\n');
    

    %chnage 4 to 5 parameters
    fprintf('Changing 4 to 5 parameters... ');    
    move4to5_Black_SimVal=sg_pts(13,i);
    move4to5_Hisp_SimVal=sg_pts(14,i);

    Parameters_Cur(inds4to5(1))=move4to5_Black_SimVal;
    Parameters_Cur(inds4to5(2))=move4to5_Hisp_SimVal;
    fprintf('done.\n');


    %change 5 to 4 parameters
    fprintf('Changing 5 to 4 parameters... ');    
    move5to4_Black_SimVal=sg_pts(15,i);
    move5to4_Hisp_SimVal=sg_pts(16,i);

    Parameters_Cur(inds5to4(1))=move5to4_Black_SimVal;
    Parameters_Cur(inds5to4(2))=move5to4_Hisp_SimVal;
    fprintf('done.\n');


    fprintf('Running HOPE... ');
    [outTest, Params]=HIVEpiModel_Sobol(Parameters_Cur, ExcelValues_Populations);
    fprintf('done.\n ');





    ann_NewInfections_Blk=sum(outTest.ann_NewInfections_Blk,2);
    ann_NewInfections_Hisp=sum(outTest.ann_NewInfections_Hisp,2);
    ann_NewInfections_Oth=sum(outTest.ann_NewInfections_Oth,2);
    
    ann_TotalNewInfections_PerPerson_Blk=ann_NewInfections_Blk./outTest.ann_popSize_Blk;
    ann_TotalNewInfections_PerPerson_Hisp=ann_NewInfections_Hisp./outTest.ann_popSize_Hisp;
    ann_TotalNewInfections_PerPerson_Oth=ann_NewInfections_Oth./outTest.ann_popSize_Oth;
    
    irr_Blk=ann_TotalNewInfections_PerPerson_Blk./ann_TotalNewInfections_PerPerson_Oth;
    irr_Hisp=ann_TotalNewInfections_PerPerson_Hisp./ann_TotalNewInfections_PerPerson_Oth;
    
    spending=outTest.ann_ARTCarePrEPTransitionAndSEPCost_Disc;

    outputMatrix=[ann_NewInfections_Blk ann_NewInfections_Hisp irr_Blk irr_Hisp spending];
       
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

