close all
clear

% initialize paramaters that we do not vary for HIV problem, as well as
% the sparse grid library.
sobolExamples_Setup
%outputRoot='./datasetInput_LongARTHiOrder/simulation_';
outputRoot='./Results_MixingAndPrEP_ExpandedOutputSpace_NoMixing/simulation_';
%outputRoot='./datasetInput_LargeScale2_Dec4/simulation_';

fprintf('setup done\n');

% Sobol indices computation

N = 4; % number of parameters...

%long ART reduction factor 
%longART_LoseVLS_ReductionFactor=.25;
%longART_ReachIncrease_AumentFactor=1.75;
%longART_LoseART_ReductionFactor=.25;

%lo mixing = .499 to .501

%Mixing level
%mixing_Range=[0 1];
hetm_PrepRange=[0 1];
hetf_PrepRange=[0 1];
msm_PrepRange=[0 1];
pwid_PrepRange=[0 1];
%pctLongART_4=[0 1];
%pctLongART_5=[0 1];

% %%PrEP_Params (multiplier, default is 1 for both):  tt_PrEPInitRateOrPctOnPrEP_Set5
% prepMult_Black=[0, 1]; %inds 1, 4, 7, 10
% prepMult_Hisp=[0, 1]; %inds 2, 5, 8, 11


num_Yrs=length(yearInds);

% CHANGE HERE! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%int_p1 = mixing_Range; % interval of variation for parameter 1 
int_p1 = hetm_PrepRange; % interval of variation for parameter 2
int_p2 = hetf_PrepRange; % interval of variation for parameter 3
int_p3 = msm_PrepRange; % interval of variation for parameter 4
int_p4 = pwid_PrepRange; % interval of variation for parameter 5
%int_p6 = move2to3_atDiag_Hisp; % interval of variation for parameter 6
%int_p7 = move2to3_afterDiag_Black; % interval of variation for parameter 7
%int_p8 = move2to3_afterDiag_Hisp; % interval of variation for parameter 8
%int_p9 = move3to4_Black; % interval of variation for parameter 9
%int_p10 = move3to4_Hisp; % interval of variation for parameter 10
%int_p11 = move4to3_Black; % interval of variation for parameter 11
%int_p12 = move4to3_Hisp; % interval of variation for parameter 12
% int_p3 = move4to5_Black; % interval of variation for parameter 13
% int_p4 = move4to5_Hisp; % interval of variation for parameter 14
% int_p5 = move5to4_Black; % interval of variation for parameter 15
% int_p6 = move5to4_Hisp; % interval of variation for parameter 16


domain = [int_p1',int_p2',int_p3', int_p4'];%, int_p5'];%, int_p6'];% int_p7' int_p8',...
%          int_p9' int_p10' int_p11' int_p12' int_p13' int_p14' int_p15' int_p16' ]; % define the N-dim. domain - we need this in the function that computes the Sobol indices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% type of knots to be used for the construction of the sparse grid
% (we use for symmetric Leja points for the moment, this should be ok!)
knots_1 = @(n) knots_leja(n,int_p1(1),int_p1(2),'sym_line');
knots_2 = @(n) knots_leja(n,int_p2(1),int_p2(2),'sym_line');
knots_3 = @(n) knots_leja(n,int_p3(1),int_p3(2),'sym_line');
knots_4 = @(n) knots_leja(n,int_p4(1),int_p4(2),'sym_line');
%knots_5 = @(n) knots_leja(n,int_p5(1),int_p5(2),'sym_line');
% knots_6 = @(n) knots_leja(n,int_p6(1),int_p6(2),'sym_line');
%knots_7 = @(n) knots_leja(n,int_p7(1),int_p7(2),'sym_line');
%knots_8 = @(n) knots_leja(n,int_p8(1),int_p8(2),'sym_line');
%knots_9 = @(n) knots_leja(n,int_p9(1),int_p9(2),'sym_line');
%knots_10 = @(n) knots_leja(n,int_p10(1),int_p10(2),'sym_line');
%knots_11 = @(n) knots_leja(n,int_p11(1),int_p11(2),'sym_line');
%knots_12 = @(n) knots_leja(n,int_p12(1),int_p12(2),'sym_line');
%knots_13 = @(n) knots_leja(n,int_p13(1),int_p13(2),'sym_line');
%knots_14 = @(n) knots_leja(n,int_p14(1),int_p14(2),'sym_line');
%knots_15 = @(n) knots_leja(n,int_p15(1),int_p15(2),'sym_line');
%knots_16 = @(n) knots_leja(n,int_p16(1),int_p16(2),'sym_line');


% combine the knots for the three parameters in one - we need this as input to the function that creates the sparse grid
knots = {knots_1,knots_2,knots_3, knots_4};%, knots_5};%, knots_6};%, knots_7, knots_8...
%         knots_9 knots_10 knots_11 knots_12 knots_13 knots_14 knots_15 knots_16}; 

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

%outputRoot='./datasetInput/simulation_';
%outputRoot='./dataset20Years/simulation_';
Alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
fprintf(strcat('Num points: ',num2str(nb_pts),'.\n'));

%outcomePP=HIVEpiModel_Sobol('HOPE Model V10_05_LB20230224_2_20231005_Fresh.xlsm', sg_pts);

%load('ExcelValues_AllParameters.mat');
%load('ExcelValues_Populations.mat')


% % CHANGE HERE! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%num_Yrs = 5; 
% read from your files and arrange the results in a convenient way
% we store the results for each group in one table 
values_g1 = 0*ones(num_Yrs,nb_pts); % incidence black
values_g2 = 0*ones(num_Yrs,nb_pts); % incidence hisp
values_g3 = 0*ones(num_Yrs,nb_pts); % irr black
values_g4 = 0*ones(num_Yrs,nb_pts); % irr hisp
values_g5 = 0*ones(num_Yrs,nb_pts); % spending
values_g6 = 0*ones(num_Yrs,nb_pts); % spending
values_g7 = 0*ones(num_Yrs,nb_pts); % spending
values_g8 = 0*ones(num_Yrs,nb_pts); % incidence black
values_g9 = 0*ones(num_Yrs,nb_pts); % incidence hisp
values_g10 = 0*ones(num_Yrs,nb_pts); % irr black
values_g11 = 0*ones(num_Yrs,nb_pts); % irr hisp
values_g12 = 0*ones(num_Yrs,nb_pts); % spending
values_g13 = 0*ones(num_Yrs,nb_pts); % spending

% read from file and copy the results in the tables
fprintf('Reading table... ');
%for k=1:nb_pts
for k=1:nb_pts
    
    tableIn=table2array(readtable(strcat(outputRoot,num2str(k),'.xls')),...
        'Range',strcat(Alphabet(1),num2str(1),':',Alphabet(1),num2str(length(yearInds)+1)));

    values_g1(:,k) = tableIn(:,1); % if you give another input to xlsread you can select columns in the file 
    values_g2(:,k) = tableIn(:,2);
    values_g3(:,k) = tableIn(:,3);
    values_g4(:,k) = tableIn(:,4);
    values_g5(:,k) = tableIn(:,5);    
    values_g6(:,k) = tableIn(:,6);    
    values_g7(:,k) = tableIn(:,7);    
    values_g8(:,k) = tableIn(:,8); % if you give another input to xlsread you can select columns in the file 
    values_g9(:,k) = tableIn(:,9);
    values_g10(:,k) = tableIn(:,10);
    values_g11(:,k) = tableIn(:,11);
    values_g12(:,k) = tableIn(:,12);    
    values_g13(:,k) = tableIn(:,13);    

%     for kk=1:size(incidenceOut,2) 
% 
%         values_(:,k,kk) = tableIn(:,kk); % if you give another input to xlsread you can select columns in the file 
%         if(kk==9)
%             values_(:,k,kk) =values_(:,k,kk)/1e9;
%         end
%         %values_(:,k) = tableIn(:,2);
%         %values_(:,k) = tableIn(:,3);
%         %values_(:,k) = tableIn(:,4);
%         %values_(:,k) = tableIn(:,5);    
%         %values_(:,k) = tableIn(:,6);    
%         %values_(:,k) = tableIn(:,7);    
% 
%     end

end
fprintf('done.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % compute Sobol indices at each year for each group
    
Sob_g1_All=[];
Sob_g2_All=[];
Sob_g3_All=[];
Sob_g4_All=[];
Sob_g5_All=[];
Sob_g6_All=[];
Sob_g7_All=[];
Sob_g8_All=[];
Sob_g9_All=[];
Sob_g10_All=[];
Sob_g11_All=[];
Sob_g12_All=[];
Sob_g13_All=[];

Tot_Sob_g1_All=[];
Tot_Sob_g2_All=[];
Tot_Sob_g3_All=[];
Tot_Sob_g4_All=[];
Tot_Sob_g5_All=[];
Tot_Sob_g6_All=[];
Tot_Sob_g7_All=[];
Tot_Sob_g8_All=[];
Tot_Sob_g9_All=[];
Tot_Sob_g10_All=[];
Tot_Sob_g11_All=[];
Tot_Sob_g12_All=[];
Tot_Sob_g13_All=[];

m1_All=[];
m2_All=[];
m3_All=[];
m4_All=[];
m5_All=[];
m6_All=[];
m7_All=[];
m8_All=[];
m9_All=[];
m10_All=[];
m11_All=[];
m12_All=[];
m13_All=[];

v1_All=[];
v2_All=[];
v3_All=[];
v4_All=[];
v5_All=[];
v6_All=[];
v7_All=[];
v8_All=[];
v9_All=[];
v10_All=[];
v11_All=[];
v12_All=[];
v13_All=[];

v1_new=[];
v2_new=[];
v3_new=[];
v4_new=[];
v5_new=[];
v6_new=[];
v7_new=[];
v8_new=[];
v9_new=[];
v10_new=[];
v11_new=[];
v12_new=[];
v13_new=[];

VTi1_All=[];
VTi2_All=[];
VTi3_All=[];
VTi4_All=[];
VTi5_All=[];
VTi6_All=[];
VTi7_All=[];
VTi8_All=[];
VTi9_All=[];
VTi10_All=[];
VTi11_All=[];
VTi12_All=[];
VTi13_All=[];

Vi1_All=[];
Vi2_All=[];
Vi3_All=[];
Vi4_All=[];
Vi5_All=[];
Vi6_All=[];
Vi7_All=[];
Vi8_All=[];
Vi9_All=[];
Vi10_All=[];
Vi11_All=[];
Vi12_All=[];
Vi13_All=[];

vti1=[];
vti2=[];
vti3=[];
vti4=[];
vti5=[];
vti6=[];
vti7=[];
vti8=[];
vti9=[];
vti10=[];
vti11=[];
vti12=[];
vti13=[];

vi1=[];vi2=[];vi3=[];vi4=[];vi5=[];vi6=[];vi7=[];vi8=[];vi9=[];vi10=[];vi11=[];vi12=[];vi13=[];

fprintf('Computing variance decomposition... ');
for j = 1:num_Yrs
      %quasi static
  [Sob_g1,Tot_Sob_g1,m1,v1] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g1(j,:),domain,'legendre'); 
  [Sob_g2,Tot_Sob_g2,m2,v2] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g2(j,:),domain,'legendre'); 
  [Sob_g3,Tot_Sob_g3,m3,v3] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g3(j,:),domain,'legendre'); 
  [Sob_g4,Tot_Sob_g4,m4,v4] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g4(j,:),domain,'legendre'); 
  [Sob_g5,Tot_Sob_g5,m5,v5] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g5(j,:),domain,'legendre'); 
  [Sob_g6,Tot_Sob_g6,m6,v6] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g6(j,:),domain,'legendre'); 
  [Sob_g7,Tot_Sob_g7,m7,v7] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g7(j,:),domain,'legendre'); 
  [Sob_g8,Tot_Sob_g8,m8,v8] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g8(j,:),domain,'legendre'); 
  [Sob_g9,Tot_Sob_g9,m9,v9] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g9(j,:),domain,'legendre'); 
  [Sob_g10,Tot_Sob_g10,m10,v10] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g10(j,:),domain,'legendre'); 
  [Sob_g11,Tot_Sob_g11,m11,v11] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g11(j,:),domain,'legendre'); 
  [Sob_g12,Tot_Sob_g12,m12,v12] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g12(j,:),domain,'legendre'); 
  [Sob_g13,Tot_Sob_g13,m13,v13] = compute_sobol_indices_from_sparse_grid(S,Sr,values_g13(j,:),domain,'legendre'); 


%    [Sob_g1,Tot_Sob_g1,m1,v1,vi1,vti1] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi1_All, VTi1_All, v1_new, values_g1(j,:),domain,'legendre'); 
%    [Sob_g2,Tot_Sob_g2,m2,v2,vi2,vti2] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi2_All, VTi2_All, v2_new, values_g2(j,:),domain,'legendre'); 
%    [Sob_g3,Tot_Sob_g3,m3,v3,vi3,vti3] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi3_All, VTi3_All, v3_new, values_g3(j,:),domain,'legendre'); 
%   [Sob_g4,Tot_Sob_g4,m4,v4,vi4,vti4] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi4_All, VTi4_All, v4_new, values_g4(j,:),domain,'legendre'); 
%   [Sob_g5,Tot_Sob_g5,m5,v5,vi5,vti5] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi5_All, VTi5_All, v5_new, values_g5(j,:),domain,'legendre'); 
%   [Sob_g6,Tot_Sob_g6,m6,v6,vi6,vti6] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi6_All, VTi6_All, v6_new, values_g6(j,:),domain,'legendre'); 
%   [Sob_g7,Tot_Sob_g7,m7,v7,vi7,vti7] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi7_All, VTi7_All, v7_new, values_g7(j,:),domain,'legendre'); 
%    [Sob_g8,Tot_Sob_g8,m8,v8,vi8,vti8] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi8_All, VTi8_All, v8_new, values_g8(j,:),domain,'legendre'); 
%   [Sob_g9,Tot_Sob_g9,m9,v9,vi9,vti9] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi9_All, VTi9_All, v9_new, values_g9(j,:),domain,'legendre'); 
%   [Sob_g10,Tot_Sob_g10,m10,v10,vi10,vti10] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi10_All, VTi10_All, v10_new, values_g10(j,:),domain,'legendre'); 
%   [Sob_g11,Tot_Sob_g11,m11,v11,vi11,vti11] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi11_All, VTi11_All, v11_new, values_g11(j,:),domain,'legendre'); 
%   [Sob_g12,Tot_Sob_g12,m12,v12,vi12,vti12] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi12_All, VTi12_All, v12_new, values_g12(j,:),domain,'legendre'); 
%   [Sob_g13,Tot_Sob_g13,m13,v13,vi13,vti13] = compute_sobol_indices_from_sparse_grid_with_history(S,Sr, Vi13_All, VTi13_All, v13_new, values_g13(j,:),domain,'legendre'); 

%     [Sob_g1,Tot_Sob_g1,m1,v1] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g1(j,:),domain,'legendre'); 
%     [Sob_g2,Tot_Sob_g2,m2,v2] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g2(j,:),domain,'legendre'); 
%     [Sob_g3,Tot_Sob_g3,m3,v3] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g3(j,:),domain,'legendre'); 
%     [Sob_g4,Tot_Sob_g4,m4,v4] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g4(j,:),domain,'legendre'); 
%     [Sob_g5,Tot_Sob_g5,m5,v5] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g5(j,:),domain,'legendre'); 
%     [Sob_g6,Tot_Sob_g6,m6,v6] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g6(j,:),domain,'legendre'); 
%     [Sob_g7,Tot_Sob_g7,m7,v7] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g7(j,:),domain,'legendre'); 
%     [Sob_g8,Tot_Sob_g8,m8,v8] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g8(j,:),domain,'legendre'); 
%     [Sob_g9,Tot_Sob_g9,m9,v9] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g9(j,:),domain,'legendre'); 
%     [Sob_g10,Tot_Sob_g10,m10,v10] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g10(j,:),domain,'legendre'); 
%     [Sob_g11,Tot_Sob_g11,m11,v11] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g11(j,:),domain,'legendre'); 
%     [Sob_g12,Tot_Sob_g12,m12,v12] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g12(j,:),domain,'legendre'); 
%     [Sob_g13,Tot_Sob_g13,m13,v13] = compute_variance_decomposition_from_sparse_grid(S,Sr,values_g13(j,:),domain,'legendre'); 

    Sob_g1_All=[Sob_g1_All Sob_g1];
    Sob_g2_All=[Sob_g2_All Sob_g2];
    Sob_g3_All=[Sob_g3_All Sob_g3];
    Sob_g4_All=[Sob_g4_All Sob_g4];
    Sob_g5_All=[Sob_g5_All Sob_g5];
    Sob_g6_All=[Sob_g6_All Sob_g6];
    Sob_g7_All=[Sob_g7_All Sob_g7];    
    Sob_g8_All=[Sob_g8_All Sob_g8];
    Sob_g9_All=[Sob_g9_All Sob_g9];
    Sob_g10_All=[Sob_g10_All Sob_g10];
    Sob_g11_All=[Sob_g11_All Sob_g11];
    Sob_g12_All=[Sob_g12_All Sob_g12];    
    Sob_g13_All=[Sob_g13_All Sob_g13];    
        
    Tot_Sob_g1_All=[Tot_Sob_g1_All Tot_Sob_g1];
    Tot_Sob_g2_All=[Tot_Sob_g2_All Tot_Sob_g2];
    Tot_Sob_g3_All=[Tot_Sob_g3_All Tot_Sob_g3];
    Tot_Sob_g4_All=[Tot_Sob_g4_All Tot_Sob_g4];    
    Tot_Sob_g5_All=[Tot_Sob_g5_All Tot_Sob_g5];   
    Tot_Sob_g6_All=[Tot_Sob_g6_All Tot_Sob_g6];   
    Tot_Sob_g7_All=[Tot_Sob_g7_All Tot_Sob_g7];   
    Tot_Sob_g8_All=[Tot_Sob_g8_All Tot_Sob_g8];
    Tot_Sob_g9_All=[Tot_Sob_g9_All Tot_Sob_g9];    
    Tot_Sob_g10_All=[Tot_Sob_g10_All Tot_Sob_g10];   
    Tot_Sob_g11_All=[Tot_Sob_g11_All Tot_Sob_g11];   
    Tot_Sob_g12_All=[Tot_Sob_g12_All Tot_Sob_g12];   
    Tot_Sob_g13_All=[Tot_Sob_g13_All Tot_Sob_g13];   

    m1_All=[m1_All m1];
    m2_All=[m2_All m2];
    m3_All=[m3_All m3];
    m4_All=[m4_All m4];
    m5_All=[m5_All m5];
    m6_All=[m6_All m6];
    m7_All=[m7_All m7];
    m8_All=[m8_All m8];
    m9_All=[m9_All m9];
    m10_All=[m10_All m10];
    m11_All=[m11_All m11];
    m12_All=[m12_All m12];
    m13_All=[m13_All m13];
    
    v1_new=v1;
    v2_new=v2;
    v3_new=v3;
    v4_new=v4;
    v5_new=v5;
    v6_new=v6;
    v7_new=v7;
    v8_new=v8;
    v9_new=v9;
    v10_new=v10;
    v11_new=v11;
    v12_new=v12;
    v13_new=v13;

    v1_All=[v1_All v1];
    v2_All=[v2_All v2];
    v3_All=[v3_All v3];
    v4_All=[v4_All v4];
    v5_All=[v5_All v5];
    v6_All=[v6_All v6];
    v7_All=[v7_All v7];
    v8_All=[v8_All v8];
    v9_All=[v9_All v9];
    v10_All=[v10_All v10];
    v11_All=[v11_All v11];
    v12_All=[v12_All v12];
    v13_All=[v13_All v13];

    VTi1_All=vti1;%[VTi1_All vti1];
    VTi2_All=vti2;%[VTi2_All vti2];
    VTi3_All=vti3;%[VTi3_All vti3];
    VTi4_All=vti4;%[VTi3_All vti3];
    VTi5_All=vti5;%[VTi3_All vti3];
    VTi6_All=vti6;%[VTi3_All vti3];
    VTi7_All=vti7;%[VTi3_All vti3];
    VTi8_All=vti8;%[VTi8_All vti8];
    VTi9_All=vti9;%[VTi8_All vti8];
    VTi10_All=vti10;%[VTi8_All vti8];
    VTi11_All=vti11;%[VTi8_All vti8];
    VTi12_All=vti12;%[VTi8_All vti8];
    VTi13_All=vti13;%[VTi8_All vti8];

    Vi1_All=vi1;
    Vi2_All=vi2;%[Vi2_All vi2];
    Vi3_All=vi3;%[Vi3_All vi3];
    Vi4_All=vi4;%[Vi3_All vi3];
    Vi5_All=vi5;%[Vi3_All vi3];
    Vi6_All=vi6;%[Vi3_All vi3];
    Vi7_All=vi7;%[Vi3_All vi3];
    Vi8_All=vi8;%[Vi8_All vi8];
    Vi9_All=vi9;%[Vi8_All vi8];
    Vi10_All=vi10;%[Vi8_All vi8];
    Vi11_All=vi11;%[Vi8_All vi8];
    Vi12_All=vi12;%[Vi8_All vi8];   
    Vi13_All=vi13;%[Vi8_All vi8];   
    
end

save('Sob_g1_All','Sob_g1_All')
save('Sob_g2_All','Sob_g2_All')
save('Sob_g3_All','Sob_g3_All')
save('Sob_g4_All','Sob_g4_All')
save('Sob_g5_All','Sob_g5_All')
save('Sob_g6_All','Sob_g6_All')
save('Sob_g7_All','Sob_g7_All')
save('Sob_g8_All','Sob_g8_All')
save('Sob_g9_All','Sob_g9_All')
save('Sob_g10_All','Sob_g10_All')
save('Sob_g11_All','Sob_g11_All')
save('Sob_g12_All','Sob_g12_All')
save('Sob_g13_All','Sob_g13_All')

save('Tot_Sob_g1_All','Tot_Sob_g1_All')
save('Tot_Sob_g2_All','Tot_Sob_g2_All')
save('Tot_Sob_g3_All','Tot_Sob_g3_All')
save('Tot_Sob_g4_All','Tot_Sob_g4_All')
save('Tot_Sob_g5_All','Tot_Sob_g5_All')
save('Tot_Sob_g6_All','Tot_Sob_g6_All')
save('Tot_Sob_g7_All','Tot_Sob_g7_All')
save('Tot_Sob_g8_All','Tot_Sob_g8_All')
save('Tot_Sob_g9_All','Tot_Sob_g9_All')
save('Tot_Sob_g10_All','Tot_Sob_g10_All')
save('Tot_Sob_g11_All','Tot_Sob_g11_All')
save('Tot_Sob_g13_All','Tot_Sob_g13_All')



fprintf('done.\n');