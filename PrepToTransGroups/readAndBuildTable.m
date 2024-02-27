close all
%%clear

% initialize paramaters that we do not vary for HIV problem, as well as
% the sparse grid library.
sobolExamples_Setup
outputRoot='./Results_MixingAndPrEP_ExpandedOutputSpace_LoMixing/simulation_';
%outputRoot='./datasetInput_LargeScale2_Dec4/simulation_';

normVar=1;
fprintf('setup done\n');

% Sobol indices computation


N = 5; % number of parameters...



%%PrEP_Params (multiplier, default is 1 for both):  tt_PrEPInitRateOrPctOnPrEP_Set5
mixing_Range=[0 1];
hetm_PrepRange=[0 1];
hetf_PrepRange=[0 1];
msm_PrepRange=[0 1];
pwid_PrepRange=[0 1];


num_Yrs=length(yearInds);
consideredYears=yearRange(yearInds);

start_YrSum=2023;
end_YrSum=2042;

startIndex=find(consideredYears==start_YrSum);
endIndex=find(consideredYears==end_YrSum);

if(normVar==1)
% CHANGE HERE! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int_p1 = [0 1]; % interval of variation for parameter 1 
    int_p2 = [0 1]; % interval of variation for parameter 2
    int_p3 = [0 1]; % interval of variation for parameter 3
    int_p4 = [0 1]; % interval of variation for parameter 4
    int_p5 = [0 1]; % interval of variation for parameter 5
%     int_p6 = [0 1]; % interval of variation for parameter 6
%     int_p7 = [0 1]; % interval of variation for parameter 7
%     int_p8 = [0 1]; % interval of variation for parameter 8
%     int_p9 = [0 1]; % interval of variation for parameter 9
%     int_p10 = [0 1]; % interval of variation for parameter 10
%     int_p11 = [0 1]; % interval of variation for parameter 11
%     int_p12 = [0 1]; % interval of variation for parameter 12
%     int_p13 = [0 1]; % interval of variation for parameter 13
%     int_p14 = [0 1]; % interval of variation for parameter 14
%     int_p15 = [0 1]; % interval of variation for parameter 15
%     int_p16 = [0 1]; % interval of variation for parameter 16    
else    
    % CHANGE HERE! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int_p1 = mixing_Range; % interval of variation for parameter 1 
    int_p2 = hetm_PrepRange; % interval of variation for parameter 2
    int_p3 = hetf_PrepRange; % interval of variation for parameter 3
    int_p4 = msm_PrepRange; % interval of variation for parameter 4
    int_p5 = pwid_PrepRange; % interval of variation for parameter 5
%     int_p6 = move2to3_atDiag_Hisp; % interval of variation for parameter 6
%     int_p7 = move2to3_afterDiag_Black; % interval of variation for parameter 7
%     int_p8 = move2to3_afterDiag_Hisp; % interval of variation for parameter 8
%     int_p9 = move3to4_Black; % interval of variation for parameter 9
%     int_p10 = move3to4_Hisp; % interval of variation for parameter 10
%     int_p11 = move4to3_Black; % interval of variation for parameter 11
%     int_p12 = move4to3_Hisp; % interval of variation for parameter 12
%     int_p13 = move4to5_Black; % interval of variation for parameter 13
%     int_p14 = move4to5_Hisp; % interval of variation for parameter 14
%     int_p15 = move5to4_Black; % interval of variation for parameter 15
%     int_p16 = move5to4_Hisp; % interval of variation for parameter 16
end

domain = [int_p1',int_p2',int_p3', int_p4', int_p5'];%, int_p6' int_p7' int_p8',...
%          int_p9' int_p10' int_p11' int_p12' int_p13' int_p14' int_p15' int_p16' ]; % define the N-dim. domain - we need this in the function that computes the Sobol indices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% type of knots to be used for the construction of the sparse grid
% (we use for symmetric Leja points for the moment, this should be ok!)
knots_1 = @(n) knots_leja(n,int_p1(1),int_p1(2),'sym_line');
knots_2 = @(n) knots_leja(n,int_p2(1),int_p2(2),'sym_line');
knots_3 = @(n) knots_leja(n,int_p3(1),int_p3(2),'sym_line');
knots_4 = @(n) knots_leja(n,int_p4(1),int_p4(2),'sym_line');
knots_5 = @(n) knots_leja(n,int_p5(1),int_p5(2),'sym_line');
% knots_6 = @(n) knots_leja(n,int_p6(1),int_p6(2),'sym_line');
% knots_7 = @(n) knots_leja(n,int_p7(1),int_p7(2),'sym_line');
% knots_8 = @(n) knots_leja(n,int_p8(1),int_p8(2),'sym_line');
% knots_9 = @(n) knots_leja(n,int_p9(1),int_p9(2),'sym_line');
% knots_10 = @(n) knots_leja(n,int_p10(1),int_p10(2),'sym_line');
% knots_11 = @(n) knots_leja(n,int_p11(1),int_p11(2),'sym_line');
% knots_12 = @(n) knots_leja(n,int_p12(1),int_p12(2),'sym_line');
% knots_13 = @(n) knots_leja(n,int_p13(1),int_p13(2),'sym_line');
% knots_14 = @(n) knots_leja(n,int_p14(1),int_p14(2),'sym_line');
% knots_15 = @(n) knots_leja(n,int_p15(1),int_p15(2),'sym_line');
% knots_16 = @(n) knots_leja(n,int_p16(1),int_p16(2),'sym_line');


% combine the knots for the three parameters in one - we need this as input to the function that creates the sparse grid
knots = {knots_1,knots_2,knots_3, knots_4, knots_5};%, knots_6, knots_7, knots_8...
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

% % CHANGE HERE! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%num_Yrs = 5; 
% read from your files and arrange the results in a convenient way
% we store the results for each group in one table 
values_g1 = 0*ones(num_Yrs,nb_pts); % incidence msm
values_g2 = 0*ones(num_Yrs,nb_pts); % incidence hetf
values_g3 = 0*ones(num_Yrs,nb_pts); % incidence hetm
values_g4 = 0*ones(num_Yrs,nb_pts); % incidence pwid
values_g5 = 0*ones(num_Yrs,nb_pts); % incidence total
%values_g6 = 0*ones(num_Yrs,nb_pts); % spending
%values_g7 = 0*ones(num_Yrs,nb_pts); % spending


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
    %values_g6(:,k) = tableIn(:,6);    
    %values_g7(:,k) = tableIn(:,7);    

end

fprintf('done.\n');

modelOutputs=[sum(values_g1(startIndex:endIndex,:)); sum(values_g2(startIndex:endIndex,:)); sum(values_g3(startIndex:endIndex,:));...
              sum(values_g4(startIndex:endIndex,:)); sum(values_g5(startIndex:endIndex,:))];%values_g6(end,:);values_g7(end,:)];



