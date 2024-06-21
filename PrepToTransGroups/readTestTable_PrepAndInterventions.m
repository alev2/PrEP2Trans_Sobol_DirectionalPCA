sobolExamples_Setup
outputRoot='./testResult_PrepAndInterventions_3/simulation_';

num_Yrs=length(yearInds);
nb_pts=12;


num_Yrs=length(yearInds);
consideredYears=yearRange(yearInds);

start_YrSum=2023;
end_YrSum=2040;

startIndex=find(consideredYears==start_YrSum);
endIndex=find(consideredYears==end_YrSum);

values_g1_Test = 0*ones(num_Yrs,nb_pts); % incidence black
values_g2_Test = 0*ones(num_Yrs,nb_pts); % incidence hisp
values_g3_Test = 0*ones(num_Yrs,nb_pts); % irr black
values_g4_Test = 0*ones(num_Yrs,nb_pts); % irr hisp
values_g5_Test = 0*ones(num_Yrs,nb_pts); % spending
values_g6_Test = 0*ones(num_Yrs,nb_pts); % spending
%values_g6_Test = 0*ones(num_Yrs,nb_pts); % spending
%values_g7_Test = 0*ones(num_Yrs,nb_pts); % spending


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

    values_g1_Test(:,k) = tableIn(:,1); % if you give another input to xlsread you can select columns in the file 
    values_g2_Test(:,k) = tableIn(:,2);
    values_g3_Test(:,k) = tableIn(:,3);
    values_g4_Test(:,k) = tableIn(:,4);
    values_g5_Test(:,k) = tableIn(:,5);    
    values_g6_Test(:,k) = tableIn(:,6);    
    %values_g7_Test(:,k) = tableIn(:,7);    

end

fprintf('done.\n');
testOutputs=[sum(values_g1_Test(startIndex:endIndex,:)); sum(values_g2_Test(startIndex:endIndex,:)); sum(values_g3_Test(startIndex:endIndex,:));...
             sum(values_g4_Test(startIndex:endIndex,:)); sum(values_g5_Test(startIndex:endIndex,:));sum(values_g6_Test(startIndex:endIndex,:));];
% modelOutputs=[(values_g1_Test(endIndex,:)); (values_g2_Test(endIndex,:)); (values_g3_Test(endIndex,:));...
%               (values_g4_Test(endIndex,:))];%; (values_g5(startIndex:endIndex,:))];%values_g6(end,:);values_g7(end,:)];

%note that you gotta first run reducedOrderSobol for this to work
testOutputsNorm=stdMat\(testOutputs-rowMeans(:,1:nb_pts));

