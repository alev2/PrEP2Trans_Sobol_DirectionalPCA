close all
%clear

% initialize paramaters that we do not vary for HIV problem, as well as
% the sparse grid library.
sobolExamples_Setup

rescaleIntervals=0;


fprintf('setup done\n');

% Sobol indices computation

N = 6; % number of parameters...



%%PrEP_Params (multiplier, default is 1 for both):  tt_PrEPInitRateOrPctOnPrEP_Set5
mixing_Range=[0 1];
hetm_PrepRange=[0 1];
hetf_PrepRange=[0 1];
msm_PrepRange=[0 1];
pwid_PrepRange=[0 1];

het_AdhRange=[0 1];
msm_AdhRange=[0 1];
pwid_AdhRange=[0 1];


het_TestRange=[0 1];
msm_TestRange=[0 1];
pwid_TestRange=[0 1];


outputRoot='./PrepToTransGroups/testResult_PrepAndInterventions_4/simulation_';
outputRootPrEP='./PrepToTransGroups/testResult_PrepAndInterventions_4/simulationNumPrEP_';
%outputRoot='./Results_NoVLS_Dec27/Verification_HopeSim/simulation_';
%outputRoot='./datasetInput/simulation_';
%outputRoot='./dataset20Years/simulation_';
Alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
fprintf(strcat('Num points: ',num2str(nb_pts),'.\n'));



min_MixMat_Base=...
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

max_MixMat_Base=...
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

mix_Delta_Base=max_MixMat_Base-min_MixMat_Base;


min_MixMat=min_MixMat_Base+.499*mix_Delta_Base;
max_MixMat=min_MixMat_Base+.501*mix_Delta_Base;
mix_Delta=max_MixMat-min_MixMat;

min_PrepHETM=[
    0.00246799716000218;
    0.0034783715630692;
    0.00299225642682675;
];
max_PrepHETM=8*[
    0.0051;
    0.00715;
    0.00615;
];
prepHETM_Delta=max_PrepHETM-min_PrepHETM;



min_PrepHETF=[
    0.00455228184738825;
    0.00351426269634512;
    0.00237119734569618;
];
max_PrepHETF=8*[
    0.0094;
    0.0072;
    0.00488;
];
prepHETF_Delta=max_PrepHETF-min_PrepHETF;



min_PrepMSM=[
    0.647381303980681;
    0.271113642540548;
    0.319019747544778;
];
max_PrepMSM=[
    3.62;
    0.75;
    0.97;
];
prepMSM_Delta=max_PrepMSM-min_PrepMSM;





min_PrepPWID=[
    0.00513814485068991;
    0.00713683105469181;
    0.0137714047406796;
];
max_PrepPWID=8*[
    0.0107;
    0.0149;
    0.029;
];
prepPWID_Delta=max_PrepPWID-min_PrepPWID;




max_AdhHET=1.0;
min_AdhHET=0.5;
adhHET_Delta=max_AdhHET-min_AdhHET;

max_AdhMSM=2.724113407;
min_AdhMSM=2.724113407*.5;
adhMSM_Delta=max_AdhMSM-min_AdhMSM;

max_AdhPWID=2.268609035;
min_AdhPWID=2.268609035*.5;
adhPWID_Delta=max_AdhPWID-min_AdhPWID;


max_TestHET=2.0;
min_TestHET=1.0;
testHET_Delta=max_TestHET-min_TestHET;

max_TestMSM=4.030993484*2;
min_TestMSM=4.030993484;
testMSM_Delta=max_TestMSM-min_TestMSM;

max_TestPWID=4.635738002*2;
min_TestPWID=4.635738002;
testPWID_Delta=max_TestPWID-min_TestPWID;



%outcomePP=HIVEpiModel_Sobol('HOPE Model V10_05_LB20230224_2_20231005_Fresh.xlsm', sg_pts);

load('ExcelValues_AllParameters.mat');
load('ExcelValues_Populations.mat');
Parameters_Cur=ExcelValues_AllParameters;
Parameters_Cur(indsDiscRate)=0.;



kFactor=1;
%evalDomain=[0 0 1 .864 .79 .83 .15 .16 1 1 .246 .252 .531 .473 .2502 .137];
evalDomain=ones(1,10)*.5;


%load('evalDomain_22Jan.mat');
%load('principalDirections_LoMix2.mat');
load('principalDirectionsIncreasedRange.mat');
%evalDomain=.5*ones(10,1);
%evalDomain=[0. 0. 0. 0. 0.];


%load('principalDirections_Rescaled_22Jan.mat');

baselineDirection=[0;0;0;0;1;1;1;0;0;0];

principalDirections=[zeros(10,1) principalDirections];
interventionDirections=  kFactor*principalDirections+evalDomain';
interventionDirections = [baselineDirection interventionDirections];
%interventionDirections= [interventionDirections [.2;.2; 1.2; 1.1; .9;.9;.15;.15; 1.3;1.3;.2; .2; .6; .6; .15;.15]];
%interventionDirections=max(interventionDirections,0);

% if rescaleIntervals==1
%     load('transformationIntervals.mat');
%     for j=1:size(interventionDirections,2)
%         interventionDirections(:,j)=valInNewInterval(interventionDirections(:,j),zeros(5,1),ones(5,1),transformationIntervals(:,1),transformationIntervals(:,2));
%     end
% 
% end

for i=1:size(interventionDirections,2)
    tic


    Parameters_Cur=ExcelValues_AllParameters;
    mixing_Matrix=min_MixMat+.5*mix_Delta;    
    Parameters_Cur(mixingInds)=mixing_Matrix;



    Parameters_Cur(indsDiscRate)=0.;
    
    fprintf('Changing HETM PrEP... ');
    prep_HETM_Factor=interventionDirections(1,i);    
    prep_HETM=min_PrepHETM+prep_HETM_Factor*prepHETM_Delta;    
    Parameters_Cur(indsPrepHETM2023)=prep_HETM;
    Parameters_Cur(indsPrepHETM2024)=prep_HETM;

    fprintf('Changing HETF PrEP... ');
    prep_HETF_Factor=interventionDirections(2,i);    
    prep_HETF=min_PrepHETF+prep_HETF_Factor*prepHETF_Delta;    
    Parameters_Cur(indsPrepHETF2023)=prep_HETF;
    Parameters_Cur(indsPrepHETF2024)=prep_HETF;

    fprintf('Changing MSM PrEP... ');
    prep_MSM_Factor=interventionDirections(3,i);    
    prep_MSM=min_PrepMSM+prep_MSM_Factor*prepMSM_Delta;    
    Parameters_Cur(indsPrepMSM2023)=prep_MSM;
    Parameters_Cur(indsPrepMSM2024)=prep_MSM;

    fprintf('Changing PWID PrEP... ');
    prep_PWID_Factor=interventionDirections(4,i);    
    prep_PWID=min_PrepPWID+prep_PWID_Factor*prepPWID_Delta;    
    Parameters_Cur(indsPrepPWID2023)=prep_PWID;
    Parameters_Cur(indsPrepPWID2024)=prep_PWID;


    fprintf('Changing HET Adherence... ');
    adh_HET_Factor=interventionDirections(5,i);    
    adh_HET=min_AdhHET+adh_HET_Factor*adhHET_Delta;    
    Parameters_Cur(ind5to4HET)=adh_HET;

    fprintf('Changing MSM Adherence... ');
    adh_MSM_Factor=interventionDirections(6,i);    
    adh_MSM=min_AdhMSM+adh_MSM_Factor*adhMSM_Delta;   
    Parameters_Cur(ind5to4MSM)=adh_MSM;

    fprintf('Changing PWID Adherence... ');
    adh_PWID_Factor=interventionDirections(7,i);    
    adh_PWID=min_AdhPWID+adh_PWID_Factor*adhPWID_Delta;
    Parameters_Cur(ind5to4PWID)=adh_PWID;



    fprintf('Changing HET Testing... ');
    test_HET_Factor=interventionDirections(8,i);    
    test_HET=min_TestHET+test_HET_Factor*testHET_Delta;    
    Parameters_Cur(ind1to2HET)=test_HET;

    fprintf('Changing MSM Testing... ');
    test_MSM_Factor=interventionDirections(9,i);    
    test_MSM=min_TestMSM+test_MSM_Factor*testMSM_Delta;   
    Parameters_Cur(ind1to2MSM)=test_MSM;

    fprintf('Changing PWID Testing... ');
    test_PWID_Factor=interventionDirections(10,i);    
    test_PWID=min_AdhPWID+test_PWID_Factor*testPWID_Delta;
    Parameters_Cur(ind1to2PWID)=test_PWID;
        

    fprintf('Running HOPE... ');
    [outTest, Params]=HIVEpiModel_Sobol(Parameters_Cur, ExcelValues_Populations);
    fprintf('done.\n ');
    

    ann_NewInfections_MSM=sum(outTest.ann_NewInfectionsMSM,2);
    ann_NewInfections_PWID=sum(outTest.ann_NewInfectionsIDU,2);
    ann_NewInfections_HETM=sum(outTest.ann_NewInfectionsHETM,2);
    ann_NewInfections_HETF=sum(outTest.ann_NewInfectionsHETF,2);
    ann_NewInfections_Total=sum(outTest.ann_TotalNewInfections,2);

    NumberOnPrEPMSM=sum(outTest.ann_NumberOnPrEP_HighRiskMSM,2);
    NumberOnPrEPPWID=sum(outTest.ann_NumberOnPrEP_HighRiskIDUs,2);
    NumberOnPrEPHETM=sum(outTest.ann_NumberOnPrEP_HighRiskHETs_M,2);
    NumberOnPrEPHETF=sum(outTest.ann_NumberOnPrEP_HighRiskHETs_F,2);
    
    spending=outTest.ann_ARTCarePrEPTransitionAndSEPCost_Disc;

%    outputMatrix=[ann_NewInfections_MSM ann_NewInfections_HETF ann_NewInfections_HETM ann_NewInfections_PWID ann_NewInfections_Total];       
    outputMatrix=[ann_NewInfections_MSM ann_NewInfections_HETF ann_NewInfections_HETM ann_NewInfections_PWID ann_NewInfections_Total spending/1e9];       
    outputMatrixNumPrEP= [NumberOnPrEPMSM NumberOnPrEPHETF NumberOnPrEPHETM NumberOnPrEPPWID];
       
    incidenceOut=array2table([outputMatrix(yearInds,:)]);% annualDeaths(1:end,:) pctAware(1:end,:)]);


    incidenceOutNumPrEP=array2table([outputMatrixNumPrEP(yearInds,:)]);% annualDeaths(1:end,:) pctAware(1:end,:)]);
    
    fprintf('Writing table... ');    
    writetable(incidenceOut,strcat(outputRoot,num2str(i),'.xls'),'Range',...
       strcat(Alphabet(1),num2str(1),':',Alphabet(N),num2str(size(incidenceOut,1)+1)));%,'WriteVariableNames',false);
    fprintf('done.\n ');
    toc

    fprintf('Writing PrEP table... ');    
    writetable(incidenceOutNumPrEP,strcat(outputRootPrEP,num2str(i),'.xls'),'Range',...
       strcat(Alphabet(1),num2str(1),':',Alphabet(N),num2str(size(incidenceOut,1)+1)));%,'WriteVariableNames',false);
    fprintf('done.\n ');
    toc


end



