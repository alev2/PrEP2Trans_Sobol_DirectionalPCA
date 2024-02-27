function [CalibTargets] = Calib_initializeTargets(ExcelFileName)
%% Purpose: A function which defines the acceptable bounds for the target 
% calibration outcomes.

% Called from: HIVEpiModel



%% 1. Import target ranges from Excel

    %Initialize variables
    importIndex = 1;
    LB = 1;
    UB = 2;    
    
    % Path
    [folder, name, ext] = fileparts(which(mfilename));
    pathName = [folder  '\'  ExcelFileName];

    % Import Values
    % The parameters will then be organized into the specific fields
    
    ExcelValues_CalibTargetRanges = xlsread(pathName, 'Import_Calibration','matlab_CalibTargets');

    % Order
    % 2006 outcomes
        % HIV Prevalence
        % AIDS Deaths
        % PLWH Aware Deaths
        % TT Distribution
        % New Infections
    % 2010 outcomes
        % New Infections
        % HIV Prevalence
        % AIDS Deaths
        % PLWH Aware Deaths
        % TT Distribution
    
    % Note: to allow the calibration functionality to run just to 2006 (not
    % continuing onto 2010), the 2006 outcomes must be listed first
    % and all together.
        % This is due to code in Calib_getDecision.m that sequententially
        % compares each model calibration outcome to the target ranges. If a 2010 target 
        % is listed earlier than one of the 2006 targets (and only 2006 outcomes
        % were collected), the code wouldn't recognize the 2010 outcome was
        % missing and it would compare with the wrong targets.
        
   % Also note: these outcomes must be in the same order as the outcomes in
   % CalibOutputStruct (defined in Calib_collectResults.m).

%% 2. Set TargetValues for Target Outcomes

    % Target Values associated with each target outcome
% New targets start here
% Deleted 24 2010 targets. Bates 6/11/19
% Allaire: Changed calib_HIVPrevalence2015 to calib_HIVPrevalence2016.
% Bates: Updated targets to new years, added 13 targets.
% Clinkscales: Updated targets to 2019. 6/29/2021
% Bates: Added in-care targets 11/4/2022
%Bates: Updated additional targets to 2019. 11/22/22
    CalibTargets.calib_HIVPrevalence2019.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_HIVPrevalence2019_B.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_HIVPrevalence2019_H.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_HIVPrevalence2019_O.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
      
    CalibTargets.calib_HIVPrevalence2019_HETM.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;    
      
    CalibTargets.calib_HIVPrevalence2019_HETF.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
      
    CalibTargets.calib_HIVPrevalence2019_MSM.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
                 
    CalibTargets.calib_HIVPrevalence2019_IDUM.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
      
    CalibTargets.calib_HIVPrevalence2019_IDUF.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;            
            
    CalibTargets.calib_HIVPrevalence2016_13_24.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_HIVPrevalence2016_25_34.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1; 
            
    CalibTargets.calib_HIVPrevalence2016_35_44.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_HIVPrevalence2016_45_54.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_HIVPrevalence2016_55_64.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_HIVPrevalence2016_65.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1; 
            
    CalibTargets.calib_HIVPrevalence2016_B_13_24.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_HIVPrevalence2016_H_13_24.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_HIVPrevalence2016_O_13_24.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;    
            
    CalibTargets.calib_HIVPrevalence2016_B_25_34.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_HIVPrevalence2016_H_25_34.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_HIVPrevalence2016_O_25_34.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;            

    CalibTargets.calib_TTdist2019_Diagnosed_B.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2019_Diagnosed_H.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
    
    CalibTargets.calib_TTdist2019_Diagnosed_O.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2016_Diagnosed_13_24.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
    
    CalibTargets.calib_TTdist2016_Diagnosed_25_34.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2016_Diagnosed_35_44.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_TTdist2016_Diagnosed_45_54.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_TTdist2016_Diagnosed_55_64.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1; 
            
    CalibTargets.calib_TTdist2016_Diagnosed_65.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_TTdist2016_Diagnosed_MSM_B_13_24.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2016_Diagnosed_MSM_H_13_24.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2016_Diagnosed_MSM_O_13_24.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2016_Diagnosed_MSM_B_25_34.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2016_Diagnosed_MSM_H_25_34.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2016_Diagnosed_MSM_O_25_34.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2019_Diagnosed_HET.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2019_Diagnosed_MSM.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2019_Diagnosed_PWID.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2019_Diagnosed_Total.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1; 

    %Added in-care targets. Laurel Bates 11/4/22

    CalibTargets.calib_TTdist2019_InCareAmongDiag_B.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2019_InCareAmongDiag_H.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
    
    CalibTargets.calib_TTdist2019_InCareAmongDiag_O.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;

    CalibTargets.calib_TTdist2019_InCareAmongDiag_Total.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        

    CalibTargets.calib_TTdist2019_VLSamongdiag_B.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2019_VLSamongdiag_H.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
    
    CalibTargets.calib_TTdist2019_VLSamongdiag_O.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
    
    CalibTargets.calib_TTdist2018_VLSamongdiag_HET.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
    
    CalibTargets.calib_TTdist2018_VLSamongdiag_MSM.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
    
    CalibTargets.calib_TTdist2018_VLSamongdiag_PWID.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;            
            
    CalibTargets.calib_TTdist2015_VLSamongdiag_13_24.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
    
    CalibTargets.calib_TTdist2015_VLSamongdiag_25_34.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2015_VLSamongdiag_35_44.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2015_VLSamongdiag_45_54.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_TTdist2015_VLSamongdiag_55_64.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_TTdist2015_VLSamongdiag_65.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_TTdist2016_VLSamongdiag_13_24_MSM_B.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2016_VLSamongdiag_13_24_MSM_H.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2016_VLSamongdiag_13_24_MSM_O.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2016_VLSamongdiag_25_34_MSM_B.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2016_VLSamongdiag_25_34_MSM_H.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_TTdist2016_VLSamongdiag_25_34_MSM_O.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_TTdist2019_VLSamongdiag_Total.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_AwarePLWHDeaths2019.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_AwarePLWHDeaths2016_45_54.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_AwarePLWHDeaths2016_55_64.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_AwarePLWHDeaths2016_65.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;

    CalibTargets.calib_AwarePLWHAIDSDeaths2019.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;

    CalibTargets.calib_UnawarePLWHAIDSDeaths2019.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
                                                
    CalibTargets.calib_NewDiagnoses2016_HETM.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NewDiagnoses2016_HETF.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1; 
            
    CalibTargets.calib_NewDiagnoses2016_MSM.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
                        
    CalibTargets.calib_NewDiagnoses2016_IDUM.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_NewDiagnoses2016_IDUF.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;    
            
    CalibTargets.calib_NewDiagnoses2016_Total.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;  
            
    %Clinkscales 04/21/2022: Added NewDiagnoses2019 calibration targets import below        
            
    CalibTargets.calib_NewDiagnoses2019_PctB1C1.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1; 
            
    CalibTargets.calib_NewDiagnoses2019_Total.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;         
            
    CalibTargets.calib_NewInfections2019_B.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_NewInfections2019_H.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NewInfections2019_O.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;   
            
    CalibTargets.calib_NewInfections2019_HETM.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NewInfections2019_HETF.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NewInfections2019_MSM.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NewInfections2019_IDUM.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NewInfections2019_IDUF.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;            
            
    CalibTargets.calib_NewInfections_2016_MSM_13_24.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NewInfections_2016_MSM_13_24_B.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1; 
            
    CalibTargets.calib_NewInfections_2016_MSM_13_24_H.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1; 
            
    CalibTargets.calib_NewInfections_2016_MSM_13_24_O.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1; 
            
    CalibTargets.calib_NewInfections_2016_MSM_25_34.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1; 
            
    CalibTargets.calib_NewInfections_2016_MSM_25_34_B.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1; 
            
     CalibTargets.calib_NewInfections_2016_MSM_25_34_H.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;  
            
    CalibTargets.calib_NewInfections_2016_MSM_25_34_O.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;              
            
    CalibTargets.calib_NewInfections2019_Total.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NewInfections2016_45_54.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_NewInfections2016_55_64.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;        
            
    CalibTargets.calib_NewInfections2016_65.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;

    %Revised number on PrEP targets 02/10/23
    CalibTargets.calib_NumOnPrEP2019_Male.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NumOnPrEP2019_Female.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NumOnPrEP2019_Black.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NumOnPrEP2019_Hispanic.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NumOnPrEP2019_Other.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NumOnPrEP2019_Total.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NumOnPrEP2021_Male.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NumOnPrEP2021_Female.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NumOnPrEP2021_Black.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NumOnPrEP2021_Hispanic.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NumOnPrEP2021_Other.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_NumOnPrEP2021_Total.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_OverallPrev2019v2010.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;           

    CalibTargets.calib_HETPrev2019v2010_LR.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
        
    CalibTargets.calib_HETPrev2019v2010_HR.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_populationMaleBlk2019.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_populationMaleHisp2019.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_populationMaleOth2019.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_populationFemaleBlk2019.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_populationFemaleHisp2019.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;            
            
    CalibTargets.calib_populationFemaleOth2019.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;

    CalibTargets.calib_population_1334_2019.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_population_3564_2019.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
    CalibTargets.calib_population_65_2019.targetValue = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
            
% Collect names of fields (for collecting other values)
TargetFields = fieldnames(CalibTargets);
% Count number of fields
nTargetOutcomes = length(TargetFields);

%% 3. Apply Weights for Target Outcomes

        for i = 1:nTargetOutcomes
            CalibTargets.(TargetFields{i}).weight = ExcelValues_CalibTargetRanges(importIndex);
            importIndex = importIndex + 1;
        end
            
%% 4. Apply Ranges for Target Outcomes
    
    % Lower bounds
    for i = 1:nTargetOutcomes
        CalibTargets.(TargetFields{i}).range(LB) = ExcelValues_CalibTargetRanges(importIndex);
        importIndex = importIndex + 1;
    end
           
    % Upper bounds
    for i = 1:nTargetOutcomes
        CalibTargets.(TargetFields{i}).range(UB) = ExcelValues_CalibTargetRanges(importIndex);
        importIndex = importIndex + 1;
    end
    
    %    Test Vector to see if bounds match
%         counter = 0;
%         for i = 1:nTargetOutcomes
%             counter = counter + 1;
%             bound_test(counter,1) = CalibTargets.(TargetFields{i}).targetValue;
%             bound_test(counter,2) = CalibTargets.(TargetFields{i}).range(LB);
%             bound_test(counter,3) = CalibTargets.(TargetFields{i}).range(UB);
%             if bound_test(counter,1) < bound_test(counter,2) | bound_test(counter,1) > bound_test(counter,3)
%                 bound_test(counter,4) = 1
%             end
%     
%         end

end