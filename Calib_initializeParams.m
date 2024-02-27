function [CalibParams, nSamples, calibMethod, calibObjValueType,initValueSource, calibStopCondition, calibTimeToRun, nRandomStartRuns, calibOVStop, calibNumBelowStop, CalibParams_reduced] = Calib_initializeParams(ExcelFileName, Params)
%% Purpose: A function which sets the initial values of the parameters that
% are calibrated

% Called from: Calib_getDecision

%% 1. Import parameters from Excel File

    % File path
    [folder, name, ext] = fileparts(which(mfilename));
    pathName = [folder  '\'  ExcelFileName];

% Import the parameters to arrays that store the values
    % The parameters will then be organized into the specific arrays
    ExcelValues_CalibParams = xlsread(pathName, 'Import_Calibration','matlab_CalibParams');

%% 2. Apply Calibration settings

% Initialize
importIndex = 1;
LB = 1;
UB = 2;
        
% If using LHS, number of samples to generate
nSamples = ExcelValues_CalibParams(importIndex);
importIndex = importIndex + 1;
    
% Method to use: 1=LHS, 2=OptFmincon
calibMethod = ExcelValues_CalibParams(importIndex);
importIndex = importIndex + 1;
if calibMethod > 1 
    nSamples=1; %nSamples only applicable if running LHS
end

% Objective value to use: 1 = OOB Penalty,  2 = Target error
calibObjValueType = ExcelValues_CalibParams(importIndex);
importIndex = importIndex + 1;

% Initial solution as base (1) or random (2)
initValueSource = ExcelValues_CalibParams(importIndex);
importIndex = importIndex + 1;

% Number of random start runs
nRandomStartRuns = ExcelValues_CalibParams(importIndex);
importIndex = importIndex + 1;

% Stopping condition 
calibStopCondition =  ExcelValues_CalibParams(importIndex);
importIndex = importIndex + 1;

% Time to run calibration if stop condition is time (in seconds) 
calibTimeToRun =  ExcelValues_CalibParams(importIndex);
importIndex = importIndex + 1;

% Objective value to stop at if OV is stop condition
calibOVStop =  ExcelValues_CalibParams(importIndex);
importIndex = importIndex + 1;

% Number of runs to collect below objective value threshold before stopping (if stop condition is OV)
calibNumBelowStop =  ExcelValues_CalibParams(importIndex);
importIndex = importIndex + 1;

%% 3. Apply Base (current Run) Values for Input Parameters
 
% Test rates and relative risks

       % Reference cases: Black, HET, CD4 > 500
        CalibParams.tt_testRefCase_1.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_testRefCase_2to4.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;

        % Period 1 relative risks (edited by JC on 11/13/2017)
        CalibParams.tt_relRiskRace_1_H.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskRace_1_O.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskPop_1_MSM.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskPop_1_IDU.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskHIVstage_1_B.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            % Note: CD4>500 (aka disease stage "C") is the ref case so it's
            % not calibrated
        CalibParams.tt_relRiskHIVstage_1_D.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskHIVstage_1_E.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskHIVstage_1_F.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_1_1824.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_1_2534.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_1_3544.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_1_4554.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_1_5564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1; 
        CalibParams.tt_relRiskAge_1_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;

        % Period 2 relative risks (edited by JC on 11/13/2017)
        CalibParams.tt_relRiskRace_2to4_H.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskRace_2to4_O.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskPop_2to4_MSM.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskPop_2to4_IDU.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskHIVstage_2to4_B.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            % Note: CD4>500 (aka disease stage "C") is the ref case so it's
            % not calibrated
        CalibParams.tt_relRiskHIVstage_2to4_D.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskHIVstage_2to4_E.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskHIVstage_2to4_F.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_2to4_1824.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_2to4_2534.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_2to4_3544.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_2to4_4554.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_2to4_5564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1; 
        CalibParams.tt_relRiskAge_2to4_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
    
% Other Continuum of care parameters
    
    % Linkage first (removed on 2018/01/25 by JC)
        %CalibParams.tt_linkageFirst_r_1_B.baseValue = ExcelValues_CalibParams(importIndex);
            %importIndex = importIndex + 1;
        %CalibParams.tt_linkageFirst_r_1_H.baseValue = ExcelValues_CalibParams(importIndex);
            %importIndex = importIndex + 1;
        %CalibParams.tt_linkageFirst_r_1_O.baseValue = ExcelValues_CalibParams(importIndex);
            %importIndex = importIndex + 1;
       % CalibParams.tt_linkageFirst_r_2to4_B.baseValue = ExcelValues_CalibParams(importIndex);
            %importIndex = importIndex + 1;
        %CalibParams.tt_linkageFirst_r_2to4_H.baseValue = ExcelValues_CalibParams(importIndex);
           % importIndex = importIndex + 1;
       % CalibParams.tt_linkageFirst_r_2to4_O.baseValue = ExcelValues_CalibParams(importIndex);
         %   importIndex = importIndex + 1;

    % Linkage after first
        CalibParams.tt_linkage_r_1_B.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_linkage_r_1_H.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_linkage_r_1_O.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_linkage_r_2to4_B.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_linkage_r_2to4_H.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_linkage_r_2to4_O.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
    % Relative risk of LTC After, by disease stage        
        CalibParams.tt_RelRiskLTC_1_E.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_RelRiskLTC_1_F.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_RelRiskLTC_2to4_E.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_RelRiskLTC_2to4_F.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
    % Drop out of care
        CalibParams.tt_dropOutProb_CareToAware_1_B.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_dropOutProb_CareToAware_1_H.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_dropOutProb_CareToAware_1_O.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_dropOutProb_CareToAware_2to4_B.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_dropOutProb_CareToAware_2to4_H.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_dropOutProb_CareToAware_2to4_O.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
    % Drop out of ANV (only ANV not in care) to Aware
        %CalibParams.tt_dropOutProb_ANVToAware_1_B.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.tt_dropOutProb_ANVToAware_1_H.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.tt_dropOutProb_ANVToAware_1_O.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.tt_dropOutProb_ANVToAware_2to4_B.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.tt_dropOutProb_ANVToAware_2to4_H.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.tt_dropOutProb_ANVToAware_2to4_O.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;

    % Drop of of ANV to Care
        CalibParams.tt_dropOutProb_ANVToCare_1_B.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_dropOutProb_ANVToCare_1_H.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_dropOutProb_ANVToCare_1_O.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_dropOutProb_ANVToCare_2to4_B.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_dropOutProb_ANVToCare_2to4_H.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_dropOutProb_ANVToCare_2to4_O.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;

    % Drop out of VLS
        CalibParams.tt_dropOutProb_VLSToANV_1_B.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_dropOutProb_VLSToANV_1_H.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_dropOutProb_VLSToANV_1_O.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_dropOutProb_VLSToANV_2to4_B.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_dropOutProb_VLSToANV_2to4_H.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_dropOutProb_VLSToANV_2to4_O.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
     % RR of dropping out of VLS by transmission group (added by JC 11/13/2017)
        CalibParams.tt_relRiskPop_VLSToANV_1_MSM.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskPop_VLSToANV_1_PWID.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskPop_VLSToANV_2to4_MSM.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskPop_VLSToANV_2to4_PWID.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
     % RR of dropping out of VLS by age group (added by JC 11/13/2017)
        CalibParams.tt_relRiskAge_VLSToANV_1_1824.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_VLSToANV_1_2534.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_VLSToANV_1_3544.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_VLSToANV_1_4554.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_VLSToANV_1_5564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_VLSToANV_1_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_VLSToANV_2to4_1824.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_VLSToANV_2to4_2534.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_VLSToANV_2to4_3544.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_VLSToANV_2to4_4554.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_VLSToANV_2to4_5564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_VLSToANV_2to4_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
    
     % ART initiation (removed on 2018/01/25 by JC)
        %CalibParams.tt_ARTInitiation_d_1_B.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.tt_ARTInitiation_d_1_C.baseValue = ExcelValues_CalibParams(importIndex);
        %        importIndex = importIndex + 1;
        %CalibParams.tt_ARTInitiation_d_1_D.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.tt_ARTInitiation_d_1_E.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.tt_ARTInitiation_d_1_F.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.tt_ARTInitiation_d_2to4_B.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.tt_ARTInitiation_d_2to4_C.baseValue = ExcelValues_CalibParams(importIndex);
        %        importIndex = importIndex + 1;
        %CalibParams.tt_ARTInitiation_d_2to4_D.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.tt_ARTInitiation_d_2to4_E.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.tt_ARTInitiation_d_2to4_F.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
            
        % Relative risk of initiating ART by race/ethnicity
        CalibParams.tt_RelRiskARTInit_r_1_B.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_RelRiskARTInit_r_1_H.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_RelRiskARTInit_r_1_O.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_RelRiskARTInit_r_2to4_B.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_RelRiskARTInit_r_2to4_H.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_RelRiskARTInit_r_2to4_O.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;

        % Probability of transitioning from ANV to VLS
        CalibParams.tt_BecomeVLSfromANV_r_1_B.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_BecomeVLSfromANV_r_1_H.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_BecomeVLSfromANV_r_1_O.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;    
        CalibParams.tt_BecomeVLSfromANV_r_2to4_B.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;    
        CalibParams.tt_BecomeVLSfromANV_r_2to4_H.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;    
        CalibParams.tt_BecomeVLSfromANV_r_2to4_O.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        
        % RR of becoming VLS by transmission group (added by JC 11/13/2017)
        CalibParams.tt_relRiskPop_ANVToVLS_1_MSM.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskPop_ANVToVLS_1_PWID.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskPop_ANVToVLS_2to4_MSM.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskPop_ANVToVLS_2to4_PWID.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
        % RR of becoming VLS by age group (added by JC 11/13/2017)
        CalibParams.tt_relRiskAge_ANVToVLS_1_1824.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_ANVToVLS_1_2534.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_ANVToVLS_1_3544.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_ANVToVLS_1_4554.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_ANVToVLS_1_5564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_ANVToVLS_1_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_ANVToVLS_2to4_1824.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_ANVToVLS_2to4_2534.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_ANVToVLS_2to4_3544.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_ANVToVLS_2to4_4554.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_ANVToVLS_2to4_5564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.tt_relRiskAge_ANVToVLS_2to4_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
    % Distribution of ANV      (removed on 2018/01/25 by JC)  
        
        % Percentage of ANV in care (ART + not on ART)
        %CalibParams.tt_PctANVWhoAreInCareOrART_1.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.tt_PctANVWhoAreInCareOrART_2to4.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;   
            
        % Percentage of in-care ANV who are on ART
        %CalibParams.tt_PctInCareANVWhoAreOnART_1.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.tt_PctInCareANVWhoAreOnART_2to4.baseValue = ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;

% Disease progression

        % Length of time in ANV
        CalibParams.hiv_durStage_ANV_LatentA_1344.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.hiv_durStage_ANV_LatentA_4564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.hiv_durStage_ANV_LatentA_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.hiv_durStage_ANV_LatentB_1344.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;       
        CalibParams.hiv_durStage_ANV_LatentB_4564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.hiv_durStage_ANV_LatentB_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.hiv_durStage_ANV_Late_1344.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.hiv_durStage_ANV_Late_4564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.hiv_durStage_ANV_Late_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;

        % CD4 decline while VLS
       CalibParams.hiv_rateVLSCD4decr_LatentA.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
       CalibParams.hiv_rateVLSCD4decr_LatentB.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
       CalibParams.hiv_rateVLSCD4decr_Late.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;

       % CD4 increase while VLS    
       CalibParams.hiv_rateVLSCD4incr_LatentB_1344.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
       CalibParams.hiv_rateVLSCD4incr_LatentB_4564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
       CalibParams.hiv_rateVLSCD4incr_LatentB_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
       CalibParams.hiv_rateVLSCD4incr_Late_1344.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
       CalibParams.hiv_rateVLSCD4incr_Late_4564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
       CalibParams.hiv_rateVLSCD4incr_Late_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
       CalibParams.hiv_rateVLSCD4incr_AIDS_1344.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
       CalibParams.hiv_rateVLSCD4incr_AIDS_4564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
       CalibParams.hiv_rateVLSCD4incr_AIDS_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
           
% Sexual mixing matrix
        % Transmission group/sex
            % KH updated HETM and HETF mixing to vary by HR and LR on 20Sep2023
        CalibParams.behav_SexualMixing_LRHETM_HETF.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixing_LRHETF_HETM.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixing_LRHETF_IDUM.baseValue =ExcelValues_CalibParams(importIndex);
             importIndex = importIndex + 1;
        CalibParams.behav_SexualMixing_HRHETM_HETF.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixing_HRHETF_HETM.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixing_HRHETF_IDUM.baseValue =ExcelValues_CalibParams(importIndex);
             importIndex = importIndex + 1;
        CalibParams.behav_SexualMixing_MSM_HETF.baseValue = ExcelValues_CalibParams(importIndex);
                importIndex = importIndex + 1;
        CalibParams.behav_SexualMixing_MSM_IDUF.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixing_IDUM_IDUF.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixing_IDUF_IDUM.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
    
    % Risk level
        % HET / IDU
        CalibParams.behav_SexualMixing_HET_LOW_LOW.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixing_HET_HIGH_HIGH.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        % MSM
        CalibParams.behav_SexualMixing_MSM_LOW_LOW.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixing_MSM_HIGH_HIGH.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
     % Race/ethnicity - HET/PWID
        CalibParams.behav_SexualMixingHET_BLK_HISP.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingHET_BLK_OTH.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingHET_HISP_BLK.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingHET_HISP_OTH.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingHET_OTH_BLK.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingHET_OTH_HISP.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
     % Race/ethnicity - MSM (added by JC on 11/08/2017)
        CalibParams.behav_SexualMixingMSM_BLK_HISP.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_BLK_OTH.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_HISP_BLK.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_HISP_OTH.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_OTH_BLK.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_OTH_HISP.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
     % Age (MSM) (edited by JC on 12/08/2017)
        % 13-17
        CalibParams.behav_SexualMixingMSM_1317_1824.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_1317_2534.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_1317_3544.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_1317_4554.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_1317_5564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_1317_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        
        % 18-24
        CalibParams.behav_SexualMixingMSM_1824_1317.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;        
        CalibParams.behav_SexualMixingMSM_1824_2534.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_1824_3544.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_1824_4554.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_1824_5564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_1824_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
        % 25-34
        CalibParams.behav_SexualMixingMSM_2534_1317.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_2534_1824.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;        
        CalibParams.behav_SexualMixingMSM_2534_3544.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_2534_4554.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_2534_5564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_2534_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        
        % 35-44
        CalibParams.behav_SexualMixingMSM_3544_1317.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_3544_1824.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_3544_2534.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;        
        CalibParams.behav_SexualMixingMSM_3544_4554.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_3544_5564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_3544_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
        % 45-54
        CalibParams.behav_SexualMixingMSM_4554_1317.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_4554_1824.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_4554_2534.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_4554_3544.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_4554_5564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_4554_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        
       % 55-64
        CalibParams.behav_SexualMixingMSM_5564_1317.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_5564_1824.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_5564_2534.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_5564_3544.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_5564_4554.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_5564_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
        % 65+
        CalibParams.behav_SexualMixingMSM_65_1317.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_65_1824.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_65_2534.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_65_3544.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_65_4554.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_SexualMixingMSM_65_5564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;

% Base transmission risks
        CalibParams.inf_vaginalInsRisk.baseValue =ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.inf_vaginalRecRisk.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.inf_analInsRisk.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.inf_analRecRisk.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.inf_sharedNeedleRisk.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
           
% Needle reduction
        CalibParams.inf_vlsNeedleReduct.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;   
            
% Percent of injections that are shared
        CalibParams.inf_pctInjectionsShared.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;   
            
% Needle mixing matrix
        CalibParams.behav_NeedleMixing_IDUM_IDUF.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_NeedleMixing_IDUF_IDUM.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;

% Percent of MSM contacts that are insertive, by age group and race (edited
% by JC on 12/05/2017)
    % Black
        CalibParams.behav_HRMSMNumContacts_Black_13.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Black_18.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Black_25.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Black_35.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Black_45.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Black_55.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Black_65.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
     % Hispanic
        CalibParams.behav_HRMSMNumContacts_Hispanic_13.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Hispanic_18.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Hispanic_25.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Hispanic_35.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Hispanic_45.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Hispanic_55.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Hispanic_65.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
     % Other
        CalibParams.behav_HRMSMNumContacts_Other_13.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Other_18.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Other_25.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Other_35.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Other_45.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Other_55.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_HRMSMNumContacts_Other_65.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
% Percent of MSM contacts that are insertive, by age group and race (edited
% by JC on 11/06/2017)
    % Black
        CalibParams.behav_MSMpctContacts_Insertive_Black_13.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Black_18.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Black_25.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Black_35.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Black_45.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Black_55.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Black_65.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
     % Hispanic
        CalibParams.behav_MSMpctContacts_Insertive_Hispanic_13.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Hispanic_18.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Hispanic_25.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Hispanic_35.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Hispanic_45.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Hispanic_55.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Hispanic_65.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
     % Other
        CalibParams.behav_MSMpctContacts_Insertive_Other_13.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Other_18.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Other_25.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Other_35.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Other_45.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Other_55.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctContacts_Insertive_Other_65.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
          
% HET AI Parameters
        %CalibParams.factor_pctContactsA_VandA_M_13.baseValue=ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.factor_pctContactsA_VandA_M_18.baseValue=ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.factor_pctContactsA_VandA_M_25.baseValue=ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.factor_pctContactsA_VandA_M_35.baseValue=ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.factor_pctContactsA_VandA_M_45.baseValue=ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.factor_pctContactsA_VandA_F_13.baseValue=ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.factor_pctContactsA_VandA_F_18.baseValue=ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.factor_pctContactsA_VandA_F_25.baseValue=ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.factor_pctContactsA_VandA_F_35.baseValue=ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;
        %CalibParams.factor_pctContactsA_VandA_F_45.baseValue=ExcelValues_CalibParams(importIndex);
        %    importIndex = importIndex + 1;

        % Percent of M-F partnerships that include AI
        CalibParams.pctPartnershipsVandA_1_13.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.pctPartnershipsVandA_1_18.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.pctPartnershipsVandA_1_25.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.pctPartnershipsVandA_1_35.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.pctPartnershipsVandA_1_45.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.pctPartnershipsVandA_1_55.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.pctPartnershipsVandA_1_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.pctPartnershipsVandA_2to4_13.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.pctPartnershipsVandA_2to4_18.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.pctPartnershipsVandA_2to4_25.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.pctPartnershipsVandA_2to4_35.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.pctPartnershipsVandA_2to4_45.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.pctPartnershipsVandA_2to4_55.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.pctPartnershipsVandA_2to4_65.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
% MSM Condom Use (added by JC on 11/07/2017)
    % Black
        CalibParams.behav_MSMpctCondomUse_Black_13.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Black_18.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Black_25.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Black_35.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Black_45.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Black_55.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Black_65.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
     % Hispanic
        CalibParams.behav_MSMpctCondomUse_Hispanic_13.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Hispanic_18.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Hispanic_25.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Hispanic_35.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Hispanic_45.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Hispanic_55.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Hispanic_65.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
     % Other
        CalibParams.behav_MSMpctCondomUse_Other_13.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Other_18.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Other_25.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Other_35.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Other_45.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Other_55.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.behav_MSMpctCondomUse_Other_65.baseValue=ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
            
    % duration in age group multiplier
        CalibParams.durInAgeMultiplier_1317.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.durInAgeMultiplier_1824.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.durInAgeMultiplier_2534.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.durInAgeMultiplier_3544.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.durInAgeMultiplier_4554.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        CalibParams.durInAgeMultiplier_5564.baseValue = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;  
        % No adjustment factor for 65+
                        
        % Collect names of fields (for collecting other values)
        ParamFields = fieldnames(CalibParams);
        % Number of fields
        nCalibratedParameters = length(ParamFields);

%% 4. Apply Lower and Upper Bounds for Input Parameters
        
    % Lower bounds
        for i = 1:nCalibratedParameters
            CalibParams.(ParamFields{i}).range(LB) = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        end
            
    % Upper bounds   
        for i = 1:nCalibratedParameters
            CalibParams.(ParamFields{i}).range(UB) = ExcelValues_CalibParams(importIndex);
            importIndex = importIndex + 1;
        end
        
    % Set paramValue for use later in Calib_updateParams (added 11/29/2018)
    for i = 1:nCalibratedParameters
        CalibParams.(ParamFields{i}).paramValue = CalibParams.(ParamFields{i}).baseValue;
    end
        
    % Remove zero range params
        CalibParams_reduced = CalibParams;
        for i = 1:nCalibratedParameters
            if CalibParams.(ParamFields{i}).range(LB) == CalibParams.(ParamFields{i}).range(UB)
                CalibParams_reduced = rmfield(CalibParams_reduced, ParamFields{i});
            end
        end
          
%    Test Vector to see if bounds match
%             counter = 0;
%         for i = 1:nCalibratedParameters
%             counter = counter + 1;
%             bound_test(counter,1) = CalibParams.(ParamFields{i}).baseValue;
%             bound_test(counter,2) = CalibParams.(ParamFields{i}).range(LB);
%             bound_test(counter,3) = CalibParams.(ParamFields{i}).range(UB);
%             if bound_test(counter,1) < bound_test(counter,2) | bound_test(counter,1) > bound_test(counter,3)
%                 bound_test(counter,4) = 1
%             end

end
            