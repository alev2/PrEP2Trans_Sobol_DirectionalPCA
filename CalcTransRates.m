function [TransRates, TTProg, Params, InfRate, InfRate_PrEP_Oral_High, InfRate_PrEP_Oral_Low, InfRate_PrEP_Inject_High, InfRate_PrEP_Inject_Low] = CalcTransRates(Params, TransRates, Year, Compartments, TTProg, ContCounterOrDiscStep)
%% Purpose: Calculates rates of each type of transition

%% Process
   
% Process (numbers do not correspond to the section numbers in the m-file)
% 1. Update TransRates struct
%   1a. New continuum of care values based on year (3 time periods)
%   1b. New infection rates based on current state of the population (every time step)

% 2. Generate 9 matrices, 1 for each transition type. The values in these 
% matrices are the rates of transitioning between each set of
% compartments. The text in brackets indicates how often the values of
% these rates need to be updated.
%   HIV progression [model start]
%   Get infected [every time step]
%   Become aware [every year]
%   Link to care [every year]
%   Drop out [every year]
%   Start ART [every year]
%   Normal death [model start]
%   AIDS death [model start]
%   Stay aged out or dead [model start]

% 3. Combine all matrices into one main TransRates matrix (discrete
% version only; in continuous this is calculated elsewhere)
%   3a. Check to be sure that no two matrices have non-zero values in the
%   same cell.
%   3b. Add matrices' values in corresponding cells together (by just
%   adding the matrices.

% Define TransRates structs
% TransRates is a 30x30x273 matrix of transition rates between all
%   compartments. The values in the matrix applied here are base values for
%   the start of the model. Many of the values in these matrices will be
%   updated every time period. Infection rates will be updated every
%   timestep.
% Rows = source states: 30 compartments based on disease status and 
%   continuum of care status + NormDeath + AIDSDeath + Aged Out
% Columns = destination states: same as rows
% Layers = stratifications: 273 different stratifications based on 
%   3 risk populations
%   2 risk levels (for HET and MSM; IDUs are all high risk) 
%   5 age groups
%   3 races
%   2 circumcision statuses (males only)
%   2 sexes

%% 1. Initialize TransRates Fields
 
    % Beginning of model
    if Year == Params.tt_modelStartYear

        % Pre-allocate key fields
        TRInitValue = 0.000;
        TransRates.BecomeAware(Params.numComparts,Params.numComparts,Params.numStrats) = TRInitValue;
        TransRates.GetInfected(Params.numComparts,Params.numComparts,Params.numStrats) = TRInitValue;
        TransRates.ProgressHIV(Params.numComparts,Params.numComparts,Params.numStrats) = TRInitValue;
        TransRates.LinkToCare(Params.numComparts,Params.numComparts,Params.numStrats) = TRInitValue;
        TransRates.DropOut(Params.numComparts,Params.numComparts,Params.numStrats) = TRInitValue;
        TransRates.StartART(Params.numComparts,Params.numComparts,Params.numStrats) = TRInitValue; %Include PrEP initiation
        TransRates.DieNormal(Params.numComparts,Params.numComparts,Params.numStrats) = TRInitValue;
        TransRates.DieAIDS(Params.numComparts,Params.numComparts,Params.numStrats) = TRInitValue;
        TransRates.StayDeadAged(Params.numComparts,Params.numComparts,Params.numStrats) = TRInitValue;
        TransRates.All(Params.numComparts,Params.numComparts,Params.numStrats) = TRInitValue;
    end
    
    % Initialize period of COVID effects to apply in T2-T4
    if Year >= Params.tt_periodTwoStartYear 
        if Year < Params.tt_periodThreeStartYear   
            COVIDperiod = Params.COVIDperiod_1;                               
        elseif  Year < Params.tt_periodFourStartYear
            COVIDperiod = Params.COVIDperiod_2;                
        else
            COVIDperiod = Params.COVIDperiod_3;                                
        end
    end
    
% Beginning of Each Year (Including Model Start)
if (Year-floor(Year)) == 0

    % % rapid, test sensitivity, and test notification for appropriate testing period
    if Year < Params.tt_YrTestSensNotifChange
        % % rapid
        TTProg.PctTestsRapid = Params.tt_pctRapid_r_TestPeriod1;
        pctRapid = TTProg.PctTestsRapid;

        % Cost per test 
        costPP_TestNegRapid = Params.costPP_TestNegRapid_TestPeriod1;
        costPP_TestPosRapid = Params.costPP_TestPosRapid_TestPeriod1;
        costPP_TestNegConv = Params.costPP_TestNegConv_TestPeriod1;
        costPP_TestPosConv = Params.costPP_TestPosConv_TestPeriod1;
            
       % Test sensitivity
        TestSens_Rapid_Acute = Params.tt_testSensAcute_Rapid_TestPeriod1;
        TestSens_Conv_Acute = Params.tt_testSensAcute_Conv_TestPeriod1;
        TestSens_Rapid_Chronic = Params.tt_testSensChronic_Rapid_TestPeriod1;
        TestSens_Conv_Chronic = Params.tt_testSensChronic_Conv_TestPeriod1;
        TestSens_Confirm_Acute = Params.tt_testSensAcute_Confirm_TestPeriod1;
        TestSens_Confirm_Chronic = Params.tt_testSensChronic_Confirm_TestPeriod1;
            
        % Notification
        NotifyRapidProb = Params.tt_probNotify_rapidPos_TestPeriod1;
        NotifyConvProb = Params.tt_probNotify_convPos_TestPeriod1;
        NotifyNegRapidProb = Params.tt_probNotify_rapidNeg_TestPeriod1;
        NotifyNegConvProb = Params.tt_probNotify_convNeg_TestPeriod1;

    else
        
        % Percent of tests rapid
        TTProg.PctTestsRapid = Params.tt_pctRapid_r_TestPeriod2;
        pctRapid = TTProg.PctTestsRapid;

        % Cost per test 
        costPP_TestNegRapid = Params.costPP_TestNegRapid_TestPeriod2;
        costPP_TestPosRapid = Params.costPP_TestPosRapid_TestPeriod2;
        costPP_TestNegConv = Params.costPP_TestNegConv_TestPeriod2;
        costPP_TestPosConv = Params.costPP_TestPosConv_TestPeriod2;
        
        % Sensitivity
        TestSens_Rapid_Acute = Params.tt_testSensAcute_Rapid_TestPeriod2;
        TestSens_Conv_Acute = Params.tt_testSensAcute_Conv_TestPeriod2;
        TestSens_Rapid_Chronic = Params.tt_testSensChronic_Rapid_TestPeriod2;
        TestSens_Conv_Chronic = Params.tt_testSensChronic_Conv_TestPeriod2;
        TestSens_Confirm_Acute = Params.tt_testSensAcute_Confirm_TestPeriod2;
        TestSens_Confirm_Chronic = Params.tt_testSensChronic_Confirm_TestPeriod2;

        % Notification
        NotifyRapidProb = Params.tt_probNotify_rapidPos_TestPeriod2;
        NotifyConvProb = Params.tt_probNotify_convPos_TestPeriod2;
        NotifyNegRapidProb = Params.tt_probNotify_rapidNeg_TestPeriod2;
        NotifyNegConvProb = Params.tt_probNotify_convNeg_TestPeriod2;
        
    end
    
    % # served by SSP for appropriate set period
    if Year < Params.tt_SEP_YrSet1Begins
        NumPWIDServedbySEP(1:Params.numRace,1) = 0; 
    elseif Year < Params.tt_SEP_YrSet2Begins
        NumPWIDServedbySEP(1:Params.numRace,1) = Params.tt_numPWIDServedbySEP_Set1;        
    elseif Year < Params.tt_SEP_YrSet3Begins
        NumPWIDServedbySEP(1:Params.numRace,1) = Params.tt_numPWIDServedbySEP_Set2;
    elseif Year < Params.tt_SEP_YrSet4Begins
        NumPWIDServedbySEP(1:Params.numRace,1) = Params.tt_numPWIDServedbySEP_Set3;
    elseif Year < Params.tt_SEP_YrSet5Begins
        NumPWIDServedbySEP(1:Params.numRace,1) = Params.tt_numPWIDServedbySEP_Set4;
    else
        NumPWIDServedbySEP(1:Params.numRace,1) = Params.tt_numPWIDServedbySEP_Set5;
    end
    
    % Save last number on SSP before last time period
    % This is to use in allocation-based programming
    if Year < Params.tt_periodFiveStartYear 
        TTProg.NumOnSEP_T1 = NumPWIDServedbySEP(1:Params.numRace,1);
    end    
    
    % PrEP initiation rates and % initiators on injectable PrEP for
    % appropriate set period
    if Year < Params.tt_PrEPInit_YrSet1Begins
        PrEPInitiationRateOrPctOnPrEP = 0;
        PctInjectPrEP = 0;
    elseif Year < Params.tt_PrEPInit_YrSet2Begins
        PrEPInitiationRateOrPctOnPrEP = Params.tt_PrEPInitRateOrPctOnPrEP_Set1;
        PctInjectPrEP = Params.tt_PctonInjectPrEP_Set1;
    elseif Year < Params.tt_PrEPInit_YrSet3Begins
        PrEPInitiationRateOrPctOnPrEP = Params.tt_PrEPInitRateOrPctOnPrEP_Set2;
        PctInjectPrEP = Params.tt_PctonInjectPrEP_Set2;
    elseif Year < Params.tt_PrEPInit_YrSet4Begins
        PrEPInitiationRateOrPctOnPrEP = Params.tt_PrEPInitRateOrPctOnPrEP_Set3;
        PctInjectPrEP = Params.tt_PctonInjectPrEP_Set3;
    elseif Year < Params.tt_PrEPInit_YrSet5Begins
        PrEPInitiationRateOrPctOnPrEP = Params.tt_PrEPInitRateOrPctOnPrEP_Set4;
        PctInjectPrEP = Params.tt_PctonInjectPrEP_Set4;
    else
        PrEPInitiationRateOrPctOnPrEP = Params.tt_PrEPInitRateOrPctOnPrEP_Set5;
        PctInjectPrEP = Params.tt_PctonInjectPrEP_Set5;
    end    
    
    %% 2. Assign progression rates that are period-specific
    for sectionPeriodSpecific = 1:1
       
        % If before time period 5, or if not using allocation-based
        % progression, set % served by SSP and PrEP initiation rates, % initiating
        % injectable (versus oral), and dropout rates
        if Year < Params.tt_periodFiveStartYear || Params.tt_progressionSource ~= 3 
            
            % Percent of active PWID served by SSP
            PctActivePWIDServedbySEP(1:Params.numRace,1) = 0;
            for nRace = 1:Params.numRace
                PctActivePWIDServedbySEP(nRace,1) = min(NumPWIDServedbySEP(nRace,1) / ...
                    max((sum(sum(Compartments(Params.UninfectedComparts,:))' .* Params.raceIndicator(:,nRace) .* Params.popIndicator(:,Params.pop_IDU)) * Params.behav_pctActivePWID),1), 1);
            end
            TTProg.PctActivePWIDServedbySEP = Params.raceIndicator * PctActivePWIDServedbySEP;
            TTProg.NumPWIDServedbySEP = NumPWIDServedbySEP;
            
            % Adjust initiation rate so that only a portion of the
            % high-risk population is eligible for PrEP. This
            % method factors percentage already on PrEP into
            % eligibility among those not on PrEP
             TotalEligPrEP = sum(Compartments(Params.UninfectedComparts,:),1)'.*Params.tt_prepPctEligible;
             EligNotOnPrEP = zeros(Params.numStrats,1);
             for counter_strats = 1:Params.numStrats 
                EligNotOnPrEP(counter_strats) = max((TotalEligPrEP(counter_strats) - sum(Compartments(Params.PrEPComparts,counter_strats))'),0);
             end
             PctOfA1PrEPElig = min((EligNotOnPrEP./(Compartments(Params.A1,:)')),1); % JCPrEPUpdate: adjusted to be no greater than 1
             PrEPInitiationRate = PrEPInitiationRateOrPctOnPrEP.*PctOfA1PrEPElig;
             
             % If in time period 2-4, apply COVID adjustment (decrease) to PrEP initiation rates
             if Year >= Params.tt_periodTwoStartYear && Year < Params.tt_periodFiveStartYear
                 
                 PrEPInitiationRate = PrEPInitiationRate .* (1 - Params.COVID_percentreduc_PrEPinit(:,COVIDperiod));
                 
             end

            % PrEP Drop off

            DropOutRate_PrEP_Oral_HighAdherence = Params.tt_rateDropOffPrEP_Oral_HighAdherence;
            DropOutRate_PrEP_Oral_LowAdherence = Params.tt_rateDropOffPrEP_Oral_LowAdherence;
            DropOutRate_PrEP_Inject_HighAdherence = Params.tt_rateDropOffPrEP_Inject_HighAdherence;
            DropOutRate_PrEP_Inject_LowAdherence = Params.tt_rateDropOffPrEP_Inject_LowAdherence;           
            
        end        
        
        if Year < Params.tt_periodTwoStartYear 
            
            %% 2.i. Period 1 TTProgression Values
            for sectionPeriod1 = 1:1
                                     
            % Testing

                TestRate(:,1) = Params.tt_testRateUninfected_rp_1;
                TestRate(:,2) = Params.tt_testRateAcute_rp_1;
                TestRate(:,3) =  Params.tt_testRateLatentA_rp_1;
                TestRate(:,4) =  Params.tt_testRateLatentB_rp_1;
                TestRate(:,5) = Params.tt_testRateLate_rp_1;
                TestRate(:,6) =  Params.tt_testRateAIDS_rp_1;
         
            % Linkage
                LinkageFirstProb = Params.tt_linkageFirst_r_1;

                LinkageAfterRate(:,1) = Params.tt_linkageAfterRate_Acute_1;
                LinkageAfterRate(:,2) = Params.tt_linkageAfterRate_LatentA_1; 
                LinkageAfterRate(:,3) = Params.tt_linkageAfterRate_LatentB_1;
                LinkageAfterRate(:,4) = Params.tt_linkageAfterRate_Late_1;
                LinkageAfterRate(:,5) = Params.tt_linkageAfterRate_AIDS_1;
         
            % Drop-out
                DropOutRate_CareToAware = Params.tt_dropOutRate_CareToAware_1;
                DropOutRate_ANVToAware = Params.tt_dropOutRate_ANVToAware_1;
                DropOutRate_ANVToCare = Params.tt_dropOutRate_ANVToCare_1;
                DropOutRate_VLSToANV = Params.tt_dropOutRate_VLSToANV_1;
                DropOutRate_VLSToLTC = Params.tt_dropOutRate_VLSToLTC_1;
                DropOutRate_VLSToAware = Params.tt_dropOutRate_VLSToAware_1;
         

            % ART initiation from LTC, no ART effects
                InitiateARTRate(:,1) = Params.tt_ARTInitRateAcute_1; 
                InitiateARTRate(:,2) = Params.tt_ARTInitRateLatentA_1; 
                InitiateARTRate(:,3) = Params.tt_ARTInitRateLatentB_1;
                InitiateARTRate(:,4) = Params.tt_ARTInitRateLate_1;
                InitiateARTRate(:,5) = Params.tt_ARTInitRateAIDS_1;

            % ART initiation from ANV
                ANVToVLSRate = Params.tt_ARTInitRateFromANV_1;
                
            % VLS
            PctWhoBecomeVLS = Params.tt_PctInitiateARTWhoBecomeVLS;

            % In Care ANV
            PctANVWhoAreInCareInclART = Params.tt_PctANVWhoAreInCareOrART_1;
            
            % ART Eligibility
            TTProg.ARTElig = Params.tt_ARTElig_1;
            
            % All-cause mortality
            pop_dyingRate = Params.pop_dyingRate;
            
            %Not Affected by ART Death
            NoARTEffectsDeathRate_NotAIDS = Params.rate_deathnoARTEffectsNotAIDS_1; 

           
            % Affected by ART Death %change to ANV
            
            ANVDeathRate_NotAIDS = Params.rate_deathANVNotAIDS_1;
            ANVDeathRate_AIDS = Params.rate_deathANVAIDS_1;
            
            
            %VLS Death
            VLSDeathRate_NotAIDS = Params.rate_deathVLSNotAIDS_1;
            VLSDeathRate_AIDS = Params.rate_deathVLSAIDS_1;
            

            end
            
        elseif Year < Params.tt_periodFiveStartYear
            %% 2.ii. Period 2-4 TTProgression Values
            for sectionPeriod2to4 = 1:1
                                       
            % Testing (with COVID adjustments)
                
                % Uninfected
                TestRate(:,1) = Params.tt_testRateUninfected_rp_2to4 .* (1 - Params.COVID_percentreduc_testrate_Uninf(:,COVIDperiod));
                % Undiagnosed PWH
                TestRate(:,2) = Params.tt_testRateAcute_rp_2to4 .* (1 - Params.COVID_percentreduc_testrate_PWH(:,COVIDperiod));
                TestRate(:,3) =  Params.tt_testRateLatentA_rp_2to4 .* (1 - Params.COVID_percentreduc_testrate_PWH(:,COVIDperiod));
                TestRate(:,4) =  Params.tt_testRateLatentB_rp_2to4 .* (1 - Params.COVID_percentreduc_testrate_PWH(:,COVIDperiod));
                TestRate(:,5) = Params.tt_testRateLate_rp_2to4 .* (1 - Params.COVID_percentreduc_testrate_PWH(:,COVIDperiod));
                TestRate(:,6) =  Params.tt_testRateAIDS_rp_2to4 .* (1 - Params.COVID_percentreduc_testrate_PWH(:,COVIDperiod));
                            
            % Linkage
                LinkageFirstProb = Params.tt_linkageFirst_r_2to4;

                LinkageAfterRate(:,1) = Params.tt_linkageAfterRate_Acute_2to4;
                LinkageAfterRate(:,2) = Params.tt_linkageAfterRate_LatentA_2to4; 
                LinkageAfterRate(:,3) = Params.tt_linkageAfterRate_LatentB_2to4;
                LinkageAfterRate(:,4) = Params.tt_linkageAfterRate_Late_2to4;
                LinkageAfterRate(:,5) = Params.tt_linkageAfterRate_AIDS_2to4;

            % Drop-out
                DropOutRate_CareToAware = Params.tt_dropOutRate_CareToAware_2to4;
                DropOutRate_ANVToAware = Params.tt_dropOutRate_ANVToAware_2to4;
                DropOutRate_ANVToCare = Params.tt_dropOutRate_ANVToCare_2to4;
                DropOutRate_VLSToANV = Params.tt_dropOutRate_VLSToANV_2to4;
                DropOutRate_VLSToLTC = Params.tt_dropOutRate_VLSToLTC_2to4;
                DropOutRate_VLSToAware = Params.tt_dropOutRate_VLSToAware_2to4;
                
            % COVID adjustment (increase) to drop-out rates
                DropOutRate_ANVToCare = DropOutRate_ANVToCare .* (1 + Params.COVID_percentincr_DropOutARTfromANV(:,COVIDperiod)); 
                DropOutRate_VLSToANV = DropOutRate_VLSToANV .* (1 + Params.COVID_percentincr_LoseVLStoANV(:,COVIDperiod));
            
            % ART initiation from LTC, no ART effects
                InitiateARTRate(:,1) = Params.tt_ARTInitRateAcute_2to4; 
                InitiateARTRate(:,2) = Params.tt_ARTInitRateLatentA_2to4; 
                InitiateARTRate(:,3) = Params.tt_ARTInitRateLatentB_2to4;
                InitiateARTRate(:,4) = Params.tt_ARTInitRateLate_2to4;
                InitiateARTRate(:,5) = Params.tt_ARTInitRateAIDS_2to4;
                
            % COVID adjustment (decrease) to ART initiation rates
                InitiateARTRate = InitiateARTRate .* (1 - Params.COVID_percentreduc_initART(:,COVIDperiod));

            % ART initiation from ANV
                ANVToVLSRate = Params.tt_ARTInitRateFromANV_2to4;
                
            % VLS
                PctWhoBecomeVLS = Params.tt_PctInitiateARTWhoBecomeVLS;
            
            % In Care ANV
            PctANVWhoAreInCareInclART = Params.tt_PctANVWhoAreInCareOrART_2to4;
            
            % ART Eligibility
            TTProg.ARTElig = Params.tt_ARTElig_2to4;
            
            % All-cause mortality
            pop_dyingRate = Params.pop_dyingRate;
            
            %Not Affected by ART Death
            NoARTEffectsDeathRate_NotAIDS = Params.rate_deathnoARTEffectsNotAIDS_2to5;


            % ART Death prob

            ANVDeathRate_NotAIDS = Params.rate_deathANVNotAIDS_2to5;
            ANVDeathRate_AIDS = Params.rate_deathANVAIDS_2to5;

            %VLS Death
            VLSDeathRate_NotAIDS = Params.rate_deathVLSNotAIDS_2to5;
            VLSDeathRate_AIDS = Params.rate_deathVLSAIDS_2to5;
            
            % COVID adjustment (increase) to mortality rates
                NoARTEffectsDeathRate_NotAIDS = NoARTEffectsDeathRate_NotAIDS .* (1 + Params.COVID_percentincr_Mortality(:,COVIDperiod));
                ANVDeathRate_NotAIDS = ANVDeathRate_NotAIDS .* (1 + Params.COVID_percentincr_Mortality(:,COVIDperiod));
                ANVDeathRate_AIDS = ANVDeathRate_AIDS .* (1 + Params.COVID_percentincr_Mortality(:,COVIDperiod));
                VLSDeathRate_NotAIDS = VLSDeathRate_NotAIDS .* (1 + Params.COVID_percentincr_Mortality(:,COVIDperiod));
                VLSDeathRate_AIDS = VLSDeathRate_AIDS .* (1 + Params.COVID_percentincr_Mortality(:,COVIDperiod));
                % Also apply adjustment to all-cause mortality
                pop_dyingRate = pop_dyingRate .* (1 + Params.COVID_percentincr_Mortality(:,COVIDperiod));
            
            end

        else
            %% 2.iii. Period 5 TTProgression Values
            for sectionPeriod5 = 1:1
                
            if Params.tt_progressionSource == 3 % allocation-based
                % Entire section revised by KH in June 2019
                
                % 2.iii.1 Allocation-Based Progression
                for sectionAllocation = 1:1
    
                %Multiplicative adjustment to apply to allocation
                    alloc_mult = 1000000;
                    
                % Population indicators
                Indicators_HET = Params.popIndicator(:,Params.pop_HET);
                Indicators_HET_Low = Params.LowRiskHETIndicator;
                Indicators_HET_High = Params.HRHIndicator;
                Indicators_HETM_High = Params.popIndicator_withHETbySex(:,Params.popHETbySex_HETM) .* Params.riskLevelIndicator(:,Params.risk_Casual);
                Indicators_HETF_High = Params.popIndicator_withHETbySex(:,Params.popHETbySex_HETF) .* Params.riskLevelIndicator(:,Params.risk_Casual);
                Indicators_MSM = Params.popIndicator(:,Params.pop_MSM);
                Indicators_MSM_Low = Params.riskLevelIndicator(:,Params.risk_Main).*Params.popIndicator(:,Params.pop_MSM);    
                Indicators_MSM_High = Params.riskLevelIndicator(:,Params.risk_Casual).*Params.popIndicator(:,Params.pop_MSM);     
                Indicators_IDU = Params.popIndicator(:,Params.pop_IDU);   
                if Params.TargetIntnsToYMSM == 0
                    Indicators_LimitPops = ones(Params.numStrats,1);
                else
                    Indicators_LimitPops = Params.YoungMSMIndicator;
                end
                Indicators_PrEPPops = Indicators_HET_High + Indicators_MSM_High + Indicators_IDU;
                Indicator_SEP = 0; % indicator set to 0 until SEP intervention
                Indicator_IsPrEP = 0; % indicator set to 0 until PrEP intervention
                
                % Set allocation period
                if Year < Params.Intn_Allocation_StartYr(2)
                    allocationPeriod = 1;
                elseif Year < Params.Intn_Allocation_StartYr(3)
                    allocationPeriod = 2;
                else
                    allocationPeriod = 3;
                end

%                 % Current interventions 
%                     intn_Testing_HET_Investment_Low = Params.intn_Testing_Investment_HET_Low;
%                     intn_Testing_HET_Investment_High = Params.intn_Testing_HET_Investment_High;
%                     intn_Testing_MSM_Investment_Low = Params.intn_Testing_MSM_Investment_Low;
%                     intn_Testing_MSM_Investment_High = Params.intn_Testing_MSM_Investment_High;
%                     intn_Testing_IDU_Investment = Params.intn_Testing_IDU_Investment;
%                     intn_LTCatDiag_Investment = Params.intn_LTCatDiag_Investment;
%                     intn_ARTInitiation_Investment = Params.intn_ARTInitiation_Investment;
%                     intn_TxAdherence_Investment = Params.intn_TxAdherence_Investment;
%                     intn_PrEP_HET_Investment = Params.intn_PrEP_HET_Investment;
%                     intn_PrEP_MSM_Investment = Params.intn_PrEP_MSM_Investment;
%                     intn_PrEP_IDU_Investment = Params.intn_PrEP_IDU_Investment;
                
               % Set up parameters that change in period 5 but aren't
               % affected by allocation
               for sectionOthProgression = 1:1

                   % InitializeCompartments
                   Comparts_CCStage4 = Params.ANVComparts; 
                   Comparts_CCStage5 = Params.VLSComparts;
                   
                    % In-care ANV
                        % Note: if this is ever updated by the alloc, update how the parameter is used in CollectResults
                    PctANVWhoAreInCareInclART = Params.tt_PctANVWhoAreInCareOrART_5; 
                    PctANVInCareWhoAreOnART = Params.tt_PctInCareANVWhoAreOnART_5;
                    
                    % ART-eligibility
                    TTProg.ARTElig= Params.tt_ARTElig_5;

                    % Drop-out rates 
                    DropOutRate_CareToAware = Params.tt_dropOutRate_CareToAware_5;
                    DropOutRate_ANVToAware = Params.tt_dropOutRate_ANVToAware_5;
                    DropOutRate_ANVToCare = Params.tt_dropOutRate_ANVToCare_5;
                    DropOutRate_VLSToANV = Params.tt_dropOutRate_VLSToANV_5;
                    DropOutRate_PrEP_Oral_HighAdherence = Params.tt_rateDropOffPrEP_Oral_HighAdherence;
                    DropOutRate_PrEP_Oral_LowAdherence = Params.tt_rateDropOffPrEP_Oral_LowAdherence;
                    DropOutRate_PrEP_Inject_HighAdherence = Params.tt_rateDropOffPrEP_Inject_HighAdherence;
                    DropOutRate_PrEP_Inject_LowAdherence = Params.tt_rateDropOffPrEP_Inject_LowAdherence;
                    DropOutRate_VLSToLTC = Params.tt_dropOutRate_VLSToLTC_5;
                    DropOutRate_VLSToAware = Params.tt_dropOutRate_VLSToAware_5;

                   % Rate of linkage to Care After Diagnosis                
                    LinkageAfterRate(:,1) = Params.tt_linkageAfterRate_Acute_5;
                    LinkageAfterRate(:,2) = Params.tt_linkageAfterRate_LatentA_5; 
                    LinkageAfterRate(:,3) = Params.tt_linkageAfterRate_LatentB_5;
                    LinkageAfterRate(:,4) = Params.tt_linkageAfterRate_Late_5;
                    LinkageAfterRate(:,5) = Params.tt_linkageAfterRate_AIDS_5;
                    
                    % All-cause mortality
                    pop_dyingRate = Params.pop_dyingRate;
                    
                    %Not Affected by ART Death
                    NoARTEffectsDeathRate_NotAIDS = Params.rate_deathnoARTEffectsNotAIDS_2to5;


                    % ART Death prob
                    
                    ANVDeathRate_NotAIDS = Params.rate_deathANVNotAIDS_2to5;
                    ANVDeathRate_AIDS = Params.rate_deathANVAIDS_2to5;
                    
                    %VLS Death
                    VLSDeathRate_NotAIDS = Params.rate_deathVLSNotAIDS_2to5;
                    VLSDeathRate_AIDS = Params.rate_deathVLSAIDS_2to5;
                    
                    % Transition to VLS from ANV 
                    ANVToVLSRate = Params.tt_ARTInitRateFromANV_5;
                    
                    % VLS
                    PctWhoBecomeVLS = Params.tt_PctInitiateARTWhoBecomeVLS;
               end

               % INTERVENTION: TESTING
               for sectionIntnTest = 1:1

                   % Define baseline variables 
                       
                    % Testing rates without allocation
                    baseTestRate(:,1) = Params.tt_testRateUninfected_rp_5;
                    baseTestRate(:,2) = Params.tt_testRateAcute_rp_5;
                    baseTestRate(:,3) = Params.tt_testRateLatentA_rp_5;
                    baseTestRate(:,4) = Params.tt_testRateLatentB_rp_5;
                    baseTestRate(:,5) = Params.tt_testRateLate_rp_5;
                    baseTestRate(:,6) = Params.tt_testRateAIDS_rp_5;
                    
                    % Testing rates in T1, used in allocation-based
                    % progression to allocate funds using the same relative
                    % differences as are reflected in these rates
                    T1TestRate(:,1) = Params.tt_testRateUninfected_rp_1;
                    T1TestRate(:,2) = Params.tt_testRateAcute_rp_1;
                    T1TestRate(:,3) =  Params.tt_testRateLatentA_rp_1;
                    T1TestRate(:,4) =  Params.tt_testRateLatentB_rp_1;
                    T1TestRate(:,5) = Params.tt_testRateLate_rp_1;
                    T1TestRate(:,6) =  Params.tt_testRateAIDS_rp_1;


                    % Number eligible for testing 
                    TestEligComparts = [Params.A1, Params.UnawareComparts];
                    numEligTest(1:Params.numStrats,1:Params.numTestEligStates) = 0;
                    for TestEligStateNum = 1:Params.numTestEligStates
                        numEligTest(:,TestEligStateNum) = ...
                            Compartments(TestEligComparts(TestEligStateNum),:)' ...
                            .* Indicators_LimitPops;
                    end

                    % Cost - Notify positive results
                    costPP_NotifyPosRapid = Params.costPP_Notify_PosRapid;
                    costPP_NotifyPosConv = Params.costPP_Notify_PosConv;

                    % Cost - notify negative results
                    costPP_NotifyNegRapid = Params.costPP_Notify_NegRapid;
                    costPP_NotifyNegConv = Params.costPP_Notify_NegConv;

                    % Number eligible for testing 
                    costPP_Test(1:Params.numStrats,1:Params.numTestEligStates)=0;
                    for TestEligStateNum = 1:Params.numTestEligStates
                        switch TestEligStateNum
                            case 1 % Cost per uninfected person tested (apply to numUninf that are tested)
                                costPP_Test(:,TestEligStateNum) = ...
                                    pctRapid .* (costPP_TestNegRapid + Params.costPP_TestOutreach'+costPP_NotifyNegRapid) ...
                                        + (1 - pctRapid) .* (costPP_TestNegConv + Params.costPP_TestOutreach' + costPP_NotifyNegConv);
                            case 2 % Cost per unaware acute person tested (apply to numUnawareAcute that are tested)
                                costPP_Test(:,TestEligStateNum) = ...
                                    (pctRapid * TestSens_Rapid_Acute * (1 - TestSens_Confirm_Acute) .*  ...
                                        (costPP_TestPosRapid + Params.costPP_TestOutreach' + Params.costPP_NATTest + NotifyRapidProb*costPP_NotifyPosRapid)) + ...
                                    (pctRapid * TestSens_Rapid_Acute * TestSens_Confirm_Acute .*  ...
                                        (costPP_TestPosRapid + Params.costPP_TestOutreach' + NotifyRapidProb*costPP_NotifyPosRapid)) + ...
                                    ((1 - pctRapid) * TestSens_Conv_Acute * (1 - TestSens_Confirm_Acute) .* ...
                                        (costPP_TestPosConv + Params.costPP_TestOutreach' + Params.costPP_NATTest + NotifyConvProb*costPP_NotifyPosConv)) + ...
                                    ((1 - pctRapid) * TestSens_Conv_Acute * TestSens_Confirm_Acute .* ...
                                        (costPP_TestPosConv + Params.costPP_TestOutreach' + NotifyConvProb*costPP_NotifyPosConv)) + ...
                                    (pctRapid * (1-TestSens_Rapid_Acute)  .*  ...
                                        (costPP_TestNegRapid + Params.costPP_TestOutreach' + NotifyRapidProb*costPP_NotifyNegRapid)) + ...
                                    ((1 - pctRapid) * (1-TestSens_Conv_Acute)  .* ...
                                        (costPP_TestNegConv + Params.costPP_TestOutreach' + NotifyConvProb*costPP_NotifyNegConv));
                            case 3 % Cost per unaware chronic person tested (apply to numUnawareChronic that are tested)
                                costPP_Test(:,TestEligStateNum) = ...
                                    (pctRapid * TestSens_Rapid_Chronic * (1 - TestSens_Confirm_Chronic) .*  ...
                                        (costPP_TestPosRapid + Params.costPP_TestOutreach' + Params.costPP_NATTest + NotifyRapidProb*costPP_NotifyPosRapid)) + ...
                                    (pctRapid * TestSens_Rapid_Chronic * TestSens_Confirm_Chronic .*  ...
                                        (costPP_TestPosRapid + Params.costPP_TestOutreach' + NotifyRapidProb*costPP_NotifyPosRapid)) + ...
                                    ((1 - pctRapid) * TestSens_Conv_Chronic * (1 - TestSens_Confirm_Chronic) .* ...
                                        (costPP_TestPosConv + Params.costPP_TestOutreach' + Params.costPP_NATTest + NotifyConvProb*costPP_NotifyPosConv)) + ...
                                    ((1 - pctRapid) * TestSens_Conv_Chronic * TestSens_Confirm_Chronic .* ...
                                        (costPP_TestPosConv + Params.costPP_TestOutreach' + NotifyConvProb*costPP_NotifyPosConv)) + ...
                                    (pctRapid * (1-TestSens_Rapid_Chronic)  .*  ...
                                        (costPP_TestNegRapid + Params.costPP_TestOutreach' + NotifyRapidProb*costPP_NotifyNegRapid)) + ...
                                    ((1 - pctRapid) * (1-TestSens_Conv_Chronic)  .* ...
                                        (costPP_TestNegConv + Params.costPP_TestOutreach' + NotifyConvProb*costPP_NotifyNegConv));
                            otherwise
                                costPP_Test(:,TestEligStateNum) = costPP_Test(:,3);
                        end
                    end
                        
                    InputType_Test = 1; % Input is rate
                    
                    % Calculate NumGetIntn for all reach rates
                    % (T1, R1, R2, R3) - result = 273 x 4

                    NumGetIntnAtT1andReachLvlRate_Test(1:Params.numStrats,1:Params.numReachLvls,1:Params.numTestEligStates)=0;
                    T1andReachLvlRate_Test(1:Params.numStrats,1:Params.numReachLvls,1:Params.numTestEligStates)=0;
                    CostT1andReachLvlRate_Test(1:Params.numStrats,1:Params.numReachLvls,1:Params.numTestEligStates)=0;
                    
                    % Calculate T1 and reach level rates and costs of each
                    % rate for each disease stage (including uninfected)
                    for TestEligStateNum = 1:Params.numTestEligStates
                        % Calculate number getting intervention at T1 and
                        % reach level rates
                        NumGetIntnAtT1andReachLvlRate_Test(:,:,TestEligStateNum) = ...
                            CalcNumIntnAtT1andReachLvlRates(numEligTest(:,TestEligStateNum),  ...
                                T1TestRate(:,TestEligStateNum), Params.intn_Testing_ReachLvls, Params.intn_Testing_MaxReach, Indicators_LimitPops, InputType_Test);
                        % Calculate T1 and reach level rates
                        T1andReachLvlRate_Test(:,:,TestEligStateNum) = ...
                            CalcT1andReachLvlIntnsPP(numEligTest(:,TestEligStateNum),  ...
                                NumGetIntnAtT1andReachLvlRate_Test(:,:,TestEligStateNum), Indicators_LimitPops, Indicator_SEP);
                        % Calculate cost of achieving T1 and reach level rates    
                        CostT1andReachLvlRate_Test(:,:,TestEligStateNum) = ...
                            CalcCostT1andReachLvl(costPP_Test(:,TestEligStateNum), Params.intn_Testing_PctCostatReachLvl, ...
                                numEligTest(:,TestEligStateNum), ...
                                baseTestRate(:,TestEligStateNum), ...
                                NumGetIntnAtT1andReachLvlRate_Test(:,:,TestEligStateNum), ...
                                InputType_Test, Indicator_IsPrEP, 0);
                    end
                    
                    PctT1andReachLvlRatesAffordedByAlloc_Test(1:Params.numStrats,1:Params.numReachLvls) = 0;
                    TestRate(Params.numStrats,Params.numTestEligStates)=0;

                    %Calculate percentage of T1 rates and each reach level
                    %afforded by allocations to the testing interventions
                    for TestSubPopNum = 1:5

                        switch TestSubPopNum %We have testing interventions targeted to each of these subpops                           
                            case 1                                
                                subPopIndicator = Indicators_HET_Low;
                                TestingInvestment = Params.intn_Testing_Investment_HET_Low(allocationPeriod) * alloc_mult;
                            case 2                                
                                subPopIndicator = Indicators_HET_High;
                                TestingInvestment = Params.intn_Testing_Investment_HET_High(allocationPeriod) * alloc_mult;
                            case 3                                
                                subPopIndicator = Indicators_MSM_Low;
                                TestingInvestment = Params.intn_Testing_Investment_MSM_Low(allocationPeriod) * alloc_mult;
                            case 4                                
                                subPopIndicator = Indicators_MSM_High;
                                TestingInvestment = Params.intn_Testing_Investment_MSM_High(allocationPeriod) * alloc_mult;
                            case 5                                
                                subPopIndicator = Indicators_IDU;
                                TestingInvestment = Params.intn_Testing_Investment_IDU(allocationPeriod) * alloc_mult;
                        end
                            PctT1andReachLvlRatesAffordedByAlloc_Test = ...
                                min(PctT1andReachLvlRatesAffordedByAlloc_Test + ...
                                    CalcAddlPctT1andReachLvlRatesAffordedByAlloc(TestingInvestment, ...
                                    CostT1andReachLvlRate_Test, subPopIndicator .* Indicators_LimitPops, Indicator_SEP),1);
                    end %TestSubPopNum
        
                    % Calculate final test rate as sum of base rates and
                    % additional rates due to allocation. 
                    for TestEligStateNum = 1:Params.numTestEligStates
                         TestRate(:,TestEligStateNum) = baseTestRate(:,TestEligStateNum) + ...
                            sum(PctT1andReachLvlRatesAffordedByAlloc_Test(:,:) .* ...
                            (T1andReachLvlRate_Test(:,:,TestEligStateNum) - ...
                            [zeros(Params.numStrats,1) T1andReachLvlRate_Test(:,1:Params.numReachLvls-1,TestEligStateNum)]),2);
                    end
                    
                    % If targeting to YMSM only, add non-YMSM rates to
                    % final testing rates
                    if and(Params.TargetIntnsToYMSM == 1, Params.tt_progressionSource == 3)
                        TestRate(:,1) = TestRate(:,1) + Params.tt_testRateUninfected_rp_5_nonYMSM;
                        TestRate(:,2) = TestRate(:,2) + Params.tt_testRateAcute_rp_5_nonYMSM;
                        TestRate(:,3) = TestRate(:,3) + Params.tt_testRateLatentA_rp_5_nonYMSM;
                        TestRate(:,4) = TestRate(:,4) + Params.tt_testRateLatentB_rp_5_nonYMSM;
                        TestRate(:,5) = TestRate(:,5) + Params.tt_testRateLate_rp_5_nonYMSM;
                        TestRate(:,6) = TestRate(:,6) + Params.tt_testRateAIDS_rp_5_nonYMSM;
                    end
               end  
               
               % Collect $$$ needed for each reach level based on
               % CostT1andReachLvlRate_Test - this is used to determine the
               % starting funding for the cost minimization optimization
               TTProg.CostT1andReachLvlRate_Test_LowRiskHETs = sum(sum(CostT1andReachLvlRate_Test .* Indicators_HET_Low,3),1);
               TTProg.CostT1andReachLvlRate_Test_HighRiskHETs = sum(sum(CostT1andReachLvlRate_Test .* Indicators_HET_High,3),1);
               TTProg.CostT1andReachLvlRate_Test_LowRiskMSM = sum(sum(CostT1andReachLvlRate_Test .* Indicators_MSM_Low,3),1);
               TTProg.CostT1andReachLvlRate_Test_HighRiskMSM = sum(sum(CostT1andReachLvlRate_Test .* Indicators_MSM_High,3),1);
               TTProg.CostT1andReachLvlRate_Test_IDU = sum(sum(CostT1andReachLvlRate_Test .* Indicators_IDU,3),1);              

                % INTERVENTION: LTC FIRST
               for sectionIntnLTCFirst = 1:1

                   % Formula:
                        % [Prob of LTCFirst without CDC funds] + 
                        % [CDC investment in LTCFirst]/[CostPP for LTC at Diag]/
                        % ([Num HIV+ people tested]* [Test Sens] * [ProbNotify])
                                
                    % Calculate number eligible for LTC First 
                    % Eligible = just notified of diagnosis 
                    % Note: Test sens captured here
                    numExpectedPosNotified_Rapid = ...
                        pctRapid .* (...
                            TestSens_Rapid_Acute .* Compartments(Params.B1,:)' .* TestRate(:,2) + ...
                            TestSens_Rapid_Chronic .* (...
                            Compartments(Params.C1,:)' .* TestRate(:,3) + ...
                            Compartments(Params.D1,:)' .* TestRate(:,4) + ...
                            Compartments(Params.E1,:)' .* TestRate(:,5) + ...
                            Compartments(Params.F1,:)' .* TestRate(:,6))).* NotifyRapidProb';
                    numExpectedPosNotified_Conv = ...
                        (1-pctRapid) .* (...
                            TestSens_Conv_Acute .* Compartments(Params.B1,:)' .* TestRate(:,2) + ...
                            TestSens_Conv_Chronic .* (...
                            Compartments(Params.C1,:)' .* TestRate(:,3) + ...
                            Compartments(Params.D1,:)' .* TestRate(:,4) + ...
                            Compartments(Params.E1,:)' .* TestRate(:,5) + ...
                            Compartments(Params.F1,:)' .* TestRate(:,6))).* NotifyConvProb'; 
                        
                    % Expected total number notified
                    numExpectedPosNotifiedTotal = ...
                          numExpectedPosNotified_Rapid + numExpectedPosNotified_Conv;

                    % Calculate probability resulting from base probability 
                    % (without allocation) and allocation 
                    % Applies code using functions below (KH 26Apr2019)

                    % Base probability
                    LinkageFirstProb_base = Params.tt_linkageFirst_r_5;
                    
                    % Calculate total probability based on base prob and allocation
                    InputType_LTCatDiag = 2;
                    
                    % Calculate number getting intervention at T1 and reach level probs                                       
                    NumGetIntnAtT1andReachLvlProbs_LTCatDiag = ...
                        CalcNumIntnAtT1andReachLvlRates(numExpectedPosNotifiedTotal, ...
                            Params.tt_linkageFirst_r_1, Params.intn_LTCatDiag_ReachLvls, Params.intn_LTCatDiag_MaxReach, Indicators_LimitPops, InputType_LTCatDiag);
                    % Calculate T1 and reach level probs    
                    T1andReachLvlProb_LTCatDiag = ...
                        CalcT1andReachLvlIntnsPP(numExpectedPosNotifiedTotal,  ...
                            NumGetIntnAtT1andReachLvlProbs_LTCatDiag, Indicators_LimitPops, Indicator_SEP);
                    % Calculate cost of achieving T1 and reach level probs    
                    CostT1andReachLvlProb_LTCatDiag = ...
                        CalcCostT1andReachLvl(Params.costPP_LTCFirst, Params.intn_LTCatDiag_PctCostatReachLvl,...
                            numExpectedPosNotifiedTotal, LinkageFirstProb_base, ...
                            NumGetIntnAtT1andReachLvlProbs_LTCatDiag, ...
                            InputType_LTCatDiag, Indicator_IsPrEP, 0);
                        
                    Alloc_LTCatDiag = (Params.intn_LTCatDiag_Investment(allocationPeriod) * alloc_mult);  
                                                                
                    %Calculate percentage of T1 rates and each reach level
                    %afforded by allocations to the LTC at diagnosis intervention
                    PctT1andReachLvlRateAffordedByAlloc_LTCatDiag = ...
                        min(CalcAddlPctT1andReachLvlRatesAffordedByAlloc(Alloc_LTCatDiag, ...
                        CostT1andReachLvlProb_LTCatDiag, Indicators_LimitPops, Indicator_SEP),1);
                    
                    % Calculate final LTC at diagnsosis rate as sum of base probs and
                    % additional probs due to allocation.
                    LinkageFirstProb = LinkageFirstProb_base + ...
                            sum(PctT1andReachLvlRateAffordedByAlloc_LTCatDiag .* ...
                            (T1andReachLvlProb_LTCatDiag - ...
                            [zeros(Params.numStrats,1) T1andReachLvlProb_LTCatDiag(:,1:Params.numReachLvls-1)]),2);
 
                    % If targeting to YMSM only, add non-YMSM rates to
                    % final linkage rates
                    if and(Params.TargetIntnsToYMSM == 1, Params.tt_progressionSource == 3)
                        LinkageFirstProb = LinkageFirstProb + Params.tt_linkageFirst_r_5_nonYMSM;
                    end
                    
                    % Collect $$$ needed for each reach level based on
                    % CostT1andReachLvlProb_LTCatDiag - this is used to determine the
                    % starting funding for the cost minimization optimization
                    TTProg.CostT1andReachLvlProb_LTCatDiag = sum(CostT1andReachLvlProb_LTCatDiag,1);
               end
               
               % INTERVENTION: LTC AFTER
               for sectionIntnLTCAfter = 1:1

                   % Formula:
                        % [Prob of LTCAfter without CDC funds] + 
                        % [CDC investment in LTCAfter]/[CostPP for LTC after Diag]/
                        % [Num HIV+ people aware and not LTC]
                                
                    % Eligible = aware and not linked to care
                    
                    % Calculate rate resulting from base rate 
                    % (without allocation) and allocation for each disease stage 
                    % Applies code using functions below (KH 2May2019)
                    
                    InputType_LTCafterDiag = 1; % Input is rate
                    
                    % Number eligible for LTC after diag 
                    LTCAfterEligComparts = [Params.B2, Params.C2, Params.D2, Params.E2, Params.F2];
                    numEligLTCAfter(1:Params.numStrats,1:Params.numHIVstages) = 0;
                    for HIVStageNum = 1:Params.numHIVstages
                        numEligLTCAfter(:,HIVStageNum) = ...
                            Compartments(LTCAfterEligComparts(HIVStageNum),:)' ...
                            .* Indicators_LimitPops;
                    end
                    
                   % Define baseline variables 
                       
                    % Linkage after rates without allocation
                    baseLinkageAfterRate(:,1) = Params.tt_linkageAfterRate_Acute_5;
                    baseLinkageAfterRate(:,2) = Params.tt_linkageAfterRate_LatentA_5;
                    baseLinkageAfterRate(:,3) =  Params.tt_linkageAfterRate_LatentB_5;
                    baseLinkageAfterRate(:,4) =  Params.tt_linkageAfterRate_Late_5;
                    baseLinkageAfterRate(:,5) = Params.tt_linkageAfterRate_AIDS_5;

                    % Linkage rates after diagnosis without allocation
                    T1LinkageAfterRate(:,1) = Params.tt_linkageAfterRate_Acute_1;
                    T1LinkageAfterRate(:,2) = Params.tt_linkageAfterRate_LatentA_1;
                    T1LinkageAfterRate(:,3) =  Params.tt_linkageAfterRate_LatentB_1;
                    T1LinkageAfterRate(:,4) =  Params.tt_linkageAfterRate_Late_1;
                    T1LinkageAfterRate(:,5) = Params.tt_linkageAfterRate_AIDS_1;
                    
                    NumGetIntnAtT1andReachLvlRate_LTCafterDiag(1:Params.numStrats,1:Params.numReachLvls,1:Params.numHIVstages)=0;
                    T1andReachLvlRate_LTCafterDiag(1:Params.numStrats,1:Params.numReachLvls,1:Params.numHIVstages)=0;
                    CostT1andReachLvlRate_LTCafterDiag(1:Params.numStrats,1:Params.numReachLvls,1:Params.numHIVstages)=0;
                    
                    % Calculate T1 and reach level rates and costs of each
                    % rate for each disease stage
                    for HIVStageNum = 1:Params.numHIVstages
                        % Calculate number getting intervention at T1 and
                        % reach level rates
                        NumGetIntnAtT1andReachLvlRate_LTCafterDiag(:,:,HIVStageNum) = ...
                            CalcNumIntnAtT1andReachLvlRates(numEligLTCAfter(:,HIVStageNum),  ...
                                T1LinkageAfterRate(:,HIVStageNum), ...
                                Params.intn_LTCafterDiag_ReachLvls, ...
                                Params.intn_LTCafterDiag_MaxReach, ...
                                Indicators_LimitPops, InputType_LTCafterDiag);
                        % Calculate T1 and reach level rates    
                        T1andReachLvlRate_LTCafterDiag(:,:,HIVStageNum) = ...
                            CalcT1andReachLvlIntnsPP(numEligLTCAfter(:,HIVStageNum),  ...
                                NumGetIntnAtT1andReachLvlRate_LTCafterDiag(:,:,HIVStageNum), Indicators_LimitPops, Indicator_SEP);
                        % Calculate costs of achieving T1 and reach level rates    
                        CostT1andReachLvlRate_LTCafterDiag(:,:,HIVStageNum) = ...
                            CalcCostT1andReachLvl(Params.costPP_LTCAfter, ...
                                Params.intn_LTCafterDiag_PctCostatReachLvl, ...
                                numEligLTCAfter(:,HIVStageNum), ...
                                baseLinkageAfterRate(:,HIVStageNum), ...
                                NumGetIntnAtT1andReachLvlRate_LTCafterDiag(:,:,HIVStageNum), ...
                                InputType_LTCafterDiag, Indicator_IsPrEP, 0);
                    end
                    
                    Alloc_LTCafterDiag = (Params.intn_LTCafterDiag_Investment(allocationPeriod) * alloc_mult);    
                                        
                    LinkageAfterRate(Params.numStrats,Params.numHIVstages)=0;

                    %Calculate percentage of T1 rates and each reach level
                    %afforded by allocations to the LTC after diagnosis intervention
                    PctT1andReachLvlRatesAffordedByAlloc_LTCafterDiag = ...
                        min(CalcAddlPctT1andReachLvlRatesAffordedByAlloc(Alloc_LTCafterDiag, ...
                            CostT1andReachLvlRate_LTCafterDiag, Indicators_LimitPops, Indicator_SEP),1);
        
                    % Calculate final LTC after rate as sum of base rates and
                    % additional rates due to allocation. 
                    for HIVStageNum = 1:Params.numHIVstages
                         LinkageAfterRate(:,HIVStageNum) = baseLinkageAfterRate(:,HIVStageNum) + ...
                            sum(PctT1andReachLvlRatesAffordedByAlloc_LTCafterDiag(:,:) .* ...
                            (T1andReachLvlRate_LTCafterDiag(:,:,HIVStageNum) - ...
                            [zeros(Params.numStrats,1) T1andReachLvlRate_LTCafterDiag(:,1:Params.numReachLvls-1,HIVStageNum)]),2);
                    end

                    
                    % If targeting to YMSM only, add non-YMSM rates to
                    % final linkage rates
                    if and(Params.TargetIntnsToYMSM == 1, Params.tt_progressionSource == 3)
                                                
                        LinkageAfterRate(:,1) = LinkageAfterRate(:,1) + Params.tt_linkageAfterRate_Acute_5_nonYMSM;
                        LinkageAfterRate(:,2) = LinkageAfterRate(:,2) + Params.tt_linkageAfterRate_LatentA_5_nonYMSM;
                        LinkageAfterRate(:,3) = LinkageAfterRate(:,3) + Params.tt_linkageAfterRate_LatentB_5_nonYMSM;
                        LinkageAfterRate(:,4) = LinkageAfterRate(:,4) + Params.tt_linkageAfterRate_Late_5_nonYMSM;
                        LinkageAfterRate(:,5) = LinkageAfterRate(:,5) + Params.tt_linkageAfterRate_AIDS_5_nonYMSM;
                        
                    end
                    
                    % Collect $$$ needed for each reach level based on
                    % CostT1andReachLvlRate_LTCafterDiag - this is used to determine the
                    % starting funding for the cost minimization optimization
                    TTProg.CostT1andReachLvlRate_LTCafterDiag = sum(sum(CostT1andReachLvlRate_LTCafterDiag,3),1);
               end
               
               % INTERVENTION: ART INITIATION
               for sectionIntnARTInit = 1:1

                % ART Initiation
                
                    InputType_ARTInit = 1; % Input is rate
                    
                    % Number eligible for ART initiation
                    ARTInitEligComparts = [Params.B3, Params.C3, Params.D3, Params.E3, Params.F3];
                    numEligARTInit(1:Params.numStrats,1:Params.numHIVstages) = 0;
                    for HIVStageNum = 1:Params.numHIVstages
                        numEligARTInit(:,HIVStageNum) = ...
                            Compartments(ARTInitEligComparts(HIVStageNum),:)' ...
                            * TTProg.ARTElig(HIVStageNum) ...
                            .* Indicators_LimitPops;
                    end
                    
                   % Define baseline variables 
                       
                    % ART init rates without allocation
                    baseARTInitRate(:,1) = Params.tt_ARTInitRateAcute_5;
                    baseARTInitRate(:,2) = Params.tt_ARTInitRateLatentA_5;
                    baseARTInitRate(:,3) =  Params.tt_ARTInitRateLatentB_5;
                    baseARTInitRate(:,4) =  Params.tt_ARTInitRateLate_5;
                    baseARTInitRate(:,5) = Params.tt_ARTInitRateAIDS_5;

                    % ART init rates in T1
                    T1ARTInitRate(:,1) = Params.tt_ARTInitRateAcute_1;
                    T1ARTInitRate(:,2) = Params.tt_ARTInitRateLatentA_1;
                    T1ARTInitRate(:,3) =  Params.tt_ARTInitRateLatentB_1;
                    T1ARTInitRate(:,4) =  Params.tt_ARTInitRateLate_1;
                    T1ARTInitRate(:,5) = Params.tt_ARTInitRateAIDS_1;
                                        
                    NumGetIntnAtT1andReachLvlRate_ARTInit(1:Params.numStrats,1:Params.numReachLvls,1:Params.numHIVstages)=0;
                    T1andReachLvlRate_ARTInit(1:Params.numStrats,1:Params.numReachLvls,1:Params.numHIVstages)=0;
                    CostT1andReachLvlRate_ARTInit(1:Params.numStrats,1:Params.numReachLvls,1:Params.numHIVstages)=0;
                    
                    % Calculate T1 and reach level rates and costs of each
                    % rate for each disease stage
                    for HIVStageNum = 1:Params.numHIVstages
                        % Calculate number getting intervention at T1 and
                        % reach level rates
                        NumGetIntnAtT1andReachLvlRate_ARTInit(:,:,HIVStageNum) = ...
                            CalcNumIntnAtT1andReachLvlRates(numEligARTInit(:,HIVStageNum),  ...
                                T1ARTInitRate(:,HIVStageNum), ...
                                Params.intn_ARTInitiation_ReachLvls, ...
                                Params.intn_ARTInitiation_MaxReach, ...
                                Indicators_LimitPops, InputType_ARTInit);
                        % Calculate T1 and reach level rates        
                        T1andReachLvlRate_ARTInit(:,:,HIVStageNum) = ...
                            CalcT1andReachLvlIntnsPP(numEligARTInit(:,HIVStageNum),  ...
                                NumGetIntnAtT1andReachLvlRate_ARTInit(:,:,HIVStageNum), Indicators_LimitPops, Indicator_SEP);
                        % Calculate costs of achieving T1 and reach level rates    
                        CostT1andReachLvlRate_ARTInit(:,:,HIVStageNum) = ...
                            CalcCostT1andReachLvl(Params.costPP_ARTInitiation, ...
                                Params.intn_ARTInitiation_PctCostatReachLvl, ...
                                numEligARTInit(:,HIVStageNum), ...
                                baseARTInitRate(:,HIVStageNum), ...
                                NumGetIntnAtT1andReachLvlRate_ARTInit(:,:,HIVStageNum), ...
                                InputType_ARTInit, Indicator_IsPrEP, 0);
                    end
                                                                                
                    
                    Alloc_ARTInit = (Params.intn_ARTInitiation_Investment(allocationPeriod) * alloc_mult);                                           

                    %Calculate percentage of T1 rates and each reach level
                    %afforded by allocations to the ART initiation intervention
                    PctT1andReachLvlRatesAffordedByAlloc_ARTInit = ...
                        min(CalcAddlPctT1andReachLvlRatesAffordedByAlloc(Alloc_ARTInit, ...
                            CostT1andReachLvlRate_ARTInit, Indicators_LimitPops, Indicator_SEP),1);
        
                    % Calculate final ART initiation rate as sum of base rates and
                    % additional rates due to allocation. Additional rates
                    % due to allocation are 
                    for HIVStageNum = 1:Params.numHIVstages
                         InitiateARTRate(:,HIVStageNum) = baseARTInitRate(:,HIVStageNum) + ...
                            sum(PctT1andReachLvlRatesAffordedByAlloc_ARTInit(:,:) .* ...
                            (T1andReachLvlRate_ARTInit(:,:,HIVStageNum) - ...
                            [zeros(Params.numStrats,1) T1andReachLvlRate_ARTInit(:,1:Params.numReachLvls-1,HIVStageNum)]),2);
                    end
                                                                          
                    % If targeting to YMSM only, add non-YMSM rates to
                    % final ART initiation rates
                    if and(Params.TargetIntnsToYMSM == 1, Params.tt_progressionSource == 3)
                        
                        InitiateARTRate(:,1) = InitiateARTRate(:,1) + Params.tt_ARTInitRateAcute_5_nonYMSM;
                        InitiateARTRate(:,2) = InitiateARTRate(:,2) + Params.tt_ARTInitRateLatentA_5_nonYMSM;
                        InitiateARTRate(:,3) = InitiateARTRate(:,3) + Params.tt_ARTInitRateLatentB_5_nonYMSM;
                        InitiateARTRate(:,4) = InitiateARTRate(:,4) + Params.tt_ARTInitRateLate_5_nonYMSM;
                        InitiateARTRate(:,5) = InitiateARTRate(:,5) + Params.tt_ARTInitRateAIDS_5_nonYMSM;
                        
                    end
             
                % Number of people not affected by ART and eligible
                    TTProg.numberInCareEligForART = ...
                          Compartments(Params.B3,:)' * TTProg.ARTElig(Params.stage_Acute)...
                        + Compartments(Params.C3,:)' * TTProg.ARTElig(Params.stage_LatentA)...
                        + Compartments(Params.D3,:)' * TTProg.ARTElig(Params.stage_LatentB)...
                        + Compartments(Params.E3,:)' * TTProg.ARTElig(Params.stage_Late)...
                        + Compartments(Params.F3,:)' * TTProg.ARTElig(Params.stage_AIDS);
                    
               % Collect $$$ needed for each reach level based on
               % CostT1andReachLvlRate_ARTInit - this is used to determine the
               % starting funding for the cost minimization optimization
               TTProg.CostT1andReachLvlRate_ARTInit = sum(sum(CostT1andReachLvlRate_ARTInit,3),1);
                
               end

               % INTERVENTION: ART ADHERENCE: BECOME VLS
               for sectionARTAdhere_BecomeVLS = 1:1

                    % Number eligible for intervention
                    numberEligARTAdh4to5 = (sum(Compartments(Comparts_CCStage4,:)))';
                     
                    % Rates without intervention
                    ANVToVLSRate_Base = ANVToVLSRate;
                    
                    
                    % Calculate total probability based on base prob and allocation
                    InputType_ARTAdhere_BecomeVLS = 1; % input is rate
                    
                    % Calculate number getting intervention at T1 and reach level rates
                    NumGetIntnAtT1andReachLvlRates_ARTAdhere_BecomeVLS = ...
                        CalcNumIntnAtT1andReachLvlRates(numberEligARTAdh4to5, ...
                            Params.tt_ARTInitRateFromANV_1, Params.intn_TxAdherence_ReachLvls, Params.intn_TxAdherence_MaxReach, Indicators_LimitPops, InputType_ARTAdhere_BecomeVLS);
                    % Calculate T1 and reach level rates
                    T1andReachLvlRates_ARTAdhere_BecomeVLS = ...
                        CalcT1andReachLvlIntnsPP(numberEligARTAdh4to5, ...
                            NumGetIntnAtT1andReachLvlRates_ARTAdhere_BecomeVLS, Indicators_LimitPops, Indicator_SEP);
                    % Calculate costs of achieving T1 and reach level rates    
                    CostT1andReachLvlRate_ARTAdhere_BecomeVLS = ...
                        CalcCostT1andReachLvl(Params.costPP_TxAdherence_BecomeVLS, ...
                            Params.intn_TxAdherence_PctCostatReachLvl, ...
                            numberEligARTAdh4to5, ANVToVLSRate_Base, ...
                            NumGetIntnAtT1andReachLvlRates_ARTAdhere_BecomeVLS, ...
                            InputType_ARTAdhere_BecomeVLS, Indicator_IsPrEP, 0);
                    
                    Alloc_ARTAdhere_BecomeVLS = (Params.intn_ARTAdher4to5_Investment(allocationPeriod) * alloc_mult);    
                                        
                    %Calculate percentage of T1 rates and each reach level
                    %afforded by allocations to the ART adherence to become VLS intervention
                    PctT1andReachLvlRateAffordedByAlloc_ARTAdhere_BecomeVLS = ...
                        min(CalcAddlPctT1andReachLvlRatesAffordedByAlloc(Alloc_ARTAdhere_BecomeVLS, ...
                        CostT1andReachLvlRate_ARTAdhere_BecomeVLS, Indicators_LimitPops, Indicator_SEP),1);
                                       
                    % Calculate final ANV to VLS rate as sum of base rates and
                    % additional rates due to allocation.
                    ANVToVLSRate = ANVToVLSRate_Base + ...
                            sum(PctT1andReachLvlRateAffordedByAlloc_ARTAdhere_BecomeVLS .* ...
                            (T1andReachLvlRates_ARTAdhere_BecomeVLS - ...
                            [zeros(Params.numStrats,1) T1andReachLvlRates_ARTAdhere_BecomeVLS(:,1:Params.numReachLvls-1)]),2);
                                                           
                    % If targeting to YMSM only, add non-YMSM rates to
                    % final become VLS rates
                    if and(Params.TargetIntnsToYMSM == 1, Params.tt_progressionSource == 3)
                        
                        ANVToVLSRate = ANVToVLSRate + Params.tt_ARTInitRateFromANV_5_nonYMSM;
                        
                    end
                    
                    % Collect $$$ needed for each reach level based on
                    % CostT1andReachLvlRate_ARTAdhere_BecomeVLS - this is used to determine the
                    % starting funding for the cost minimization optimization
                    TTProg.CostT1andReachLvlRate_ARTAdhere_BecomeVLS = sum(CostT1andReachLvlRate_ARTAdhere_BecomeVLS,1);
                    
               end
               
               
               % INTERVENTION: ART ADHERENCE: REMAIN VLS
               for sectionARTAdhere_RemainVLS = 1:1
                   
                    % Number eligible for intervention
                     numberEligARTAdh5to4 = (sum(Compartments(Comparts_CCStage5,:)))';

                    % Number losing VLS and transitioning to ANV based on period 1 rates
                    probDropOutVLS_p1 = RateToProb(Params.tt_dropOutRate_VLSToANV_1);
                    probRemainVLS_p1 = 1 - probDropOutVLS_p1;
                    probRemainVLS_ReachLvls = Params.intn_TxAdherence_ReachLvls;
                     
                    % Rates without intervention
                    DropOutRate_VLSToANV_Base = DropOutRate_VLSToANV;
                    pctDropOutVLS_Base = min(RateToProb(DropOutRate_VLSToANV_Base),0.9999);
                    pctRemainVLS_Base = 1 - pctDropOutVLS_Base;
                  
                    % Calculate total probability based on base prob and allocation
                    InputType_ARTAdhere_RemainVLS = 2; % input is prob 
                    
                    % Calculate number getting intervention at T1 and reach
                    % level probs
                    NumGetIntnAtT1andReachLvlRates_ARTAdhere_RemainVLS = ...
                        CalcNumIntnAtT1andReachLvlRates(numberEligARTAdh5to4, ...
                            probRemainVLS_p1, probRemainVLS_ReachLvls, Params.intn_TxAdherence_MaxReach, Indicators_LimitPops, InputType_ARTAdhere_RemainVLS);
                    % Calculate T1 and reach level probs
                    T1andReachLvlProbs_ARTAdhere_RemainVLS = ...
                        CalcT1andReachLvlIntnsPP(numberEligARTAdh5to4, ...
                            NumGetIntnAtT1andReachLvlRates_ARTAdhere_RemainVLS, Indicators_LimitPops, Indicator_SEP);
                    % Calculate costs of achieving T1 and reach level probs
                    CostT1andReachLvlRate_ARTAdhere_RemainVLS = ...
                        CalcCostT1andReachLvl(Params.costPP_TxAdherence_RemainVLS, ...
                            Params.intn_TxAdherence_PctCostatReachLvl,...
                            numberEligARTAdh5to4, pctRemainVLS_Base, ...
                            NumGetIntnAtT1andReachLvlRates_ARTAdhere_RemainVLS, ...
                            InputType_ARTAdhere_RemainVLS, Indicator_IsPrEP, 0);
                    
                    Alloc_ARTAdhere_RemainVLS = (Params.intn_ARTAdher5to4_Investment(allocationPeriod) * alloc_mult);
                    
                    %Calculate percentage of T1 probs and each reach level
                    %afforded by allocations to the ART adherence to remain VLS intervention
                    PctT1andReachLvlRateAffordedByAlloc_ARTAdhere_RemainVLS = ...
                        min(CalcAddlPctT1andReachLvlRatesAffordedByAlloc(Alloc_ARTAdhere_RemainVLS, ...
                        CostT1andReachLvlRate_ARTAdhere_RemainVLS, Indicators_LimitPops, Indicator_SEP),1);                                     
                    
                    % Calculate final percent remaining VLS as sum of base probs and
                    % additional probs due to allocation.
                    PctRemainVLS = pctRemainVLS_Base + ...
                            sum(PctT1andReachLvlRateAffordedByAlloc_ARTAdhere_RemainVLS .* ...
                            (T1andReachLvlProbs_ARTAdhere_RemainVLS - ...
                            [zeros(Params.numStrats,1) T1andReachLvlProbs_ARTAdhere_RemainVLS(:,1:Params.numReachLvls-1)]),2);
                                        
                    % Create 1's for non-eligible pops to prevent infinite
                    % rates - check if this is still needed
                    PctRemainVLS = min(PctRemainVLS + (1- Indicators_LimitPops),1);
                    DropOutRate_VLSToANV = ProbToRate(1 - PctRemainVLS);                    
                   
                    % If targeting to YMSM only, add non-YMSM rates to
                    % final dropout rates
                    if and(Params.TargetIntnsToYMSM == 1, Params.tt_progressionSource == 3)
                        
                        DropOutRate_VLSToANV = DropOutRate_VLSToANV + Params.tt_dropOutRate_VLSToANV_5_nonYMSM;
                        
                    end
                    
                    % Collect $$$ needed for each reach level based on
                    % CostT1andReachLvlRate_ARTAdhere_RemainVLS - this is used to determine the
                    % starting funding for the cost minimization optimization
                    TTProg.CostT1andReachLvlRate_ARTAdhere_RemainVLS = sum(CostT1andReachLvlRate_ARTAdhere_RemainVLS,1);

                end
                   
               % INTERVENTION: SYRINGE EXCHANGE PROGRAM (SEP)
               for sectionSyringeExchange = 1:1

                   PctActivePWIDServedbySEP(1:Params.numRace,1) = 0;
                   TTProg.CostT1andReachLvlPctServed_SEP(1:Params.numRace,1:Params.numReachLvls) = 0;
                   
                   % # non-YMSM served by SSP for appropriate set period
                   % (if needed)
                   if Params.TargetIntnsToYMSM == 1
                        if Year < Params.tt_SEP_YrSet1Begins
                            NumPWIDServedbySEP_NonYMSM(1:Params.numRace,1) = 0; 
                        elseif Year < Params.tt_SEP_YrSet2Begins
                            NumPWIDServedbySEP_NonYMSM(1:Params.numRace,1) = Params.tt_numPWIDServedbySEP_Set1_nonYMSM;        
                        elseif Year < Params.tt_SEP_YrSet3Begins
                            NumPWIDServedbySEP_NonYMSM(1:Params.numRace,1) = Params.tt_numPWIDServedbySEP_Set2_nonYMSM;
                        elseif Year < Params.tt_SEP_YrSet4Begins
                            NumPWIDServedbySEP_NonYMSM(1:Params.numRace,1) = Params.tt_numPWIDServedbySEP_Set3_nonYMSM;
                        elseif Year < Params.tt_SEP_YrSet5Begins
                            NumPWIDServedbySEP_NonYMSM(1:Params.numRace,1) = Params.tt_numPWIDServedbySEP_Set4_nonYMSM;
                        else
                            NumPWIDServedbySEP_NonYMSM(1:Params.numRace,1) = Params.tt_numPWIDServedbySEP_Set5_nonYMSM;
                        end
                   end
                   
                   % Loop through each race/ethnicity
                   for nRace = 1:Params.numRace
                        NumPWIDServedbySEP_Base = NumPWIDServedbySEP(nRace,1);                    
    
                        % Number eligible for intervention
                        numberElig_SEP = (sum(Compartments(Params.UninfectedComparts,:)))' .* Params.raceIndicator(:,nRace) .* Params.popIndicator(:,Params.pop_IDU) * ...
                               Params.behav_pctActivePWID;
                        
                        PctActivePWIDServedbySEP_Base = NumPWIDServedbySEP_Base ./ max(sum(numberElig_SEP),1);
                                              
                        InputType_SEP = 2; % input is % reached
                        Indicator_SEP = 1; % indicator to let function 8b know not to stratify by subpop
                        
                        % Apply only if not targeting interventions only to MSM
                        if Params.TargetIntnsToYMSM == 0
                            
                            % If first year of time period 5, collect total number
                            % eligible for SEP (on SEP or not) for calculating T1 reach 
                            % PctActivePWIDServedbySEP(nRace,1) = min(NumPWIDServedbySEP(nRace,1) / ...
                            % max((sum(sum(Compartments(Params.UninfectedComparts,:))' .* Params.raceIndicator(:,nRace) .* Params.popIndicator(:,Params.pop_IDU)) * Params.behav_pctActivePWID),1), 1);
                            if Year == Params.tt_periodFiveStartYear 
                                TTProg.PctOnSEP_T1 = TTProg.NumOnSEP_T1(nRace,1) ./  ...
                                    max(sum(numberElig_SEP),0.000000000000000000001);
                            end
                        
                            % Calculate number getting intervention at T1 and reach
                            % level % reached
                            NumGetIntnAtT1andReachLvlPctServed_SEP = ...
                                CalcNumIntnAtT1andReachLvlRates(sum(numberElig_SEP), ...
                                    TTProg.PctOnSEP_T1, Params.intn_SEP_ReachLvls, Params.intn_SEP_MaxReach, 1, InputType_SEP);
                            % Calculate T1 and reach level % reached
                            T1andReachLvlPctServed_SEP = ...
                                CalcT1andReachLvlIntnsPP(sum(numberElig_SEP), ...
                                    NumGetIntnAtT1andReachLvlPctServed_SEP, 1, Indicator_SEP);
                            % Calculate costs of achieving T1 and reach level %
                            % reached
                            CostT1andReachLvlPctServed_SEP = ...
                                CalcCostT1andReachLvl(Params.costPP_SyringeExchange, ...
                                    Params.intn_SEP_PctCostatReachLvl,...
                                    sum(numberElig_SEP), NumPWIDServedbySEP_Base, ...
                                    NumGetIntnAtT1andReachLvlPctServed_SEP, ...
                                    InputType_SEP, Indicator_IsPrEP, 0);   
                            
                            % Collect allocation for correct race
                            switch nRace
                                case 1
                                    Alloc_SEP = Params.intn_SEP_Investment_B(allocationPeriod) * alloc_mult; 
                                case 2
                                    Alloc_SEP = Params.intn_SEP_Investment_H(allocationPeriod) * alloc_mult; 
                                case 3
                                    Alloc_SEP = Params.intn_SEP_Investment_O(allocationPeriod) * alloc_mult; 
                            end
    
                            %Calculate percentage of T1 % and each reach level
                            %afforded by allocations to the SEP intervention
                            PctT1andReachLvlPctServedAffordedByAlloc_SEP = ...
                                min(CalcAddlPctT1andReachLvlRatesAffordedByAlloc(Alloc_SEP, ...
                                CostT1andReachLvlPctServed_SEP, 1, Indicator_SEP),1);
    
                            % Calculate final percent reached with SEP as sum of base probs and
                            % additional probs due to allocation.
                            PctActivePWIDServedbySEP(nRace,1) = PctActivePWIDServedbySEP_Base + ...
                                    sum(PctT1andReachLvlPctServedAffordedByAlloc_SEP .* ...
                                    (T1andReachLvlPctServed_SEP - ...
                                    [zeros(1,1) T1andReachLvlPctServed_SEP(:,1:Params.numReachLvls-1)]),2);
    
                            % Collect $$$ needed for each reach level based on
                            % CostT1andReachLvlPctServed_SEP - this is used to determine the
                            % starting funding for the cost minimization optimization                           
                            TTProg.CostT1andReachLvlPctServed_SEP(nRace,:) = sum(CostT1andReachLvlPctServed_SEP,1);                                    
                                
                        else                                                
                            PctActivePWIDServedbySEP_noMaxReach = (NumPWIDServedbySEP_Base + NumPWIDServedbySEP_NonYMSM(nRace,1)) / ...
                                max(1,sum(numberElig_SEP));
    
                            PctActivePWIDServedbySEP(nRace,1) = min(PctActivePWIDServedbySEP_noMaxReach, Params.intn_SEP_MaxReach);
                        end                         
                                                 
                   end
                       
                    % Apply pct served by race to all 273 subpopulations
                    TTProg.PctActivePWIDServedbySEP = (Params.raceIndicator * PctActivePWIDServedbySEP) .* ...
                        Params.popIndicator(:,Params.pop_IDU);

                    % Total number eligible for intervention
                    numberElig_SEP_allraces = (sum(Compartments(Params.UninfectedComparts,:)))' .* ...
                        Params.popIndicator(:,Params.pop_IDU) * ...
                          Params.behav_pctActivePWID;

                    % Calculate # served by subpopulation
                    TTProg.NumPWIDServedbySEP = ((TTProg.PctActivePWIDServedbySEP .* numberElig_SEP_allraces)' * Params.raceIndicator)';                         
                       

                    Indicator_SEP = 0; % indicator to let function 8b for next interventions know to stratify by subpop (not stratified for SEP only)
                    
               end
                
               % INTERVENTION: PrEP INITIATION 
               for sectionPrEPIntn = 1:1
                   
                    % Pull base PrEP init rate or % on PrEP w/o funding for appropriate year
                    %  In years when DE progression used: Annual rate of initiating PrEP per eligible person
                    %  In years when AB progression used: Percentage of eligible people on PrEP, given no funding
                    PctOnPrEP_Oral_base = PrEPInitiationRateOrPctOnPrEP .* (1 - PctInjectPrEP);
                    PctOnPrEP_Inject_base = PrEPInitiationRateOrPctOnPrEP .* PctInjectPrEP;
                    
                    % numA1_PrEPElig: Everyone who is PrEP-eligible 
                    % (limited by % eligible) and NOT on PrEP already 
                    numA1_PrEPElig = sum(Compartments(Params.A1,:),1)' .* ...
                        Indicators_PrEPPops  .* Indicators_LimitPops .* Params.tt_prepPctEligible;
                    % numA1andPrEP_PrEPElig: Everyone who is PrEP-eligible 
                    % (limited by % eligible), regardless of PrEP status  
                    numA1andPrEP_PrEPElig = (sum(Compartments([Params.A1, Params.PrEPComparts],:),1)' .* ...
                        Indicators_PrEPPops .* Indicators_LimitPops) .* Params.tt_prepPctEligible;
                                                                      
                                       
                    InputType_PrEP = 2; % allocation covers % of pop, so assuming a %
                    Indicator_IsPrEP = 1; % indicator to let function 8c know to apply different intervention costs PP for those on PrEP and those initating PrEP
                                
                    % numA1andPrEPPops: 
                    % PrEP-eligible subpops: uninfected not on PrEP and people on oral PrEP or injectable PrEP 
                    
                    % If first year of time period 5, collect total number
                    % eligible for PrEP (on PrEP or not), on oral PrEP, and on injectable PrEP
                    % for calculating T1 reach 
                    if Year == Params.tt_periodFiveStartYear 
                        numA1andPrEP_PrEPEligT1 = numA1andPrEP_PrEPElig;
                        numOnOralPrEPT1 = sum(Compartments(Params.OralPrEPComparts,:),1)';
                        numOnInjectPrEPT1 = sum(Compartments(Params.InjectPrEPComparts,:),1)';
                        TTProg.PctOnOralPrEPT1 = numOnOralPrEPT1 ./  max(numA1andPrEP_PrEPEligT1,0.000000000000000000001);
                        TTProg.PctOnInjectPrEPT1 = numOnInjectPrEPT1 ./  max(numA1andPrEP_PrEPEligT1,0.000000000000000000001);                       
                    end
                    
                    Alloc_PrEP(1:Params.numIntns_PrEPType,1:Params.numPop_plusHETbySex,1:Params.numRace) = 0;
                    % Calculate initiation rate among A1 based on only a portion of the
                    % high-risk population eligible for PrEP. This
                    % method factors percentage already on PrEP into
                    % eligibility among those not on PrEP. Also considers
                    % max reach. Initiation and dropout rates are calculated
                    % each year (previously only the first year of each
                    % allocation period)                        
                    for PrEPType = 1:2 % 1 = Oral; 2 = Injectable

                        switch PrEPType
                            case 1 % Oral
                                PctOnPrEP_base = PctOnPrEP_Oral_base;
                                % numOnPrEP_byPrEPType: Everyone who is currently on oral PrEP 
                                numOnPrEP_byPrEPType = sum(Compartments(Params.OralPrEPComparts,:),1)';
                                                               
                                PctonPrEPT1 = TTProg.PctOnOralPrEPT1;

                                intn_PrEP_ReachLvls = Params.intn_PrEP_Oral_ReachLvls;
                                intn_PrEP_PctCostatReachLvl = Params.intn_PrEP_Oral_PctCostatReachLvl; 
                                intn_PrEP_MaxReach = Params.intn_PrEP_Oral_MaxReach;
                                % Annual PrEP costs per person
                                % initiating PrEP, calculated as a
                                % weighted average based on the % of
                                % oral PrEP initiators with high and
                                % low adherence
                                AnnualPrEPCostsPP_Initiate = ((Params.PctHighAdherence_OralPrEP * Params.hiv_annualCostPerCompart_5(Params.A6)) ...
                                + ((1 - Params.PctHighAdherence_OralPrEP) * Params.hiv_annualCostPerCompart_5(Params.A7))) .* Indicators_PrEPPops;
                                % Annual PrEP costs per person
                                % currently on PrEP, calculated as a
                                % weighted average based on the number
                                % of people on high adherence and low
                                % adherence oral PrEP
                                AnnualPrEPCostsPP_OnPrEP = AnnualPrEPCostsPP_Initiate;
                                for subpop = 1:Params.numStrats
                                    if sum(Compartments(Params.OralPrEPComparts,subpop))>0
                                        AnnualPrEPCostsPP_OnPrEP(subpop,1) = ((Compartments(Params.A6,subpop) * Params.hiv_annualCostPerCompart_5(Params.A6)) ...
                                        + (Compartments(Params.A7,subpop) * Params.hiv_annualCostPerCompart_5(Params.A7))) ...
                                        /  sum(Compartments(Params.OralPrEPComparts,subpop));                                  
                                    end
                                end
                                
                                % Alloc_PrEP(# PrEP type, # transmission
                                % group with HET by sex, # race)
                                Alloc_PrEP(PrEPType,1,1) = Params.intn_PrEP_Oral_Investment_HETM_B(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,1,2) = Params.intn_PrEP_Oral_Investment_HETM_H(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,1,3) = Params.intn_PrEP_Oral_Investment_HETM_O(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,2,1) = Params.intn_PrEP_Oral_Investment_HETF_B(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,2,2) = Params.intn_PrEP_Oral_Investment_HETF_H(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,2,3) = Params.intn_PrEP_Oral_Investment_HETF_O(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,3,1) = Params.intn_PrEP_Oral_Investment_MSM_B(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,3,2) = Params.intn_PrEP_Oral_Investment_MSM_H(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,3,3) = Params.intn_PrEP_Oral_Investment_MSM_O(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,4,1) = Params.intn_PrEP_Oral_Investment_IDU_B(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,4,2) = Params.intn_PrEP_Oral_Investment_IDU_H(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,4,3) = Params.intn_PrEP_Oral_Investment_IDU_O(allocationPeriod) * alloc_mult;
                            case 2 % Injectable
                                PctOnPrEP_base = PctOnPrEP_Inject_base;
                                % numOnPrEP: Everyone who currently on PrEP 
                                numOnPrEP_byPrEPType = sum(Compartments(Params.InjectPrEPComparts,:),1)';                                
                                
                                PctonPrEPT1 = TTProg.PctOnInjectPrEPT1;

                                intn_PrEP_ReachLvls = Params.intn_PrEP_Inject_ReachLvls;
                                intn_PrEP_PctCostatReachLvl = Params.intn_PrEP_Inject_PctCostatReachLvl; 
                                intn_PrEP_MaxReach = Params.intn_PrEP_Inject_MaxReach;
                                % Annual PrEP costs per person
                                % initiating PrEP, calculated as a
                                % weighted average based on the % of
                                % oral PrEP initiators with high and
                                % low adherence
                                AnnualPrEPCostsPP_Initiate = ((Params.PctHighAdherence_InjectPrEP * Params.hiv_annualCostPerCompart_5(Params.A8)) ...
                                + ((1 - Params.PctHighAdherence_InjectPrEP) * Params.hiv_annualCostPerCompart_5(Params.A9))) .* Indicators_PrEPPops;                                  
                                % Annual PrEP costs per person
                                % currently on PrEP, calculated as a
                                % weighted average based on the number
                                % of people on high adherence and low
                                % adherence oral PrEP
                                AnnualPrEPCostsPP_OnPrEP = AnnualPrEPCostsPP_Initiate;
                                for subpop = 1:Params.numStrats
                                    if sum(Compartments(Params.InjectPrEPComparts,subpop))>0
                                        AnnualPrEPCostsPP_OnPrEP(subpop,1) = ((Compartments(Params.A8,subpop) * Params.hiv_annualCostPerCompart_5(Params.A8)) ...
                                        + (Compartments(Params.A9,subpop) * Params.hiv_annualCostPerCompart_5(Params.A9))) ...
                                        /  sum(Compartments(Params.InjectPrEPComparts,subpop));
    %                                     AnnualPrEPCostsPP_OnPrEP = ((Compartments(Params.A8,:)' * Params.hiv_annualCostPerCompart_5(Params.A8)) ...
    %                                     + (Compartments(Params.A9,:)' * Params.hiv_annualCostPerCompart_5(Params.A9))) ...
    %                                     ./  max(sum(Compartments(Params.InjectPrEPComparts,:),1)',0.000000000000000000001) .* Indicators_PrEPPops;                                   
                                    end
                                end    

                                % Alloc_PrEP(# PrEP type, # transmission
                                % group with HET by sex, # race)
                                Alloc_PrEP(PrEPType,1,1) = Params.intn_PrEP_Inject_Investment_HETM_B(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,1,2) = Params.intn_PrEP_Inject_Investment_HETM_H(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,1,3) = Params.intn_PrEP_Inject_Investment_HETM_O(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,2,1) = Params.intn_PrEP_Inject_Investment_HETF_B(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,2,2) = Params.intn_PrEP_Inject_Investment_HETF_H(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,2,3) = Params.intn_PrEP_Inject_Investment_HETF_O(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,3,1) = Params.intn_PrEP_Inject_Investment_MSM_B(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,3,2) = Params.intn_PrEP_Inject_Investment_MSM_H(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,3,3) = Params.intn_PrEP_Inject_Investment_MSM_O(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,4,1) = Params.intn_PrEP_Inject_Investment_IDU_B(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,4,2) = Params.intn_PrEP_Inject_Investment_IDU_H(allocationPeriod) * alloc_mult;
                                Alloc_PrEP(PrEPType,4,3) = Params.intn_PrEP_Inject_Investment_IDU_O(allocationPeriod) * alloc_mult;
                        end
                        
                        % Calculate number getting intervention at current
                        % and reach level % on PrEP
                        NumGetIntnAtT1andReachLvlPctOnPrEP = ...
                            CalcNumIntnAtT1andReachLvlRates(numA1andPrEP_PrEPElig, ...
                                PctonPrEPT1, intn_PrEP_ReachLvls, intn_PrEP_MaxReach, Indicators_LimitPops, InputType_PrEP); 
                        % Calculate current and reach level % on PrEP
                        T1andReachLvlPctOnPrEP = ...
                            CalcT1andReachLvlIntnsPP(numA1andPrEP_PrEPElig, ...
                                NumGetIntnAtT1andReachLvlPctOnPrEP, Indicators_LimitPops, Indicator_SEP);


                        % Calculate costs of achieving current and reach level % on PrEP    
                        CostT1andReachLvlPctOnPrEP = ...
                            CalcCostT1andReachLvl(AnnualPrEPCostsPP_OnPrEP, ...
                                intn_PrEP_PctCostatReachLvl, ...
                                numA1andPrEP_PrEPElig, PctOnPrEP_base, ...
                                NumGetIntnAtT1andReachLvlPctOnPrEP, ...
                                InputType_PrEP, Indicator_IsPrEP, AnnualPrEPCostsPP_Initiate);

                        PctT1andReachLvlPctOnPrEPAffordedByAlloc_PrEP = ...
                            zeros(Params.numStrats,Params.numReachLvls);

                        %Calculate percentage of T1 % on PrEP and each
                        %reach levels
                        %afforded by allocations to the oral and injectable PrEP interventions

                        TTProg.CostT1andReachLvlRate_PrEP(1:Params.numIntns_PrEPType,1:Params.numPop_plusHETbySex,1:Params.numRace,1:Params.numReachLvls) = 0;

                        for PrEPSubPopNum = 1:Params.numPop_plusHETbySex

                            for nRace = 1:Params.numRace

                                subPopIndicator = Params.popIndicator_withHETbySex(:,PrEPSubPopNum) .* Params.raceIndicator(:,nRace) .* Params.riskLevelIndicator(:,Params.risk_Casual);

                                % Calculate $$$ needed for each reach level 
                                % CostT1andReachLvlPctServed_PrEP - this is used to determine the
                                % starting funding for the cost minimization optimization
                                TTProg.CostT1andReachLvlRate_PrEP(PrEPType,PrEPSubPopNum,nRace,:) = sum(CostT1andReachLvlPctOnPrEP .* subPopIndicator,1);

                                PrEPInvestment = Alloc_PrEP(PrEPType,PrEPSubPopNum,nRace);
                                
                                    PctT1andReachLvlPctOnPrEPAffordedByAlloc_PrEP = ...
                                        min(PctT1andReachLvlPctOnPrEPAffordedByAlloc_PrEP + ...
                                            CalcAddlPctT1andReachLvlRatesAffordedByAlloc(PrEPInvestment, ...
                                            CostT1andReachLvlPctOnPrEP, subPopIndicator .* Indicators_LimitPops, Indicator_SEP),1);

                            end
                        end %PrEPSubPopNum

                        % Calculate total % on PrEP among eligible and reachable afforded as sum of base rates and
                        % additional % on PrEP due to allocation.                      
                        PctOnPrEP_Among_Elig_ReachableAfforded_Alloc = PctOnPrEP_base + ...
                                sum(PctT1andReachLvlPctOnPrEPAffordedByAlloc_PrEP .* ...
                                (T1andReachLvlPctOnPrEP - ...
                                [zeros(Params.numStrats,1) T1andReachLvlPctOnPrEP(:,1:Params.numReachLvls-1)]),2);                        

                        NumIntns_Alloc_PrEP = PctOnPrEP_Among_Elig_ReachableAfforded_Alloc .* numA1andPrEP_PrEPElig;    

                        % Each year, start over so that 
                        % everyone previously on PrEP drops off and 
                        % allocation restarts. Do that by calculating
                        % dropout and initiation rates to get net number on
                        % PrEP from allocation correct.
                        %if allocationPeriod == 1                           
                        %    PrEPInitiationRate_AmongElig = ...
                        %        Rate_Among_Elig_ReachableAfforded_Alloc;
                        %else
                        PrEPInitiationRate_AmongElig(1:Params.numStrats,1)=0;
                        DropOutRate_PrEP(1:Params.numStrats,1) = 0;
                        for i = 1:Params.numStrats
                            if NumIntns_Alloc_PrEP(i)<numOnPrEP_byPrEPType(i)
                                if numOnPrEP_byPrEPType(i)>0
                                    DropOutRate_PrEP(i) = ...
                                        DropOutRate_PrEP(i) + ...
                                        ProbToRate((numOnPrEP_byPrEPType(i) - NumIntns_Alloc_PrEP(i)) / numOnPrEP_byPrEPType(i));
%                                            (numOnPrEP(i) - NumIntns_Alloc_PrEP(i)) / numOnPrEP(i);
                                else
                                    DropOutRate_PrEP(i) = 0; 
                                end
                                PrEPInitiationRate_AmongElig(i) = 0;
                            else 
                                DropOutRate_PrEP(i) = 0;
                                if numA1_PrEPElig(i)>0
                                    PrEPInitiationRate_AmongElig(i) = ...
                                        ProbToRate((NumIntns_Alloc_PrEP(i) - numOnPrEP_byPrEPType(i)) / numA1_PrEPElig(i));
%                                         (NumIntns_Alloc_PrEP(i) - numOnPrEP_byPrEPType(i)) / numA1_PrEPElig(i); 
                                else
                                    PrEPInitiationRate_AmongElig(i) = 0;
                                end
                            end
                        end
                        %end


                        % Adjust overall initiation rate for % of pop eligible PrEP, as specified by user 
                        PrEPInitiationRate = ...
                            PrEPInitiationRate_AmongElig .* Params.tt_prepPctEligible;

                        % Collect initiation rates, dropout rates, and $$$ needed for each reach level 
                        % CostT1andReachLvlPctServed_PrEP - this is used to determine the
                        % starting funding for the cost minimization optimization
                        switch PrEPType
                            case 1                                 
                                PrEPInitiationRate_A6 = PrEPInitiationRate .* Params.PctHighAdherence_OralPrEP;
                                PrEPInitiationRate_A7 = PrEPInitiationRate .* (1 - Params.PctHighAdherence_OralPrEP);

                                DropOutRate_PrEP_Oral_HighAdherence = DropOutRate_PrEP;
                                DropOutRate_PrEP_Oral_LowAdherence = DropOutRate_PrEP;                                
                            case 2
                                PrEPInitiationRate_A8 = PrEPInitiationRate .* Params.PctHighAdherence_InjectPrEP;
                                PrEPInitiationRate_A9 = PrEPInitiationRate .* (1 - Params.PctHighAdherence_InjectPrEP);

                                DropOutRate_PrEP_Inject_HighAdherence = DropOutRate_PrEP;
                                DropOutRate_PrEP_Inject_LowAdherence = DropOutRate_PrEP;                                
                        end

                    end % PrEPType
                        
                        
                        
                    %else
                       
                        %PrEPInitiationRate_AmongElig = PrEPInitiationRate_base;
                        
                        % Collect $$$ needed for each reach level based on
                        % CostT1andReachLvlPctServed_PrEP - this is used to determine the
                        % starting funding for the cost minimization
                        % optimization (just collecting $0 here as it is
                        % not an initiation year)
                        %TTProg.CostT1andReachLvlRate_PrEP_HET = 0;
                        %TTProg.CostT1andReachLvlRate_PrEP_MSM = 0;
                        %TTProg.CostT1andReachLvlRate_PrEP_IDU = 0;
                        
                    %end
                    
                    
                    
                    % If targeting to YMSM only, add non-YMSM rates to
                    % final PrEP initiation rates 
                    % JC note (11/09/2021): still need to modify this once
                    % we confirm the above approach is working correctly.
                    if and(Params.TargetIntnsToYMSM == 1, Params.tt_progressionSource == 3)
                        
                        % Set PrEP initiation rates for all non-YMSM
                        if Year < Params.tt_PrEPInit_YrSet1Begins
                            PrEPInitiationRate_nonYMSM = 0;
                        elseif Year < Params.tt_PrEPInit_YrSet2Begins
                            PrEPInitiationRate_nonYMSM = Params.tt_PrEPInitRate_Set1_nonYMSM;
                        elseif Year < Params.tt_PrEPInit_YrSet3Begins
                            PrEPInitiationRate_nonYMSM = Params.tt_PrEPInitRate_Set2_nonYMSM;
                        elseif Year < Params.tt_PrEPInit_YrSet4Begins
                            PrEPInitiationRate_nonYMSM = Params.tt_PrEPInitRate_Set3_nonYMSM;
                        elseif Year < Params.tt_PrEPInit_YrSet5Begins
                            PrEPInitiationRate_nonYMSM = Params.tt_PrEPInitRate_Set4_nonYMSM;
                        else
                            PrEPInitiationRate_nonYMSM = Params.tt_PrEPInitRate_Set5_nonYMSM;
                        end
                        
                        % Adjust initiation rate so that only a portion of the
                        % high-risk population is eligible for PrEP. This
                        % method factors percentage already on PrEP into
                        % eligibility among those not on PrEP
                         TotalEligPrEP = sum(Compartments(Params.UninfectedComparts,:),1)'.*Params.tt_prepPctEligible;
                         EligNotOnPrEP = zeros(Params.numStrats,1);
                         for counter_strats = 1:Params.numStrats 
                            EligNotOnPrEP(counter_strats) = max((TotalEligPrEP(counter_strats) - Compartments(Params.PrEPComparts,counter_strats)'),0);
                         end
                         PctOfA1PrEPElig = EligNotOnPrEP./(Compartments(Params.A1,:)');
                         PrEPInitiationRate_nonYMSM = PrEPInitiationRate_nonYMSM.*PctOfA1PrEPElig;
                         
                         PrEPInitiationRate = PrEPInitiationRate + PrEPInitiationRate_nonYMSM;
                        
                    end

               end
                 
                end

                 
                
            else %not allocation-based (direct-entry or HET-testing-frequency in third time period)

                % 2.iii.2 Non-Allocation-Based Progression
                for sectionNonAlloc = 1:1
                                                            
              % Test Rates
                TestRate(:,1) = Params.tt_testRateUninfected_rp_5;
                TestRate(:,2) = Params.tt_testRateAcute_rp_5;
                TestRate(:,3) =  Params.tt_testRateLatentA_rp_5;
                TestRate(:,4) =  Params.tt_testRateLatentB_rp_5;
                TestRate(:,5) = Params.tt_testRateLate_rp_5;
                TestRate(:,6) =  Params.tt_testRateAIDS_rp_5;
                
                LinkageFirstProb = Params.tt_linkageFirst_r_5;

                LinkageAfterRate(:,1) = Params.tt_linkageAfterRate_Acute_5;
                LinkageAfterRate(:,2) = Params.tt_linkageAfterRate_LatentA_5; 
                LinkageAfterRate(:,3) = Params.tt_linkageAfterRate_LatentB_5;
                LinkageAfterRate(:,4) = Params.tt_linkageAfterRate_Late_5;
                LinkageAfterRate(:,5) = Params.tt_linkageAfterRate_AIDS_5;

                DropOutRate_CareToAware = Params.tt_dropOutRate_CareToAware_5;
                DropOutRate_ANVToAware = Params.tt_dropOutRate_ANVToAware_5;
                DropOutRate_ANVToCare = Params.tt_dropOutRate_ANVToCare_5;
                DropOutRate_VLSToANV = Params.tt_dropOutRate_VLSToANV_5;
                DropOutRate_VLSToLTC = Params.tt_dropOutRate_VLSToLTC_5;
                DropOutRate_VLSToAware = Params.tt_dropOutRate_VLSToAware_5;
                
                % ART initiation from LTC, no ART effects
                InitiateARTRate(:,1) = Params.tt_ARTInitRateAcute_5; 
                InitiateARTRate(:,2) = Params.tt_ARTInitRateLatentA_5; 
                InitiateARTRate(:,3) = Params.tt_ARTInitRateLatentB_5;
                InitiateARTRate(:,4) = Params.tt_ARTInitRateLate_5;
                InitiateARTRate(:,5) = Params.tt_ARTInitRateAIDS_5;
                
                % ART init from ANV
                ANVToVLSRate = Params.tt_ARTInitRateFromANV_5;
                
                PctWhoBecomeVLS = Params.tt_PctInitiateARTWhoBecomeVLS;
                
                TTProg.ARTElig = Params.tt_ARTElig_5; 
                
                % All-cause mortality
                pop_dyingRate = Params.pop_dyingRate;
                
                %Not Affected by ART Death
                NoARTEffectsDeathRate_NotAIDS = Params.rate_deathnoARTEffectsNotAIDS_2to5;

                % ART Death prob
                ANVDeathRate_NotAIDS = Params.rate_deathANVNotAIDS_2to5;
                ANVDeathRate_AIDS = Params.rate_deathANVAIDS_2to5;
                    
                %VLS Death
                VLSDeathRate_NotAIDS = Params.rate_deathVLSNotAIDS_2to5;
                VLSDeathRate_AIDS = Params.rate_deathVLSAIDS_2to5;

                PctANVWhoAreInCareInclART = Params.tt_PctANVWhoAreInCareOrART_5;
                
                % HET frequency-based method
                if Params.tt_progressionSource == 2
                   
                % Replace HET testing rates with those calc'ed from HET frequency inputs
                    % Code can apply intervals to both low
                    % and high risk HET (i.e., All-HET) or just HRH
                    % depending on user entry
                    
                    TestRate(:,2) = ...
                        (TestRate(:,2) .* (1-Params.popIndicator(:,Params.pop_HET))) + ...
                        Params.HRHIndicator .* Params.tt_HETFreq_TestRate_HRH + ...
                        Params.LowRiskHETIndicator .* (...
                            Params.tt_HETFreq_ApplyLRH * Params.tt_HETFreq_TestRate_LRH + ...
                            (1-Params.tt_HETFreq_ApplyLRH) * TestRate(:,2));
                                  
                    TestRate(:,3) = ...
                       (TestRate(:,3) .* (1-Params.popIndicator(:,Params.pop_HET))) + ...
                        Params.HRHIndicator .* Params.tt_HETFreq_TestRate_HRH + ...
                        Params.LowRiskHETIndicator .* (...
                            Params.tt_HETFreq_ApplyLRH * Params.tt_HETFreq_TestRate_LRH + ...
                            (1-Params.tt_HETFreq_ApplyLRH) * TestRate(:,3));
                        
                    TestRate(:,4) = ...
                        (TestRate(:,4) .* (1-Params.popIndicator(:,Params.pop_HET))) + ...
                        Params.HRHIndicator .* Params.tt_HETFreq_TestRate_HRH + ...
                        Params.LowRiskHETIndicator .* (...
                            Params.tt_HETFreq_ApplyLRH * Params.tt_HETFreq_TestRate_LRH + ...
                            (1-Params.tt_HETFreq_ApplyLRH) * TestRate(:,4));
                        
                    TestRate(:,5) = ...
                        (TestRate(:,5) .* (1-Params.popIndicator(:,Params.pop_HET))) + ...
                        Params.HRHIndicator .* Params.tt_HETFreq_TestRate_HRH + ...
                        Params.LowRiskHETIndicator .* (...
                            Params.tt_HETFreq_ApplyLRH * Params.tt_HETFreq_TestRate_LRH + ...
                            (1-Params.tt_HETFreq_ApplyLRH) * TestRate(:,5));
                        
                    TestRate_AIDS_DE = TestRate(:,6) ;
                        
                    TestRate_AIDS_Int = ...
                       (TestRate(:,6) .* (1-Params.popIndicator(:,Params.pop_HET))) + ...
                        Params.HRHIndicator .* Params.tt_HETFreq_TestRate_HRH + ...
                        Params.LowRiskHETIndicator .* (...
                            Params.tt_HETFreq_ApplyLRH * Params.tt_HETFreq_TestRate_LRH + ...
                            (1-Params.tt_HETFreq_ApplyLRH) * TestRate_AIDS);

                    TestRate_AIDS(1:Params.numStrats) = 0;
                    % Code for running TF with symptomatic testing for PLWH with AIDS
                    for findMaxIndex = 1:Params.numStrats
                        TestRate_AIDS(findMaxIndex)=max(TestRate_AIDS_Int(findMaxIndex),TestRate_AIDS_DE(findMaxIndex));
                    end                        
                        
                    TestRate(:,1) = ...
                         (TestRate(:,1) .* (1-Params.popIndicator(:,Params.pop_HET))) + ...
                        Params.HRHIndicator .* Params.tt_HETFreq_TestRate_HRH + ...
                        Params.LowRiskHETIndicator .* (...
                            Params.tt_HETFreq_ApplyLRH * Params.tt_HETFreq_TestRate_LRH + ...
                            (1-Params.tt_HETFreq_ApplyLRH) * TestRate(:,1));
                                        

                end
                end
            end
            end
        end
    end

%% 3. Record current values of key parameters 
    for sectionNonCDC = 1:1
                
        % Record test sensitivities for use in Collect Results
        TTProg.TestSens_Rapid_Acute = TestSens_Rapid_Acute;
        TTProg.TestSens_Conv_Acute = TestSens_Conv_Acute;
        TTProg.TestSens_Rapid_Chronic = TestSens_Rapid_Chronic;
        TTProg.TestSens_Conv_Chronic = TestSens_Conv_Chronic;
        TTProg.TestSens_Confirm_Acute = TestSens_Confirm_Acute;
        TTProg.TestSens_Confirm_Chronic = TestSens_Confirm_Chronic;
        
        % Record cost per test by result and type of test
            % for use in Collect results
            TTProg.costPP_TestNegRapid = costPP_TestNegRapid;
            TTProg.costPP_TestPosRapid = costPP_TestPosRapid;                   
            TTProg.costPP_TestNegConv = costPP_TestNegConv;
            TTProg.costPP_TestPosConv = costPP_TestPosConv;            
    
            % Non-allocation progression
        if Params.tt_progressionSource ~= 3 || Year < Params.tt_periodFiveStartYear
      
            % Probability of notification
           
                TTProg.NotifyPosConvProb = NotifyConvProb;
                TTProg.NotifyNegConvProb =  NotifyNegConvProb;
                TTProg.NotifyPosRapidProb =  NotifyRapidProb;
                TTProg.NotifyNegRapidProb =  NotifyNegRapidProb;

       end
    end
    
%% 4. Calculate transition rates defined at the model start  and 
    % do not change over the model's time horizon
    for sectionTPModelStart = 1:1
        if Year == Params.tt_modelStartYear % Only at model start
      
        % 4.i. HIV progression and death from AIDS
        for HIVProgressionSection = 1:1
            
                % Eligible: % All PLWH
                    % on ART/VLS and have AIDS
                        % This population gets the ART death rates
 
                % Calculate rate of HIV progression given years in each stage
                    % [1/(number of years in HIV stage)]
                    HIVprogRate_NatHistory = 1./Params.hiv_durHIVStage_NaturalHistory;
                    HIVprogRate_ANV_LatentA = 1./Params.hiv_durHIVStage_ANV_LatentA;
                    HIVprogRate_ANV_LatentB = 1./Params.hiv_durHIVStage_ANV_LatentB;
                    HIVprogRate_ANV_Late = 1./Params.hiv_durHIVStage_ANV_Late;
                    
                % Transitions while VLS are determined by rates read from
                % InitParams.m
                 
                    % Initialize index variables
                        Acute = Params.stage_Acute;
                        LatentA = Params.stage_LatentA;
                        LatentB = Params.stage_LatentB;
                        Late = Params.stage_Late;
                        AIDS = Params.stage_AIDS;
                
                 % Unaware
                 TransRates.ProgressHIV(Params.B1,Params.C1,:) = HIVprogRate_NatHistory(Acute,:);
                 TransRates.ProgressHIV(Params.C1,Params.D1,:) = HIVprogRate_NatHistory(LatentA,:);
                 TransRates.ProgressHIV(Params.D1,Params.E1,:) = HIVprogRate_NatHistory(LatentB,:);
                 TransRates.ProgressHIV(Params.E1,Params.F1,:) = HIVprogRate_NatHistory(Late,:);
                 TransRates.ProgressHIV(Params.F1,Params.Compart_AIDSDeath,:) = HIVprogRate_NatHistory(AIDS,:);

                 % Aware
                 TransRates.ProgressHIV(Params.B2,Params.C2,:) = HIVprogRate_NatHistory(Acute,:);
                 TransRates.ProgressHIV(Params.C2,Params.D2,:) = HIVprogRate_NatHistory(LatentA,:);
                 TransRates.ProgressHIV(Params.D2,Params.E2,:) = HIVprogRate_NatHistory(LatentB,:);
                 TransRates.ProgressHIV(Params.E2,Params.F2,:) = HIVprogRate_NatHistory(Late,:);
                 TransRates.ProgressHIV(Params.F2,Params.Compart_AIDSDeath,:) = HIVprogRate_NatHistory(AIDS,:);

                 % In Care
                 TransRates.ProgressHIV(Params.B3,Params.C3,:) = HIVprogRate_NatHistory(Acute,:);
                 TransRates.ProgressHIV(Params.C3,Params.D3,:) = HIVprogRate_NatHistory(LatentA,:);
                 TransRates.ProgressHIV(Params.D3,Params.E3,:) = HIVprogRate_NatHistory(LatentB,:);
                 TransRates.ProgressHIV(Params.E3,Params.F3,:) = HIVprogRate_NatHistory(Late,:);
                 TransRates.ProgressHIV(Params.F3,Params.Compart_AIDSDeath,:) = HIVprogRate_NatHistory(AIDS,:);

                 % On ART-not-VLS
                 TransRates.ProgressHIV(Params.C4,Params.D4,:) = HIVprogRate_ANV_LatentA;
                 TransRates.ProgressHIV(Params.D4,Params.E4,:) = HIVprogRate_ANV_LatentB;
                 TransRates.ProgressHIV(Params.E4,Params.F4,:) = HIVprogRate_ANV_Late;
                 % AIDS/ART-not-VLS: receive the ART death rate
                 
                 % VLS (CD4 count drop)
                 TransRates.ProgressHIV(Params.C5,Params.D5,:) = Params.hiv_rateVLSCD4decr_LatentA;
                 TransRates.ProgressHIV(Params.D5,Params.E5,:) = Params.hiv_rateVLSCD4decr_LatentB;
                 TransRates.ProgressHIV(Params.E5,Params.F5,:) = Params.hiv_rateVLSCD4decr_Late;
                 % AIDS/VLS: receive the ART death rate
        
                 % VLS (CD4 count increase)
                 % LatentA: no increase
                 TransRates.ProgressHIV(Params.D5,Params.C5,:) = Params.hiv_rateVLSCD4incr_LatentB;
                 TransRates.ProgressHIV(Params.E5,Params.D5,:) = Params.hiv_rateVLSCD4incr_Late;
                 TransRates.ProgressHIV(Params.F5,Params.E5,:) = Params.hiv_rateVLSCD4incr_AIDS;
        end

        % 4.ii. Normal death for [HIV-] and [PLWH not affected by ART]
        %This code should now only be used for HIV-
        for normDeathSection = 1:1
                
                % Calculate only at the beginning (t = 0)
                
                % Eligible:
                    % HIV-
                    % PLWH, not on ART or VLS
                    
                % Note: 
                    % PLWH, not AIDS, on ART/VLS do have "normal death" due to
                    % all-cause mortality, but it is is another section
                    % because it is updated each time period
                    
                    % PLWH, not on ART or VLS, but with AIDS can also have
                    % "AIDS death" included in other sections
                
                % Do not apply to absorbing states 
                    %Logical operator for people who have left the model
                    gone = zeros(Params.numComparts,1);

                    % Apply to absorbing states
                    for i = Params.AbsorbingComparts 
                        gone(i,1) = 1;
                    end
                               
                % Do not apply to CC stages
                    % Create matrix where CC compartments are
                        %Initialize
                        indic_ContinuumComparts = zeros(Params.numComparts,1);

                        % Apply 1s to the CC compartments; 0s
                        % elsewhere
                        indic_ContinuumComparts(Params.CCStages345,1) = 1;
                        indic_ContinuumComparts(Params.AwareComparts,1) = 1;
                        indic_ContinuumComparts(Params.UnawareComparts,1) = 1;

                % Apply all-cause mortality to TransRates struct
                    % Applied to all compartments except those removed
                TransRates.DieNormal(:, Params.Compart_NormDeath, :) = ((ones(Params.numComparts,1)-gone-indic_ContinuumComparts)*pop_dyingRate');
                
                % Clear unneeded variables
                clear gone
                 
        end
               
                   
        end

    end
    
%% 5. Calculate transition rates defined at the beginning of each year
    for sectionTRYearStart = 1:1
       
        % 5.i. Become aware and not immediately linked to care
        for awareSection = 1:1 
                        
            % Supporting calculations
                       
                    % Calculation of rate that an HIV+ individual will receive a
                    % positive test results and be notified
                    
                        % Adjusted 12Nov2014 to tag notification as pos
                        % only (previously the same value applied to both
                        % pos/neg notification - when alloc is applied, intn 
                        % dollars only go towards pos notification)
                        
                        WtdSensAndNotify_Acute = ...
                            TTProg.PctTestsRapid * TestSens_Rapid_Acute .* NotifyRapidProb' + ...
                            (1-TTProg.PctTestsRapid) * TestSens_Conv_Acute .* NotifyConvProb' ;

                        TTProg.WtdSensAndNotify_Acute = WtdSensAndNotify_Acute;

                        WtdSensAndNotify_Chronic = ...
                            TTProg.PctTestsRapid * TestSens_Rapid_Chronic .* NotifyRapidProb' + ...
                            (1-TTProg.PctTestsRapid) * TestSens_Conv_Chronic .* NotifyConvProb';

                        TTProg.WtdSensAndNotify_Chronic = WtdSensAndNotify_Chronic;
            
                   
            % Test Rates for non-CDC tests (used to calc the number of
            % tests)
            
                    % Assigned here to be read through to CalcNumIntnAffected.m
                    % file
                    
                    % <-- possibly update when allocation outcomes are
                    % finalized
            
                    TTProg.UninfectedTestRate = TestRate(:,1);
                    TTProg.AcuteTestRate = TestRate(:,2);
                    TTProg.LatentATestRate = TestRate(:,3);
                    TTProg.LatentBTestRate = TestRate(:,4);
                    TTProg.LateTestRate = TestRate(:,5);
                    TTProg.AIDSTestRate = TestRate(:,6);
            
                        
         % Fill in TransRates.Become aware struct:
              
            %  Eligible: Unaware
            % Calc = 
            %   [Test Rate] x  
            %   ([Percent Rapid] x [Test Sens Rapid] x [Prob Notify Rapid] +
            %   (1 - [Percent Rapid]) x [Test Sens Conv] x [Prob Notify Conv])
            %   x (1- [Prob LTC First])
             
              % Note: LinkageFirstProb is intentionally left as a
              % probability. 
            
                TransRates.BecomeAware(Params.B1, Params.B2, :) = TestRate(:,2) .* WtdSensAndNotify_Acute .* (1- LinkageFirstProb);
                TransRates.BecomeAware(Params.C1, Params.C2, :) = TestRate(:,3) .* WtdSensAndNotify_Chronic .* (1- LinkageFirstProb);
                TransRates.BecomeAware(Params.D1, Params.D2, :) = TestRate(:,4) .* WtdSensAndNotify_Chronic  .* (1- LinkageFirstProb);
                TransRates.BecomeAware(Params.E1, Params.E2, :) = TestRate(:,5) .* WtdSensAndNotify_Chronic .* (1- LinkageFirstProb);
                TransRates.BecomeAware(Params.F1, Params.F2, :) = TestRate(:,6) .* WtdSensAndNotify_Chronic .* (1- LinkageFirstProb);
         
        end
        
        
        % 5.ii. Link to Care
        
            % 5.ii.1. Become aware and immediately link to care
            for LTCFirstSection = 1:1
                
            %   Eligible: Unaware
            % Calc = 
            %   [Test Rate] x  
            %   ([Percent Rapid] x [Test Sens Rapid] x [Prob Notify Rapid] +
            %   (1 - [Percent Rapid]) x [Test Sens Conv] x [Prob Notify Conv])
            %   x [Prob LTC First]
            
            
            % Note: LinkageFirstProb is intentionally left as a probability. 
     
                TransRates.LinkToCare(Params.B1, Params.B3, :) = TestRate(:,2) .* WtdSensAndNotify_Acute .* LinkageFirstProb;
                TransRates.LinkToCare(Params.C1, Params.C3, :) = TestRate(:,3) .* WtdSensAndNotify_Chronic .* LinkageFirstProb;
                TransRates.LinkToCare(Params.D1, Params.D3, :) = TestRate(:,4) .* WtdSensAndNotify_Chronic .* LinkageFirstProb;
                TransRates.LinkToCare(Params.E1, Params.E3, :) = TestRate(:,5) .* WtdSensAndNotify_Chronic .* LinkageFirstProb;
                TransRates.LinkToCare(Params.F1, Params.F3, :) = TestRate(:,6) .* WtdSensAndNotify_Chronic .* LinkageFirstProb;
        
        
            end
            
            % 5.ii.2. Link to care after diagnosis
            for LTCAfterFirstSection = 1:1
            % Link to care in a period after when they were newly diagnosed
            %   Eligible: Aware, not in care
            %   Calc:(Pct of aware who are linked to care after diagnosis) 
                
            
                % Applies to aware, not in care rows
                TransRates.LinkToCare(Params.B2, Params.B3, :) = LinkageAfterRate(:,1); 
                TransRates.LinkToCare(Params.C2, Params.C3, :) = LinkageAfterRate(:,2);
                TransRates.LinkToCare(Params.D2, Params.D3, :) = LinkageAfterRate(:,3);
                TransRates.LinkToCare(Params.E2, Params.E3, :) = LinkageAfterRate(:,4);
                TransRates.LinkToCare(Params.F2, Params.F3, :) = LinkageAfterRate(:,5);
            end

            
        % 5.iii. Drop out of care, ART-not-VLS, VLS, or PrEP
            for dropOutSection = 1:1
                
            %   Eligible: [In care-not-on-ART], [ART-not-VLS], [VLS]
            %   Calc: (Pct of people in a stage who drop out of the stage)
                %   Previous assumption: model does not allow 
                %   people at AIDS stage to drop
                %   Revised 13Mar2015 under advisement of CDC SMEs
                

                % Drop out of [in care]
                    % Drop out of care from in care (In Care -> Aware)
                TransRates.DropOut(Params.B3, Params.B2, :) = DropOutRate_CareToAware;
                TransRates.DropOut(Params.C3, Params.C2, :) = DropOutRate_CareToAware;
                TransRates.DropOut(Params.D3, Params.D2, :) = DropOutRate_CareToAware;
                TransRates.DropOut(Params.E3, Params.E2, :) = DropOutRate_CareToAware;
                TransRates.DropOut(Params.F3, Params.F2, :) = DropOutRate_CareToAware;
                
                
                % Drop off of [ART not VLS]
                    
                    % Move to Aware (drop entirely from ART/Care)
                    % [ART not VLS] to [Aware]
                    
                        % Eligible: [ANV] who are not in care or ART
                    
                    TransRates.DropOut(Params.C4, Params.C2, :) = DropOutRate_ANVToAware * (1 - PctANVWhoAreInCareInclART);
                    TransRates.DropOut(Params.D4, Params.D2, :) = DropOutRate_ANVToAware * (1 - PctANVWhoAreInCareInclART);
                    TransRates.DropOut(Params.E4, Params.E2, :) = DropOutRate_ANVToAware * (1 - PctANVWhoAreInCareInclART);
                    TransRates.DropOut(Params.F4, Params.F2, :) = DropOutRate_ANVToAware * (1 - PctANVWhoAreInCareInclART);
                
                    
                    % Move to In Care 
                    % [ART not VLS] to [In Care]
                    
                    TransRates.DropOut(Params.C4, Params.C3, :) = DropOutRate_ANVToCare;
                    TransRates.DropOut(Params.D4, Params.D3, :) = DropOutRate_ANVToCare;
                    TransRates.DropOut(Params.E4, Params.E3, :) = DropOutRate_ANVToCare;
                    TransRates.DropOut(Params.F4, Params.F3, :) = DropOutRate_ANVToCare;
               
                    
               % Drop off [VLS]
                   % [VLS] to [ART not VLS]
                   
                    TransRates.DropOut(Params.C5, Params.C4, :) = DropOutRate_VLSToANV;
                    TransRates.DropOut(Params.D5, Params.D4, :) = DropOutRate_VLSToANV;
                    TransRates.DropOut(Params.E5, Params.E4, :) = DropOutRate_VLSToANV;
                    TransRates.DropOut(Params.F5, Params.F4, :) = DropOutRate_VLSToANV;
                    
                   % [VLS] to [LTC]
                   
                    TransRates.DropOut(Params.C5, Params.C3, :) = DropOutRate_VLSToLTC;
                    TransRates.DropOut(Params.D5, Params.D3, :) = DropOutRate_VLSToLTC;
                    TransRates.DropOut(Params.E5, Params.E3, :) = DropOutRate_VLSToLTC;
                    TransRates.DropOut(Params.F5, Params.F3, :) = DropOutRate_VLSToLTC;
                    
                   % [VLS] to [Aware]
                   
                    TransRates.DropOut(Params.C5, Params.C2, :) = DropOutRate_VLSToAware;
                    TransRates.DropOut(Params.D5, Params.D2, :) = DropOutRate_VLSToAware;
                    TransRates.DropOut(Params.E5, Params.E2, :) = DropOutRate_VLSToAware;
                    TransRates.DropOut(Params.F5, Params.F2, :) = DropOutRate_VLSToAware;
                
                   % Drop off [PrEP] % JCPrEPUpdate: modified to add dropout
                   % rates for each PrEP state
                    % [PrEP] to [Sus-no-PrEP]
                    TransRates.DropOut(Params.A6, Params.A1, :) = DropOutRate_PrEP_Oral_HighAdherence;
                    TransRates.DropOut(Params.A7, Params.A1, :) = DropOutRate_PrEP_Oral_LowAdherence;
                    TransRates.DropOut(Params.A8, Params.A1, :) = DropOutRate_PrEP_Inject_HighAdherence;
                    TransRates.DropOut(Params.A9, Params.A1, :) = DropOutRate_PrEP_Inject_LowAdherence;
 
                    
            end


        % 5.iv. Start ART
            for startArtSection = 1:1
             
            % Start ART Transitions
                % LTC->ANV
                % LTC->VLS
                % ANV->VLS (Only for NHAS analysis)
                % Sus-no-PrEP->On-PrEP 
                    
 
                
            % From LTC, no ART effects    
                % Eligible: In care, not affected by ART                
                % Note: Acute/In Care flows to LatentA/ANV or LatentA/VLS

                    % Move to ART-not-VLS
                        % Calc: [Rate of ART prescription] * (1 - [Pct become VLS])
                        TransRates.StartART(Params.B3, Params.C4, :) = InitiateARTRate(:,1) * (1-PctWhoBecomeVLS);
                        TransRates.StartART(Params.C3, Params.C4, :) = InitiateARTRate(:,2) * (1-PctWhoBecomeVLS);
                        TransRates.StartART(Params.D3, Params.D4, :) = InitiateARTRate(:,3) * (1-PctWhoBecomeVLS);
                        TransRates.StartART(Params.E3, Params.E4, :) = InitiateARTRate(:,4) * (1-PctWhoBecomeVLS);
                        TransRates.StartART(Params.F3, Params.F4, :) = InitiateARTRate(:,5) * (1-PctWhoBecomeVLS);


                   % Move to VLS
                        % Calc: [Rate of ART prescription] * [Pct become VLS]
                        TransRates.StartART(Params.B3, Params.C5, :) = InitiateARTRate(:,1) * PctWhoBecomeVLS;
                        TransRates.StartART(Params.C3, Params.C5, :) = InitiateARTRate(:,2) * PctWhoBecomeVLS;
                        TransRates.StartART(Params.D3, Params.D5, :) = InitiateARTRate(:,3) * PctWhoBecomeVLS;
                        TransRates.StartART(Params.E3, Params.E5, :) = InitiateARTRate(:,4) * PctWhoBecomeVLS;
                        TransRates.StartART(Params.F3, Params.F5, :) = InitiateARTRate(:,5) * PctWhoBecomeVLS;
                
                        
            % From ANV 
                
                % [ANV->VLS]
                
                    % Calc: [Rate of ART prescription] * [Pct become VLS]
                    TransRates.StartART(Params.C4, Params.C5, :) = ANVToVLSRate;
                    TransRates.StartART(Params.D4, Params.D5, :) = ANVToVLSRate;
                    TransRates.StartART(Params.E4, Params.E5, :) = ANVToVLSRate;
                    TransRates.StartART(Params.F4, Params.F5, :) = ANVToVLSRate;
                       
            % Start PrEP
                % [Sus-no-PrEP]->[PrEP]
                    if Year < Params.tt_periodFiveStartYear || Params.tt_progressionSource ~= 3 
                        PrEPInitiationRate_A6 = PrEPInitiationRate .* (1 - PctInjectPrEP) .* Params.PctHighAdherence_OralPrEP;
                        PrEPInitiationRate_A7 = PrEPInitiationRate .* (1 - PctInjectPrEP) .* (1 - Params.PctHighAdherence_OralPrEP);
                        PrEPInitiationRate_A8 = PrEPInitiationRate .* PctInjectPrEP .* Params.PctHighAdherence_InjectPrEP;
                        PrEPInitiationRate_A9 = PrEPInitiationRate .* PctInjectPrEP .* (1 - Params.PctHighAdherence_InjectPrEP);
                    else %KH added this else statement to make sure variables have values if the 
                        % if statement condition isn't met 
                        % (I think only time period 5 with AB)
                        % This variable will be calculated above in the AB
                        % code.
                    end
                    
                    TransRates.StartART(Params.A1, Params.A6, :) = PrEPInitiationRate_A6;
                    TransRates.StartART(Params.A1, Params.A7, :) = PrEPInitiationRate_A7;
                    TransRates.StartART(Params.A1, Params.A8, :) = PrEPInitiationRate_A8;
                    TransRates.StartART(Params.A1, Params.A9, :) = PrEPInitiationRate_A9;
             
            end
    
            
        % 5.v. Deaths of individuals with No ART, ANV, and VLS
            for deathAIDSARTSection = 1:1
            
            % Death of individuals who are on ART or are VLS
                % Considered normal death for PLWH non-AIDS
                % Considered AIDS death for PLWH with AIDS
                
                % Note: not stratified by ART-not-VLS and VLS; same rate is
                % applied for both continuum stages

                % Calc: death rates already set and apply to all 273 strata 
                % ensure dimensionality is right here - (30 x 30 x 273)
            
            % Convert to rates

                
                %NoARTEffectsDeathRate_NotAIDS = ProbToRate(NoARTEffectsDeathProbNotAIDS);
                %ANVDeathRate_NotAIDS = ProbToRate(ANVDeathProbNotAIDS);
                %ANVDeathRate_AIDS = ProbToRate(ANVDeathProbAIDS);
                %VLSDeathRate_NotAIDS = ProbToRate(VLSDeathProbNotAIDS);
                %VLSDeathRate_AIDS = ProbToRate(VLSDeathProbAIDS);             
            
                
            % NORMAL DEATH    
               
                % Eligible: all PLWH on ART and VLS except those with AIDS
                
                %No ART Effects
                %How can I condense code here?
                TransRates.DieNormal(Params.B1,Params.Compart_NormDeath,:) = NoARTEffectsDeathRate_NotAIDS;
                TransRates.DieNormal(Params.C1,Params.Compart_NormDeath,:) = NoARTEffectsDeathRate_NotAIDS;
                TransRates.DieNormal(Params.D1,Params.Compart_NormDeath,:) = NoARTEffectsDeathRate_NotAIDS;
                TransRates.DieNormal(Params.E1,Params.Compart_NormDeath,:) = NoARTEffectsDeathRate_NotAIDS;
                
                TransRates.DieNormal(Params.B2,Params.Compart_NormDeath,:) = NoARTEffectsDeathRate_NotAIDS;
                TransRates.DieNormal(Params.C2,Params.Compart_NormDeath,:) = NoARTEffectsDeathRate_NotAIDS;
                TransRates.DieNormal(Params.D2,Params.Compart_NormDeath,:) = NoARTEffectsDeathRate_NotAIDS;
                TransRates.DieNormal(Params.E2,Params.Compart_NormDeath,:) = NoARTEffectsDeathRate_NotAIDS;
                
                TransRates.DieNormal(Params.B3,Params.Compart_NormDeath,:) = NoARTEffectsDeathRate_NotAIDS;
                TransRates.DieNormal(Params.C3,Params.Compart_NormDeath,:) = NoARTEffectsDeathRate_NotAIDS;
                TransRates.DieNormal(Params.D3,Params.Compart_NormDeath,:) = NoARTEffectsDeathRate_NotAIDS;
                TransRates.DieNormal(Params.E3,Params.Compart_NormDeath,:) = NoARTEffectsDeathRate_NotAIDS;
                
                
                
                % ART-not-VLS
                
                TransRates.DieNormal(Params.C4,Params.Compart_NormDeath,:) = ANVDeathRate_NotAIDS;
                TransRates.DieNormal(Params.D4,Params.Compart_NormDeath,:) = ANVDeathRate_NotAIDS;
                TransRates.DieNormal(Params.E4,Params.Compart_NormDeath,:) = ANVDeathRate_NotAIDS;
                
                % VLS
                
                TransRates.DieNormal(Params.C5,Params.Compart_NormDeath,:) = VLSDeathRate_NotAIDS;
                TransRates.DieNormal(Params.D5,Params.Compart_NormDeath,:) = VLSDeathRate_NotAIDS;
                TransRates.DieNormal(Params.E5,Params.Compart_NormDeath,:) = VLSDeathRate_NotAIDS;
                
               
            % AIDS DEATH
                
                % Eligible: PLWH on ART or VLS with AIDS
                
                % ART-not-VLS
                TransRates.DieAIDS(Params.F4,Params.Compart_AIDSDeath,:) = ANVDeathRate_AIDS;
                % VLS
                TransRates.DieAIDS(Params.F5,Params.Compart_AIDSDeath,:) = VLSDeathRate_AIDS;
                
            end
    
    end

end  
%% 6. Transition Rates Calculated Every Time Step
for sectionTimeStepTRs = 1:1
    
    % 6.i. Assign Period-Specific Beta Values 
    for sectionPeriodBetas = 1:1
        
        % InfRate initialized to 0 so it can be listed as an output
        % argument
        if Params.SolveCont == 1
            InfRate = 0;
            InfRate_PrEP_Oral_High = 0;
            InfRate_PrEP_Inject_High = 0;
        end
    
        % First time period
        if Year < Params.tt_periodTwoStartYear
            TTProg.oneMinusBetaVonly_NoPrEP = Params.oneMinusBetaVonly_NoPrEP_period1;
            TTProg.oneMinusBetaAonly_NoPrEP = Params.oneMinusBetaAonly_NoPrEP_period1;
            TTProg.oneMinusBetaAsomeVI_NoPrEP = Params.oneMinusBetaAsomeVI_NoPrEP_period1;
            TTProg.oneMinusBetaVsomeAI_NoPrEP = Params.oneMinusBetaVsomeAI_NoPrEP_period1;
            if Params.CalcPrEPSpecificBetas == Params.ind_PrEP
                
                TTProg.oneMinusBetaVonly_PrEP = Params.oneMinusBetaVonly_PrEP_period1;
                TTProg.oneMinusBetaAonly_PrEP = Params.oneMinusBetaAonly_PrEP_period1;
                TTProg.oneMinusBetaAsomeVI_PrEP = Params.oneMinusBetaAsomeVI_PrEP_period1;
                TTProg.oneMinusBetaVsomeAI_PrEP = Params.oneMinusBetaVsomeAI_PrEP_period1;
                TTProg.oneMinusBetaN_PrEP = Params.oneMinusBetaN_PrEP_1to4;
            end

            TTProg.oneMinusBetaN_NoPrEP = Params.oneMinusBetaN_NoPrEP_1to4;            
            
            TTProg.pctMFPartnersA_IfSomeAInMF = Params.pctMFPartnersA_IfSomeAInMF_1;
            TTProg.pctPeopleWithA_InMF = Params.pctPeopleA_InMF_1;


        % Second time period    
        elseif Year < Params.tt_periodThreeStartYear
            TTProg.oneMinusBetaVonly_NoPrEP = Params.oneMinusBetaVonly_NoPrEP_period2;
            TTProg.oneMinusBetaAonly_NoPrEP = Params.oneMinusBetaAonly_NoPrEP_period2;
            TTProg.oneMinusBetaAsomeVI_NoPrEP = Params.oneMinusBetaAsomeVI_NoPrEP_period2;
            TTProg.oneMinusBetaVsomeAI_NoPrEP = Params.oneMinusBetaVsomeAI_NoPrEP_period2;
            if Params.CalcPrEPSpecificBetas == Params.ind_PrEP
            
                TTProg.oneMinusBetaVonly_PrEP = Params.oneMinusBetaVonly_PrEP_period2;
                TTProg.oneMinusBetaAonly_PrEP = Params.oneMinusBetaAonly_PrEP_period2;
                TTProg.oneMinusBetaAsomeVI_PrEP = Params.oneMinusBetaAsomeVI_PrEP_period2;
                TTProg.oneMinusBetaVsomeAI_PrEP = Params.oneMinusBetaVsomeAI_PrEP_period2;
                TTProg.oneMinusBetaN_PrEP = Params.oneMinusBetaN_PrEP_1to4;
            end
            
            TTProg.oneMinusBetaN_NoPrEP = Params.oneMinusBetaN_NoPrEP_1to4;            
            
            TTProg.pctMFPartnersA_IfSomeAInMF = Params.pctMFPartnersA_IfSomeAInMF_2to4;
            TTProg.pctPeopleWithA_InMF = Params.pctPeopleA_InMF_2to4;
        
        % Third time period    
        elseif Year < Params.tt_periodFourStartYear
            TTProg.oneMinusBetaVonly_NoPrEP = Params.oneMinusBetaVonly_NoPrEP_period3;
            TTProg.oneMinusBetaAonly_NoPrEP = Params.oneMinusBetaAonly_NoPrEP_period3;
            TTProg.oneMinusBetaAsomeVI_NoPrEP = Params.oneMinusBetaAsomeVI_NoPrEP_period3;
            TTProg.oneMinusBetaVsomeAI_NoPrEP = Params.oneMinusBetaVsomeAI_NoPrEP_period3;
            if Params.CalcPrEPSpecificBetas == Params.ind_PrEP
                
                TTProg.oneMinusBetaVonly_PrEP = Params.oneMinusBetaVonly_PrEP_period3;
                TTProg.oneMinusBetaAonly_PrEP = Params.oneMinusBetaAonly_PrEP_period3;
                TTProg.oneMinusBetaAsomeVI_PrEP = Params.oneMinusBetaAsomeVI_PrEP_period3;
                TTProg.oneMinusBetaVsomeAI_PrEP = Params.oneMinusBetaVsomeAI_PrEP_period3;
                TTProg.oneMinusBetaN_PrEP = Params.oneMinusBetaN_PrEP_1to4;
            end

            TTProg.oneMinusBetaN_NoPrEP = Params.oneMinusBetaN_NoPrEP_1to4;            
            
            TTProg.pctMFPartnersA_IfSomeAInMF = Params.pctMFPartnersA_IfSomeAInMF_2to4;
            TTProg.pctPeopleWithA_InMF = Params.pctPeopleA_InMF_2to4;
         
        % Fourth time period    
        elseif Year < Params.tt_periodFiveStartYear
            TTProg.oneMinusBetaVonly_NoPrEP = Params.oneMinusBetaVonly_NoPrEP_period4;
            TTProg.oneMinusBetaAonly_NoPrEP = Params.oneMinusBetaAonly_NoPrEP_period4;
            TTProg.oneMinusBetaAsomeVI_NoPrEP = Params.oneMinusBetaAsomeVI_NoPrEP_period4;
            TTProg.oneMinusBetaVsomeAI_NoPrEP = Params.oneMinusBetaVsomeAI_NoPrEP_period4;
            if Params.CalcPrEPSpecificBetas == Params.ind_PrEP
                
                TTProg.oneMinusBetaVonly_PrEP = Params.oneMinusBetaVonly_PrEP_period4;
                TTProg.oneMinusBetaAonly_PrEP = Params.oneMinusBetaAonly_PrEP_period4;
                TTProg.oneMinusBetaAsomeVI_PrEP = Params.oneMinusBetaAsomeVI_PrEP_period4;
                TTProg.oneMinusBetaVsomeAI_PrEP = Params.oneMinusBetaVsomeAI_PrEP_period4;
                TTProg.oneMinusBetaN_PrEP = Params.oneMinusBetaN_PrEP_1to4;
            end
            
            TTProg.oneMinusBetaN_NoPrEP = Params.oneMinusBetaN_NoPrEP_1to4;            

            TTProg.pctMFPartnersA_IfSomeAInMF = Params.pctMFPartnersA_IfSomeAInMF_2to4;
            TTProg.pctPeopleWithA_InMF = Params.pctPeopleA_InMF_2to4;

        % Fifth time period
        else
            TTProg.oneMinusBetaVonly_NoPrEP = Params.oneMinusBetaVonly_NoPrEP_period5;
            TTProg.oneMinusBetaAonly_NoPrEP = Params.oneMinusBetaAonly_NoPrEP_period5;
            TTProg.oneMinusBetaAsomeVI_NoPrEP = Params.oneMinusBetaAsomeVI_NoPrEP_period5;
            TTProg.oneMinusBetaVsomeAI_NoPrEP = Params.oneMinusBetaVsomeAI_NoPrEP_period5;
            if Params.CalcPrEPSpecificBetas == Params.ind_PrEP
                
                TTProg.oneMinusBetaVonly_PrEP = Params.oneMinusBetaVonly_PrEP_period5;
                TTProg.oneMinusBetaAonly_PrEP = Params.oneMinusBetaAonly_PrEP_period5;
                TTProg.oneMinusBetaAsomeVI_PrEP = Params.oneMinusBetaAsomeVI_PrEP_period5;
                TTProg.oneMinusBetaVsomeAI_PrEP = Params.oneMinusBetaVsomeAI_PrEP_period5;
                TTProg.oneMinusBetaN_PrEP = Params.oneMinusBetaN_PrEP_5;
            end
            
            TTProg.oneMinusBetaN_NoPrEP = Params.oneMinusBetaN_NoPrEP_5;           

            TTProg.pctMFPartnersA_IfSomeAInMF = Params.pctMFPartnersA_IfSomeAInMF_5;
            TTProg.pctPeopleWithA_InMF = Params.pctPeopleA_InMF_5;       
        end     
    end 
    
    % Discrete version only
    if Params.SolveCont ~= 1   

        % 6.ii. Discrete Version: Calculate Infection Rates
        for sectionDiscNewInfections = 1:1       

        % Determine infection rate
            % Note: this code is duplicated in RatesODE.m
        
            % Indicators    
            lambdasNoPrEP = Params.ind_NoPrEP; 
            lambdasPrEP = Params.ind_PrEP;
            
            % Normal infection rate
            InfRate = CalcInfectionRates(Params, Compartments, TTProg, lambdasNoPrEP, Params.ind_OralPrEP_High); 

            % PrEP infection rate

                % Not calculating PrEP-specific lambdas
                if Params.CalcPrEPSpecificBetas == Params.ind_NoPrEP

                    % Reduce normal lambdas by PrEP reduction %
                    % JCPrEPUpdate: modified to calculated reduced lambdas
                    % for oral and injectable PrEP
                        InfRate_PrEP_Oral_High.totalInfRate = InfRate.totalInfRate .* (1 - Params.inf_prepIncidenceReduct_Oral_HighAdherence);
                        InfRate_PrEP_Oral_High.lambdaVonly = InfRate.lambdaVonly .* (1 - Params.inf_prepIncidenceReduct_Oral_HighAdherence);
                        InfRate_PrEP_Oral_High.lambdaAonly = InfRate.lambdaAonly .* (1 - Params.inf_prepIncidenceReduct_Oral_HighAdherence);
                        InfRate_PrEP_Oral_High.lambdaVsomeAI = InfRate.lambdaVsomeAI .* (1 - Params.inf_prepIncidenceReduct_Oral_HighAdherence);
                        InfRate_PrEP_Oral_High.lambdaAsomeVI = InfRate.lambdaAsomeVI .* (1 - Params.inf_prepIncidenceReduct_Oral_HighAdherence);
                        InfRate_PrEP_Oral_High.lambdaN = InfRate.lambdaN .* (1 - Params.inf_prepIncidenceReduct_Oral_HighAdherence);
                        
                        InfRate_PrEP_Oral_Low.totalInfRate = InfRate.totalInfRate .* (1 - Params.inf_prepIncidenceReduct_Oral_LowAdherence);
                        InfRate_PrEP_Oral_Low.lambdaVonly = InfRate.lambdaVonly .* (1 - Params.inf_prepIncidenceReduct_Oral_LowAdherence);
                        InfRate_PrEP_Oral_Low.lambdaAonly = InfRate.lambdaAonly .* (1 - Params.inf_prepIncidenceReduct_Oral_LowAdherence);
                        InfRate_PrEP_Oral_Low.lambdaVsomeAI = InfRate.lambdaVsomeAI .* (1 - Params.inf_prepIncidenceReduct_Oral_LowAdherence);
                        InfRate_PrEP_Oral_Low.lambdaAsomeVI = InfRate.lambdaAsomeVI .* (1 - Params.inf_prepIncidenceReduct_Oral_LowAdherence);
                        InfRate_PrEP_Oral_Low.lambdaN = InfRate.lambdaN .* (1 - Params.inf_prepIncidenceReduct_Oral_LowAdherence);
                        
                        InfRate_PrEP_Inject_High.totalInfRate = InfRate.totalInfRate .* (1 - Params.inf_prepIncidenceReduct_Inject_HighAdherence);
                        InfRate_PrEP_Inject_High.lambdaVonly = InfRate.lambdaVonly .* (1 - Params.inf_prepIncidenceReduct_Inject_HighAdherence);
                        InfRate_PrEP_Inject_High.lambdaAonly = InfRate.lambdaAonly .* (1 - Params.inf_prepIncidenceReduct_Inject_HighAdherence);
                        InfRate_PrEP_Inject_High.lambdaVsomeAI = InfRate.lambdaVsomeAI .* (1 - Params.inf_prepIncidenceReduct_Inject_HighAdherence);
                        InfRate_PrEP_Inject_High.lambdaAsomeVI = InfRate.lambdaAsomeVI .* (1 - Params.inf_prepIncidenceReduct_Inject_HighAdherence);
                        InfRate_PrEP_Inject_High.lambdaN = InfRate.lambdaN .* (1 - Params.inf_prepIncidenceReduct_Inject_HighAdherence);

                        InfRate_PrEP_Inject_Low.totalInfRate = InfRate.totalInfRate .* (1 - Params.inf_prepIncidenceReduct_Inject_LowAdherence);
                        InfRate_PrEP_Inject_Low.lambdaVonly = InfRate.lambdaVonly .* (1 - Params.inf_prepIncidenceReduct_Inject_LowAdherence);
                        InfRate_PrEP_Inject_Low.lambdaAonly = InfRate.lambdaAonly .* (1 - Params.inf_prepIncidenceReduct_Inject_LowAdherence);
                        InfRate_PrEP_Inject_Low.lambdaVsomeAI = InfRate.lambdaVsomeAI .* (1 - Params.inf_prepIncidenceReduct_Inject_LowAdherence);
                        InfRate_PrEP_Inject_Low.lambdaAsomeVI = InfRate.lambdaAsomeVI .* (1 - Params.inf_prepIncidenceReduct_Inject_LowAdherence);
                        InfRate_PrEP_Inject_Low.lambdaN = InfRate.lambdaN .* (1 - Params.inf_prepIncidenceReduct_Inject_LowAdherence);
                        
                    % Relative rate
                        InfRate_PrEP_Oral_High.RelativeInfRates = InfRate.RelativeInfRates;
                        InfRate_PrEP_Oral_Low.RelativeInfRates = InfRate.RelativeInfRates;
                        InfRate_PrEP_Inject_High.RelativeInfRates = InfRate.RelativeInfRates;
                        InfRate_PrEP_Inject_Low.RelativeInfRates = InfRate.RelativeInfRates;

                else % Calculate PrEP-specific lambdas (incorporate behavior change) % JCPrEPUpdate: modified to call CalcInfectionRates separately for oral and injectable PrEP
                   
                    InfRate_PrEP_Oral_High = CalcInfectionRates(Params, Compartments, TTProg, lambdasPrEP, Params.ind_OralPrEP_High);
                    InfRate_PrEP_Oral_Low = CalcInfectionRates(Params, Compartments, TTProg, lambdasPrEP, Params.ind_OralPrEP_Low);
                    InfRate_PrEP_Inject_High = CalcInfectionRates(Params, Compartments, TTProg, lambdasPrEP, Params.ind_InjectPrEP_High);
                    InfRate_PrEP_Inject_Low = CalcInfectionRates(Params, Compartments, TTProg, lambdasPrEP, Params.ind_InjectPrEP_Low);
                    
                end
            

        % Apply to TransRates % JCPrEPUpdate: modified to include
        % transitions from all PrEP states
        	TransRates.GetInfected(Params.A1,Params.B1,:) = InfRate.totalInfRate;
            TransRates.GetInfected(Params.A6,Params.B3,:) = InfRate_PrEP_Oral_High.totalInfRate;
            TransRates.GetInfected(Params.A7,Params.B3,:) = InfRate_PrEP_Oral_Low.totalInfRate;
            TransRates.GetInfected(Params.A8,Params.B3,:) = InfRate_PrEP_Inject_High.totalInfRate;
            TransRates.GetInfected(Params.A9,Params.B3,:) = InfRate_PrEP_Inject_Low.totalInfRate;

        end

        % 6.iii. Discrete Version: Calculate number of people tested
            % Value 1 is an input for the discrete version.
        TTProg = CalcNumIntnAffected(Compartments,Params,TTProg,Year,1);

        % 6.iv. Discrete Version: Calculate Full Transition Rate Matrix
        for sectionDiscFullTP = 1:1
            
        % Check that no two matrices have non-zero values at the same spot.          
            chknonzero = (TransRates.BecomeAware>0)   + (TransRates.GetInfected>0)   ...
                       + (TransRates.ProgressHIV>0)   + (TransRates.LinkToCare>0)     ...
                       + (TransRates.DropOut>0)   + (TransRates.StartART>0)       ...
                       + (TransRates.DieNormal>0)      ...
                       + (TransRates.DieAIDS>0)    + (TransRates.StayDeadAged>0);            

            if max(max(max(chknonzero)))>1
                error('ERROR! Elements of TransRates matrix conflict! Nonzero elems in more than 1 Trans')
            end
               clear chknonzero

              
        % Determine full transition rate matrix
        TransRates.AnnualRates = ...
                         TransRates.BecomeAware   + TransRates.GetInfected   ...
                       + TransRates.ProgressHIV   + TransRates.LinkToCare     ...
                       + TransRates.DropOut   + TransRates.StartART       ...
                       + TransRates.DieNormal     ...
                       + TransRates.DieAIDS    ;

        % Make Diagonal equal to the negative sum of the rows
        for i = 1:Params.numStrats
            TransRates.AnnualRates(:,:,i) = TransRates.AnnualRates(:,:,i)+...
                diag(-sum(TransRates.AnnualRates(:,:,i),2));
        end
        
        
        
        % Convert from annual rates to timestep rates        
        TransRates.TSRates =  real(TransRates.AnnualRates * Params.tt_tstep);
       
        % Clear
        clear diagsum diagones chkdiag Transitions

        end
    
    end
end


%% 7. Conversion Functions
    function OutRate = ProbToRate(InAnnualProb)
        if (InAnnualProb >= 1) 
            InAnnualProb = 0.99999;
        end
            
        OutRate = -log(1 - InAnnualProb);
    end

    function OutProb = RateToProb(InAnnualRate)
        OutProb = 1 - exp(-InAnnualRate);
    end

%% 8. Functions for estimating allocation-based progression rates

%% 8a. Calculate number interventions given if applying time period 1 (T1) rates / probs and intervention cost reach threshold rates / probs
    function NumIntnAtT1andReachLvlRates = ...
            CalcNumIntnAtT1andReachLvlRates(nNumEligible, ...
                sngT1RatesProbs, sngReachLvls, sngMaxReach, nIndicators_LimitPops, nInputType)
            
        %nInputType: 1 = transition rate, 2 = % of eligible
        
        % Calculate number of interventions per eligible person given T1
        if nInputType == 1 % Input is rate (Test, LTC after, ART init, BecomeVLS)
            NumIntnPPAtT1Rates = sngT1RatesProbs;
            sngT1Prob = RateToProb(sngT1RatesProbs);
            sngReachLvl(:,1) = ProbToRate(sngT1Prob + max(sngMaxReach-sngT1Prob,0)*sngReachLvls(1)); 
            sngReachLvl(:,2) = ProbToRate(sngT1Prob + max(sngMaxReach-sngT1Prob,0)*(sngReachLvls(1)+sngReachLvls(2))); 
            sngMaxReach = ProbToRate(sngMaxReach);
        else % Input is % of eligible (LTC at diag, RemainVLS, SEP, PrEP)
            NumIntnPPAtT1Rates = sngT1RatesProbs;
            sngReachLvl(:,1) = sngT1RatesProbs + max(sngMaxReach-sngT1RatesProbs,0)*sngReachLvls(1);
            sngReachLvl(:,2) = sngT1RatesProbs + max(sngMaxReach-sngT1RatesProbs,0)*(sngReachLvls(1)+sngReachLvls(2));
        %else % Input is prob transition
            %NumIntnPPAtT1Rates = ProbToRate(sngT1RatesProbs);          
            %sngReachLvl(:,1) = sngT1RatesProbs + max(sngMaxReach-sngT1RatesProbs,0)*sngReachLvls(1);
            %sngReachLvl(:,2) = sngT1RatesProbs + max(sngMaxReach-sngT1RatesProbs,0)*(sngReachLvls(1)+sngReachLvls(2));
        end
        
        % Calculate number of interventions given in each subpop if applying time 
        % period 1 (T1) rates / probs, applying max reach constraints
        NumIntnAtT1andReachLvlRates(:,1) = nNumEligible .* NumIntnPPAtT1Rates .* nIndicators_LimitPops;
        NumIntnAtT1andReachLvlRates(:,2) = nNumEligible .* sngReachLvl(:,1) .* nIndicators_LimitPops;
        NumIntnAtT1andReachLvlRates(:,3) = nNumEligible .* sngReachLvl(:,2) .* nIndicators_LimitPops;
        NumIntnAtT1andReachLvlRates(:,4) = nNumEligible .* max(NumIntnPPAtT1Rates,sngMaxReach) .* nIndicators_LimitPops;
                       
    end

%% 8b. Calculate rates / probs for intervention cost thresholds
    function T1andReachLvlIntnsPP = ...
            CalcT1andReachLvlIntnsPP(nNumEligible, ...
                sngNumIntnAtT1andRchLvlRates, nIndicators_LimitPops, nIndicator_SEP)
        
        if nIndicator_SEP == 1
            T1andReachLvlIntnsPP (1,1:4) = 0;
        else
            T1andReachLvlIntnsPP (1:Params.numStrats,1:4) = 0;
        end
    
        for Lvl = 1:4 
            
            T1andReachLvlIntnsPP (:,Lvl) = (sngNumIntnAtT1andRchLvlRates(:,Lvl) ./ max(nNumEligible,0.000000000000000000001)) .* nIndicators_LimitPops;
        
        end
        
    end

%% 8c. Calculate cost to reach time period 1 rates / probs and other 3 reach levels, given base rates get you part way
    function CostT1andReachLvl = ...
            CalcCostT1andReachLvl(sngIntnCostPP, sngPctIntnCostByReachLvl, nNumEligible, sngBaseRatesProbs, ...
            nNumIntnAtT1andReachLvlRates, nInputType, IsPrEP, sngIntnCostPP_BeyondCurrent)
        
        %nInputType: 1 = transition rate, 2 = percent of eligible

        if nInputType == 1 % Input is rate (Test, LTC after, ART init, BecomeVLS)
            NumIntnPPAtBaseRates = sngBaseRatesProbs;
        else  % Input is % of eligible (LTC at diag, RemainVLS, SEP, PrEP)
            NumIntnPPAtBaseRates = sngBaseRatesProbs;
        %else % Input is probability of transition
            %NumIntnPPAtBaseRates = ProbToRate(sngBaseRatesProbs);
        end
        
        if IsPrEP == 0 % Intervention is not PrEP
            sngIntnCostPP_BeyondCurrent = sngIntnCostPP; %PrEP has a different cost PP based 
            % on if costs are for those already on PrEP or if costs are for
            % those initiating PrEP. If intervention is not PrEP, there is
            % only one intervention cost PP
        end
        
        % Calculate number in each subpop that get intn at base rates                    
        NumIntnAtBaseRates = nNumEligible .* NumIntnPPAtBaseRates;
        
        % Calculate cost of reaching time period 1 rates, given base rates
        % may get you part way or all the way to time period 1 rates
        CostT1andReachLvl(:,1) = max(nNumIntnAtT1andReachLvlRates(:,1) - NumIntnAtBaseRates,0) .* ...
            sngIntnCostPP;
        
        for ReachLevel = 1:3
            
            CostT1andReachLvl(:,ReachLevel + 1) = ...
                CostT1andReachLvl(:,ReachLevel) + ...
                max(nNumIntnAtT1andReachLvlRates(:,ReachLevel + 1) ...
                - max(nNumIntnAtT1andReachLvlRates(:,ReachLevel),NumIntnAtBaseRates),0) .* ...
            (sngIntnCostPP_BeyondCurrent * sngPctIntnCostByReachLvl(ReachLevel,1));
            
        end
    
    end

%% 8d. Calculate percentage of T1 and reach level rates in a subpopulation afforded by each intervention

    function AddlPctT1andReachLvlRatesAffordedByAlloc = ...
            CalcAddlPctT1andReachLvlRatesAffordedByAlloc(nAlloc, ...
            nCostT1andReachLvlRate, nIndicators_LimitPops, nIndicator_SEP)
        
        if nIndicator_SEP == 1
            AddlPctT1andReachLvlRatesAffordedByAlloc(1,Params.numReachLvls)=0;
        else
            AddlPctT1andReachLvlRatesAffordedByAlloc(Params.numStrats,Params.numReachLvls)=0;
        end
%         for ReachLevel = 1:Params.numReachLvls
%             
%             AddlPctT1andReachLvlRatesAffordedByAlloc(:,ReachLevel) = ...
%                 (min(nAlloc / ...
%                 max(sum(nCostT1andReachLvlRate(:,ReachLevel,:),3)' * nIndicators_LimitPops,1),1)) ...
%                 .* nIndicators_LimitPops;

        nCostT1andReachLvlRate_acrossTestEligStages = sum(nCostT1andReachLvlRate(:,:,:),3);
        % Identify reach level thresholds around the reach obtained by the allocation
        if nAlloc >= max(nCostT1andReachLvlRate_acrossTestEligStages(:,1)'*nIndicators_LimitPops,1)
            if nAlloc >= max(nCostT1andReachLvlRate_acrossTestEligStages(:,4)'*nIndicators_LimitPops,1)
                nHighestReachLvl = 4;
                nAlloc = ...
                    nCostT1andReachLvlRate_acrossTestEligStages(:,4)'*nIndicators_LimitPops;
            elseif nAlloc >= max(nCostT1andReachLvlRate_acrossTestEligStages(:,3)'*nIndicators_LimitPops,1)
                nHighestReachLvl = 4;
            elseif nAlloc >= max(nCostT1andReachLvlRate_acrossTestEligStages(:,2)'*nIndicators_LimitPops,1)
                nHighestReachLvl = 3;
            else %nAlloc >= max(nCostT1andReachLvlRate_acrossTestEligStages(:,1)'*nIndicators_LimitPops,1)
                nHighestReachLvl = 2;
            end
            nCostPassedReachLvl = nCostT1andReachLvlRate_acrossTestEligStages(:,nHighestReachLvl-1) .* nIndicators_LimitPops;
            nCostHighestReachLvl = nCostT1andReachLvlRate_acrossTestEligStages(:,nHighestReachLvl) .* nIndicators_LimitPops;
            AddlPctT1andReachLvlRatesAffordedByAlloc(:,1:nHighestReachLvl-1)=nIndicators_LimitPops*ones(1,nHighestReachLvl-1);
        else
            nHighestReachLvl = 1;
            nCostPassedReachLvl = 0;
            nCostHighestReachLvl = nCostT1andReachLvlRate_acrossTestEligStages(:,1) .* nIndicators_LimitPops;
%             AddlPctT1andReachLvlRatesAffordedByAlloc(:,1) = ...
%                 (min(nAlloc / ...
%                 max(nCostT1andReachLvlRate_acrossTestEligStages(:,1)' * nIndicators_LimitPops,1),1)) ...
%                 .* nIndicators_LimitPops;            
        end
        

        if sum(max(nCostHighestReachLvl - nCostPassedReachLvl,0)) > 0 %only calculate if someone is reachable
            % If allocation is not more than the cost of reaching current rates,
            % then the rate afforded over current rates is 0. Otherwise, calculate
            % rate afforded by allocation over current rates.

            % Calculate rate that can be purchased with remaining allocation
            % (beyond cost of applying current rates), given number of people
            % eligible and reachable and increased cost of intervention
            % when rate is beyond current
            AddlPctT1andReachLvlRatesAffordedByAlloc(:,nHighestReachLvl) = ...
                (nAlloc - sum(nCostPassedReachLvl)) / ...
                sum(nCostHighestReachLvl - nCostPassedReachLvl) .* nIndicators_LimitPops;
            

        end
    end

end

