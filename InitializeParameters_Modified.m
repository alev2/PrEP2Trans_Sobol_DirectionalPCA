function [Params]=InitializeParameters_Modified(ExcelValues_AllParameters,ExcelValues_Populations)
%% Purpose: Import model inputs from Excel file and initialize them to be in a form usable to Matlab

%% 1. Nomenclature

% Struct: all parameters are being initialized within the Params struct
%   Params:

% Excel model sheet: variables are grouped according to the relevant sheet
%  in the Excel model
%   pop_        Population-related variables
%   tt_         Test and treat progression
%   behav_      Behaviors
%   hiv_        HIV Progression
%   inf_        Infectivity

% Parameter name: name of parameter

% stratification: if the parameter is stratified, the name is followed by the 
% stratification codes in the order of the dimensions
%   a           age (13-17, 18-24, 25-34, 35-44, 45-54, 55-64, 65+)
%   c           circumcision status (uncircumcised, circumcised)
%   d           disease progression (Acute, CD4>500, CD4 350-500, CD4 200-
%               350, CD4<200) 
%                   equivalent to: Acute, LatentA, LatentB, Late, AIDS
%   l           risk level (Main, Casual) 
%                   equivalent to: Low, High
%   p           risk group/population (HET, MSM, IDU)
%   r           race (Black, Hispanic, Other)
%   s           sex (Male, Female)
%   t           test and treat continuum 
%                   (Unaware, Diagnosed, In care, ANV, VLS, PrEP) 
%                   equivalent to (1,2,3,4,5, 6)
%   y           time periods the rates correspond to 
%                   STELLA nomenclature: (init, intermediate, interv) or (init, pre, post)

%% 2. Define variables that describe model structure

    % 2.i. Define the Compartments reference variables
    for compartmentsSection = 1:1
   
    % INDIVIDUAL COMPARTMENTS
    Params.A1=1; % A1 Sus, no-PrEP
    Params.A6=2; % A6 Sus, on-oral-PrEP, high adherence
    Params.A7=3; % A7 Sus, on-oral-PrEP, low adherence JCPrEPUpdate: new PrEP states starting here - other compartment nums shifted by 3
    Params.A8=4; % A8 Sus, on-injectable-PrEP, high adherence
    Params.A9=5; % A9 Sus, on-injectable-PrEP, low adherence
    Params.B1=6; % B1 Acute UA
    Params.B2=7; % B2 Acute Aware
    Params.B3=8; % B3 Acute In Care
    Params.C1=9; % C1 CD4 500 UA
    Params.C2=10; % C2 CD4 500 Aware
    Params.C3=11; % C3 CD4 500 In Care
    Params.C4=12;% C4 CD4 500 ART-not-VLS
    Params.C5=13;% C5 CD4 500 VLS
    Params.D1=14;% D1 CD4 350 UA
    Params.D2=15;% D2 CD4 350 Aware
    Params.D3=16;% D3 CD4 350 In Care
    Params.D4=17;% D4 CD4 350 ART-not-VLS
    Params.D5=18;% D5 CD4 350 VLS
    Params.E1=19;% E1 CD4 200 UA
    Params.E2=20;% E2 CD4 200 Aware
    Params.E3=21;% E3 CD4 200 In Care
    Params.E4=22;% E4 CD4 200 ART-not-VLS
    Params.E5=23;% E5 CD4 200 VLS
    Params.F1=24;% F1 AIDS UA
    Params.F2=25;% F2 AIDS Aware
    Params.F3=26;% F3 AIDS In Care
    Params.F4=27;% F4 AIDS ART-not-VLS
    Params.F5=28;% F5 AIDS VLS
    Params.Compart_NormDeath =29;% Normal Death
    Params.Compart_AIDSDeath =30;% AIDS/ART Death

        
% TT CONTINUUM COMPARTMENTS

    % Unaware
    UnawareComparts = [Params.B1, Params.C1,Params.D1, Params.E1, Params.F1];
    
    % Aware
    AwareComparts = [Params.B2:Params.B3,Params.C2:Params.C5,Params.D2...
        :Params.D5,Params.E2:Params.E5,Params.F2:Params.F5];
    % In Care
    CCStages345 = [Params.B3,Params.C3:Params.C5,Params.D3...
        :Params.D5,Params.E3:Params.E5,Params.F3:Params.F5];
    
    % ART only
    ANVComparts = [Params.C4, Params.D4,Params.E4,Params.F4]; 
        
    % ART and VLS
   ARTandVLSComparts = [Params.C4:Params.C5,Params.D4:Params.D5,...
       Params.E4:Params.E5,Params.F4:Params.F5]; 
    % VLS
   VLSComparts = [Params.C5, Params.D5, Params.E5, Params.F5];
   % PrEP
   PrEPComparts = [Params.A6:Params.A9];
   % Oral PrEP
   OralPrEPComparts = [Params.A6:Params.A7];
   % Injectable PrEP
   InjectPrEPComparts = [Params.A8:Params.A9];
        
    % Apply to Params
    Params.UnawareComparts = UnawareComparts;
    Params.AwareComparts = AwareComparts;
    Params.CCStages345 = CCStages345;
    Params.ANVComparts = ANVComparts;
    Params.ANVandVLSComparts = ARTandVLSComparts;
    Params.VLSComparts = VLSComparts;
    Params.PrEPComparts = PrEPComparts;
    Params.OralPrEPComparts = OralPrEPComparts;
    Params.InjectPrEPComparts = InjectPrEPComparts;
        
    % DISEASE PROGRESSION COMPARTMENTS
    UninfectedComparts = [Params.A1 PrEPComparts]; % JCPrEPUpdate: other PrEP states included in uninf states
    AcuteComparts = [Params.B1:Params.B3];
    LatentAComparts = [Params.C1:Params.C5];
    LatentBComparts = [Params.D1:Params.D5];
    LateComparts = [Params.E1:Params.E5];
    AIDSComparts = [Params.F1:Params.F5];
    HIVComparts = [Params.B1:Params.F5];
    
    Params.NonAbsorbingComparts = [Params.A1:Params.F5];
    Params.AbsorbingComparts = [Params.Compart_NormDeath:Params.Compart_AIDSDeath];
    
    Params.UninfectedComparts = UninfectedComparts;
    Params.AcuteComparts = AcuteComparts;
    Params.LatentAComparts = LatentAComparts;
    Params.LatentBComparts = LatentBComparts;
    Params.LateComparts = LateComparts;
    Params.AIDSComparts = AIDSComparts;
    Params.HIVComparts = HIVComparts;
    
    
    end
    
    % 2.ii. Define stratification variables
    for stratificationsSection = 1:1
        Params.numStrats = 273; % all subpopulations
        Params.numComparts = numel([Params.NonAbsorbingComparts Params.AbsorbingComparts]); % JCPrEPUpdate
        Params.numAbsorbingStates = numel(Params.NonAbsorbingComparts);
        Params.numHIVNegStates = numel(UninfectedComparts); % JCPrEPUpdate
        Params.numPrEPStates = numel(PrEPComparts); % JCPrEPUpdate: added
        Params.numTestEligStates = numel([Params.A1 UnawareComparts]); % JCPrEPUpdate: PrEP states are not included
        Params.numReachLvls = 4; %Used for intervention costs that vary by reach
        Params.numAge = 7;
        numPop = 3;
        Params.numPop=numPop;
        numPop_plusHETbySex = 4;
        Params.numPop_plusHETbySex = numPop_plusHETbySex;
        numCirc=2;
        numSex=2;
        numRace=3;
        Params.numRace = numRace;
        numRiskLevel = 2;
        numPeriod = 5; %Number of rate periods. Updated by MClinkscales on 5/24/2022
        Params.numRateInputPeriod = 3; %Number of periods in which progression rates are specificed. Added by MClinkscales on 5/25/2022
        Params.numPeriod = numPeriod;
        Params.numSEPTimePeriods = 5;
        Params.numPrEPTimePeriods = 5;
        numHIVstages = 5;
        Params.numHIVstages = numHIVstages;
        numChronicDiseaseStages = 4;
        numDistinctContStagesAcute = 3; % Unaware,Aware,LTC
        numDistinctContStagesWithCosts = 4; % See below for notes:
            % Based on the "annual care and treatment costs input in Excel"
            % The 4 refers to unaware, aware, LTC, ART
 
        Params.numModelsContStagesIfChronic = 5; % Unaware, Aware no ART, LTC no ART, ANV, VLS
        
        numTransGroupSex_SexualCombos = 18; % KH updated from 13 to 18 to add in HET mixing by sex 20Sep2023
        numTransGroupSex_NeedleCombos = 4;
        Params.numIntnReachCategories = 12; % Used for min / max reach inputs
        Params.numIntnIncrCostCategories = 8;
        % Number of each type of intervention
        Params.numIntns_Testing = 5; %Testing
        Params.numIntns_Testing_YMSM = 2; %Testing (YMSM only)
        Params.numIntns_LTCatDiag = 1; % LTC at diagnosis 
        Params.numIntns_LTCafterDiag = 1; % LTC after diagnosis
        Params.numIntns_ARTPrescription = 1;% ART prescription
        Params.numIntns_ARTAdherence = 2;% ART adherence
        Params.numIntns_SEP = 3; % Syringe services by race/ethnicity
        Params.numIntns_PrEPType = 2; % Oral and injectable
        Params.numIntns_PrEP = 24;    % PrEP % JCPrEPUpdate: 12 for oral PrEP and 12 for injectable PrEP
        Params.numIntns_PrEP_YMSM = 6;    % PrEP (YMSM only) % JCPrEPUpdate: 3 for oral PrEP and 3 for injectable PrEP
        Params.numIntns = Params.numIntns_Testing+...
            Params.numIntns_LTCatDiag + ...
            Params.numIntns_LTCafterDiag + ...
            Params.numIntns_ARTPrescription+ ...
            Params.numIntns_ARTAdherence+ ...
            Params.numIntns_SEP + ...
            Params.numIntns_PrEP;
        Params.numIntns_YMSM = Params.numIntns_Testing_YMSM+...
            Params.numIntns_LTCatDiag + ...
            Params.numIntns_LTCafterDiag + ...
            Params.numIntns_ARTPrescription+ ...
            Params.numIntns_ARTAdherence+ ...            
            Params.numIntns_PrEP_YMSM;
        Params.YMSMint = [3:4, 6:10, 20:22, 32:34]; % edit if list of interventions change % JCPrEPUpdate: updated to all PrEP MSM int (20:22 = oral PrEP; 32:34 = inject PrEP)
        % Number of reach levels beyond T2 rates/probabilities
        Params.numReachLevels = 3;
        % Number of allocation periods
        Params.numMaxAllocationPeriods = 3;
    end

    % 2.iii. Define indicator variables
    for indicatorSection = 1:1
        
    % Preallocate indicator arrays
        Params.popIndicator = zeros(Params.numStrats,numPop);
        Params.popIndicator_withHETbySex = zeros(Params.numStrats,numPop+1); % JCPrEPUpdate: new indicator to align with transmission/sex groups for PrEP initiation
        Params.circIndicator = zeros(Params.numStrats,numCirc);
        Params.riskLevelIndicator = zeros(Params.numStrats,numRiskLevel);
        Params.sexIndicator= zeros(Params.numStrats,numSex);
        Params.raceIndicator = zeros(Params.numStrats,numRace);
        Params.popSexIndicator = zeros(Params.numStrats,numSex*numPop - 1);

                
% INDICATOR: RACE/ETHNICITY

    % Fill indicator
        Params.raceIndicator(1:91,1) = 1; % Black
        Params.raceIndicator(92:182,2) = 1; % Hispanic
        Params.raceIndicator(183:273,3) = 1; % Other

    % Corresponding column numbers
        Params.race_B = 1;
        Params.race_H = 2;
        Params.race_O = 3;
        
 % INDICATOR: AGE
    
    % Fill indicator
        Params.ageIndicator = repmat(eye(Params.numAge),[Params.numStrats/Params.numAge,1]);
    
    % Corresponding column numbers
        Params.age_1 = 1;
        Params.age_2 = 2;
        Params.age_3 = 3;
        Params.age_4 = 4;
        Params.age_5 = 5;
        Params.age_6 = 6;
        Params.age_7 = 7;

        
% INDICATOR: CIRCUMCISION STATUS

    % Fill indicator
        % Uncirc
        Params.circIndicator(1:35,1)=1;
        Params.circIndicator(71:126,1)=1;
        Params.circIndicator(162:217,1)=1;
        Params.circIndicator(253:273,1)=1;
        
        % Circ
        Params.circIndicator(36:70,2)=1;
        Params.circIndicator(127:161,2)=1;
        Params.circIndicator(218:252,2)=1;

    % Corresponding column number
        Params.circ_U = 1;
        Params.circ_C = 2;
        
        
% INDICATOR: TRANSISSION GROUP AKA POPULATION
    % Note: labeled as "pop" throughout because of the terminology that was
    % being used when it was initially coded
        
    % Column numbers
        Params.pop_HET = 1;
        Params.pop_MSM = 2;
        Params.pop_IDU = 3;
    
    % Fill indicator
        %HET Column
            storeHET = repmat([21;21;7],[3,1]);
            count = 1;
            for a = 1:9
                Params.popIndicator(count:count+13,Params.pop_HET)=1;
                count = count+storeHET(a)+14;
            end

        %MSM Column
            storeMSM = repmat([21;42],[3,1]);
            count = 15;
            for a = 1:6
               Params.popIndicator(count:count+13,Params.pop_MSM)=1;
               count = count+storeMSM(a)+14;
            end

        %IDU Column
            storeIDU = repmat([28;14;28],[3,1]);
            count = 29;
            for a = 1:9
               Params.popIndicator(count:count+6,Params.pop_IDU)=1;
               count = count+storeIDU(a)+7;
            end

% INDICATOR: SEX

    % Fill indicator
        % Male
        count = 1;
        for a = 1:3
            Params.sexIndicator(count:count+69,1)=1;
            count = count+21+70;
        end
        
        % Female
        Params.sexIndicator(:,2) = abs(1-Params.sexIndicator(:,1));

    % Corresponding column numbers for sex indicator
        Params.sex_Male = 1;
        Params.sex_Female = 2;
                  
            
% INDICATOR: POP (HET BY SEX) % JCPrEPUpdate: new indicator to align with transmission/sex groups for PrEP initiation
           
    Params.popHETbySex_HETM = 1;
    Params.popHETbySex_HETF = 2;
    Params.popHETbySex_MSM = 3;
    Params.popHETbySex_IDU = 4;

    Params.popIndicator_withHETbySex(:,Params.popHETbySex_HETM) = Params.popIndicator(:,Params.pop_HET).*Params.sexIndicator(:,Params.sex_Male);        
    Params.popIndicator_withHETbySex(:,Params.popHETbySex_HETF) = Params.popIndicator(:,Params.pop_HET).*Params.sexIndicator(:,Params.sex_Female);
    Params.popIndicator_withHETbySex(:,Params.popHETbySex_MSM) = Params.popIndicator(:,Params.pop_MSM);
    Params.popIndicator_withHETbySex(:,Params.popHETbySex_IDU) = Params.popIndicator(:,Params.pop_IDU);
    
% INDICATOR: RISK LEVEL

    % Fill indicator
        % Main
        storeMain= [7;14;7;14;14;7;14;7;14;14;7;14;7;14;1];
        count=1;
        for a = 1:15
            Params.riskLevelIndicator(count:count+6,1)=1;
            count = count+storeMain(a)+7;
        end
        
        % Casual
        Params.riskLevelIndicator(:,2)=abs(1-Params.riskLevelIndicator(:,1));

   % Corresponding column numbers
        Params.risk_Main = 1;
        Params.risk_Casual = 2;

        

        
% INDICATOR: CONTINUUM OF CARE RATE PERIODS
        % Doesn't have an indicator array because typically periods are
        % applied using separate variables
        Params.period_1 = 1;
        Params.period_2to4 = 2;       
        Params.period_5 = 3;
        
% INDICATOR: COVID PERIOD
        Params.COVIDperiod_1 = 1; % model time period 2
        Params.COVIDperiod_2 = 2; % model time period 3
        Params.COVIDperiod_3 = 3; % model time period 4
        
% INDICATOR: POP / SEX
    
    % Column numbers
        Params.popSex_HETM = 1;
        Params.popSex_HETF = 2;
        Params.popSex_MSM = 3;
        Params.popSex_IDUM = 4;
        Params.popSex_IDUF = 5;

    % Fill indicator
        Params.popSexIndicator(:,Params.popSex_HETM) = Params.popIndicator(:,Params.pop_HET).*Params.sexIndicator(:,Params.sex_Male);
        Params.popSexIndicator(:,Params.popSex_HETF) = Params.popIndicator(:,Params.pop_HET).*Params.sexIndicator(:,Params.sex_Female);
        Params.popSexIndicator(:,Params.popSex_MSM) = Params.popIndicator(:,Params.pop_MSM).*Params.sexIndicator(:,Params.sex_Male);
        Params.popSexIndicator(:,Params.popSex_IDUM) = Params.popIndicator(:,Params.pop_IDU).*Params.sexIndicator(:,Params.sex_Male);
        Params.popSexIndicator(:,Params.popSex_IDUF) = Params.popIndicator(:,Params.pop_IDU).*Params.sexIndicator(:,Params.sex_Female);
    

    
% INDICATOR: HET BY RISK LEVEL
    Params.HRHIndicator = Params.riskLevelIndicator(:,Params.risk_Casual) .* Params.popIndicator(:,Params.pop_HET);
    Params.LowRiskHETIndicator = Params.riskLevelIndicator(:,Params.risk_Main) .* Params.popIndicator(:,Params.pop_HET);
    Params.HighRiskMSMIndicator = Params.riskLevelIndicator(:,Params.risk_Casual) .* Params.popIndicator(:,Params.pop_MSM);
    
% INDICATOR: YOUNG MSM
    Params.YoungMSMIndicator = sum(Params.ageIndicator(:,[Params.age_1, Params.age_2,Params.age_3]),2) .* Params.popIndicator(:,Params.pop_MSM);
    Params.nonYMSMIndicator = 1 - Params.YoungMSMIndicator;
    
% PrEP
    % Note: these are defined as 1 and 2 (rather than 0 and 1) to support
    % the PrEP for loop in CalcBetas.m
    
    Params.ind_NoPrEP = 1;
    Params.ind_PrEP = 2;
    
% ORAL AND INJECTABLE PrEP  % JCPrEPUpdate: new indicator used for
% CalcInfectionRates

    Params.ind_OralPrEP_High = 1;
    Params.ind_OralPrEP_Low = 2;
    Params.ind_InjectPrEP_High = 3;
    Params.ind_InjectPrEP_Low = 4;

% SEP
    Params.ind_NoSEP = 1;
    Params.ind_SEP = 2;
    
% DISEASE STAGES
    Params.stage_Acute = 1; % B
    Params.stage_LatentA = 2; % C
    Params.stage_LatentB = 3; % D
    Params.stage_Late = 4; % E
    Params.stage_AIDS = 5; % F

% Model's Continuum of care stages if chronic HIV
        % Note: Chronic is designated because Acute only has 3 continuum stages
    Params.chronicContStage_Unaware = 1;
    Params.chronicContStage_AwareNoART = 2;
    Params.chronicContStage_LTCNoART = 3;
    Params.chronicContStage_ANV = 4;
    Params.chronicContStage_VLS = 5;

% Distinct Continuum stages
        % Used for costing
    Params.distinctContStage_Unaware = 1;
    Params.distinctContStage_Aware = 2;
    Params.distinctContStage_LTC = 3;
    Params.distinctContStage_ART = 4;
    
% INDICATOR: FORCE OF INFECTION    
    Params.force_Vonly = 1;
    Params.force_Aonly = 2;
    Params.force_VsomeA = 3;
    Params.force_AsomeV = 4;
    Params.force_N = 5;


% INDICATOR: AWARE COMPARTMENTS
    % Whether the infected person knows he is infected
    % Binary: yes is 1, no is 0
    % [273x30]
    
    inf_awareComparts_273x273x30=zeros(Params.numStrats,Params.numStrats,Params.numComparts);
    inf_awareComparts_273x273x30(:,:,AwareComparts)=1;
    
    
% INDICATOR: VLS COMPARTMENTS
   
    Params.inf_vls = zeros(Params.numStrats,Params.numStrats,Params.numComparts);
    Params.inf_vls(:,:,Params.VLSComparts) = 1; 
    
% INDICATOR: Transmission group/sex combination for sexual mixing
    % Used to bring in sexual mixing parameters by trans group/sex
    % Updated by KH on 20Sep2023 to incorporate HET mixing varying by risk
    % level
    Params.sexMix_LRHETM_HETF = 1;
    Params.sexMix_LRHETM_IDUF = 2;
    Params.sexMix_LRHETF_HETM = 3;
    Params.sexMix_LRHETF_MSM = 4;
    Params.sexMix_LRHETF_IDUM = 5;
    Params.sexMix_HRHETM_HETF = 6;
    Params.sexMix_HRHETM_IDUF = 7;
    Params.sexMix_HRHETF_HETM = 8;
    Params.sexMix_HRHETF_MSM = 9;
    Params.sexMix_HRHETF_IDUM = 10;
    Params.sexMix_MSM_HETF = 11;
    Params.sexMix_MSM_MSM = 12;
    Params.sexMix_MSM_IDUF = 13;
    Params.sexMix_IDUM_HETF = 14;
    Params.sexMix_IDUM_IDUF = 15;
    Params.sexMix_IDUF_HETM = 16;
    Params.sexMix_IDUF_MSM = 17;
    Params.sexMix_IDUF_IDUM = 18;

    
% INDICATOR: Transmission group/sex combination for needle mixing
    % Used to import needle mixing parameters
    
    Params.needleMix_IDUM_IDUM = 1;
    Params.needleMix_IDUM_IDUF = 2;
    Params.needleMix_IDUF_IDUM = 3;
    Params.needleMix_IDUF_IDUF = 4;

    end

    % 2.iv. Define mixing indices
    for mixingIndexSection = 1:1
        
        for indicesPopSex = 1:1
            
       % Determine indices for PopSex mixing
        % Revised by KH on 20Sep2023 to incorporate HET mixing by risk
        % level to replace Params.idx_HETM_[] variables with 
        % Params.idx_LRHETM_[] and Params.idx_HRHETM_[] variables  

          
        % Index for HETM-HETF
            Params.idx_LRHETM_HETF = ones(Params.numStrats,Params.numStrats);
            Params.idx_HRHETM_HETF = ones(Params.numStrats,Params.numStrats);

            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_LRHETM_HETF(:,j) = ...
                   Params.idx_LRHETM_HETF(:,j) .* Params.popSexIndicator(:,Params.popSex_HETM) .* Params.riskLevelIndicator(:,Params.risk_Main);
               Params.idx_HRHETM_HETF(:,j) = ...
                   Params.idx_HRHETM_HETF(:,j) .* Params.popSexIndicator(:,Params.popSex_HETM) .* Params.riskLevelIndicator(:,Params.risk_Casual);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_LRHETM_HETF(i,:) = ...
                   Params.idx_LRHETM_HETF(i,:) .* Params.popSexIndicator(:,Params.popSex_HETF)'; 
               Params.idx_HRHETM_HETF(i,:) = ...
                   Params.idx_HRHETM_HETF(i,:) .* Params.popSexIndicator(:,Params.popSex_HETF)'; 
            end

         % Index for HETM-IDUF
            Params.idx_LRHETM_IDUF = ones(Params.numStrats,Params.numStrats);
            Params.idx_HRHETM_IDUF = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_LRHETM_IDUF(:,j) = ...
                   Params.idx_LRHETM_IDUF(:,j) .* Params.popSexIndicator(:,Params.popSex_HETM) .* Params.riskLevelIndicator(:,Params.risk_Main);
               Params.idx_HRHETM_IDUF(:,j) = ...
                   Params.idx_HRHETM_IDUF(:,j) .* Params.popSexIndicator(:,Params.popSex_HETM) .* Params.riskLevelIndicator(:,Params.risk_Casual);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_LRHETM_IDUF(i,:) = ...
                   Params.idx_LRHETM_IDUF(i,:) .* Params.popSexIndicator(:,Params.popSex_IDUF)'; 
               Params.idx_HRHETM_IDUF(i,:) = ...
                   Params.idx_HRHETM_IDUF(i,:) .* Params.popSexIndicator(:,Params.popSex_IDUF)'; 
            end

        % Index for HETF-HETM
            Params.idx_LRHETF_HETM = ones(Params.numStrats,Params.numStrats);
            Params.idx_HRHETF_HETM = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_LRHETF_HETM(:,j) = ...
                   Params.idx_LRHETF_HETM(:,j) .* Params.popSexIndicator(:,Params.popSex_HETF) .* Params.riskLevelIndicator(:,Params.risk_Main);
               Params.idx_HRHETF_HETM(:,j) = ...
                   Params.idx_HRHETF_HETM(:,j) .* Params.popSexIndicator(:,Params.popSex_HETF) .* Params.riskLevelIndicator(:,Params.risk_Casual);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_LRHETF_HETM(i,:) = ...
                   Params.idx_LRHETF_HETM(i,:) .* Params.popSexIndicator(:,Params.popSex_HETM)'; 
               Params.idx_HRHETF_HETM(i,:) = ...
                   Params.idx_HRHETF_HETM(i,:) .* Params.popSexIndicator(:,Params.popSex_HETM)'; 
            end

            % Index for HETF-MSM
            Params.idx_LRHETF_MSM = ones(Params.numStrats,Params.numStrats);
            Params.idx_HRHETF_MSM = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_LRHETF_MSM(:,j) = ...
                   Params.idx_LRHETF_MSM(:,j) .* Params.popSexIndicator(:,Params.popSex_HETF) .* Params.riskLevelIndicator(:,Params.risk_Main);
               Params.idx_HRHETF_MSM(:,j) = ...
                   Params.idx_HRHETF_MSM(:,j) .* Params.popSexIndicator(:,Params.popSex_HETF) .* Params.riskLevelIndicator(:,Params.risk_Casual);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_LRHETF_MSM(i,:) = Params.idx_LRHETF_MSM(i,:) .* Params.popSexIndicator(:,Params.popSex_MSM)'; 
               Params.idx_HRHETF_MSM(i,:) = Params.idx_HRHETF_MSM(i,:) .* Params.popSexIndicator(:,Params.popSex_MSM)'; 
            end

            % Index for HETF-IDUM
            Params.idx_LRHETF_IDUM = ones(Params.numStrats,Params.numStrats);
            Params.idx_HRHETF_IDUM = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_LRHETF_IDUM(:,j) = ...
                   Params.idx_LRHETF_IDUM(:,j) .* Params.popSexIndicator(:,Params.popSex_HETF) .* Params.riskLevelIndicator(:,Params.risk_Main);
               Params.idx_HRHETF_IDUM(:,j) = ...
                   Params.idx_HRHETF_IDUM(:,j) .* Params.popSexIndicator(:,Params.popSex_HETF) .* Params.riskLevelIndicator(:,Params.risk_Casual);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_LRHETF_IDUM(i,:) = Params.idx_LRHETF_IDUM(i,:) .* Params.popSexIndicator(:,Params.popSex_IDUM)'; 
               Params.idx_HRHETF_IDUM(i,:) = Params.idx_HRHETF_IDUM(i,:) .* Params.popSexIndicator(:,Params.popSex_IDUM)'; 
            end

            % Index for MSM-HETF
            Params.idx_MSM_HETF = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_MSM_HETF(:,j) = Params.idx_MSM_HETF(:,j) .* Params.popSexIndicator(:,3);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_MSM_HETF(i,:) = Params.idx_MSM_HETF(i,:) .* Params.popSexIndicator(:,2)'; 
            end

            % Index for MSM-MSM
            Params.idx_MSM_MSM = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_MSM_MSM(:,j) = Params.idx_MSM_MSM(:,j) .* Params.popSexIndicator(:,3);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_MSM_MSM(i,:) = Params.idx_MSM_MSM(i,:) .* Params.popSexIndicator(:,3)'; 
            end

            % Index for MSM-IDUF
            Params.idx_MSM_IDUF = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_MSM_IDUF(:,j) = Params.idx_MSM_IDUF(:,j) .* Params.popSexIndicator(:,3);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_MSM_IDUF(i,:) = Params.idx_MSM_IDUF(i,:) .* Params.popSexIndicator(:,5)'; 
            end

            % Index for IDUM-HETF
            Params.idx_IDUM_HETF = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_IDUM_HETF(:,j) = Params.idx_IDUM_HETF(:,j) .* Params.popSexIndicator(:,4);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_IDUM_HETF(i,:) = Params.idx_IDUM_HETF(i,:) .* Params.popSexIndicator(:,2)'; 
            end

            % Index for IDUM-IDUM
            Params.idx_IDUM_IDUM = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_IDUM_IDUM(:,j) = Params.idx_IDUM_IDUM(:,j) .* Params.popSexIndicator(:,4);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_IDUM_IDUM(i,:) = Params.idx_IDUM_IDUM(i,:) .* Params.popSexIndicator(:,4)'; 
            end  

              % Index for IDUM-IDUF
            Params.idx_IDUM_IDUF = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_IDUM_IDUF(:,j) = Params.idx_IDUM_IDUF(:,j) .* Params.popSexIndicator(:,4);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_IDUM_IDUF(i,:) = Params.idx_IDUM_IDUF(i,:) .* Params.popSexIndicator(:,5)'; 
            end

            % Index for IDUF-HETM
            Params.idx_IDUF_HETM = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_IDUF_HETM(:,j) = Params.idx_IDUF_HETM(:,j) .* Params.popSexIndicator(:,5);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_IDUF_HETM(i,:) = Params.idx_IDUF_HETM(i,:) .* Params.popSexIndicator(:,1)'; 
            end

            % Index for IDUF-MSM
            Params.idx_IDUF_MSM = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_IDUF_MSM(:,j) = Params.idx_IDUF_MSM(:,j) .* Params.popSexIndicator(:,5);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_IDUF_MSM(i,:) = Params.idx_IDUF_MSM(i,:) .* Params.popSexIndicator(:,3)'; 
            end

            % Index for IDUF-IDUM
            Params.idx_IDUF_IDUM = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_IDUF_IDUM(:,j) = Params.idx_IDUF_IDUM(:,j) .* Params.popSexIndicator(:,5);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_IDUF_IDUM(i,:) = Params.idx_IDUF_IDUM(i,:) .* Params.popSexIndicator(:,4)'; 
            end  

            % Index for IDUF-IDUF
            Params.idx_IDUF_IDUF = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_IDUF_IDUF(:,j) = Params.idx_IDUF_IDUF(:,j) .* Params.popSexIndicator(:,5);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_IDUF_IDUF(i,:) = Params.idx_IDUF_IDUF(i,:) .* Params.popSexIndicator(:,5)'; 
            end  


            % The index below is used for needle mixing

                    % Index for IDU-IDU
                    Params.idx_IDU_IDU = ones(Params.numStrats,Params.numStrats);
                    % 1st partner
                    for j = 1:Params.numStrats
                       Params.idx_IDU_IDU(:,j) = Params.idx_IDU_IDU(:,j) .* Params.popIndicator(:,Params.pop_IDU);
                    end
                    % 2nd partner
                    for i = 1:Params.numStrats
                       Params.idx_IDU_IDU(i,:) = Params.idx_IDU_IDU(i,:) .* Params.popIndicator(:,Params.pop_IDU)'; 
                    end  

            % Male-Female
                Params.idx_Male_Female = ones(Params.numStrats,Params.numStrats);
                % 1st partner
                for j = 1:Params.numStrats
                   Params.idx_Male_Female(:,j) = Params.idx_Male_Female(:,j) .* Params.sexIndicator(:,Params.sex_Male);
                end
                % 2nd partner
                for i = 1:Params.numStrats
                   Params.idx_Male_Female(i,:) = Params.idx_Male_Female(i,:) .* Params.sexIndicator(:,Params.sex_Female)'; 
                end     
                
            % Female-Male
                Params.idx_Female_Male = ones(Params.numStrats,Params.numStrats);
                % 1st partner
                for j = 1:Params.numStrats
                   Params.idx_Female_Male(:,j) = Params.idx_Female_Male(:,j) .* Params.sexIndicator(:,Params.sex_Female);
                end
                % 2nd partner
                for i = 1:Params.numStrats
                   Params.idx_Female_Male(i,:) = Params.idx_Female_Male(i,:) .* Params.sexIndicator(:,Params.sex_Male)'; 
                end      


            % Params.idx_MSMandMSM and Params.idx_MSM_MSM look exactly the
            % same, but they are different variable types
            % MSMandMSM is (similar to?) a boolean yes/no type array so
            % it can be applied to referencing and MSM_MSM is a normal double type
            % so it can be multiplied with other parameters
            
            % They're both used in this m file and Calib_UpdateParams.m
            Params.idx_MSMandMSM = Params.idx_MSM_MSM == 1;
           
            
      
            
        end

        for indicesRiskLevel = 1:1

            % in the riskLevelIndicator - main is column1

            Params.idx_low_low = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_low_low(:,j) = Params.idx_low_low(:,j) .* Params.riskLevelIndicator(:,1);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_low_low(i,:) = Params.idx_low_low(i,:) .* Params.riskLevelIndicator(:,1)'; 
            end

            Params.idx_low_high = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_low_high(:,j) = Params.idx_low_high(:,j) .* Params.riskLevelIndicator(:,1);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_low_high(i,:) = Params.idx_low_high(i,:) .* Params.riskLevelIndicator(:,2)'; 
            end

            Params.idx_high_low = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_high_low(:,j) = Params.idx_high_low(:,j) .* Params.riskLevelIndicator(:,2);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_high_low(i,:) = Params.idx_high_low(i,:) .* Params.riskLevelIndicator(:,1)'; 
            end

            Params.idx_high_high = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_high_high(:,j) = Params.idx_high_high(:,j) .* Params.riskLevelIndicator(:,2);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_high_high(i,:) = Params.idx_high_high(i,:) .* Params.riskLevelIndicator(:,2)'; 
            end       

        end

        for indicesRaceMixing = 1:1

            % Index for Black-Black
            Params.idx_black_black = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_black_black(:,j) = Params.idx_black_black(:,j) .* Params.raceIndicator(:,1);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_black_black(i,:) = Params.idx_black_black(i,:) .* Params.raceIndicator(:,1)'; 
            end

            % Index for Black-Hispanic
            Params.idx_black_hispanic = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_black_hispanic(:,j) = Params.idx_black_hispanic(:,j) .* Params.raceIndicator(:,1);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_black_hispanic(i,:) = Params.idx_black_hispanic(i,:) .* Params.raceIndicator(:,2)'; 
            end

            % Index for Black-Other
            Params.idx_black_other = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_black_other(:,j) = Params.idx_black_other(:,j) .* Params.raceIndicator(:,1);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_black_other(i,:) = Params.idx_black_other(i,:) .* Params.raceIndicator(:,3)'; 
            end

            % Index for Hispanic-black
            Params.idx_hispanic_black = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_hispanic_black(:,j) = Params.idx_hispanic_black(:,j) .* Params.raceIndicator(:,2);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_hispanic_black(i,:) = Params.idx_hispanic_black(i,:) .* Params.raceIndicator(:,1)'; 
            end

            % Index for Hispanic-hispanic
            Params.idx_hispanic_hispanic = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_hispanic_hispanic(:,j) = Params.idx_hispanic_hispanic(:,j) .* Params.raceIndicator(:,2);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_hispanic_hispanic(i,:) = Params.idx_hispanic_hispanic(i,:) .* Params.raceIndicator(:,2)'; 
            end

            % Index for Hispanic-other
            Params.idx_hispanic_other = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_hispanic_other(:,j) = Params.idx_hispanic_other(:,j) .* Params.raceIndicator(:,2);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_hispanic_other(i,:) = Params.idx_hispanic_other(i,:) .* Params.raceIndicator(:,3)'; 
            end

            % Index for Other-black
            Params.idx_other_black = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_other_black(:,j) = Params.idx_other_black(:,j) .* Params.raceIndicator(:,3);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_other_black(i,:) = Params.idx_other_black(i,:) .* Params.raceIndicator(:,1)'; 
            end

            % Index for Other-hispanic
            Params.idx_other_hispanic = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_other_hispanic(:,j) = Params.idx_other_hispanic(:,j) .* Params.raceIndicator(:,3);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_other_hispanic(i,:) = Params.idx_other_hispanic(i,:) .* Params.raceIndicator(:,2)'; 
            end

            % Index for Other-other
            Params.idx_other_other = ones(Params.numStrats,Params.numStrats);
            % 1st partner
            for j = 1:Params.numStrats
               Params.idx_other_other(:,j) = Params.idx_other_other(:,j) .* Params.raceIndicator(:,3);
            end
            % 2nd partner
            for i = 1:Params.numStrats
               Params.idx_other_other(i,:) = Params.idx_other_other(i,:) .* Params.raceIndicator(:,3)'; 
            end

        end
    end
    
    % 2.v. Misc
    for miscSection = 1:1
       
        % for pulling the best solution vector after a calibration by opt
        Params.calib_num_outputVectorErrorValues = 5;
        
    end

%% 3. Import data from Excel to Matlab

    % Define path to Excel file
    %[folder, name, ext] = fileparts(which(mfilename));
    %pathName = [folder  '\'  ExcelFileName];

    %Import the parameters
    %ExcelValues_AllParameters = xlsread(pathName, 'ParameterList','matlab_ParameterList');
    %ExcelValues_Populations = xlsread(pathName, 'Import_Populations','matlab_InitAndNewPops');
    %save('ExcelValues_AllParameters','ExcelValues_AllParameters');
    %save('ExcelValues_Populations','ExcelValues_Populations');

%     ExcelValues_AllParameters = xlsread(pathName, 'ParameterList','L388:L1787');
%     ExcelValues_Populations = xlsread(pathName, 'Import_Populations','C48:C6872');
%     ExcelValues_AllParametersStruct = load('ExcelValues_AllParameters.mat');
%     ExcelValues_PopulationsStruct = load('ExcelValues_Populations.mat');
%     ExcelValues_AllParameters = ExcelValues_AllParametersStruct.ExcelValues_AllParameters;
%     ExcelValues_Populations = ExcelValues_PopulationsStruct.ExcelValues_Populations;
    
%% 4. Define Model Settings

    % Parameter Index keeps track of the variables from the parameter list
    parameter_Index = 1;

    
    % ONE-WAY SA OUTCOME
    % These are defined on the 1-way SA results page, model mechanics worksheet
    % and in the HIVEpi model main m file.

    Params.outcomeNumber = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.outcomeNumber);
    
    
    % PROGRESSION SOURCE
    % Determines how progression is applied
    
        % 1 - [normal] - direct entry of probabilities
        % 2 - [test frequency] - testing frequency values for HETs; direct
        %       entry for all other progression and transmission groups
        % 3 - [allocation-based progression] - allocation-based progression
        
    Params.tt_progressionSource = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_progressionSource);

    % TYPE OF SOLVER
    % Discrete or continuous
        % Discrete = 0
        % Continuous = 1
    Params.SolveCont = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.SolveCont);

    % TYPE OF MODEL RUN
    % Epi Model
    % Calibration
    % Optimization
    % One-way SA
    
    Params.ModelRunType = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.ModelRunType);

    % NEW LINES - OUTCOMES TO COLLECT

    Params.OutcomesToCollect = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.ModelRunType);   
    
    % YEAR MODEL STARTS
    % [1x1]
    
    Params.tt_modelStartYear = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_modelStartYear);
    
    % TIME HORIZON
    % in years
    
    Params.tt_timeHorizon = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_timeHorizon);

    % TIME STEP / TIME SPAN
    % Model time step in years if running discrete DEs
    % Time step = 1 if running continuous DEs
    if Params.SolveCont == 0
        Params.tt_tstep=ExcelValues_AllParameters(parameter_Index);
    else
        Params.tt_tstep=1;
    end
    parameter_Index=parameter_Index+numel(Params.tt_tstep);

    % PERIOD 2 START YEAR    
    % Value is the year at which the second set of rates (e.g., testing,
    % ART initiation) begin
    
    Params.tt_periodTwoStartYear = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_periodTwoStartYear);

    % PERIOD 3 START YEAR    
    % Value is the year at which the third set of rates (e.g., testing,
    % ART initiation) begin   
    
    Params.tt_periodThreeStartYear = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_periodThreeStartYear);
    
    % PERIOD 4 START YEAR - Added by MClinkscales on 5/12/2022    
    % Value is the year at which the third set of rates (e.g., testing,
    % ART initiation) begin   
    
    Params.tt_periodFourStartYear = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_periodFourStartYear);
    
    % PERIOD 5 START YEAR - Added by MClinkscales on 5/12/2022      
    % Value is the year at which the third set of rates (e.g., testing,
    % ART initiation) begin   
    
    Params.tt_periodFiveStartYear = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_periodFiveStartYear);

    % OUTCOME COLLECTION: START YEAR
    Params.outcomeCollectionStartYr = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;  
    
    % OUTCOME COLLECTION: END YEAR
    Params.outcomeCollectionEndYr = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;  

    % First and last outcome years
    FirstOutcomeYr = max(Params.tt_modelStartYear,Params.outcomeCollectionStartYr);
    LastOutcomeYr = min(Params.tt_modelStartYear+Params.tt_timeHorizon-1,Params.outcomeCollectionEndYr);

    % CALCULATED: INDEX OF KEY MODEL YEARS 
    % Value is how many years from start of the model to particular year

    %Eg: given start year of 2000, the 2006 year index will be:
        % 2007 - 2000 = index is 7
    % Given start year of 2006, the year 2006 will be from index #1
        % 2007 - 2006 = 1
    
    
    % 2006
    Params.year_2006 = max(2007 - FirstOutcomeYr,0);
    
    % 2009
    Params.year_2009 = max(2010 - FirstOutcomeYr,0);

    % 2010
    Params.year_2010 = max(2011 - FirstOutcomeYr,0);
    
    % 2012
    Params.year_2012 = max(2013 - FirstOutcomeYr,0);
    
    % 2013
    Params.year_2013 = max(2014 - FirstOutcomeYr,0);
    
    %2014
    Params.year_2014 = max(2015 - FirstOutcomeYr,0); 

    % 2015
    Params.year_2015 = max(2016 - FirstOutcomeYr,0);
    
    % 2016
    Params.year_2016 = max(2017 - FirstOutcomeYr,0);
    
    % 2017
    Params.year_2017 = max(2018 - FirstOutcomeYr,0);
    
    % 2018
    Params.year_2018 = max(2019 - FirstOutcomeYr,0);
    
    %2019 - Added 6/30/2021 Clinkscales
    Params.year_2019 = max(2020 - FirstOutcomeYr,0);
    
    % 2020
    Params.year_2020 = max(2021 - FirstOutcomeYr,0);

    % 2021 - Added 02/10/2023. Bates
    Params.year_2021 = max(2022 - FirstOutcomeYr,0);

    % YEAR TO COLLECT TEST AND DIAG RATES
    Params.diagRateYr = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
    % MINIMUM AGE INCLUDED
    Params.MinAgeIncluded = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;

    % DOLLAR YEAR
    Params.DollarYear = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;

    % CALCULATE PREP-SPECIFIC BETAS
    % 1 = no
    % 2 = yes
    
    Params.CalcPrEPSpecificBetas = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
   
    % DISCOUNT FACTORS
 
        % QALYS
        Params.Disc_QALYs = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(Params.Disc_QALYs);

        % Life years
        Params.Disc_LifeYears = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(Params.Disc_LifeYears);

        % Costs
        Params.Disc_Costs = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(Params.Disc_Costs);

        % New infections
        Params.alloc_discRateNewInfections = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(Params.alloc_discRateNewInfections); 

    % Target interventions only to young MSM?
    % 0 = no
    % 1 = yes
    
    Params.TargetIntnsToYMSM = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
    % Apply COVID-19 effects?
    % 0 = no
    % 1 = yes
    
    Params.ApplyCOVIDEffects_ContBeh = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;

    Params.ApplyCOVIDEffects_Mort = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;

    % Collect outcomes for infections by source of infection?
    % 0 = no
    % 1 = yes

    Params.CollectInfbySource = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;

    % Are input values for the second time period being calibrated?
    % 0 = no
    % 1 = yes
    
    Params.Calib2ndPeriod = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;

%% 5. Define NHAS Optimization Settings

    % OUTCOME TO TARGET
    Params.nhas_TargetOutcome = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.nhas_TargetOutcome);

    % TARGET VALUE FOR OUTCOME
    Params.nhas_TargetOutcomeVal = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.nhas_TargetOutcomeVal);

    % INPUT TO BE VARIED
    Params.nhas_InputToVary = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.nhas_InputToVary);

    % CURRENT VALUE OF SELECTED INPUT
    Params.nhas_InputToVary_Val = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.nhas_InputToVary_Val);

    % RANGE OF VALUES TO CONSIDER FOR SELECTED INPUT
    % Lower Bound
    Params.nhas_InputToVary_LB = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.nhas_InputToVary_LB);
    
    % Upper Bound
    Params.nhas_InputToVary_UB = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.nhas_InputToVary_UB);        
    
%% 6. Define Population inputs
    
    % Percent of MSM that are high-risk (for calibration) (added by JC on
    %12/05/2017)
    % Read in from Excel
    Params.MSM_PercentHR = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.MSM_PercentHR);

    %MSM population sizes by race (for calibration) (added by JC on
    %12/05/2017)
    % Read in from Excel
    Params.MSMInitPopSize_Race(:,1)=ExcelValues_AllParameters(parameter_Index:parameter_Index+Params.numRace-1);
    parameter_Index = parameter_Index + numel(Params.MSMInitPopSize_Race);
    
    % Initial HIV prevalence by race (for calibration) (added by JC on
    %08/06/2018)
    % Read in from Excel
    Params.InitHIVPrevalence_Race(:,1)=ExcelValues_AllParameters(parameter_Index:parameter_Index+Params.numRace-1);
    parameter_Index = parameter_Index + numel(Params.InitHIVPrevalence_Race);

    %Growth Rate by Transmission Group / Sex / Race
    Params.growthRate = zeros(numPop,2,numRace);
    
    % triple loop creates array that stores growth rate.
    counter_loop = 0;
    counter_index = 0;
    for i = 1:Params.numPop
        for j = 1:2
            for k = 1:Params.numRace
                if i == Params.pop_MSM && j == Params.sex_Female %this captures the non-existent Female MSMs
                    counter_loop = counter_loop + 1;
                    Params.growthRate(counter_loop) = 0;
                else
                    counter_loop = counter_loop + 1;
                    Params.growthRate(counter_loop) = ExcelValues_AllParameters(parameter_Index + counter_index);
                    counter_index = counter_index + 1;
                end
            end
        end
    end
    
    parameter_Index = parameter_Index + counter_index;
    
    % Percent of entry rate applied to each age group
    Params.entryRateApplied = zeros(Params.numAge,1);
    
    for i = 1:Params.numAge
        Params.entryRateApplied(i) = ExcelValues_AllParameters(parameter_Index+i-1);
    end
    
    parameter_Index = parameter_Index + numel(Params.entryRateApplied);
    
    % PERCENT CIRCUMCISED
    % by race [3x1]
    Params.pop_pctCirc_r = ExcelValues_AllParameters(parameter_Index:parameter_Index+numRace-1);
    parameter_Index = parameter_Index + numel(Params.pop_pctCirc_r);

    
    % ALL-CAUSE MORTALITY
    % Annual risk of non-HIV related death

    % Read in from Excel
        pop_dyingRate_psa(numPop,numSex,Params.numAge) = 0; 

        %Loop through the pops (only HET and MSM)
        for p = 1:numPop - 1
            %Loop through the sexes
            for s = 1:numSex
            pop_dyingRate_psa(p,s,:)=ExcelValues_AllParameters(parameter_Index + Params.numAge*(s-1)+Params.numAge*numSex*(p-1): ...
                (parameter_Index + Params.numAge-1)+Params.numAge*(s-1)+(p-1)*Params.numAge*numSex);

            end
        end
        parameter_Index = parameter_Index + Params.numAge*numSex*(numPop-1);

        %IDU only
        for s = 1:numSex
           pop_dyingRate_psa(numPop,s,:)=ExcelValues_AllParameters(parameter_Index+Params.numAge*(s-1):...
               (parameter_Index+Params.numAge-1)+Params.numAge*(s-1));
        end
        parameter_Index = parameter_Index + Params.numAge*numSex;

    % Adjust stratification to be compatible with model's array dimensions
    
        %Pre-allocate
        store = zeros(Params.numStrats,numSex,Params.numAge);
        storParams.E2 = zeros(Params.numStrats,Params.numStrats,Params.numAge);
        storParams.E3= zeros(Params.numStrats,Params.numAge);

        for a = 1:Params.numAge
            store(:,:,a) = Params.popIndicator*pop_dyingRate_psa(:,:,a);   
        end
        for a = 1:Params.numAge
            storParams.E2(:,:,a)=store(:,:,a)*Params.sexIndicator';
        end
        for a = 1:Params.numAge
            storParams.E3(:,a)= diag(storParams.E2(:,:,a));
        end
        transAge = transpose(Params.ageIndicator);

        pop_dyingRate = diag(storParams.E3*transAge);
        Params.pop_dyingRate = ProbToRate(pop_dyingRate);

%% 7. Define Continuum-of-care inputs

    % HET FREQUENCY: POPULATION (IF RELEVANT) 
    % Population to apply intervals to
    applyIntervalAllHET = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(applyIntervalAllHET);
    
   if applyIntervalAllHET == 1
        Params.tt_HETFreq_ApplyLRH = 1;
    else
        Params.tt_HETFreq_ApplyLRH = 0;
    end
      
    % Average intervals
    Params.tt_HETFreq_Intervals_l(:,1)=ExcelValues_AllParameters(parameter_Index:parameter_Index+1);
    parameter_Index = parameter_Index + numel(Params.tt_HETFreq_Intervals_l);
    
    
    % HET FREQUENCY: TESTING RATES

    % The values from Excel are coming in as rates

    % Low Risk HET
        % Read in from Excel
        storeLowRiskRate(:,1)=ExcelValues_AllParameters(parameter_Index:parameter_Index+Params.numAge-1);
        parameter_Index = parameter_Index + numel(storeLowRiskRate);

        % Apply to Params struct
        Params.tt_HETFreq_TestRate_LRH = Params.popIndicator(:,Params.pop_HET) .* ...
            (Params.ageIndicator * storeLowRiskRate .* Params.riskLevelIndicator(:,Params.risk_Main));
    
    % High risk HET
        % Read in from Excel
        storeHighRiskRate(:,1)=ExcelValues_AllParameters(parameter_Index:parameter_Index+Params.numAge-1);
        parameter_Index = parameter_Index + numel(storeHighRiskRate);
        
        % Apply to Params struct
        Params.tt_HETFreq_TestRate_HRH = Params.popIndicator(:,Params.pop_HET) .* ...
            (Params.ageIndicator * storeHighRiskRate .* Params.riskLevelIndicator(:,Params.risk_Casual));

        
    % RELATIVE RISKS OF BEING TESTED
    % Period 1 relative risks

        % Note: the factors are also included under the "Params" struct because the reference
        % cases for each are used in the "Calib_updateParams" m file.
            % E.g., Params.tt_relRiskPop_1 is used in Calib_updateParams 
                % The first element (Parmams.tt_relRiskPop_1(Params.pop_HET)
                % is used but the values for MSM and IDU relative risks are 
                % replaced with the values from the calibration set when the 
                % test rates are calculated for that calibration run.

        % Reference case
            % Black, HET, CD4 >500
        Params.tt_testRefCase_1= ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(Params.tt_testRefCase_1);

        % Relative risk by transmission group
            %HET, MSM, IDU
        Params.tt_relRiskPop_1(numPop,1)=0;
        Params.tt_relRiskPop_1(:,1)=ExcelValues_AllParameters(parameter_Index:...
            parameter_Index+numPop-1);
        parameter_Index = parameter_Index + numel(Params.tt_relRiskPop_1);

        % Relative risk by race/eth
            %Black, Hispanic, Other
        Params.tt_relRiskRace_1(numRace,1)=0;
        Params.tt_relRiskRace_1(:,1)=ExcelValues_AllParameters(parameter_Index:...
            parameter_Index+numRace-1);
        parameter_Index = parameter_Index + numel(Params.tt_relRiskRace_1);

        % Relative risk by disease stage
            % Acute, LatentA, LatentB, Late, AIDS, Uninfected
                %Note: reference case is LatentA. It is not the first element
                %so that it remains consistent with the order of the
                %previously-defined indicators (e.g., Params.stage_LatentA = 2)         
        Params.tt_relRiskHIVstage_1(numHIVstages+1,1)=0;
        Params.tt_relRiskHIVstage_1(:,1)=ExcelValues_AllParameters(parameter_Index:...
            parameter_Index+numHIVstages+1-1);
        parameter_Index = parameter_Index + numel(Params.tt_relRiskHIVstage_1);
        
        % Relative risk by age group (added by JC on 11/13/2017)
            % 13-17, 18-24, 25-34, 35-44, 45-64         
        Params.tt_relRiskAge_1(Params.numAge,1)=0;
        Params.tt_relRiskAge_1(:,1)=ExcelValues_AllParameters(parameter_Index:...
            parameter_Index+Params.numAge-1);
        parameter_Index = parameter_Index + numel(Params.tt_relRiskAge_1);

    % Period 2 to 4 Relative Risks

        % Reference case
            % Black, HET, CD4 >500
        Params.tt_testRefCase_2to4= ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(Params.tt_testRefCase_2to4);

        % Relative risk by transmission group
            %HET, MSM, IDU
        Params.tt_relRiskPop_2to4(numPop,1)=0;
        Params.tt_relRiskPop_2to4(:,1)=ExcelValues_AllParameters(parameter_Index:...
            parameter_Index+numPop-1);
        parameter_Index = parameter_Index + numel(Params.tt_relRiskPop_2to4);

        % Relative risk by race/eth
            %Black,Hispanic,Other
        Params.tt_relRiskRace_2to4(numRace,1)=0;
        Params.tt_relRiskRace_2to4(:,1)=ExcelValues_AllParameters(parameter_Index:...
            parameter_Index+numRace-1);
        parameter_Index = parameter_Index + numel(Params.tt_relRiskRace_2to4);

        % Relative risk by disease stage
            % Acute, LatentA, LatentB, Late, AIDS, uninfected
                %Note: reference case is LatentA
        Params.tt_relRiskHIVstage_2to4(numHIVstages+1,1)=0;
        Params.tt_relRiskHIVstage_2to4(:,1)=ExcelValues_AllParameters(parameter_Index:...
            parameter_Index+numHIVstages+1-1);
        parameter_Index = parameter_Index + numel(Params.tt_relRiskHIVstage_2to4);
        
        % Relative risk by age group (added by JC on 11/13/2017)
            % 13-17, 18-24, 25-34, 35-44, 45-54, 55-64, 65+         
        Params.tt_relRiskAge_2to4(Params.numAge,1)=0;
        Params.tt_relRiskAge_2to4(:,1)=ExcelValues_AllParameters(parameter_Index:...
            parameter_Index+Params.numAge-1);
        parameter_Index = parameter_Index + numel(Params.tt_relRiskAge_2to4);

    %Period 5 Relative Risks

        % Reference case
            % Black, HET, CD4 >500
        Params.tt_testRefCase_5= ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(Params.tt_testRefCase_5);
        
        % Reference case for non-YMSM (only relevant for societal
        % allocation targeted to YMSM only scenarios)
            % Black, HET, CD4 >500
        Params.tt_testRefCase_5_nonYMSM= ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(Params.tt_testRefCase_5_nonYMSM);

        % Relative risk by transmission group
            %HET, MSM, IDU
        Params.tt_relRiskPop_5(numPop,1)=0;
        Params.tt_relRiskPop_5(:,1)=ExcelValues_AllParameters(parameter_Index:...
            parameter_Index+numPop-1);
        parameter_Index = parameter_Index + numel(Params.tt_relRiskPop_5);

        % Relative risk by race/eth
            %Black,Hispanic,Other
        Params.tt_relRiskRace_5(numRace,1)=0;
        Params.tt_relRiskRace_5(:,1)=ExcelValues_AllParameters(parameter_Index:...
            parameter_Index+numRace-1);
        parameter_Index = parameter_Index + numel(Params.tt_relRiskRace_5);

        % Relative risk by disease stage
            % Acute, LatentA, LatentB, Late, AIDS, Uninfected
                %Note: reference case is LatentA
        Params.tt_relRiskHIVstage_5(numHIVstages+1,1)=0;
        Params.tt_relRiskHIVstage_5(:,1)=ExcelValues_AllParameters(parameter_Index:...
            parameter_Index+numHIVstages+1-1);
        parameter_Index = parameter_Index + numel(Params.tt_relRiskHIVstage_5);
        
        % Relative risk by age group (added by JC on 11/13/2017)
            % 13-17, 18-24, 25-34, 35-44, 45-54, 55-64, 65+         
        Params.tt_relRiskAge_5(Params.numAge,1)=0;
        Params.tt_relRiskAge_5(:,1)=ExcelValues_AllParameters(parameter_Index:...
            parameter_Index+Params.numAge-1);
        parameter_Index = parameter_Index + numel(Params.tt_relRiskAge_5);

    % CALCULATE: ANNUAL TESTING RATES

    %Preallocate
    Params.tt_testRateAcute_rp_1(Params.numStrats)=1;
    Params.tt_testRateLatentA_rp_1(Params.numStrats)=1;
    Params.tt_testRateLatentB_rp_1(Params.numStrats)=1;
    Params.tt_testRateLate_rp_1(Params.numStrats)=1;
    Params.tt_testRateAIDS_rp_1(Params.numStrats)=1;
    Params.tt_testRateUninfected_rp_1(Params.numStrats)=1;

    Params.tt_testRateAcute_rp_2to4(Params.numStrats)=1;
    Params.tt_testRateLatentA_rp_2to4(Params.numStrats)=1;
    Params.tt_testRateLatentB_rp_2to4(Params.numStrats)=1;
    Params.tt_testRateLate_rp_2to4(Params.numStrats)=1;
    Params.tt_testRateAIDS_rp_2to4(Params.numStrats)=1;
    Params.tt_testRateUninfected_rp_2to4(Params.numStrats)=1;
    
    Params.tt_testRateAcute_rp_5(Params.numStrats)=1;
    Params.tt_testRateLatentA_rp_5(Params.numStrats)=1;
    Params.tt_testRateLatentB_rp_5(Params.numStrats)=1;
    Params.tt_testRateLate_rp_5(Params.numStrats)=1;
    Params.tt_testRateAIDS_rp_5(Params.numStrats)=1;
    Params.tt_testRateUninfected_rp_5(Params.numStrats)=1;
    
    
    %Dummy variable
    testingStore=zeros(Params.numStrats,1);
    
    %Annual rate of being tested for period 1
    for s = 1:(numHIVstages+1)
        for r = 1:numRace
          for p = 1:numPop
            for a = 1:Params.numAge
      
                pTest = Params.raceIndicator(:,r).* ...
                    Params.popIndicator(:,p).* Params.ageIndicator(:,a)* (Params.tt_relRiskRace_1(r)...
                    * Params.tt_relRiskPop_1(p) * Params.tt_relRiskAge_1(a) * Params.tt_relRiskHIVstage_1(s) ...
                    * Params.tt_testRefCase_1);
                
                testingStore = pTest+testingStore;
                
            end
      end
    end
       switch s
           case 1
                Params.tt_testRateAcute_rp_1 = testingStore;
                clear testingStore
                testingStore = zeros(Params.numStrats,1);
           case 2
                Params.tt_testRateLatentA_rp_1 = testingStore;
                clear testingStore 
                testingStore = zeros(Params.numStrats,1);
           case 3
                Params.tt_testRateLatentB_rp_1 = testingStore;
                clear testingStore
                testingStore = zeros(Params.numStrats,1);
           case 4
                Params.tt_testRateLate_rp_1 = testingStore;
                clear testingStore 
                testingStore = zeros(Params.numStrats,1);
           case 5
                Params.tt_testRateAIDS_rp_1 = testingStore;
                clear testingStore 
                testingStore = zeros(Params.numStrats,1);
           case 6
                Params.tt_testRateUninfected_rp_1 = testingStore;
                clear testingStore 
                testingStore = zeros(Params.numStrats,1);       
       end
 
end

    %Annual rate of being tested for periods 2 to 4
    for s = 1:(numHIVstages+1)
        for r = 1:numRace
          for p = 1:numPop
            for a = 1:Params.numAge
      
                pTest = Params.raceIndicator(:,r).* ...
                    Params.popIndicator(:,p) .* Params.ageIndicator(:,a)* ...
                    (Params.tt_relRiskRace_2to4(r)* ...
                    Params.tt_relRiskPop_2to4(p) * Params.tt_relRiskAge_2to4(a) * ...
                    Params.tt_relRiskHIVstage_2to4(s) ...
                    * Params.tt_testRefCase_2to4);
                
                testingStore = pTest+testingStore;
                
            end
      end 
    end
      
       switch s
           case 1
                Params.tt_testRateAcute_rp_2to4 = testingStore;
                clear testingStore
                testingStore = zeros(Params.numStrats,1);
           case 2
                Params.tt_testRateLatentA_rp_2to4 = testingStore;
                clear testingStore 
                testingStore = zeros(Params.numStrats,1);
           case 3
                Params.tt_testRateLatentB_rp_2to4 = testingStore;
                clear testingStore 
                testingStore = zeros(Params.numStrats,1);
           case 4
                Params.tt_testRateLate_rp_2to4 = testingStore;
                clear testingStore
                testingStore = zeros(Params.numStrats,1);
           case 5
                Params.tt_testRateAIDS_rp_2to4 = testingStore;
                clear testingStore
                testingStore = zeros(Params.numStrats,1);
           case 6
                Params.tt_testRateUninfected_rp_2to4 = testingStore;
                clear testingStore
                testingStore = zeros(Params.numStrats,1);
       end
end
    

    %Annual rate of being tested for period 5
    
    % if running allocation-based progression and targeting to YMSM only,
    % loop through testing rate collection twice.
    if and(Params.TargetIntnsToYMSM == 1, Params.tt_progressionSource == 3)
        
        testloops = 2;
        
    else
        
        testloops = 1;
        
    end    
    
    for loop = 1:testloops
        
        switch loop
            case 1
                refcasetesting = Params.tt_testRefCase_5;                
            case 2
                refcasetesting = Params.tt_testRefCase_5_nonYMSM;           
        end
    
        for s = 1:(numHIVstages+1)
            for r = 1:numRace
              for p = 1:numPop
                for a = 1:Params.numAge

                    pTest = Params.raceIndicator(:,r).* ...
                        Params.popIndicator(:,p) .* Params.ageIndicator(:,a) * ...
                        (Params.tt_relRiskRace_5(r)* ...
                        Params.tt_relRiskPop_5(p) * Params.tt_relRiskAge_5(a) *...
                        Params.tt_relRiskHIVstage_5(s) ...
                        * refcasetesting);

                    testingStore = pTest+testingStore;

                end
              end 
            end

           switch s
               case 1
                    switch loop
                        case 1
                            Params.tt_testRateAcute_rp_5 = testingStore;
                        case 2
                            Params.tt_testRateAcute_rp_5_nonYMSM = testingStore .* Params.nonYMSMIndicator;
                    end
                    clear testingStore
                    testingStore = zeros(Params.numStrats,1);
               case 2
                    switch loop
                        case 1
                            Params.tt_testRateLatentA_rp_5 = testingStore;
                        case 2
                            Params.tt_testRateLatentA_rp_5_nonYMSM = testingStore .* Params.nonYMSMIndicator;
                    end
                    clear testingStore 
                    testingStore = zeros(Params.numStrats,1);
               case 3
                    switch loop
                        case 1
                            Params.tt_testRateLatentB_rp_5 = testingStore;
                        case 2
                            Params.tt_testRateLatentB_rp_5_nonYMSM = testingStore .* Params.nonYMSMIndicator;
                    end
                    clear testingStore 
                    testingStore = zeros(Params.numStrats,1);
               case 4
                    switch loop
                        case 1
                            Params.tt_testRateLate_rp_5 = testingStore;
                        case 2
                            Params.tt_testRateLate_rp_5_nonYMSM = testingStore .* Params.nonYMSMIndicator;
                    end
                    clear testingStore
                    testingStore = zeros(Params.numStrats,1);
               case 5
                    switch loop
                        case 1
                            Params.tt_testRateAIDS_rp_5 = testingStore;
                        case 2
                            Params.tt_testRateAIDS_rp_5_nonYMSM = testingStore .* Params.nonYMSMIndicator;
                    end
                    clear testingStore
                    testingStore = zeros(Params.numStrats,1);
               case 6
                    switch loop
                        case 1
                            Params.tt_testRateUninfected_rp_5 = testingStore;
                        case 2
                            Params.tt_testRateUninfected_rp_5_nonYMSM = testingStore .* Params.nonYMSMIndicator;
                    end
                    clear testingStore
                    testingStore = zeros(Params.numStrats,1);
           end
        end
    end

% PROBABILITY OF BEING NOTIFIED, GIVEN TESTED
    % By time period
    % By type of test (rapid or conventional)
    % By result of test (negative or positive)
    
        %Rapid, negative, first testing period
    Params.tt_probNotify_rapidNeg_TestPeriod1 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_probNotify_rapidNeg_TestPeriod1);
    
        %Rapid, positive, first testing period
    Params.tt_probNotify_rapidPos_TestPeriod1 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_probNotify_rapidPos_TestPeriod1);

        %Conventional, negative, first testing period
    Params.tt_probNotify_convNeg_TestPeriod1 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_probNotify_convNeg_TestPeriod1);
    
        %Conventional, positive, first testing period
    Params.tt_probNotify_convPos_TestPeriod1 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_probNotify_convPos_TestPeriod1);

        %Rapid, negative, second testing period
    Params.tt_probNotify_rapidNeg_TestPeriod2 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_probNotify_rapidNeg_TestPeriod2);
    
        %Rapid, positive, second testing period
    Params.tt_probNotify_rapidPos_TestPeriod2 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_probNotify_rapidPos_TestPeriod2);

        %Conventional, negative, second testing period
    Params.tt_probNotify_convNeg_TestPeriod2 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_probNotify_convNeg_TestPeriod2);
    
        %Conventional, positive, second testing period
    Params.tt_probNotify_convPos_TestPeriod2 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_probNotify_convPos_TestPeriod2);

% YEAR IN WHICH PERCENTAGE OF TESTS RAPID, TEST SENSITIVITY, AND TEST
% NOTIFICATION INPUTS CHANGE

    Params.tt_YrTestSensNotifChange = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_YrTestSensNotifChange);    
    
% PERCENT OF TESTS PERFORMED IN NON-CLINICAL (VS. CLINICAL) SETTINGS

    Params.tt_pctnonclinical = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_pctnonclinical);
    
    
% PERCENT OF TESTS THAT ARE RAPID

    %Populate pctRapid_clinical_ry(3,2)
        %Percent of tests that are rapid
        %Params.D1: race
        %Params.D2: time period
        tt_pctRapid_clinical_ry(numRace,2) = 0;
        tt_pctRapid_ry(numRace,2) = 0; 

        for k=1:2
            tt_pctRapid_clinical_ry(:, k) = ExcelValues_AllParameters((parameter_Index+numRace*(k-1)): ...
                ((parameter_Index + numRace-1)+numRace*(k-1)));
        end

        parameter_Index = parameter_Index + numel(tt_pctRapid_clinical_ry);
        
        tt_pctRapid_nonclinical_TestPeriod1 = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index+numel(tt_pctRapid_nonclinical_TestPeriod1);
        
        tt_pctRapid_nonclinical_TestPeriod2 = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index+numel(tt_pctRapid_nonclinical_TestPeriod2);
        
        
        for k=1:numRace
            tt_pctRapid_ry(k,1) = tt_pctRapid_clinical_ry(k,1)*(1-Params.tt_pctnonclinical) + ...
                (Params.tt_pctnonclinical * tt_pctRapid_nonclinical_TestPeriod1);
        end
        
        for k=1:numRace
            tt_pctRapid_ry(k,2) = tt_pctRapid_clinical_ry(k,2)*(1-Params.tt_pctnonclinical) + ...
                (Params.tt_pctnonclinical * tt_pctRapid_nonclinical_TestPeriod2);
        end
           

        % Preallocate
        Params.tt_pctRapid_r_TestPeriod1(Params.numStrats)=0;
        Params.tt_pctRapid_r_TestPeriod2(Params.numStrats)=0;

        Params.tt_pctRapid_r_TestPeriod1= Params.raceIndicator * tt_pctRapid_ry(:,1);
        %Params.tt_pctRapid_r_2 uses the same rates as 1 because the input has
        %the same rate for both the first and the second time periods
        Params.tt_pctRapid_r_TestPeriod2= Params.raceIndicator * tt_pctRapid_ry(:,2);

        
% TEST SENSITIVITIES

    % Acute - Rapid, first testing period
    Params.tt_testSensAcute_Rapid_TestPeriod1 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index+numel(Params.tt_testSensAcute_Rapid_TestPeriod1);

    % Acute - Conv, first testing period
    Params.tt_testSensAcute_Conv_TestPeriod1 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index+numel(Params.tt_testSensAcute_Conv_TestPeriod1);

    % Chronic - Rapid, first testing period
    Params.tt_testSensChronic_Rapid_TestPeriod1 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index+numel(Params.tt_testSensChronic_Rapid_TestPeriod1);

    % Chronic - Conv, first testing period
    Params.tt_testSensChronic_Conv_TestPeriod1 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index+numel(Params.tt_testSensChronic_Conv_TestPeriod1);

    % Acute - Rapid, second testing period
    Params.tt_testSensAcute_Rapid_TestPeriod2 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index+numel(Params.tt_testSensAcute_Rapid_TestPeriod2);

    % Acute - Conv, second testing period
    Params.tt_testSensAcute_Conv_TestPeriod2 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index+numel(Params.tt_testSensAcute_Conv_TestPeriod2);
    
    % Chronic - Rapid, second testing period
    Params.tt_testSensChronic_Rapid_TestPeriod2 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index+numel(Params.tt_testSensChronic_Rapid_TestPeriod2);
    
    % Chronic - Conv, second testing period
    Params.tt_testSensChronic_Conv_TestPeriod2 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index+numel(Params.tt_testSensChronic_Conv_TestPeriod2);

    % Confirmatory test sensitivities
        % first testing period
        Params.tt_testSensAcute_Confirm_TestPeriod1 = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index+numel(Params.tt_testSensAcute_Confirm_TestPeriod1);

        Params.tt_testSensChronic_Confirm_TestPeriod1 = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index+numel(Params.tt_testSensChronic_Confirm_TestPeriod1);

        % second testing period
        Params.tt_testSensAcute_Confirm_TestPeriod2 = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index+numel(Params.tt_testSensAcute_Confirm_TestPeriod2);

        Params.tt_testSensChronic_Confirm_TestPeriod2 = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index+numel(Params.tt_testSensChronic_Confirm_TestPeriod2);



    % LINKAGE FIRST
    % Probability of diagnosed individual linked to care within 1 year
    
    % Note: intentionally left as a probability (rather than a rate)
    % because it is applied as the percentage of individuals who transition
    % immediately into care vs. those who transition into into the aware
    % stages
    
    % Dim 1: race
    % Dim 2: years
    
    tt_linkageFirst_ry(numRace,Params.numRateInputPeriod) = 0; 
    tt_linkageFirst_r_nonYMSM(numRace,1) = 0;

    for k=1:Params.numRateInputPeriod
        tt_linkageFirst_ry(:, k) = ExcelValues_AllParameters((parameter_Index+numRace*(k-1)): ...
            ((parameter_Index + numRace-1)+numRace*(k-1)));
    end
    parameter_Index = parameter_Index + numel(tt_linkageFirst_ry);
    
    tt_linkageFirst_r_nonYMSM(:, 1) = ExcelValues_AllParameters(parameter_Index: ...
            ((parameter_Index + numRace-1)));
    parameter_Index = parameter_Index + numel(tt_linkageFirst_r_nonYMSM);    
    
    %Preallocate
    Params.tt_linkageFirst_r_1(Params.numStrats)= 0;
    Params.tt_linkageFirst_r_2to4(Params.numStrats)= 0;
    Params.tt_linkageFirst_r_5(Params.numStrats)= 0;
    Params.tt_linkageFirst_r_5_nonYMSM(Params.numStrats)= 0;
    
    %Restratify
    Params.tt_linkageFirst_r_1 = Params.raceIndicator*tt_linkageFirst_ry(:,1);
    Params.tt_linkageFirst_r_2to4 = Params.raceIndicator*tt_linkageFirst_ry(:,2);
    Params.tt_linkageFirst_r_5 = Params.raceIndicator*tt_linkageFirst_ry(:,3);
    Params.tt_linkageFirst_r_5_nonYMSM = (Params.raceIndicator*tt_linkageFirst_r_nonYMSM(:,1)).* Params.nonYMSMIndicator;
 
    % LINKAGE TO CARE AFTER DIAGNOSIS (VS. IMMEDIATE)
    % Probability of diagnosed individual linked each year after the first year

    % Pre-allocate
    store_baselinkageProb_ry(numRace,Params.numRateInputPeriod) = 0;
    store_baselinkageProb_r_nonYMSM(numRace,1) = 0;

    % Bring in probabilities by race from Excel
        % [3x3] 
        % Dim 1: by race
        % Dim 2: by time period
    
    for k=1:Params.numRateInputPeriod
        store_baselinkageProb_ry(:, k) = ExcelValues_AllParameters((parameter_Index+numRace*(k-1)): ...
            ((parameter_Index + numRace-1)+numRace*(k-1)));
    end
    parameter_Index = parameter_Index + numel(store_baselinkageProb_ry);
    
    store_baselinkageProb_r_nonYMSM(:, 1) = ExcelValues_AllParameters(parameter_Index: ...
            ((parameter_Index + numRace-1)));
    parameter_Index = parameter_Index + numel(store_baselinkageProb_r_nonYMSM); 
    
    % Convert to rates
        store_baseLinkageRates_ry = -log(1-store_baselinkageProb_ry);
        store_baselinkageRates_r_nonYMSM = -log(1-store_baselinkageProb_r_nonYMSM);
        
    % Apply the probabilities by race and period to full subpopulation
    % variables [273x1]    
            
        rate_baseLinkageAfter_1 = Params.raceIndicator * store_baseLinkageRates_ry(:,Params.period_1);
        rate_baseLinkageAfter_2to4 = Params.raceIndicator * store_baseLinkageRates_ry(:,Params.period_2to4);
        rate_baseLinkageAfter_5 = Params.raceIndicator * store_baseLinkageRates_ry(:,Params.period_5);
        rate_baseLinkageAfter_5_nonYMSM = Params.raceIndicator * store_baselinkageRates_r_nonYMSM(:,1);

    % RELATIVE RISK OF LINKAGE TO CARE AFTER, BY DISEASE STAGE
    % Three sets of [5x1] variables
    % Acute/LatentA/LatentB/Late/AIDS

    % Read in from Excel

        % Period 1 [5x1]
        Params.store_relRiskLinkageAfter_1 = ExcelValues_AllParameters(parameter_Index: ...
            parameter_Index + numHIVstages-1);
    
        parameter_Index = parameter_Index + numel(Params.store_relRiskLinkageAfter_1);
            
        % Period 2 to 4 [5x1]
        Params.store_relRiskLinkageAfter_2to4 = ExcelValues_AllParameters(parameter_Index: ...
            parameter_Index + numHIVstages-1);
            
        parameter_Index = parameter_Index + numel(Params.store_relRiskLinkageAfter_2to4);
        
        % Period 5 [5x1]
        Params.store_relRiskLinkageAfter_5 = ExcelValues_AllParameters(parameter_Index: ...
            parameter_Index + numHIVstages-1);
            
        parameter_Index = parameter_Index + numel(Params.store_relRiskLinkageAfter_5);
        
    
    % CALCULATE LINKAGE TO CARE AFTER, BY DISEASE STAGE
    % [273x1]
    % Apply relative risks of LTC After by disease state
               
        
    % Period 1

        % Loop through HIV stages
        for nHIVstage = Params.stage_Acute:Params.stage_AIDS

            % Calculation
            calculatedVariable = rate_baseLinkageAfter_1 * Params.store_relRiskLinkageAfter_1(nHIVstage);

            % Record calculated variable as a parameter 
            switch nHIVstage
                case Params.stage_Acute
                    Params.tt_linkageAfterRate_Acute_1 = calculatedVariable;

                case Params.stage_LatentA
                    Params.tt_linkageAfterRate_LatentA_1 = calculatedVariable;

                case Params.stage_LatentB
                    Params.tt_linkageAfterRate_LatentB_1 = calculatedVariable;

                case Params.stage_Late
                    Params.tt_linkageAfterRate_Late_1 = calculatedVariable;

                case Params.stage_AIDS
                    Params.tt_linkageAfterRate_AIDS_1 = calculatedVariable;
            end
        end
          
        
    % Period 2 to 4

         % Loop through HIV stages
        for nHIVstage = Params.stage_Acute:Params.stage_AIDS

            % Calculation
            calculatedVariable = rate_baseLinkageAfter_2to4 * Params.store_relRiskLinkageAfter_2to4(nHIVstage);

            % Record calculated variable as a parameter 
            switch nHIVstage
                case Params.stage_Acute
                    Params.tt_linkageAfterRate_Acute_2to4 = calculatedVariable; 

                case Params.stage_LatentA
                    Params.tt_linkageAfterRate_LatentA_2to4 = calculatedVariable;

                case Params.stage_LatentB
                    Params.tt_linkageAfterRate_LatentB_2to4 = calculatedVariable;

                case Params.stage_Late
                    Params.tt_linkageAfterRate_Late_2to4 = calculatedVariable;

                case Params.stage_AIDS
                    Params.tt_linkageAfterRate_AIDS_2to4 = calculatedVariable;
            end
        end

        
   % Period 5
   
   % if running allocation-based progression and targeting to YMSM only,
    % loop through testing rate collection twice.
    if and(Params.TargetIntnsToYMSM == 1, Params.tt_progressionSource == 3)
        
        LTCafterloops = 2;
        
    else
        
        LTCafterloops = 1;
        
    end    
    
    for loop = 1:LTCafterloops
        
        switch loop
            case 1
                refcaseLTCafter = rate_baseLinkageAfter_5;                
            case 2
                refcaseLTCafter = rate_baseLinkageAfter_5_nonYMSM;           
        end
        
        % Loop through HIV stages
        for nHIVstage = Params.stage_Acute:Params.stage_AIDS

            % Calculation
            calculatedVariable = refcaseLTCafter * Params.store_relRiskLinkageAfter_5(nHIVstage);

            % Record calculated variable as a parameter
            switch nHIVstage
                case Params.stage_Acute
                    switch loop
                        case 1
                            Params.tt_linkageAfterRate_Acute_5 = calculatedVariable; 
                        case 2
                            Params.tt_linkageAfterRate_Acute_5_nonYMSM = calculatedVariable.* Params.nonYMSMIndicator;
                    end
                    
                case Params.stage_LatentA
                    switch loop
                        case 1
                            Params.tt_linkageAfterRate_LatentA_5 = calculatedVariable;
                        case 2
                            Params.tt_linkageAfterRate_LatentA_5_nonYMSM = calculatedVariable.* Params.nonYMSMIndicator;
                    end

                case Params.stage_LatentB
                    switch loop
                        case 1
                            Params.tt_linkageAfterRate_LatentB_5 = calculatedVariable;
                        case 2
                            Params.tt_linkageAfterRate_LatentB_5_nonYMSM = calculatedVariable.* Params.nonYMSMIndicator;
                    end

                case Params.stage_Late
                    switch loop
                        case 1
                            Params.tt_linkageAfterRate_Late_5 = calculatedVariable;
                        case 2
                            Params.tt_linkageAfterRate_Late_5_nonYMSM = calculatedVariable.* Params.nonYMSMIndicator;
                    end
                    
                case Params.stage_AIDS
                    switch loop
                        case 1
                            Params.tt_linkageAfterRate_AIDS_5 = calculatedVariable;
                        case 2
                            Params.tt_linkageAfterRate_AIDS_5_nonYMSM = calculatedVariable.* Params.nonYMSMIndicator;
                    end
            end

        end

    end
    
    % DROP OUT

    % Dropping from [In Care] to [Aware]
        % Annual rate of dropping out of in care
        
        % Preallocate
        store_probDropOutCareToAware(numRace,Params.numRateInputPeriod) = 0; 

        % Read probabilities from Excel
            %Dim 1: race
            %Dim 2: period
    
        for k = 1:Params.numRateInputPeriod
            store_probDropOutCareToAware(:, k) = ExcelValues_AllParameters((parameter_Index+numRace*(k-1)): ...
                ((parameter_Index + numRace-1)+numRace*(k-1)));
        end
        parameter_Index = parameter_Index + numel(store_probDropOutCareToAware);

        % Convert to rates
        store_rateDropOutCareToAware = -log(1-store_probDropOutCareToAware);
        
        % Apply race-specific values to full 273x1 subpopulation variables            
        Params.tt_dropOutRate_CareToAware_1 = Params.raceIndicator*store_rateDropOutCareToAware(:,Params.period_1);
        Params.tt_dropOutRate_CareToAware_2to4 = Params.raceIndicator*store_rateDropOutCareToAware(:,Params.period_2to4);
        Params.tt_dropOutRate_CareToAware_5 = Params.raceIndicator*store_rateDropOutCareToAware(:,Params.period_5);

        
   % Dropping from [ART-not-VLS] to [Aware]
    
        store_probDropOutANVToAware(numRace,Params.numRateInputPeriod) = 0; 

        % Read probabilities from Excel
            %Dim 1: race
            %Dim 2: period
    
        for k = 1:Params.numRateInputPeriod
            store_probDropOutANVToAware(:, k) = ExcelValues_AllParameters((parameter_Index+numRace*(k-1)): ...
                ((parameter_Index + numRace-1)+numRace*(k-1)));
        end
        parameter_Index = parameter_Index + numel(store_probDropOutANVToAware);

        % Convert to rates
        store_rateDropOutANVToAware = -log(1-store_probDropOutANVToAware);
        
        % Apply race-specific values to full 273x1 subpopulation variables            
        Params.tt_dropOutRate_ANVToAware_1 = Params.raceIndicator*store_rateDropOutANVToAware(:,Params.period_1);
        Params.tt_dropOutRate_ANVToAware_2to4 = Params.raceIndicator*store_rateDropOutANVToAware(:,Params.period_2to4);
        Params.tt_dropOutRate_ANVToAware_5 = Params.raceIndicator*store_rateDropOutANVToAware(:,Params.period_5);
        
        
    % Dropping from [ART-not-VLS] to [In Care]
        
        store_probDropOutANVToCare(numRace,Params.numRateInputPeriod) = 0; 

        % Read probabilities from Excel
            %Dim 1: race
            %Dim 2: period
    
        for k = 1:Params.numRateInputPeriod
            store_probDropOutANVToCare(:, k) = ExcelValues_AllParameters((parameter_Index+numRace*(k-1)): ...
                ((parameter_Index + numRace-1)+numRace*(k-1)));
        end
        parameter_Index = parameter_Index + numel(store_probDropOutANVToCare);

        % Convert to rates
        store_rateDropOutANVToCare = -log(1-store_probDropOutANVToCare);
        
        % Apply race-specific values to full 273x1 subpopulation variables            
        tt_dropOutRate_ANVToCare_1 = Params.raceIndicator*store_rateDropOutANVToCare(:,Params.period_1);
        tt_dropOutRate_ANVToCare_2to4 = Params.raceIndicator*store_rateDropOutANVToCare(:,Params.period_2to4);
        tt_dropOutRate_ANVToCare_5 = Params.raceIndicator*store_rateDropOutANVToCare(:,Params.period_5);
        
        % Read in relative risks of dropping off ART by age group (added
        % 07/31/2018)
        % Read values in from Excel
            % Period 1 [7x1]    
            store_relRiskDropOutANVAge_1(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + Params.numAge - 1));
                parameter_Index = parameter_Index + Params.numAge;
            
            % Period 2 to 4 [7x1]    
            store_relRiskDropOutANVAge_2to4(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + Params.numAge - 1));
                parameter_Index = parameter_Index + Params.numAge;
            
            % Period 5 [7x1]    
            store_relRiskDropOutANVAge_5(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + Params.numAge - 1));
                parameter_Index = parameter_Index + Params.numAge;
    
         % Apply relative risks to [273x1] matrix
            tt_relRiskDropOutANVAge_1 = Params.ageIndicator * store_relRiskDropOutANVAge_1;
            tt_relRiskDropOutANVAge_2to4 = Params.ageIndicator * store_relRiskDropOutANVAge_2to4;
            tt_relRiskDropOutANVAge_5 = Params.ageIndicator * store_relRiskDropOutANVAge_5;
        
        % Apply relative risks to VLS drop-out rates    
        Params.tt_dropOutRate_ANVToCare_1 = tt_dropOutRate_ANVToCare_1 .* tt_relRiskDropOutANVAge_1;
        Params.tt_dropOutRate_ANVToCare_2to4 = tt_dropOutRate_ANVToCare_2to4 .* tt_relRiskDropOutANVAge_2to4;
        Params.tt_dropOutRate_ANVToCare_5 = tt_dropOutRate_ANVToCare_5 .* tt_relRiskDropOutANVAge_5;
        
        % Read in odds ratios for dropping off ART by age group (added
        % 07/31/2018) - for calibration calculations
        % Read values in from Excel
        Params.oddsRatio_DropOutANVAge(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + Params.numAge - 1));
                parameter_Index = parameter_Index + Params.numAge;
    
    % Dropping from [VLS] to [ART-not-VLS]
    
     % Read in probabilities from Excel
        store_probDropOutVLSToANV(numRace,Params.numRateInputPeriod) = 0; 
        store_probDropOutVLSToANV_nonYMSM(numRace,1) = 0; 

        % Read probabilities from Excel
            %Dim 1: race
            %Dim 2: period
    
        for k = 1:Params.numRateInputPeriod
            store_probDropOutVLSToANV(:, k) = ExcelValues_AllParameters((parameter_Index+numRace*(k-1)): ...
                ((parameter_Index + numRace-1)+numRace*(k-1)));
        end
        parameter_Index = parameter_Index + numel(store_probDropOutVLSToANV);
        
        store_probDropOutVLSToANV_nonYMSM(:, 1) = ExcelValues_AllParameters(parameter_Index: ...
            ((parameter_Index + numRace-1)));
        parameter_Index = parameter_Index + numel(store_probDropOutVLSToANV_nonYMSM);

        % Convert to rates
        store_rateDropOutVLSToANV = -log(1-store_probDropOutVLSToANV);
        store_rateDropOutVLSToANV_nonYMSM = -log(1-store_probDropOutVLSToANV_nonYMSM);
        
        % Apply race-specific values to full 273x1 subpopulation variables            
        tt_dropOutRate_VLSToANV_1 = Params.raceIndicator*store_rateDropOutVLSToANV(:,Params.period_1);
        tt_dropOutRate_VLSToANV_2to4 = Params.raceIndicator*store_rateDropOutVLSToANV(:,Params.period_2to4);
        tt_dropOutRate_VLSToANV_5 = Params.raceIndicator*store_rateDropOutVLSToANV(:,Params.period_5);
        tt_dropOutRate_VLSToANV_5_nonYMSM = Params.raceIndicator*store_rateDropOutVLSToANV_nonYMSM(:,1);
        
     % Read in relative risks of losing VLS by transmission group (added by
     % JC on 11/08/2017)
        % Read values in from Excel
            % Period 1 [3x1]    
            Params.store_relRiskLoseVLSPop_1(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + numPop - 1));
                parameter_Index = parameter_Index + numPop;
            
            % Period 2 to 4 [3x1]    
            Params.store_relRiskLoseVLSPop_2to4(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + numPop - 1));
                parameter_Index = parameter_Index + numPop;
            
            % Period 5 [3x1]    
            Params.store_relRiskLoseVLSPop_5(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + numPop - 1));
                parameter_Index = parameter_Index + numPop;
    
         % Apply relative risks to [273x1] matrix
            tt_relRiskLoseVLSPop_1 = Params.popIndicator * Params.store_relRiskLoseVLSPop_1;
            tt_relRiskLoseVLSPop_2to4 = Params.popIndicator * Params.store_relRiskLoseVLSPop_2to4;
            tt_relRiskLoseVLSPop_5 = Params.popIndicator * Params.store_relRiskLoseVLSPop_5;
            
     % Read in relative risks of losing VLS by age group
        % Read values in from Excel
            % Period 1 [7x1]    
            Params.store_relRiskLoseVLSAge_1(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + Params.numAge - 1));
                parameter_Index = parameter_Index + Params.numAge;
            
            % Period 2 to 4 [7x1]    
            Params.store_relRiskLoseVLSAge_2to4(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + Params.numAge - 1));
                parameter_Index = parameter_Index + Params.numAge;
            
            % Period 5 [7x1]    
            Params.store_relRiskLoseVLSAge_5(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + Params.numAge - 1));
                parameter_Index = parameter_Index + Params.numAge;
    
         % Apply relative risks to [273x1] matrix
            tt_relRiskLoseVLSAge_1 = Params.ageIndicator * Params.store_relRiskLoseVLSAge_1;
            tt_relRiskLoseVLSAge_2to4 = Params.ageIndicator * Params.store_relRiskLoseVLSAge_2to4;
            tt_relRiskLoseVLSSge_5 = Params.ageIndicator * Params.store_relRiskLoseVLSAge_5;
        
        % Apply relative risks to VLS drop-out rates    
        Params.tt_dropOutRate_VLSToANV_1 = tt_dropOutRate_VLSToANV_1 .* tt_relRiskLoseVLSPop_1 .* tt_relRiskLoseVLSAge_1;
        Params.tt_dropOutRate_VLSToANV_2to4 = tt_dropOutRate_VLSToANV_2to4 .* tt_relRiskLoseVLSPop_2to4 .* tt_relRiskLoseVLSAge_2to4;
        Params.tt_dropOutRate_VLSToANV_5 = tt_dropOutRate_VLSToANV_5 .* tt_relRiskLoseVLSPop_5 .* tt_relRiskLoseVLSSge_5;
        Params.tt_dropOutRate_VLSToANV_5_nonYMSM = (tt_dropOutRate_VLSToANV_5_nonYMSM .* tt_relRiskLoseVLSPop_5 .* tt_relRiskLoseVLSSge_5).* Params.nonYMSMIndicator;
        
    % Dropping from [VLS] to [LTC]
    
     % Read in probabilities from Excel
        store_probDropOutVLSToLTC(numRace,Params.numRateInputPeriod) = 0; 

        % Read probabilities from Excel
            %Dim 1: race
            %Dim 2: period
    
        for k = 1:Params.numRateInputPeriod
            store_probDropOutVLSToLTC(:, k) = ExcelValues_AllParameters((parameter_Index+numRace*(k-1)): ...
                ((parameter_Index + numRace-1)+numRace*(k-1)));
        end
        parameter_Index = parameter_Index + numel(store_probDropOutVLSToLTC);

        % Convert to rates
        store_rateDropOutVLSToLTC = -log(1-store_probDropOutVLSToLTC);
        
        % Apply race-specific values to full 273x1 subpopulation variables            
        Params.tt_dropOutRate_VLSToLTC_1 = Params.raceIndicator*store_rateDropOutVLSToLTC(:,Params.period_1);
        Params.tt_dropOutRate_VLSToLTC_2to4 = Params.raceIndicator*store_rateDropOutVLSToLTC(:,Params.period_2to4);
        Params.tt_dropOutRate_VLSToLTC_5 = Params.raceIndicator*store_rateDropOutVLSToLTC(:,Params.period_5);
      
        
   % Dropping from [VLS] to [Aware]
    
     % Read in probabilities from Excel
        store_probDropOutVLSToAware(numRace,Params.numRateInputPeriod) = 0; 

        % Read probabilities from Excel
            %Dim 1: race
            %Dim 2: period
    
        for k = 1:Params.numRateInputPeriod
            store_probDropOutVLSToAware(:, k) = ExcelValues_AllParameters((parameter_Index+numRace*(k-1)): ...
                ((parameter_Index + numRace-1)+numRace*(k-1)));
        end
        parameter_Index = parameter_Index + numel(store_probDropOutVLSToAware);

        % Convert to rates
        store_rateDropOutVLSToAware = -log(1-store_probDropOutVLSToAware);
        
        % Apply race-specific values to full 273x1 subpopulation variables            
        Params.tt_dropOutRate_VLSToAware_1 = Params.raceIndicator*store_rateDropOutVLSToAware(:,Params.period_1);
        Params.tt_dropOutRate_VLSToAware_2to4 = Params.raceIndicator*store_rateDropOutVLSToAware(:,Params.period_2to4);
        Params.tt_dropOutRate_VLSToAware_5 = Params.raceIndicator*store_rateDropOutVLSToAware(:,Params.period_5);         
          
          
% ART ELIGIBILITY
    % Applies a 0 or 1 to each disease stage showing whether an individual
    % is eligible
    
    % Value from Excel is the number of the disease stage that one is
    % eligible for ART
        % E.g., a value of 4 means PLWH are eligible for ART starting at the
        % 4th disease stage = CD4 200-350
        
        stageEligibleByPeriod = ExcelValues_AllParameters(parameter_Index:parameter_Index+Params.numRateInputPeriod - 1);
        parameter_Index = parameter_Index + numel(stageEligibleByPeriod);

    % Preallocate
        Params.tt_ARTElig_1=zeros(Params.numHIVstages,1);
        Params.tt_ARTElig_2to4=zeros(Params.numHIVstages,1);
        Params.tt_ARTElig_5=zeros(Params.numHIVstages,1);

    % Populate MATLAB parameter
        % Matlab variables are binary by disease stage (either eligible = 1 or
        % ineligible = 0)
        
        % This code puts a 1 at every disease stage an individual is
        % eligible
        Params.tt_ARTElig_1(stageEligibleByPeriod(Params.period_1):end) = 1;
        Params.tt_ARTElig_2to4(stageEligibleByPeriod(Params.period_2to4):end) = 1;
        Params.tt_ARTElig_5(stageEligibleByPeriod(Params.period_5):end) = 1;


% ART PRESCRIPTION, NOT ART-EFFECTED
    
    % Annual rate of initiating ART, given eligible for ART and not
    % affected by ART
    
        % Read probabilities from Excel
            % Period 1 [5x1]
                store_initARTProbGivenElig_1(:,1) =  ...
                    ExcelValues_AllParameters(parameter_Index :(parameter_Index + Params.numHIVstages - 1));
                parameter_Index = parameter_Index + Params.numHIVstages;

            % Period 2 to 4 [5x1]
                store_initARTProbGivenElig_2to4(:,1) =  ...
                    ExcelValues_AllParameters(parameter_Index :(parameter_Index + Params.numHIVstages - 1));
                parameter_Index = parameter_Index + Params.numHIVstages;

            % Period 5 [5x1]
                store_initARTProbGivenElig_5(:,1) =  ...
                    ExcelValues_AllParameters(parameter_Index :(parameter_Index + Params.numHIVstages - 1));
                parameter_Index = parameter_Index + Params.numHIVstages;
                
            % Period 5 - non-YMSM [5x1] (only relevant when running
            % allocation-based progression targeted to YMSM only)
                store_initARTProbGivenElig_5_nonYMSM(:,1) =  ...
                    ExcelValues_AllParameters(parameter_Index :(parameter_Index + Params.numHIVstages - 1));
                parameter_Index = parameter_Index + Params.numHIVstages;

        % Adjust for eligibility
            % Only individuals in eligible stages can initiate ART
            base_initARTProb_1 = Params.tt_ARTElig_1 .* store_initARTProbGivenElig_1;
            base_initARTProb_2to4 = Params.tt_ARTElig_2to4 .* store_initARTProbGivenElig_2to4;
            base_initARTProb_5 = Params.tt_ARTElig_5 .* store_initARTProbGivenElig_5;
            base_initARTProb_5_nonYMSM = Params.tt_ARTElig_5 .* store_initARTProbGivenElig_5_nonYMSM;

        % Flip to rates
            base_initARTRate_1 = -log(1-base_initARTProb_1);
            base_initARTRate_2to4 = -log(1-base_initARTProb_2to4);
            base_initARTRate_5 = -log(1-base_initARTProb_5);
            base_initARTRate_5_nonYMSM = -log(1-base_initARTProb_5_nonYMSM);
    
        % Save base rates for Calib_updateParams (added by JC on
        % 2018/01/25)
            Params.rate_ARTInitRate_B_1 = base_initARTRate_1(1,1);
            Params.rate_ARTInitRate_C_1 = base_initARTRate_1(2,1);
            Params.rate_ARTInitRate_D_1 = base_initARTRate_1(3,1);
            Params.rate_ARTInitRate_E_1 = base_initARTRate_1(4,1);
            Params.rate_ARTInitRate_F_1 = base_initARTRate_1(5,1);
            
            Params.rate_ARTInitRate_B_2to4 = base_initARTRate_2to4(1,1);
            Params.rate_ARTInitRate_C_2to4 = base_initARTRate_2to4(2,1);
            Params.rate_ARTInitRate_D_2to4 = base_initARTRate_2to4(3,1);
            Params.rate_ARTInitRate_E_2to4 = base_initARTRate_2to4(4,1);
            Params.rate_ARTInitRate_F_2to4 = base_initARTRate_2to4(5,1);
            
    % Relative risk of ART Initiation by race 
            
        % Read values in from Excel
            % Period 1 [3x1]    
            store_relRiskARTInitRace_1(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + numRace - 1));
                parameter_Index = parameter_Index + numRace;
            
            % Period 2 to 4 [3x1]    
            store_relRiskARTInitRace_2to4(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + numRace - 1));
                parameter_Index = parameter_Index + numRace;
            
            % Period 5 [3x1]    
            store_relRiskARTInitRace_5(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + numRace - 1));
                parameter_Index = parameter_Index + numRace;
    
         % Apply relative risks to [273x1] matrix
            Params.tt_relRiskARTInitRate_1 = Params.raceIndicator * store_relRiskARTInitRace_1;
            Params.tt_relRiskARTInitRate_2to4 = Params.raceIndicator * store_relRiskARTInitRace_2to4;
            Params.tt_relRiskARTInitRate_5 = Params.raceIndicator * store_relRiskARTInitRace_5;
                
     % ART initiation rates by disease stage and race [273x1]
               
        % Period 1
            % Loop through HIV stages
            for nHIVstage = Params.stage_Acute:Params.stage_AIDS
                
                % Calculation: [1 base rate] * [273x1 of relative risks]
                calculatedVariable = base_initARTRate_1(nHIVstage) * Params.tt_relRiskARTInitRate_1;

                switch nHIVstage
                    case Params.stage_Acute
                        Params.tt_ARTInitRateAcute_1 = calculatedVariable;
                        
                    case Params.stage_LatentA
                        Params.tt_ARTInitRateLatentA_1 = calculatedVariable;

                    case Params.stage_LatentB
                        Params.tt_ARTInitRateLatentB_1 = calculatedVariable;

                    case Params.stage_Late
                        Params.tt_ARTInitRateLate_1 = calculatedVariable;

                    case Params.stage_AIDS
                        Params.tt_ARTInitRateAIDS_1 = calculatedVariable;
                end
            end
        
        % Period 2 to 4
             % Loop through HIV stages
            for nHIVstage = Params.stage_Acute:Params.stage_AIDS
                
                % Calculation
                calculatedVariable = base_initARTRate_2to4(nHIVstage) * Params.tt_relRiskARTInitRate_2to4;

                switch nHIVstage
                    case Params.stage_Acute
                        Params.tt_ARTInitRateAcute_2to4 = calculatedVariable;

                    case Params.stage_LatentA
                        Params.tt_ARTInitRateLatentA_2to4 = calculatedVariable;

                    case Params.stage_LatentB
                        Params.tt_ARTInitRateLatentB_2to4 = calculatedVariable;

                    case Params.stage_Late
                        Params.tt_ARTInitRateLate_2to4 = calculatedVariable;

                    case Params.stage_AIDS
                        Params.tt_ARTInitRateAIDS_2to4 = calculatedVariable;
                end
            end
            
         % Period 5
         if and(Params.TargetIntnsToYMSM == 1, Params.tt_progressionSource == 3)
        
            ARTinitloops = 2;
        
         else
        
            ARTinitloops = 1;
        
         end
         
         for loop = 1:ARTinitloops
        
            switch loop
                case 1
                    refcaseARTinit = base_initARTRate_5;                
                case 2
                    refcaseARTinit = base_initARTRate_5_nonYMSM;           
            end
            
             % Loop through HIV stages
            for nHIVstage = Params.stage_Acute:Params.stage_AIDS
                
                % Calculation
                calculatedVariable = refcaseARTinit(nHIVstage) * Params.tt_relRiskARTInitRate_5;

                switch nHIVstage
                    case Params.stage_Acute
                        switch loop
                            case 1                             
                                Params.tt_ARTInitRateAcute_5 = calculatedVariable;
                            case 2
                                Params.tt_ARTInitRateAcute_5_nonYMSM = calculatedVariable .* Params.nonYMSMIndicator;
                        end
                    
                    case Params.stage_LatentA
                        switch loop
                            case 1
                                Params.tt_ARTInitRateLatentA_5 = calculatedVariable;
                            case 2
                                Params.tt_ARTInitRateLatentA_5_nonYMSM = calculatedVariable .* Params.nonYMSMIndicator;
                        end    

                    case Params.stage_LatentB
                        switch loop
                            case 1
                                Params.tt_ARTInitRateLatentB_5 = calculatedVariable;
                            case 2
                                Params.tt_ARTInitRateLatentB_5_nonYMSM = calculatedVariable .* Params.nonYMSMIndicator;
                        end

                    case Params.stage_Late
                        switch loop
                            case 1
                                Params.tt_ARTInitRateLate_5 = calculatedVariable;
                            case 2
                                Params.tt_ARTInitRateLate_5_nonYMSM = calculatedVariable .* Params.nonYMSMIndicator;
                        end

                    case Params.stage_AIDS
                        switch loop
                            case 1
                                Params.tt_ARTInitRateAIDS_5 = calculatedVariable;
                            case 2
                                Params.tt_ARTInitRateAIDS_5_nonYMSM = calculatedVariable .* Params.nonYMSMIndicator;
                        end
                end
            end
         end
            
% ART PRESCRIPTION, AFFECTED BY ART, NOT VLS  

    % Informs transition [ANV to VLS]
        % Note: not all prescribed from ANV transition to VLS. These rates
        % are multiplied by [% prescribed ART who become VLS]
        
        store_probPrescribedARTFromANV(numRace,Params.numRateInputPeriod) = 0; 
        store_probPrescribedARTFromANV_nonYMSM(numRace,1) = 0; 

        % Read probabilities from Excel
            %Dim 1: race
            %Dim 2: period
    
        for k = 1:Params.numRateInputPeriod
            store_probPrescribedARTFromANV(:, k) = ExcelValues_AllParameters((parameter_Index+numRace*(k-1)): ...
                ((parameter_Index + numRace-1)+numRace*(k-1)));
        end
        parameter_Index = parameter_Index + numel(store_probPrescribedARTFromANV);
        
        store_probPrescribedARTFromANV_nonYMSM(:,1) = ...
            ExcelValues_AllParameters(parameter_Index :(parameter_Index + numRace - 1));
        parameter_Index = parameter_Index + numRace;

        % Convert to rates
        store_ratePrescribedARTFromANV = -log(1-store_probPrescribedARTFromANV);
        store_ratePrescribedARTFromANV_nonYMSM = -log(1-store_probPrescribedARTFromANV_nonYMSM);
        
        % Apply race-specific values to full 273x1 subpopulation variables            
        tt_ARTInitRateFromANV_1 = Params.raceIndicator*store_ratePrescribedARTFromANV(:,Params.period_1);
        tt_ARTInitRateFromANV_2to4 = Params.raceIndicator*store_ratePrescribedARTFromANV(:,Params.period_2to4);
        tt_ARTInitRateFromANV_5 = Params.raceIndicator*store_ratePrescribedARTFromANV(:,Params.period_5);
        tt_ARTInitRateFromANV_5_nonYMSM = Params.raceIndicator*store_ratePrescribedARTFromANV_nonYMSM(:,1);
        
    % Read in relative risks of becoming VLS by transmission group (added by
     % JC on 11/08/2017)
        % Read values in from Excel
            % Period 1 [3x1]    
            Params.store_relRiskBecomeVLSPop_1(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + numPop - 1));
                parameter_Index = parameter_Index + numPop;
            
            % Period 2 to 4 [3x1]    
            Params.store_relRiskBecomeVLSPop_2to4(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + numPop - 1));
                parameter_Index = parameter_Index + numPop;
            
            % Period 5 [3x1]    
            Params.store_relRiskBecomeVLSPop_5(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + numPop - 1));
                parameter_Index = parameter_Index + numPop;
    
         % Apply relative risks to [273x1] matrix
            tt_relRiskBecomeVLSPop_1 = Params.popIndicator * Params.store_relRiskBecomeVLSPop_1;
            tt_relRiskBecomeVLSPop_2to4 = Params.popIndicator * Params.store_relRiskBecomeVLSPop_2to4;
            tt_relRiskBecomeVLSPop_5 = Params.popIndicator * Params.store_relRiskBecomeVLSPop_5;
            
     % Read in relative risks of becoming VLS by age group
        % Read values in from Excel
            % Period 1 [7x1]    
            Params.store_relRiskBecomeVLSAge_1(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + Params.numAge - 1));
                parameter_Index = parameter_Index + Params.numAge;
            
            % Period 2 to 4 [7x1]    
            Params.store_relRiskBecomeVLSAge_2to4(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + Params.numAge - 1));
                parameter_Index = parameter_Index + Params.numAge;
            
            % Period 5 [7x1]    
            Params.store_relRiskBecomeVLSAge_5(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + Params.numAge - 1));
                parameter_Index = parameter_Index + Params.numAge;
    
         % Apply relative risks to [273x1] matrix
            tt_relRiskBecomeVLSAge_1 = Params.ageIndicator * Params.store_relRiskBecomeVLSAge_1;
            tt_relRiskBecomeVLSAge_2to4 = Params.ageIndicator * Params.store_relRiskBecomeVLSAge_2to4;
            tt_relRiskBecomeVLSSge_5 = Params.ageIndicator * Params.store_relRiskBecomeVLSAge_5;
        
        % Apply relative risks to ANV-to-VLS rates    
        Params.tt_ARTInitRateFromANV_1 = tt_ARTInitRateFromANV_1 .* tt_relRiskBecomeVLSPop_1 .* tt_relRiskBecomeVLSAge_1;
        Params.tt_ARTInitRateFromANV_2to4 = tt_ARTInitRateFromANV_2to4 .* tt_relRiskBecomeVLSPop_2to4 .* tt_relRiskBecomeVLSAge_2to4;
        Params.tt_ARTInitRateFromANV_5 = tt_ARTInitRateFromANV_5 .* tt_relRiskBecomeVLSPop_5 .* tt_relRiskBecomeVLSSge_5;
        Params.tt_ARTInitRateFromANV_5_nonYMSM = (tt_ARTInitRateFromANV_5_nonYMSM .* tt_relRiskBecomeVLSPop_5 .* tt_relRiskBecomeVLSSge_5) .* Params.nonYMSMIndicator;
        
        
% PERCENT WHO BECOME VLS WHEN PRESCRIBED ART
    %[1x1]
   
    Params.tt_PctInitiateARTWhoBecomeVLS = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_PctInitiateARTWhoBecomeVLS);    
    
    
% DISTRIBUTION OF PLWH IN ANV    
   % Percent of ANV in care (ART + not on ART)
   Params.tt_PctANVWhoAreInCareOrART_1 = ExcelValues_AllParameters(parameter_Index);
   parameter_Index = parameter_Index + numel(Params.tt_PctANVWhoAreInCareOrART_1);
   
   Params.tt_PctANVWhoAreInCareOrART_2to4 = ExcelValues_AllParameters(parameter_Index);
   parameter_Index = parameter_Index + numel(Params.tt_PctANVWhoAreInCareOrART_2to4);
   
   Params.tt_PctANVWhoAreInCareOrART_5 = ExcelValues_AllParameters(parameter_Index);
   parameter_Index = parameter_Index + numel(Params.tt_PctANVWhoAreInCareOrART_5);
   
            
   % Percent of in-care ANV who are on ART
   Params.tt_PctInCareANVWhoAreOnART_1 = ExcelValues_AllParameters(parameter_Index);
   parameter_Index = parameter_Index + numel(Params.tt_PctInCareANVWhoAreOnART_1);
   
   Params.tt_PctInCareANVWhoAreOnART_2to4 = ExcelValues_AllParameters(parameter_Index);
   parameter_Index = parameter_Index + numel(Params.tt_PctInCareANVWhoAreOnART_2to4);
   
   Params.tt_PctInCareANVWhoAreOnART_5 = ExcelValues_AllParameters(parameter_Index);
   parameter_Index = parameter_Index + numel(Params.tt_PctInCareANVWhoAreOnART_5);
   
% NUMBER OF PWID SERVED BY SYRINGE SERVICES PROGRAMS (SSPS)

    % YEARS # SERVED BY SSPs  CHANGES
    Params.tt_SEP_YrSet1Begins =  ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_SEP_YrSet1Begins);
   
    Params.tt_SEP_YrSet2Begins =  ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_SEP_YrSet2Begins);
    
    Params.tt_SEP_YrSet3Begins =  ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_SEP_YrSet3Begins);
   
    Params.tt_SEP_YrSet4Begins =  ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_SEP_YrSet4Begins);
    
    Params.tt_SEP_YrSet5Begins =  ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_SEP_YrSet5Begins);   

    
    store_tt_numPWIDServedbySEP(numRace,Params.numSEPTimePeriods) = 0;
    store_tt_numPWIDServedbySEP_nonYMSMscenario(numRace,Params.numSEPTimePeriods) = 0;
    for nSet = 1:Params.numSEPTimePeriods
        % Read in from Excel [3x5]
        store_tt_numPWIDServedbySEP(:,nSet) = ExcelValues_AllParameters(parameter_Index: parameter_Index+numRace-1);
        parameter_Index = parameter_Index + numel(store_tt_numPWIDServedbySEP(:,nSet));
        
    end

    for nSet = 1:Params.numSEPTimePeriods
        % Read in from Excel for allocation-based when targeting YMSM only [3x5]
        store_tt_numPWIDServedbySEP_nonYMSMscenario(:,nSet) = ExcelValues_AllParameters(parameter_Index: parameter_Index+numRace-1);
        parameter_Index = parameter_Index + numel(store_tt_numPWIDServedbySEP_nonYMSMscenario(:,nSet));
    end
    
    if and(Params.TargetIntnsToYMSM == 1, Params.tt_progressionSource == 3)
        
        SEPinitloops = 2;

    else

        SEPinitloops = 1;

    end

    for loop = 1:SEPinitloops

        switch loop
            case 1
                refcaseSEPinit = store_tt_numPWIDServedbySEP;                
            case 2
                refcaseSEPinit = store_tt_numPWIDServedbySEP_nonYMSMscenario;           
        end

        for nSet = 1:Params.numSEPTimePeriods

            tempVariable = refcaseSEPinit(:,nSet);            

            switch nSet
                case 1
                    switch loop
                        case 1
                            Params.tt_numPWIDServedbySEP_Set1 = tempVariable;
                        case 2
                            Params.tt_numPWIDServedbySEP_Set1_nonYMSM = tempVariable;
                    end
                case 2
                    switch loop
                        case 1
                            Params.tt_numPWIDServedbySEP_Set2 = tempVariable;
                        case 2
                            Params.tt_numPWIDServedbySEP_Set2_nonYMSM = tempVariable;
                    end    
                case 3
                    switch loop
                        case 1
                            Params.tt_numPWIDServedbySEP_Set3 = tempVariable;
                        case 2
                            Params.tt_numPWIDServedbySEP_Set3_nonYMSM = tempVariable;
                    end    
                case 4
                    switch loop
                        case 1
                            Params.tt_numPWIDServedbySEP_Set4 = tempVariable;
                        case 2
                            Params.tt_numPWIDServedbySEP_Set4_nonYMSM = tempVariable;
                    end    
                case 5
                    switch loop
                        case 1
                            Params.tt_numPWIDServedbySEP_Set5 = tempVariable;
                        case 2
                            Params.tt_numPWIDServedbySEP_Set5_nonYMSM = tempVariable;
                    end    
            end

        end
    end
   

% PERCENTAGE OF HIGH-RISK THAT ARE PREP-ELIGIBLE
    
    % Read in from Excel [3x1]
    store_prepPctEligible(:,1) = ExcelValues_AllParameters(parameter_Index:parameter_Index+numPop-1);
    parameter_Index = parameter_Index + numel(store_prepPctEligible);

    % Apply to full Param [273x1], only high-risk
    Params.tt_prepPctEligible = ...
        (Params.popIndicator * store_prepPctEligible).*Params.riskLevelIndicator(:,Params.risk_Casual);
    
   
% YEARS PrEP INITIATION CHANGES
    Params.tt_PrEPInit_YrSet1Begins =  ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_PrEPInit_YrSet1Begins);
   
    Params.tt_PrEPInit_YrSet2Begins =  ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_PrEPInit_YrSet2Begins);
    
    Params.tt_PrEPInit_YrSet3Begins =  ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_PrEPInit_YrSet3Begins);
   
    Params.tt_PrEPInit_YrSet4Begins =  ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_PrEPInit_YrSet4Begins);
    
    Params.tt_PrEPInit_YrSet5Begins =  ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.tt_PrEPInit_YrSet5Begins);   
   
    
% PrEP INITIATION RATES %JCPrEPUpdate: Updated to read in HET rates by sex
    % By transmission group and by set of PrEP years    
    % Read in from Excel [12x5]
        % Dim1: transmission groups (+1 since HETs are strat by sex) and race/ethnicity, Dim2: PrEP periods
        store_prepInitRates((numPop + 1) * numRace,Params.numPrEPTimePeriods) = 0;
        for nSet = 1:Params.numPrEPTimePeriods
            store_prepInitRates(:,nSet) = ExcelValues_AllParameters(parameter_Index:parameter_Index + ((numPop + 1) * numRace) - 1);
            parameter_Index = parameter_Index + numel(store_prepInitRates(:,nSet));
        end
    
    % initiation rates for non-YMSM when using allocation-based progression for intns to YMSM only    
    % Read in from Excel [4x5]
        % Dim1: transmission groups (+1 since HETs are strat by sex) and race/ethnicity, Dim2: PrEP periods
        store_prepInitRates_nonYMSM((numPop + 1) * numRace,Params.numPrEPTimePeriods) = 0;
        for nSet = 1:Params.numPrEPTimePeriods
            store_prepInitRates_nonYMSM(:,nSet) = ExcelValues_AllParameters(parameter_Index:parameter_Index + ((numPop + 1) * numRace) - 1);
            parameter_Index = parameter_Index + numel(store_prepInitRates_nonYMSM(:,nSet));
        end    
    

    % Apply to [273x1]
        % Steps
            % Loop through each set of 4 PrEP init rates
            % Apply to full [273x1] variable
            % Only include high-risk subpopulations
            % Loop through twice if allocation-based progression targeted
            % to YMSM only
            
        if and(Params.TargetIntnsToYMSM == 1, Params.tt_progressionSource == 3)
        
            PrEPinitloops = 2;
        
        else
        
            PrEPinitloops = 1;
        
        end
        
        for loop = 1:PrEPinitloops
        
            switch loop
                case 1
                    refcasePrEPinit = store_prepInitRates;                
                case 2
                    refcasePrEPinit = store_prepInitRates_nonYMSM;           
            end
    
            for nSet = 1:Params.numPrEPTimePeriods

                tempVariable(1:Params.numStrats,1) = 0;
                for nPop_HETbySex = 1:numPop_plusHETbySex
                    for nRace = 1:numRace
                        tempVariable = tempVariable + refcasePrEPinit(numRace * (nPop_HETbySex - 1) + nRace,nSet) .* Params.popIndicator_withHETbySex(:,nPop_HETbySex) .* Params.raceIndicator(:,nRace);
                    end
                end
                tempVariable = tempVariable .* Params.riskLevelIndicator(:,Params.risk_Casual);

                %tempVariable = ...
                %    Params.popIndicator_withHETbySex * refcasePrEPinit(:,nSet) .* Params.riskLevelIndicator(:,Params.risk_Casual);

                switch nSet
                    case 1
                        switch loop
                            case 1
                                Params.tt_PrEPInitRateOrPctOnPrEP_Set1 = tempVariable;
                            case 2
                                Params.tt_PrEPInitRate_Set1_nonYMSM = tempVariable .* Params.nonYMSMIndicator;
                        end
                    case 2
                        switch loop
                            case 1
                                Params.tt_PrEPInitRateOrPctOnPrEP_Set2 = tempVariable;
                            case 2
                                Params.tt_PrEPInitRate_Set2_nonYMSM = tempVariable .* Params.nonYMSMIndicator;
                        end    
                    case 3
                        switch loop
                            case 1
                                Params.tt_PrEPInitRateOrPctOnPrEP_Set3 = tempVariable;
                            case 2
                                Params.tt_PrEPInitRate_Set3_nonYMSM = tempVariable .* Params.nonYMSMIndicator;
                        end    
                    case 4
                        switch loop
                            case 1
                                Params.tt_PrEPInitRateOrPctOnPrEP_Set4 = tempVariable;
                            case 2
                                Params.tt_PrEPInitRate_Set4_nonYMSM = tempVariable .* Params.nonYMSMIndicator;
                        end    
                    case 5
                        switch loop
                            case 1
                                Params.tt_PrEPInitRateOrPctOnPrEP_Set5 = tempVariable;
                            case 2
                                Params.tt_PrEPInitRate_Set5_nonYMSM = tempVariable .* Params.nonYMSMIndicator;
                        end    
                end

            end
        end

% PERCENT ON INJECTABLE (VERSUS ORAL) PrEP %JCPrEPUpdate: New code to read in new
% inputs
    % By transmission group and by set of PrEP years    
    % Read in from Excel [4x5]
        % Dim1: transmission groups (+1 since HETs are strat by sex), Dim2: PrEP periods
        store_preppctInj(numPop + 1,Params.numPrEPTimePeriods) = 0;
        for nSet = 1:Params.numPrEPTimePeriods
            store_preppctInj(:,nSet) = ExcelValues_AllParameters(parameter_Index:parameter_Index + numPop);
            parameter_Index = parameter_Index + numel(store_preppctInj(:,nSet));
            
            % Apply to [273x1]
            % Steps
                % Loop through each set of 4 percent on injectable PrEP values
                % Apply to full [273x1] variable
                % Only include high-risk subpopulations
            tempVariable = ...
                Params.popIndicator_withHETbySex * store_preppctInj(:,nSet) .* Params.riskLevelIndicator(:,Params.risk_Casual);
            
                switch nSet
                    case 1                        
                        Params.tt_PctonInjectPrEP_Set1 = tempVariable;                           
                    case 2
                        Params.tt_PctonInjectPrEP_Set2 = tempVariable;    
                    case 3
                        Params.tt_PctonInjectPrEP_Set3 = tempVariable;   
                    case 4
                        Params.tt_PctonInjectPrEP_Set4 = tempVariable;
                    case 5
                        Params.tt_PctonInjectPrEP_Set5 = tempVariable;
                end
            
        end                                               
        
% PERCENT WITH HIGH ADHERENCE %JCPrEPUpdate: New code to read in new
% inputs
    % By transmission group and by set of PrEP years    
    % Read in from Excel                
         % Oral [3x1]    
            store_PctHighAdherence_OralPrEP(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + numPop - 1));
                parameter_Index = parameter_Index + numPop;
    
         % Apply to [273x1] matrix
            Params.PctHighAdherence_OralPrEP = Params.popIndicator * store_PctHighAdherence_OralPrEP .* Params.riskLevelIndicator(:,Params.risk_Casual);
            
         % Injectable [3x1]    
            store_PctHighAdherence_InjectPrEP(:,1) = ...
                ExcelValues_AllParameters(parameter_Index :(parameter_Index + numPop - 1));
                parameter_Index = parameter_Index + numPop;
    
         % Apply to [273x1] matrix
            Params.PctHighAdherence_InjectPrEP = Params.popIndicator * store_PctHighAdherence_InjectPrEP .* Params.riskLevelIndicator(:,Params.risk_Casual);   
        
% DROP OUT: [SUS-ON-PrEP] -> [SUS-NO-PrEP]   %JCPrEPUpdate: Modified to read in values by PrEP delivery and adherence level      
        
     % Read in rates from Excel 
        % Oral
        % High adherence
        Params.tt_rateDropOffPrEP_Oral_HighAdherence = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(Params.tt_rateDropOffPrEP_Oral_HighAdherence);
        % Low adherence
        Params.tt_rateDropOffPrEP_Oral_LowAdherence = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(Params.tt_rateDropOffPrEP_Oral_LowAdherence);
        
        % Injectable
        % High adherence
        Params.tt_rateDropOffPrEP_Inject_HighAdherence = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(Params.tt_rateDropOffPrEP_Inject_HighAdherence);
        % Low adherence
        Params.tt_rateDropOffPrEP_Inject_LowAdherence = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(Params.tt_rateDropOffPrEP_Inject_LowAdherence);

        
% COVID Inputs - Added by MClinkscales 6/2/2022; JC updated on 8/29/2022

    % Effects on annual testing rates: uninfected people
    
        % Dim1: race/ethn, Dim2: COVID time periods (t2-t4)
        COVID_percentreduc_testrate(numRace,Params.numPeriod-2) = 0;
        for nSet = 1:Params.numPeriod-2
            COVID_percentreduc_testrate(:,nSet) = ExcelValues_AllParameters(parameter_Index:parameter_Index + Params.numPeriod - 3);
            parameter_Index = parameter_Index + numPeriod-2;
        end
        
        % Apply to 273 subpopulations using race indicator (results in 273
        % x 3) and turn on/off effect based on model settings
        Params.COVID_percentreduc_testrate_Uninf = (Params.raceIndicator * COVID_percentreduc_testrate) * Params.ApplyCOVIDEffects_ContBeh;        

    % Effects on annual testing rates: undiagnosed PWH
    
        % Dim1: race/ethn, Dim2: COVID time periods (t2-t4)
        COVID_percentreduc_testrate(numRace,Params.numPeriod-2) = 0;
        for nSet = 1:Params.numPeriod-2
            COVID_percentreduc_testrate(:,nSet) = ExcelValues_AllParameters(parameter_Index:parameter_Index + Params.numPeriod - 3);
            parameter_Index = parameter_Index + numPeriod-2;
        end
        
        % Apply to 273 subpopulations using race indicator (results in 273 x 3)
        Params.COVID_percentreduc_testrate_PWH = (Params.raceIndicator * COVID_percentreduc_testrate) * Params.ApplyCOVIDEffects_ContBeh;
        
    % Effects on ART
        COVID_percentreduc_initART(1,1:Params.numPeriod-2) = ExcelValues_AllParameters(parameter_Index:parameter_Index + Params.numPeriod - 3);
        parameter_Index = parameter_Index + numel(COVID_percentreduc_initART);
        
        % Turn on/off effects based on model settings
        Params.COVID_percentreduc_initART = COVID_percentreduc_initART * Params.ApplyCOVIDEffects_ContBeh;
                
    % Effects on backwards continuum progression        
        
        % ANV to LTC (row 4 to row 3)
        % Dim1: age groups, Dim2: COVID time periods (t2-t4)
        COVID_percentincr_DropOutARTfromANV(Params.numAge,Params.numPeriod-2) = 0;       
        
        for nSet = 1:Params.numPeriod-2
            COVID_percentincr_DropOutARTfromANV(:,nSet) = ExcelValues_AllParameters(parameter_Index:parameter_Index + Params.numAge - 1);
            parameter_Index = parameter_Index + Params.numAge;
        end
        
        % Apply to 273 subpopulations using age indicator (results in 273 x
        % 3) and turn on/off effect based on model settings
        Params.COVID_percentincr_DropOutARTfromANV = (Params.ageIndicator * COVID_percentincr_DropOutARTfromANV) * Params.ApplyCOVIDEffects_ContBeh;
        
        % VLS to ANV (row 5 to row 4)
        % Dim1: age groups, Dim2: COVID time periods (t2-t4)
        COVID_percentincr_LoseVLStoANV(Params.numAge,Params.numPeriod-2) = 0;
        for nSet = 1:Params.numPeriod-2
            COVID_percentincr_LoseVLStoANV(:,nSet) = ExcelValues_AllParameters(parameter_Index:parameter_Index + Params.numAge - 1);
            parameter_Index = parameter_Index + Params.numAge;
        end
        
        % Apply to 273 subpopulations using age indicator (results in 273 x
        % 3) and turn on/off effect based on model settings
        Params.COVID_percentincr_LoseVLStoANV = (Params.ageIndicator * COVID_percentincr_LoseVLStoANV) * Params.ApplyCOVIDEffects_ContBeh;
        
    % Effects on PrEP initiation
        COVID_percentreduc_PrEPinit(1,1:Params.numPeriod-2) = ExcelValues_AllParameters(parameter_Index:parameter_Index + Params.numPeriod - 3);
        parameter_Index = parameter_Index + numel(COVID_percentreduc_PrEPinit); 
        
        % Turn on/off effects based on model settings
        Params.COVID_percentreduc_PrEPinit = COVID_percentreduc_PrEPinit * Params.ApplyCOVIDEffects_ContBeh;
        
    % Effects on sexual contacts
    
        % Dim1: transmission groups, Dim2: COVID time periods (t2-t4)
        store_COVID_percentreduc_sexcont(numPop,Params.numPeriod-2) = 0;
        for nSet = 1:Params.numPeriod-2
            store_COVID_percentreduc_sexcont(:,nSet) = ExcelValues_AllParameters(parameter_Index:parameter_Index + Params.numPeriod - 3);
            parameter_Index = parameter_Index + numPeriod-2;
        end
        
        % Apply to 273 subpopulations using pop indicator (results in 273 x
        % 3) and turn on/off effect based on model settings
        Params.COVID_percentreduc_sexcont = (Params.popIndicator * store_COVID_percentreduc_sexcont.* Params.riskLevelIndicator(:,Params.risk_Casual)) * Params.ApplyCOVIDEffects_ContBeh;
        
       
    % Effects on mortality
    
        % Dim1: age groups, Dim2: COVID time periods (t2-t4)
        COVID_percentincr_Mortality(Params.numAge,Params.numPeriod-2) = 0;
        for nSet = 1:Params.numPeriod-2
            COVID_percentincr_Mortality(:,nSet) = ExcelValues_AllParameters(parameter_Index:parameter_Index + Params.numAge - 1);
            parameter_Index = parameter_Index + Params.numAge;
        end
        
        % Apply to 273 subpopulations using age indicator (results in 273 x
        % 3) and turn on/off effect based on model settings
        Params.COVID_percentincr_Mortality = (Params.ageIndicator * COVID_percentincr_Mortality) * Params.ApplyCOVIDEffects_Mort;
        
        
        
%% 8. Define Behavior inputs (except mixing and AI)
   
% CONDOM EFFICACY
    % Percent of contacts in which condom provides effective protection
    
    % HETs
    Params.behav_condomEfficacyHET(1) = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.behav_condomEfficacyHET);
    
    % MSMs
    Params.behav_condomEfficacyMSM(1) = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.behav_condomEfficacyMSM);

% ANNUAL NUMBER OF SEXUAL CONTACTS PER UNINFECTED PERSON (added by JC on
% 12/5/2017) - only used for Calib_UpdateParams 
    Params.store_behav_MSMOverallNumContacts(:,1) =...
        ExcelValues_AllParameters(parameter_Index: parameter_Index+ Params.numAge-1);
        parameter_Index = parameter_Index + Params.numAge;
        
% ANNUAL NUMBER OF SEXUAL CONTACTS (BASE)    
    
    %Read in subgroups
    Params.store_HETMaleLow(:,1) = ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numRace - 1);
        parameter_Index = parameter_Index + Params.numRace;
    Params.store_HETFemaleLow(:,1) = ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numRace - 1);
        parameter_Index = parameter_Index + Params.numRace;
        %adding age stratification. Laurel 4/28/17
    store_MSMLow_Black(:,1)= ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numAge - 1);
        parameter_Index = parameter_Index + Params.numAge;
    store_MSMLow_Hisp(:,1)= ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numAge - 1);
        parameter_Index = parameter_Index + Params.numAge;
    store_MSMLow_Other(:,1)= ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numAge - 1);
        parameter_Index = parameter_Index + Params.numAge;  
        %
    Params.store_HETMaleHigh(:,1) = ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numRace - 1);
        parameter_Index = parameter_Index + Params.numRace;
    Params.store_HETFemaleHigh(:,1) = ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numRace - 1);
        parameter_Index = parameter_Index + Params.numRace;
        %adding age stratification. Laurel 4/28/17
    store_MSMHigh_Black(:,1)= ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numAge - 1);
        parameter_Index = parameter_Index + Params.numAge;
    store_MSMHigh_Hisp(:,1)= ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numAge - 1);
        parameter_Index = parameter_Index + Params.numAge;
    store_MSMHigh_Other(:,1)= ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numAge - 1);
        parameter_Index = parameter_Index + Params.numAge;  
        %
    Params.store_IDUMaleHigh(:,1) = ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numRace - 1);
        parameter_Index = parameter_Index + Params.numRace;
    Params.store_IDUFemaleHigh(:,1) = ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numRace - 1);
        parameter_Index = parameter_Index + Params.numRace;
        
        
 %Combining MSM race/age groups into high and low -code below doesn't work
%MSMLow = Params.ageIndicator * store_MSMLow_Black ...
  %   + Params.ageIndicator * store_MSMLow_Hisp ...
   %  + Params.ageIndicator * store_MSMLow_Other

% MSMHigh = Params.ageIndicator * store_MSMHigh_Black ...
    % + Params.ageIndicator * store_MSMHigh_Hisp ...
    % + Params.ageIndicator * store_MSMHigh_Other
 
    % Convert subgroups to a single Parameter [273x1]
    %adding age stratification. Laurel 4/28/17
    baseSexContacts = ...
        Params.raceIndicator * Params.store_HETMaleLow ...
            .* Params.popSexIndicator(:,Params.popSex_HETM) ...
            .* Params.riskLevelIndicator(:,Params.risk_Main)...        
        + Params.raceIndicator * Params.store_HETFemaleLow ...
            .* Params.popSexIndicator(:,Params.popSex_HETF) ...
            .* Params.riskLevelIndicator(:,Params.risk_Main) ...
         + Params.ageIndicator * store_MSMLow_Black ...            
           .* Params.popSexIndicator(:,Params.popSex_MSM) ...
           .* Params.riskLevelIndicator(:,Params.risk_Main) ...
           .* Params.raceIndicator(:,Params.race_B)...
        + Params.ageIndicator * store_MSMLow_Hisp ...            
          .* Params.popSexIndicator(:,Params.popSex_MSM) ...
           .* Params.riskLevelIndicator(:,Params.risk_Main) ...
           .* Params.raceIndicator(:,Params.race_H)...
        + Params.ageIndicator * store_MSMLow_Other ...            
          .* Params.popSexIndicator(:,Params.popSex_MSM) ...
          .* Params.riskLevelIndicator(:,Params.risk_Main) ... 
           .* Params.raceIndicator(:,Params.race_O)...
        + Params.raceIndicator * Params.store_HETMaleHigh ...
            .* Params.popSexIndicator(:,Params.popSex_HETM) ...
            .* Params.riskLevelIndicator(:,Params.risk_Casual) ...
        + Params.raceIndicator * Params.store_HETFemaleHigh ...
            .* Params.popSexIndicator(:,Params.popSex_HETF) ...
            .* Params.riskLevelIndicator(:,Params.risk_Casual) ...
      + Params.ageIndicator * store_MSMHigh_Black ...            
           .* Params.popSexIndicator(:,Params.popSex_MSM) ...
           .* Params.riskLevelIndicator(:,Params.risk_Casual) ...
           .* Params.raceIndicator(:,Params.race_B)...
      + Params.ageIndicator * store_MSMHigh_Hisp ...            
           .* Params.popSexIndicator(:,Params.popSex_MSM) ...
           .* Params.riskLevelIndicator(:,Params.risk_Casual) ...
           .* Params.raceIndicator(:,Params.race_H)...
       + Params.ageIndicator * store_MSMHigh_Other ...            
          .* Params.popSexIndicator(:,Params.popSex_MSM) ...
          .* Params.riskLevelIndicator(:,Params.risk_Casual) ...
          .* Params.raceIndicator(:,Params.race_O)...
        + Params.raceIndicator * Params.store_IDUMaleHigh ...
            .* Params.popSexIndicator(:,Params.popSex_IDUM) ...
            .* Params.riskLevelIndicator(:,Params.risk_Casual) ...
        + Params.raceIndicator * Params.store_IDUFemaleHigh ...
            .* Params.popSexIndicator(:,Params.popSex_IDUF) ...
            .* Params.riskLevelIndicator(:,Params.risk_Casual);
                   
 % EFFECT OF DIAGNOSIS: REDUCTION IN SEXUAL ACTS
    %Reduction in number of sexual acts for diagnosed versus undiagnosed
    
    Params.hiv_reductSexActsDueToAwareness = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.hiv_reductSexActsDueToAwareness);

 
% EFFECT OF PrEP: INCREASE IN SEXUAL CONTACTS
    % Scalar
    % Percent increase in number of sex contactsfor individuals on PrEP
    % vs. not on PrEP
    
    Params.behav_incrSexContactsDueToPrEP = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.behav_incrSexContactsDueToPrEP);  

    
% ANNUAL NUMBER OF SEXUAL CONTACTS (ADJUSTED)
    % Adjustments due to awareness and PrEP for each model time period
    %Loops through 4 periods only. Values for period 1 are also applied to period 5.
    for Period = 1:(Params.numPeriod - 1)
        
        switch Period
            case 1
                sexContacts = baseSexContacts;
            case {2, 3, 4}
                sexContacts = baseSexContacts .* (1 - Params.COVID_percentreduc_sexcont(:, Period - 1));
        end
    
        % Apply to full matrix [273x273x30]
            factor_adjNumSexContact_NoPrEP = repmat(sexContacts, [1,Params.numStrats, Params.numComparts]);

        % Apply awareness reduction
            factor_adjNumSexContact_NoPrEP(:,:,AwareComparts) = ...
                factor_adjNumSexContact_NoPrEP(:,:,AwareComparts) * (1 - Params.hiv_reductSexActsDueToAwareness);

        % Apply PrEP reduction
            factor_adjNumSexContact_PrEP = ...
                factor_adjNumSexContact_NoPrEP * (1 + Params.behav_incrSexContactsDueToPrEP);
            
       switch Period
           case 1
                Params.factor_adjNumSexContact_NoPrEP_1and5 = factor_adjNumSexContact_NoPrEP;
                Params.factor_adjNumSexContact_PrEP_1and5 = factor_adjNumSexContact_PrEP;
           case 2
                Params.factor_adjNumSexContact_NoPrEP_2 = factor_adjNumSexContact_NoPrEP;
                Params.factor_adjNumSexContact_PrEP_2 = factor_adjNumSexContact_PrEP;
           case 3
                Params.factor_adjNumSexContact_NoPrEP_3 = factor_adjNumSexContact_NoPrEP;
                Params.factor_adjNumSexContact_PrEP_3 = factor_adjNumSexContact_PrEP;
           case 4
                Params.factor_adjNumSexContact_NoPrEP_4 = factor_adjNumSexContact_NoPrEP;
                Params.factor_adjNumSexContact_PrEP_4 = factor_adjNumSexContact_PrEP;
        end     
        
    end
    
% PERCENTAGE OF PWID ACTIVELY INJECTING
    % scalar
    Params.behav_pctActivePWID = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.behav_pctActivePWID);    
        
% NUMBER OF INJECTIONS PER YEAR PER ACTIVE PWID
    % scalar

    store_behav_numInjectionsPerYear = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(store_behav_numInjectionsPerYear);
    
    % Calculate number of injections per year for all PWID
    Params.behav_numInjectionsPerYear = store_behav_numInjectionsPerYear * Params.behav_pctActivePWID;

% EFFECT OF PREP: INCREASE IN NEEDLE SHARING
    % Scalar
    % Percent increase in number of needles shared for individuals on PrEP
    % vs. not on PrEP
    
    Params.behav_incrNeedleShareDueToPrEP = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.behav_incrNeedleShareDueToPrEP);
    
    
% PERCENTAGE OF INJECTIONS THAT ARE SHARED 
    % scalar
    behav_pctNeedlesShared = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(behav_pctNeedlesShared);

% DEFAULT NUMBER OF NEEDLE PARTNERS
    % Default number of needle-sharing partners for IDUs
    % This parameter is used to keep the number of needles shared per
    % partner constant when the number of needle-sharing partners is
    % changed.
    
    Params.behav_defaultpartnersNeedle = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.behav_defaultpartnersNeedle);        
    
% NUMBER OF NEEDLE PARTNERS
    % Annual number of needle-sharing partners for IDUs
    
    Params.behav_partnersNeedle = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.behav_partnersNeedle);
    
    
% CALCULATED: NUMBER OF NEEDLES SHARED PER PARTNER (BASE)
    % Will be adjusted later adjusted based on awareness
    % Calculated as: 
        % [Number injections per year] * [Pct inj shared] / [Number of partners]
    
    if Params.behav_defaultpartnersNeedle > 0

        behav_baseNeedlesPerPartner = ...
               Params.behav_numInjectionsPerYear ...
            .* behav_pctNeedlesShared ...
            / Params.behav_defaultpartnersNeedle;
        
    else
        Params.behav_defaultpartnersNeedle = 0;
    end

    
% PERCENT OF CONTACTS PROTECTED WITH A CONDOM

    %VI

        % Percent of vaginal contacts with a condom
        % Re-coded by JC on 7/31/2018
        
        % Pre-allocate    
        store_basePctCondomUseV_Black(Params.numAge,Params.numPop) = 0;
        store_basePctCondomUseV_Hisp(Params.numAge,Params.numPop) = 0;
        store_basePctCondomUseV_Other(Params.numAge,Params.numPop) = 0;
        behav_basePctCondomUseV(Params.numStrats, 1) = 0;
            
        % Loop through transmission groups    
        for p=1:numPop
            store_basePctCondomUseV_Black(:,p) = ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numAge - 1);
            parameter_Index = parameter_Index + Params.numAge;
            store_basePctCondomUseV_Hisp(:,p) = ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numAge - 1);
            parameter_Index = parameter_Index + Params.numAge;
            store_basePctCondomUseV_Other(:,p) = ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numAge - 1);
            parameter_Index = parameter_Index + Params.numAge;
            
            % [273x1] VI condom use for population p
            store_basePctCondomV_pop = ...
                ((Params.ageIndicator * store_basePctCondomUseV_Black(:,p) ...
                .* Params.raceIndicator(:,Params.race_B)  ...
                + Params.ageIndicator * store_basePctCondomUseV_Hisp(:,p) ...
                .* Params.raceIndicator(:,Params.race_H) ...
                + Params.ageIndicator * store_basePctCondomUseV_Other(:,p) ...
                .* Params.raceIndicator(:,Params.race_O))...
                .* Params.popIndicator(:,p));
            
            % Add to single parameter [273x1]
            behav_basePctCondomUseV = behav_basePctCondomUseV + store_basePctCondomV_pop;
            
        end
                    

       % Convert subgroups to a single Parameter [273x1]
       %behav_basePctCondomUseV = ...
       %    Params.raceIndicator * store_condomUseV_HET .* Params.popIndicator(:,Params.pop_HET) + ...
       %    Params.raceIndicator * store_condomUseV_MSM .* Params.popIndicator(:,Params.pop_MSM) + ...
       %    Params.raceIndicator * store_condomUseV_IDU .* Params.popIndicator(:,Params.pop_IDU);

    % AI

        % Import
            % Male-female
            % Female-male
            % Male-male
            %Modified 11/07/17 by JC to stratify by race for MSM
            
        store_basePctCondomA_HET(:,1) = ExcelValues_AllParameters(parameter_Index: parameter_Index + numSex - 1);
        parameter_Index = parameter_Index + numel(store_basePctCondomA_HET);
        
        store_basePctCondomA_MSM_Black(:,1) = ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numAge - 1);
        parameter_Index = parameter_Index + numel(store_basePctCondomA_MSM_Black);
        store_basePctCondomA_MSM_Hisp(:,1) = ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numAge - 1);
        parameter_Index = parameter_Index + numel(store_basePctCondomA_MSM_Hisp);
        store_basePctCondomA_MSM_Other(:,1) = ExcelValues_AllParameters(parameter_Index: parameter_Index + Params.numAge - 1);
        parameter_Index = parameter_Index + numel(store_basePctCondomA_MSM_Other);


        % [273x273] Male-Female and Female-Male AI condom use
        storeHET_basePctCondomA = ...
            Params.sexIndicator(:,Params.sex_Male)*transpose(Params.sexIndicator(:,Params.sex_Female).*store_basePctCondomA_HET(Params.sex_Male)) + ...
            Params.sexIndicator(:,Params.sex_Female)*transpose(Params.sexIndicator(:,Params.sex_Male).*store_basePctCondomA_HET(Params.sex_Female));    
            % save value for calibUpdateParams
        Params.storeHET_basePctCondomA = storeHET_basePctCondomA;
        % [273x273] Male-Male AI condom use
        storeMSM_basePctCondomA = ...
            ((Params.ageIndicator * store_basePctCondomA_MSM_Black ...
            .* Params.raceIndicator(:,Params.race_B)  ...
            + Params.ageIndicator * store_basePctCondomA_MSM_Hisp ...
            .* Params.raceIndicator(:,Params.race_H) ...
            + Params.ageIndicator * store_basePctCondomA_MSM_Other ...
            .* Params.raceIndicator(:,Params.race_O))...
            .* Params.popIndicator(:, Params.pop_MSM))...
            * transpose(Params.popIndicator(:, Params.pop_MSM));


        % Sum values for HET and MSM-MSM partnerships [273x273]
        store_basePctCondomA = storeHET_basePctCondomA + storeMSM_basePctCondomA;

        % Expand to [273x273x30]
        behav_basePctCondomA = repmat(store_basePctCondomA,[1,1,Params.numComparts]);


% PERCENT MSM RECEPTIVE                                            
    %Percent of MSM contacts receptive. 
    %Stratified by age 4/28/17. Laurel
    
    store_behav_pctMSMRec_Black(:,1) =...
        ExcelValues_AllParameters(parameter_Index: parameter_Index+ Params.numAge-1);
        parameter_Index = parameter_Index + Params.numAge;
    store_behav_pctMSMRec_Hisp(:,1) =...
        ExcelValues_AllParameters(parameter_Index: parameter_Index+ Params.numAge-1);
        parameter_Index = parameter_Index + Params.numAge;
    store_behav_pctMSMRec_Other(:,1) =...
        ExcelValues_AllParameters(parameter_Index: parameter_Index+ Params.numAge-1);
        parameter_Index = parameter_Index + Params.numAge;
    
        
    behav_pctMSMRec = ...
        Params.ageIndicator * store_behav_pctMSMRec_Black ...
            .* Params.popSexIndicator(:,Params.popSex_MSM) ...
            .* Params.raceIndicator(:,Params.race_B) ...
        + Params.ageIndicator * store_behav_pctMSMRec_Hisp ...
            .* Params.popSexIndicator(:,Params.popSex_MSM) ...
            .* Params.raceIndicator(:,Params.race_H) ...
        + Params.ageIndicator * store_behav_pctMSMRec_Other ...
            .* Params.popSexIndicator(:,Params.popSex_MSM) ...
            .* Params.raceIndicator(:,Params.race_O);


% EFFECT OF DIAGNOSIS: INCREASE IN PCT CONDOM USE
    %Increase in percent of sexual contacts protected with a condom for the
    %diagnosed versus undiagnosed
    
    behav_incrCondomUseDueToDiag_1to4= ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(behav_incrCondomUseDueToDiag_1to4);
    
    behav_incrCondomUseDueToDiag_5= ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(behav_incrCondomUseDueToDiag_5);
    
    % save value for calibUpdateParams
    Params.behav_incrCondomUseDueToDiag_1to4 = behav_incrCondomUseDueToDiag_1to4;
    Params.behav_incrCondomUseDueToDiag_5 = behav_incrCondomUseDueToDiag_5;
    
% EFFECT OF PrEP: DECREASE IN PCT CONDOM USE
    % Decrease in percent of sexual contacts protected with a condom for
    % people who use PrEP vs. don't
    behav_decrCondomUseDueToPrEP= ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(behav_decrCondomUseDueToPrEP);
    
    % save value for calibUpdateParams
    Params.behav_decrCondomUseDueToPrEP = behav_decrCondomUseDueToPrEP;
    
    
% PERCENT OF CONTACTS PROTECTED WITH A CONDOM (ADJUSTED)
    % Adjusted based on awareness
    
    % Based on increase in sexual acts protected with a condom when the 
    % infected partner has been diagnosed    
    % Updated by MClinksales for time periods 4 & 5
    
    % VI Contacts
    
        % Expand base number of contacts to be [273x273x30]
        expand_basePctCondomUseV = repmat(behav_basePctCondomUseV,[1,Params.numStrats,Params.numComparts]);

        % Apply increase in condom use based on awareness
            store_adjPctCondomUseV_NoPrEP_1to4 = expand_basePctCondomUseV;
            store_adjPctCondomUseV_NoPrEP_5 = expand_basePctCondomUseV;
            store_adjPctCondomUseV_NoPrEP_1to4(:,:,AwareComparts) = ...
                store_adjPctCondomUseV_NoPrEP_1to4(:,:,AwareComparts) * (1 + behav_incrCondomUseDueToDiag_1to4);
            store_adjPctCondomUseV_NoPrEP_5(:,:,AwareComparts) = ...
                store_adjPctCondomUseV_NoPrEP_5(:,:,AwareComparts) * (1 + behav_incrCondomUseDueToDiag_5);

        % Apply decrease due to PrEP
            store_adjPctCondomUseV_PrEP_1to4 = ...
                store_adjPctCondomUseV_NoPrEP_1to4 * (1-behav_decrCondomUseDueToPrEP);            
            store_adjPctCondomUseV_PrEP_5 = ...
                store_adjPctCondomUseV_NoPrEP_5 * (1-behav_decrCondomUseDueToPrEP);
            
        % Limit maximum pct condom usage to 1
            Params.factor_adjPctCondomUseV_NoPrEP_1to4 = min(store_adjPctCondomUseV_NoPrEP_1to4,1);
            Params.factor_adjPctCondomUseV_PrEP_1to4 = min(store_adjPctCondomUseV_PrEP_1to4,1);            
            Params.factor_adjPctCondomUseV_NoPrEP_5 = min(store_adjPctCondomUseV_NoPrEP_5,1);
            Params.factor_adjPctCondomUseV_PrEP_5 = min(store_adjPctCondomUseV_PrEP_5,1);

     % AI Contacts
       
        % Apply increase due to awareness
        store_factor_adjPctCondomUseA_NoPrEP_1to4 = behav_basePctCondomA .* ...
            (1 + inf_awareComparts_273x273x30 * behav_incrCondomUseDueToDiag_1to4);
        store_factor_adjPctCondomUseA_NoPrEP_5 = behav_basePctCondomA .* ...
            (1 + inf_awareComparts_273x273x30 * behav_incrCondomUseDueToDiag_5);

        % Apply decrease due to PrEP
        store_factor_adjPctCondomUseA_PrEP_1to4 = ...
            store_factor_adjPctCondomUseA_NoPrEP_1to4 * (1 - behav_decrCondomUseDueToPrEP);
        store_factor_adjPctCondomUseA_PrEP_5 = ...
            store_factor_adjPctCondomUseA_NoPrEP_5 * (1 - behav_decrCondomUseDueToPrEP);
        
        % Limit maximum pct condom usage to 1
        Params.factor_adjPctCondomUseA_NoPrEP_1to4 = min(store_factor_adjPctCondomUseA_NoPrEP_1to4, 1);
        Params.factor_adjPctCondomUseA_PrEP_1to4 = min(store_factor_adjPctCondomUseA_PrEP_1to4, 1);
        Params.factor_adjPctCondomUseA_NoPrEP_5 = min(store_factor_adjPctCondomUseA_NoPrEP_5, 1);
        Params.factor_adjPctCondomUseA_PrEP_5 = min(store_factor_adjPctCondomUseA_PrEP_5, 1);
        

% EFFECT OF DIAGNOSIS: REDUCTION IN SHARED NEEDLES
    % Decrease in needle sharing for diagnosed versus undiagnosed
    
    Params.behav_needleReductDueToAware_1to4 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.behav_needleReductDueToAware_1to4);
    
    Params.behav_needleReductDueToAware_5 = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.behav_needleReductDueToAware_5);
    
% PERCENT REDUCTION IN NEEDLE-SHARING PARTNERS FOR DIAGNOSED VS UNDIAGNOSED

    Params.pctReductDiagNeedlePartners = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
% ANNUAL NUMBER OF SEXUAL PARTNERS 
    %Params.D1: race
    %Params.D2: risk level
    %Params.D3: sex
    %Params.D4: population
    
    % ET note: if you ever need to adjust this parameter to include
    % different stratifications, I would scrap this code and rewrite it.
    % It's one of the last hold-overs from when we had a 
    % different dimension for each stratification type (we switched in 2013
    % to have a single dimension of 273 cohorts).
    
behav_partnersSex_rlsp(numRace,numRiskLevel,numSex,numPop) = 0;


% EFFECT OF CONDOM USE IN REDUCING TRANSMISSION 
    % This code will combine the receptive and insertive MSM
    % condom efficacy parameters and the HET parameters into
    % on large vector.
    
    Params.behav_condomEfficacy = Params.behav_condomEfficacyHET(1)* Params.popIndicator(:,1)'...
                                 +Params.behav_condomEfficacyMSM(1)* Params.popIndicator(:,2)' ;
    
    

%loop through the populations
for a = 1:numPop
    %loop through the sexes
    for b = 1:numSex
        %loop through the risk levels
        for c = 1:numRiskLevel
        behav_partnersSex_rlsp(:,c,b,a) = ExcelValues_AllParameters(( ...
            parameter_Index + numRace*(c-1)+numRace*numRiskLevel*(b-1)+...
            numRace*numRiskLevel*numSex*(a-1)): (parameter_Index + numRace-1) ...
            + numRace*(c-1)+numRace*numRiskLevel*(b-1)+numRace*...
            numRiskLevel*numSex*(a-1));
        end
    end
end
parameter_Index = parameter_Index + numel(behav_partnersSex_rlsp);

%Save # partners for MSM for calibration (added by JC on 12/05/2017)
behav_partnersSex_MSM(numRace,numRiskLevel) =  0;
behav_partnersSex_MSM(:,:) = behav_partnersSex_rlsp(:,:,1,2);
Params.behav_partnersSex_MSM = behav_partnersSex_MSM;


%Restratify behav_partnersSex_rlsp to behav_partnersSex
behav_partnersSex(Params.numStrats)=0;

storePS(Params.numStrats,numRiskLevel,numSex,numPop)=0;
storePS2(Params.numStrats,Params.numStrats,numSex,numPop)=0;
storePS3(Params.numStrats,numSex,numPop)=0;
storePS4(Params.numStrats,Params.numStrats,numPop)=0;
storePS5(Params.numStrats,numPop)=0;

for s = 1:numSex
    for p = 1:numPop
        storePS(:,:,s,p) = Params.raceIndicator*behav_partnersSex_rlsp(:,:,s,p);
    end
end
for s = 1:numSex
    for p = 1:numPop
        storePS2(:,:,s,p)= storePS(:,:,s,p)*Params.riskLevelIndicator';     
    end
end
for s = 1:numSex
    for p = 1:numPop
        storePS3(:,s,p)=diag(storePS2(:,:,s,p));
    end
end
for p = 1:numPop
   storePS4(:,:,p)= storePS3(:,:,p)*Params.sexIndicator';
end
for p = 1:numPop
    storePS5(:,p)=diag(storePS4(:,:,p));
end

behav_partnersSex = diag(storePS5*Params.popIndicator');

Params.partnersSex = behav_partnersSex;

% REDUCTION IN SEXUAL PARTNERSHIPS BY AGE
    % 13-17
    behav_pctReductAgePartnersSex(1) = ExcelValues_AllParameters(parameter_Index);
    % 45-54
    behav_pctReductAgePartnersSex(2) = ExcelValues_AllParameters(parameter_Index+1);
    % 55-64
    behav_pctReductAgePartnersSex(3) = ExcelValues_AllParameters(parameter_Index+2);
    % 65+
    behav_pctReductAgePartnersSex(4) = ExcelValues_AllParameters(parameter_Index+3);
    parameter_Index = parameter_Index + numel(behav_pctReductAgePartnersSex);

% Apply reduction for ages 13-17, 45-54, 55-64, and 65+

    idx_age1 = Params.ageIndicator(:,1) == 1;
    idx_age5 = Params.ageIndicator(:,5) == 1;
    idx_age6 = Params.ageIndicator(:,6) == 1;
    idx_age7 = Params.ageIndicator(:,7) == 1;

    Params.partnersSex(idx_age1)=Params.partnersSex(idx_age1)* ...
       (1- behav_pctReductAgePartnersSex(1));
    Params.partnersSex(idx_age5)=Params.partnersSex(idx_age5)* ...
       (1- behav_pctReductAgePartnersSex(2));
    Params.partnersSex(idx_age6)=Params.partnersSex(idx_age6)* ...
       (1- behav_pctReductAgePartnersSex(3));
    Params.partnersSex(idx_age7)=Params.partnersSex(idx_age7)* ...
       (1- behav_pctReductAgePartnersSex(4));
 
% PERCENT REDUCTION IN SEXUAL PARTNERS BY RISK GROUP FOR DIAGNOSED VS UNDIAGNOSED

    Params.pctReductDiagPartners(Params.pop_HET) = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    Params.pctReductDiagPartners(Params.pop_MSM) = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    Params.pctReductDiagPartners(Params.pop_IDU) = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
%% 9. Define Mixing matrix (Behaviors sheet)

% SEXUAL MIXING: TRANSMISSION GROUP AND SEX    
    
    % Import from Excel
    pctSexMixing_TransGroupSex = ExcelValues_AllParameters(parameter_Index:parameter_Index + numTransGroupSex_SexualCombos-1);
    parameter_Index = parameter_Index + numel(pctSexMixing_TransGroupSex);
    
    % Apply to storage parameters [273x273]
    
        % Transmission group/sex combinations [273x273]
            % 1 - Pull appropriate value
            % 2 - Use indicator to apply to appropriate cells
               
% Adapted by KH on 20Sep2023 to accomodate mixing that varies for HR and LR HETs        
        behav_sexMix_LRHETM_HETF = Params.idx_LRHETM_HETF * pctSexMixing_TransGroupSex(Params.sexMix_LRHETM_HETF);
        behav_sexMix_LRHETM_IDUF = Params.idx_LRHETM_IDUF * pctSexMixing_TransGroupSex(Params.sexMix_LRHETM_IDUF);
        behav_sexMix_LRHETF_HETM = Params.idx_LRHETF_HETM * pctSexMixing_TransGroupSex(Params.sexMix_LRHETF_HETM);
        behav_sexMix_LRHETF_MSM = Params.idx_LRHETF_MSM * pctSexMixing_TransGroupSex(Params.sexMix_LRHETF_MSM);
        behav_sexMix_LRHETF_IDUM = Params.idx_LRHETF_IDUM * pctSexMixing_TransGroupSex(Params.sexMix_LRHETF_IDUM);
        behav_sexMix_HRHETM_HETF = Params.idx_HRHETM_HETF * pctSexMixing_TransGroupSex(Params.sexMix_HRHETM_HETF);
        behav_sexMix_HRHETM_IDUF = Params.idx_HRHETM_IDUF * pctSexMixing_TransGroupSex(Params.sexMix_HRHETM_IDUF);
        behav_sexMix_HRHETF_HETM = Params.idx_HRHETF_HETM * pctSexMixing_TransGroupSex(Params.sexMix_HRHETF_HETM);
        behav_sexMix_HRHETF_MSM = Params.idx_HRHETF_MSM * pctSexMixing_TransGroupSex(Params.sexMix_HRHETF_MSM);
        behav_sexMix_HRHETF_IDUM = Params.idx_HRHETF_IDUM * pctSexMixing_TransGroupSex(Params.sexMix_HRHETF_IDUM);
        behav_sexMix_MSM_HETF = Params.idx_MSM_HETF * pctSexMixing_TransGroupSex(Params.sexMix_MSM_HETF);
        behav_sexMix_MSM_MSM = Params.idx_MSM_MSM * pctSexMixing_TransGroupSex(Params.sexMix_MSM_MSM);
        behav_sexMix_MSM_IDUF = Params.idx_MSM_IDUF * pctSexMixing_TransGroupSex(Params.sexMix_MSM_IDUF);
        behav_sexMix_IDUM_HETF = Params.idx_IDUM_HETF * pctSexMixing_TransGroupSex(Params.sexMix_IDUM_HETF);
        behav_sexMix_IDUM_IDUF = Params.idx_IDUM_IDUF * pctSexMixing_TransGroupSex(Params.sexMix_IDUM_IDUF);
        behav_sexMix_IDUF_HETM = Params.idx_IDUF_HETM * pctSexMixing_TransGroupSex(Params.sexMix_IDUF_HETM);    
        behav_sexMix_IDUF_MSM = Params.idx_IDUF_MSM * pctSexMixing_TransGroupSex(Params.sexMix_IDUF_MSM);  
        behav_sexMix_IDUF_IDUM = Params.idx_IDUF_IDUM * pctSexMixing_TransGroupSex(Params.sexMix_IDUF_IDUM);


    % Create full [273x273] matrix of each trans group/sex sexual combo
        %KH incorporated HET mixing by risk level on 20Sep2023 
        behav_sexMix_popSex = ...
              behav_sexMix_LRHETM_HETF + behav_sexMix_LRHETM_IDUF ...
            + behav_sexMix_LRHETF_HETM + behav_sexMix_LRHETF_MSM + behav_sexMix_LRHETF_IDUM ...
            + behav_sexMix_HRHETM_HETF + behav_sexMix_HRHETM_IDUF ...
            + behav_sexMix_HRHETF_HETM + behav_sexMix_HRHETF_MSM + behav_sexMix_HRHETF_IDUM ...
            + behav_sexMix_MSM_HETF + behav_sexMix_MSM_MSM + behav_sexMix_MSM_IDUF ...
            + behav_sexMix_IDUM_HETF + behav_sexMix_IDUM_IDUF...
            + behav_sexMix_IDUF_HETM + behav_sexMix_IDUF_MSM + behav_sexMix_IDUF_IDUM;
        
        % Store IDUF-MSM as a parameter (to be used in Calib_updateParams)
        Params.mix_IDUF_MSM_matrix = behav_sexMix_IDUF_MSM;
        Params.mix_IDUF_MSM_value = pctSexMixing_TransGroupSex(Params.sexMix_IDUF_MSM);
        
% NEEDLE MIXING: TRANSMISSION GROUP AND SEX            


    % Import from Excel
        pctNeedleMixing_TransGroupSex = ExcelValues_AllParameters(parameter_Index:parameter_Index + numTransGroupSex_NeedleCombos-1);
        parameter_Index = parameter_Index + numel(pctNeedleMixing_TransGroupSex);

     % Transmission group/sex combinations [273x273]
            % 1 - Pull appropriate value
            % 2 - Use indicator to apply to appropriate cells
        behav_needleMix_IDUM_IDUM = Params.idx_IDUM_IDUM * pctNeedleMixing_TransGroupSex(Params.needleMix_IDUM_IDUM);
        behav_needleMix_IDUM_IDUF = Params.idx_IDUM_IDUF * pctNeedleMixing_TransGroupSex(Params.needleMix_IDUM_IDUF);
        behav_needleMix_IDUF_IDUM = Params.idx_IDUF_IDUM * pctNeedleMixing_TransGroupSex(Params.needleMix_IDUF_IDUM);  
        behav_needleMix_IDUF_IDUF = Params.idx_IDUF_IDUF * pctNeedleMixing_TransGroupSex(Params.needleMix_IDUF_IDUF);
    
    % Create full [273x273] matrix with each needle combo
        behav_needleMix_popSex = ...
              behav_needleMix_IDUM_IDUM + behav_needleMix_IDUM_IDUF ...
            + behav_needleMix_IDUF_IDUM + behav_needleMix_IDUF_IDUF;
        
   
% SEXUAL MIXING: OTHER FACTORS


    % RISK LEVEL
    
    
        % HET subpopulations

            % Pre-allocate
                sexMix_riskLevel_HET(numRiskLevel,numRiskLevel) = 0;
                sexMix_riskLevel_MSM(numRiskLevel,numRiskLevel) = 0;
                sexMix_riskLevel_IDU(numRiskLevel,numRiskLevel) = 0;
                
            % Import from Excel [2x2]
                
                % HET
                    for l = 1:numRiskLevel
                        sexMix_riskLevel_HET(l,:)=ExcelValues_AllParameters(parameter_Index+numRiskLevel*...
                            (l-1):parameter_Index+numRiskLevel-1+numRiskLevel*(l-1));   
                    end

                    parameter_Index = parameter_Index + numel(sexMix_riskLevel_HET);
                % MSM    
                     for l = 1:numRiskLevel
                        sexMix_riskLevel_MSM(l,:)=ExcelValues_AllParameters(parameter_Index+numRiskLevel*...
                            (l-1):parameter_Index+numRiskLevel-1+numRiskLevel*(l-1));   
                    end

                    parameter_Index = parameter_Index + numel(sexMix_riskLevel_MSM);
                % IDU
                    % Assumed equal to HET <-- CHECK THIS
                    sexMix_riskLevel_IDU = sexMix_riskLevel_HET;
                    
        % Define indicators
            % Combine risk level and transmission group
            % [273x2]: dim 1 is low risk, dim 2 is high risk
            % If the subpopulation isn't specified, both dim are 0
            riskLevelIndicator_HET(Params.numStrats,numRiskLevel) = 0;
            riskLevelIndicator_MSM(Params.numStrats,numRiskLevel) = 0;
            riskLevelIndicator_IDU(Params.numStrats,numRiskLevel) = 0;
            for nRL = 1:numRiskLevel
               riskLevelIndicator_HET(:,nRL) = Params.riskLevelIndicator(:,nRL) .* Params.popIndicator(:,Params.pop_HET);
            end

            for nRL = 1:numRiskLevel
               riskLevelIndicator_MSM(:,nRL) = Params.riskLevelIndicator(:,nRL) .* Params.popIndicator(:,Params.pop_MSM);
            end

            for nRL = 1:numRiskLevel
               riskLevelIndicator_IDU(:,nRL) = Params.riskLevelIndicator(:,nRL) .* Params.popIndicator(:,Params.pop_IDU);
            end
            
        % Apply to full [273x273] matrix        
        behav_sexMix_riskLevel = ...
                riskLevelIndicator_HET * sexMix_riskLevel_HET * Params.riskLevelIndicator' ...
              + riskLevelIndicator_MSM * sexMix_riskLevel_MSM * Params.riskLevelIndicator' ...
              + riskLevelIndicator_IDU * sexMix_riskLevel_IDU * Params.riskLevelIndicator';

        % Adjust when IDU is partner #2
            % Whenever a population has IDU as a partner, the risk level
            % mix is 1.
            
            % Record index when IDU is partner 2
            idx_IDUisPartner2 = Params.popIndicator(:,Params.pop_IDU) == 1;
        
            % Adjust in full matrix
            behav_sexMix_riskLevel(:,idx_IDUisPartner2) = 1;

       % Apply to Params struct
            Params.behav_mix_riskLevel_S = behav_sexMix_riskLevel;
            
            
    % RACE/ETHNICITY (HET/PWID) (edited by JC on 11/07/2017)

        % Pre-allocate
            store_mix_race_S_HET(numRace,numRace)=0;
       
        % Import from Excel [3x3]
            for r = 1:numRace
                store_mix_race_S_HET(r,:)=ExcelValues_AllParameters(parameter_Index+numRace*...
                    (r-1):parameter_Index+numRace-1+numRace*(r-1));   
            end

            parameter_Index = parameter_Index + numel(store_mix_race_S_HET);

        % Apply to full [273x273] matrix 
            % Dim 1 is partner 1, dim 2 is partner 2
        
            Params.behav_mix_race_S_HET = ...
                  Params.raceIndicator ...
                * store_mix_race_S_HET ...
                * Params.raceIndicator';
        
        % zero out MSM strats
            HET_indicator = repmat(Params.popIndicator(:,Params.pop_HET),1,Params.numStrats);
            IDU_indicator = repmat(Params.popIndicator(:,Params.pop_IDU),1,Params.numStrats);
            Params.HET_IDU_indicator = HET_indicator + IDU_indicator;
            
            Params.behav_mix_race_S_HET = Params.behav_mix_race_S_HET .* Params.HET_IDU_indicator;
            
    % RACE/ETHNICITY (MSM) (added by JC on 11/07/2017)

        % Pre-allocate
            store_mix_race_S_MSM(numRace,numRace)=0;
       
        % Import from Excel [3x3]
            for r = 1:numRace
                store_mix_race_S_MSM(r,:)=ExcelValues_AllParameters(parameter_Index+numRace*...
                    (r-1):parameter_Index+numRace-1+numRace*(r-1));   
            end

            parameter_Index = parameter_Index + numel(store_mix_race_S_MSM);

        % Apply to full [273x273] matrix 
            % Dim 1 is partner 1, dim 2 is partner 2
        
            Params.behav_mix_race_S_MSM = ...
                  Params.raceIndicator ...
                * store_mix_race_S_MSM ...
                * Params.raceIndicator';
            
        % Zero out HET / IDU Strats
            Params.MSM_indicator = repmat(Params.popIndicator(:,Params.pop_MSM),1,Params.numStrats);
            
            Params.behav_mix_race_S_MSM = Params.behav_mix_race_S_MSM .* Params.MSM_indicator;
        
        Params.behav_mix_race_S = Params.behav_mix_race_S_HET + Params.behav_mix_race_S_MSM;
        
    % AGE (HETS)

        % Pre-allocate
        store_mix_age_S_HET(Params.numAge,Params.numAge)=0;

        % Bring in from Excel [5x5]
            for a = 1:Params.numAge
                store_mix_age_S_HET(a,:)=ExcelValues_AllParameters(parameter_Index+Params.numAge*...
                    (a-1):parameter_Index+Params.numAge-1+Params.numAge*(a-1));   
            end

            parameter_Index = parameter_Index + numel(store_mix_age_S_HET);

        % Apply to full [273x273] matrix
            Params.behav_mix_age_S_HET(Params.numStrats,Params.numStrats)=0;

            Params.behav_mix_age_S_HET = ...
                  Params.ageIndicator ...
                * store_mix_age_S_HET ...
                * Params.ageIndicator';
            
            % zero out MSM strats
            Params.behav_mix_age_S_HET = Params.behav_mix_age_S_HET .* Params.HET_IDU_indicator;
 
   % AGE (MSM)[Added by Laurel 4/24/17]
        % Pre-allocate
        store_mix_age_S_MSM(Params.numAge,Params.numAge)=0;

        % Bring in from Excel [5x5]
            for a = 1:Params.numAge
                store_mix_age_S_MSM(a,:)=ExcelValues_AllParameters(parameter_Index+Params.numAge*...
                    (a-1):parameter_Index+Params.numAge-1+Params.numAge*(a-1));   
            end

            parameter_Index = parameter_Index + numel(store_mix_age_S_MSM);

        % Apply to full [273x273] matrix
            Params.behav_mix_age_S_MSM(Params.numStrats,Params.numStrats)=0;

            Params.behav_mix_age_S_MSM = ...
                  Params.ageIndicator ...
                * store_mix_age_S_MSM ...
                * Params.ageIndicator';  
            
            % Zero out HET / IDU Strats
            Params.behav_mix_age_S_MSM = Params.behav_mix_age_S_MSM .* Params.MSM_indicator;
            
            % Combined age mixing
            Params.behav_mix_age_S = Params.behav_mix_age_S_HET + Params.behav_mix_age_S_MSM;

        
% CALCULATE FULL SEXUAL MIXING MATRIX

    % Multiply factors
        % Trans group/sex
        % Risk level
        % Age
        % Race

        behav_store_sexMixing = ...
               Params.behav_mix_riskLevel_S .* behav_sexMix_popSex ...
            .* Params.behav_mix_race_S .* Params.behav_mix_age_S;

    % Adjust for circumcision
        % Applied proportionally
        
        % Pre-allocate [273x3]
            % By race b/c that's how parametr is defined
        malesUncirc = zeros(Params.numStrats, numRace);
        malesCirc = zeros(Params.numStrats, numRace);

        % Populate men who aren't circumcised [273x3]
            for r = 1:numRace
               malesUncirc(:,r) = ...
                       Params.raceIndicator(:,r) .* Params.circIndicator(:,Params.circ_U) ...
                    .* Params.sexIndicator(:,Params.sex_Male) * (1-Params.pop_pctCirc_r(r));
            end

        %Populate men who are circumcised [273x3]
        for r = 1:numRace
           malesCirc(:,r) = ...
                    Params.raceIndicator(:,r) .* Params.circIndicator(:,Params.circ_C)...
                 .* Params.sexIndicator(:,Params.sex_Male)*(Params.pop_pctCirc_r(r));
        end

        % Sum across all races
            % Note: zeroes in this parameter must be populated with ones
            Params.maleUncircMultiplier = sum(malesUncirc, 2); 
            Params.maleCircMultiplier = sum(malesCirc, 2);


        % Replace zeros with 1's
            idx = Params.maleUncircMultiplier == 0;
            Params.maleUncircMultiplier(idx) = 1;

            idx2 = Params.maleCircMultiplier == 0;
            Params.maleCircMultiplier(idx2) = 1;

        % Apply Cicumcision factor to mixing matrix
    
            % Pre-allocate
            Params.behav_SexualMixing(Params.numStrats,Params.numStrats)=0;
        
            % Loop through all the rows of the mixing matrix
            for a = 1:Params.numStrats
                Params.behav_SexualMixing(a,:) = ...
                       Params.maleUncircMultiplier' .* Params.maleCircMultiplier' ...
                    .* behav_store_sexMixing(a,:);
            end

    % Check: sum(Params.behav_SexualMixing,2) = 1 for all rows

    
    
    % Use sexual mixing matrix to define indices for all eligible mixing
    % combinatoins
    
        % Figure out where in the mixing matrix people don't mix with each other
        Params.idx_doesnotmix = Params.behav_SexualMixing == 0;

        % Determine all of the non-MSM sexual mixing combos

        % Initially - everyone can mix together
        Params.allMixCombos = ones(Params.numStrats, Params.numStrats);

        % Remove all the people who don't mix with each other
        Params.allMixCombos(Params.idx_doesnotmix) = 0;

        %Set MSM-MSM mixing to 0 (this will apply a 0% vaginal to
        %MSM-MSM interactions)
        allNonMSM_MSM_MixCombos = Params.allMixCombos;
        allNonMSM_MSM_MixCombos(Params.idx_MSMandMSM)=0;
        Params.idx_allMaleAndFemaleMixing = (allNonMSM_MSM_MixCombos == 1);

   
    
% NEEDLE MIXING (OTHER FACTORS)

    % Needle mixing by race/eth

        store_mix_race_N(numRace,numRace)=0;
        %Where dimension 1 is race of partner 1 and dimension 2 is race of partner
        %2
        for r = 1:numRace
            store_mix_race_N(r,:)=ExcelValues_AllParameters(parameter_Index+numRace*...
                (r-1):parameter_Index+numRace-1+numRace*(r-1));   
        end
        parameter_Index = parameter_Index + numel(store_mix_race_N);

        %Restratify store_mix_race_N from 3x3 to 273x273
        %Params.D1 is partner 1, Params.D2 is partner 2
        Params.behav_mix_race_N(Params.numStrats,Params.numStrats)=0;
        Params.behav_mix_race_N = Params.raceIndicator*store_mix_race_N*...
            Params.raceIndicator';

    % Needle mixing by age

        store_mix_age_N(Params.numAge,Params.numAge)=0;
        %Where dimension 1 is age of partner 1 and dim 2 is age of partner 2
        for a = 1:Params.numAge
            store_mix_age_N(a,:)=ExcelValues_AllParameters(parameter_Index+Params.numAge*...
                (a-1):parameter_Index+Params.numAge-1+Params.numAge*(a-1));   
        end

        parameter_Index = parameter_Index + numel(store_mix_age_N);

        %Restratify store_mix_age_N
        %from 5x5 to 273x273
        Params.behav_mix_age_N(Params.numStrats,Params.numStrats)=0;

        Params.behav_mix_age_N = Params.ageIndicator*store_mix_age_N*...
            Params.ageIndicator';    

        
% CALCULATE FULL NEEDLE MIXING MATRIX

    %Factors
        %needleMix_SexRiskType
        %Params.behav_mix_race_N
        %Params.behav_mix_race_N

    behav_mixing_store_N =  behav_needleMix_popSex .* ...
        Params.behav_mix_race_N .* Params.behav_mix_age_N;

    %Use multipliers to take weighted averages of uncirc and circ populations
    %Factors
        %maleUncircMultiplier
        %maleCircMultiplier

    Params.behav_NeedleMixing(Params.numStrats,Params.numStrats)=0;
    %Loop through all the rows of the mixing matrix

    for a = 1:Params.numStrats
         Params.behav_NeedleMixing(a,:)=Params.maleUncircMultiplier'.*Params.maleCircMultiplier' ...
             .* behav_mixing_store_N(a,:);
    end

    % Check: sum(Params.behav_NeedleMixing,2) = 1 for all rows with IDUs

%% 10. Define anal intercourse behavior inputs (Behaviors sheet)

    % In male-male partnerships, the only act is AI
    % In male-female partnerships, the default act is VI, some partnerships include AI.
    
    
% PERCENTAGE OF PEOPLE WHO HAVE AI IN THEIR MALE-FEMALE PARTNERSHIPS     
    % [273x1]    
        
    % Populate for each time period
    for RatePeriod = 1:Params.numRateInputPeriod
        
        % Preallocate
        AIPrev_Overall_273 = zeros(Params.numStrats,1);
        
        % Loop through the races to pull the AI prev for each one by one
        for r = 1:Params.numRace
            
        % Pull the age- and sex- specific data for the race specified by the loop
            % [10x1] where 5 age categories and 2 sexes
        store_raceSpecificAIAge(1:Params.numAge*2,1) = ExcelValues_AllParameters(parameter_Index:parameter_Index +...
            numSex * Params.numAge - 1);
        parameter_Index = parameter_Index + numel(store_raceSpecificAIAge);
        
        % Apply the age- and sex- specific AI rates to the full 273
        % parameter. Zero out cohorts which are not applicable to the
        % current race.
            % Moving from [10x1] -> [273x1]
        AIPrev_raceSpecific_273 = (Params.ageIndicator * ...
            store_raceSpecificAIAge(1:Params.numAge)...
            .* Params.sexIndicator(:,Params.sex_Male) + ...
            Params.ageIndicator * ...
            store_raceSpecificAIAge(Params.numAge+1:Params.numAge*2)...
             .* Params.sexIndicator(:,Params.sex_Female)) .* ...
            Params.raceIndicator(:,r);
        
        % Apply the race specific AI parameter to the all race parameter.
            % Because cohorts for the races are mutually exclusive, this 
            % should never add
            % values on top of each other - it will merely plug the race
            % specific values into the race-specific cohorts of the
            % overall273 parameter
        AIPrev_Overall_273 = AIPrev_Overall_273 + AIPrev_raceSpecific_273;
        
        end
        
        % Apply the calculated array to the period-specific parameter
        switch RatePeriod
            case 1
                Params.pctPeopleA_InMF_1 = AIPrev_Overall_273;
            case 2
                Params.pctPeopleA_InMF_2to4 = AIPrev_Overall_273;
            case 3
                Params.pctPeopleA_InMF_5 = AIPrev_Overall_273;
        end
        
        clear AIPrev_Overall273
        
    end  

    
% AMONG PEOPLE WHO HAVE AI IN THEIR MALE-FEMALE PARTNERSHIPS,
% PERCENTAGE OF PARTNERSHIPS THAT INCLUDE AI (modified 07/31/2018)
    

    % Import
    store_pctPartnershipsAIfSomeA_1(:,1) = ExcelValues_AllParameters(parameter_Index : parameter_Index+ Params.numAge- 1);
    parameter_Index = parameter_Index + numel(store_pctPartnershipsAIfSomeA_1);
    store_pctPartnershipsAIfSomeA_2to4(:,1) = ExcelValues_AllParameters(parameter_Index : parameter_Index+ Params.numAge- 1);
    parameter_Index = parameter_Index + numel(store_pctPartnershipsAIfSomeA_2to4);
    store_pctPartnershipsAIfSomeA_5(:,1) = ExcelValues_AllParameters(parameter_Index : parameter_Index+ Params.numAge- 1);
    parameter_Index = parameter_Index + numel(store_pctPartnershipsAIfSomeA_5);
    
    % Apply to 273x1
    pctMFPartnersA_IfSomeAInMF_1 = Params.ageIndicator * store_pctPartnershipsAIfSomeA_1;
    pctMFPartnersA_IfSomeAInMF_2to4 = Params.ageIndicator * store_pctPartnershipsAIfSomeA_2to4;
    pctMFPartnersA_IfSomeAInMF_5 = Params.ageIndicator * store_pctPartnershipsAIfSomeA_5;
    
    pctMFPartnersA_IfSomeAInMF_1_2D(Params.numStrats,Params.numStrats) = 0;
    pctMFPartnersA_IfSomeAInMF_2to4_2D(Params.numStrats,Params.numStrats) = 0;
    pctMFPartnersA_IfSomeAInMF_5_2D(Params.numStrats,Params.numStrats) = 0;
    
    % Loop over partner cohorts
    for s = 1:Params.numStrats     
       pctMFPartnersA_IfSomeAInMF_1_2D(:,s) = allNonMSM_MSM_MixCombos(:,s) .* ...
           pctMFPartnersA_IfSomeAInMF_1;
       pctMFPartnersA_IfSomeAInMF_2to4_2D(:,s) = allNonMSM_MSM_MixCombos(:,s) .* ...
           pctMFPartnersA_IfSomeAInMF_2to4;
       pctMFPartnersA_IfSomeAInMF_5_2D(:,s) = allNonMSM_MSM_MixCombos(:,s) .* ...
           pctMFPartnersA_IfSomeAInMF_5;       
    end
    
    % Ceate 3D matrix
    pctMFPartnersA_IfSomeAInMF_1_3D = repmat(pctMFPartnersA_IfSomeAInMF_1_2D,[1,1,Params.numComparts]);
    pctMFPartnersA_IfSomeAInMF_2to4_3D = repmat(pctMFPartnersA_IfSomeAInMF_2to4_2D,[1,1,Params.numComparts]);
    pctMFPartnersA_IfSomeAInMF_5_3D = repmat(pctMFPartnersA_IfSomeAInMF_5_2D,[1,1,Params.numComparts]);
    % Re-arrange matrix for Params struct to be used in CalcInfectionRates (30x273x273)
    Params.pctMFPartnersA_IfSomeAInMF_1 = permute(pctMFPartnersA_IfSomeAInMF_1_3D,[3,1,2]);
    Params.pctMFPartnersA_IfSomeAInMF_2to4 = permute(pctMFPartnersA_IfSomeAInMF_2to4_3D,[3,1,2]);
    Params.pctMFPartnersA_IfSomeAInMF_5 = permute(pctMFPartnersA_IfSomeAInMF_5_3D,[3,1,2]);


% IN M-F PARTNERSHIPS WITH AI, PERCENTAGE OF CONTACTS VI (OR AI)
    % Parameter 1: Percent of sexual contacts which are vaginal given the
    %              partnership includes both vaginal and anal sex
    
    % Parameter 2: Percent of sexual contacts which are anal given the
    %              partnership includes both vaginal and anal sex
    
    % Note: these are defined by separate parameters rather than one being
    % the complement because there are mixing combinations that should
    % always be 0. E.g., if you took [pctAnal] = (1 - [pctVaginal]), you would 
    % see a 100% in the invalid combination HETM-IDUM.
    
    
    % Loop through each period
    for RatePeriod = 1:Params.numRateInputPeriod
        

        % Import sex and age specific data
            store_AgeMale = ExcelValues_AllParameters(parameter_Index:...
                parameter_Index + Params.numAge - 1);
            parameter_Index = parameter_Index + numel(store_AgeMale);

            store_AgeFemale = ExcelValues_AllParameters(parameter_Index:...
                parameter_Index + Params.numAge -1);
            parameter_Index = parameter_Index + numel(store_AgeFemale);

        % Percent contacts anal [273x1]
            pctContactsAnal = Params.ageIndicator * store_AgeMale .* ...
                Params.sexIndicator(:,Params.sex_Male) + ...
                Params.ageIndicator * store_AgeFemale .* ...
                Params.sexIndicator(:,Params.sex_Female);

        % Percent contacts vaginal [273x1]
            pctContactsVaginal = 1 - pctContactsAnal;


        % Apply to 273x273x30 variables

            % Pre-allocate
            pctVaginal_2D = zeros(Params.numStrats,Params.numStrats);
            pctAnal_2D = zeros(Params.numStrats,Params.numStrats);

            % Loop over partner cohorts
            for s = 1:Params.numStrats 

               pctVaginal_2D(:,s) = allNonMSM_MSM_MixCombos(:,s) .* ...
                   pctContactsVaginal;
               pctAnal_2D(:,s) = allNonMSM_MSM_MixCombos(:,s) .* ...
                   pctContactsAnal;

            end

            % Expand to a 3D matrix 
                % Apply the 2D mixing matrix to every Compartment
                % Note: So that it can be applied to the betas

            pctVaginal_3D = repmat(pctVaginal_2D,[1,1,Params.numComparts]);
            pctAnal_3D = repmat(pctAnal_2D,[1,1,Params.numComparts]);


        % Apply to Beta Factors - Updated by MClinkscales on 6/1/2022 to
        % have T2-T4 betas calculate during the second case
            % Note: used in CalcBetas.m

        switch RatePeriod
            case 1
                BetaFactor.pctContactsVaginal_InMFWithA_1 = pctVaginal_3D;
                BetaFactor.pctContactsAnal_InMFWithA_1 = pctAnal_3D;
                
                % Save values of period 1 beta factors for use in
                % Calib_updateParams
                Params.pctContactsAnal_1 = pctContactsAnal;
                Params.pctContactsVaginal_1 = pctContactsVaginal;
                %Params.pctContactsVaginal_InMFWithA_1 = BetaFactor.pctContactsVaginal_InMFWithA_1;
                %Params.pctContactsAnal_InMFWithA_1 = BetaFactor.pctContactsAnal_InMFWithA_1;
            case 2
                BetaFactor.pctContactsVaginal_InMFWithA_2to4 = pctVaginal_3D;
                BetaFactor.pctContactsAnal_InMFWithA_2to4 = pctAnal_3D;                
                % Save values of Period 2 to 4 beta factors for use in
                % Calib_updateParams
                Params.pctContactsAnal_2to4 = pctContactsAnal;
                Params.pctContactsVaginal_2to4 = pctContactsVaginal;
                                                             
            case 3
                BetaFactor.pctContactsVaginal_InMFWithA_5 = pctVaginal_3D;             
                BetaFactor.pctContactsAnal_InMFWithA_5 = pctAnal_3D;
        end
        
        
    end
        
%% 11. Define HIV progression and cost inputs

% COMPARTMENT COSTS
    % Comparised of HIV cost and PrEP costs
    % Goal is to apply the correct costs to each compartment [30x1]

    % HIV costs: applied to all PLWH
        % Stratifed by disease stage and continuum of care stages
        % Varies by period
        
    % Overview of parameter
        % Each distinct continuum of care stage (Unaware/Aware/LTC/ART) has
        % a cost which varies by severity of disease
        % These costs directly correspond to the model diagrams rows 1-3
        % and row 5, respectively.
        
        % However, row 4 of the model diagram (ANV) is comprised of
        % multiple continuum of care stages (Aware/LTC/ART), and the
        % distribution changes each time period. The applied cost is a
        % weighted average.
            
        % The goal of all the code below is to apply the appropriate costs
        % to the correct compartment.
            % First, we bring in the distinct costs (Unaware/Aware/LTC/ART)
            % Then, we calculate the costs for the ANV row
            % Finally, we apply them to the HIV cost parameters.
        
       % Acute is handled separately from Chronic because it doesn't have
       % ANV or VLS stages.
    
       % Note that "changes over time" means that the base value changes
       % each time period. Within a time period, all of the costs are fixed.
            % Discounting is handled in CollectResults.m.
       
    % Bring in Acute values
        % These don't change over time (because no ANV stage)
        costsByDistinctContStageAcute = ...
            ExcelValues_AllParameters(parameter_Index:parameter_Index+Params.B3-Params.B1); 
        parameter_Index = parameter_Index + numel(costsByDistinctContStageAcute);
        
    % Bring in Chronic Initial Matrix from Excel
        % This does not change over time
    
        % Variable is set up as [distinct cont stages with costs x chronic disease stages]
        % Eg.
                %               LatentA     LatentB     Late    AIDS
                %   Unaware      $X          $X          $X      $X
                %   Aware        $X          $X          $X      $X
                %   LTC          $X          $X          $X      $X                 
                %   ART          $X          $X          $X      $X
    
        store_ChronicTreatCostsByDistinctContStage = zeros(numDistinctContStagesWithCosts,numChronicDiseaseStages);
        
        store_ChronicTreatCostsByDistinctContStage(:) = ...
            ExcelValues_AllParameters(parameter_Index:parameter_Index+numDistinctContStagesWithCosts*numChronicDiseaseStages-1); 
       
        parameter_Index = parameter_Index + numel(store_ChronicTreatCostsByDistinctContStage);


    % Apply to compartment-compatible matrix
        % Does change over time based on the distribution of the population
        % in the ANV compartments
        
        % Variable is set up as [(model's cont stages if chronic HIV) x (chronic disease stages)]
        % Eg. 
                %                          LatentA     LatentB     Late    AIDS
                %   Unaware                 $X          $X          $X      $X
                %   Aware, no ART effects   $X          $X          $X      $X
                %   LTC, no ART effects     $X          $X          $X      $X                 
                %   ANV                     $X          $X          $X      $X
                %   VLS                     $X          $X          $X      $X
                
         % Pre-allocate
         costsByModelsContStageChronic_1 = zeros(Params.numModelsContStagesIfChronic,numChronicDiseaseStages);
         costsByModelsContStageChronic_2to4 = zeros(Params.numModelsContStagesIfChronic,numChronicDiseaseStages);
         costsByModelsContStageChronic_5 = zeros(Params.numModelsContStagesIfChronic,numChronicDiseaseStages);
         
         % First 3 rows (Unaware/Aware/LTC) don't change
            % 1st row: based on Unaware costs
            % 2nd row: based on Aware costs
            % 3rd row: based on Care costs
         costsByModelsContStageChronic_1(Params.chronicContStage_Unaware:Params.chronicContStage_LTCNoART,:) = store_ChronicTreatCostsByDistinctContStage(Params.distinctContStage_Unaware:Params.distinctContStage_LTC,:);
         costsByModelsContStageChronic_2to4(Params.chronicContStage_Unaware:Params.chronicContStage_LTCNoART,:) = store_ChronicTreatCostsByDistinctContStage(Params.distinctContStage_Unaware:Params.distinctContStage_LTC,:);
         costsByModelsContStageChronic_5(Params.chronicContStage_Unaware:Params.chronicContStage_LTCNoART,:) = store_ChronicTreatCostsByDistinctContStage(Params.distinctContStage_Unaware:Params.distinctContStage_LTC,:);
         
         % Last row (VLS) doesn't change
            % Based on ART costs
         costsByModelsContStageChronic_1(Params.chronicContStage_VLS,:) = store_ChronicTreatCostsByDistinctContStage(Params.distinctContStage_ART,:);
         costsByModelsContStageChronic_2to4(Params.chronicContStage_VLS,:) = store_ChronicTreatCostsByDistinctContStage(Params.distinctContStage_ART,:);
         costsByModelsContStageChronic_5(Params.chronicContStage_VLS,:) = store_ChronicTreatCostsByDistinctContStage(Params.distinctContStage_ART,:);
         
         
         % 4th row does change over time (because distribution of PLWH in ANV changes over time)
            % Based on Aware, Care, and ART costs proportionally to the
            % distribution of people in each of those stages in ANV
            
            % Calc (weighted average):
                % [ANV Costs] =
                %     [Aware costs] * (1-[Pct ANV in care])
                %   + [Care costs] * [Pct ANV in care] * (1 - [Pct in-care ANV on ART])
                %   + [ART costs] * [Pct ANV in care] * [Pct in-care ANV on ART]
                
         
         costsByModelsContStageChronic_1(Params.chronicContStage_ANV,:) = ...
               store_ChronicTreatCostsByDistinctContStage(Params.distinctContStage_Aware,:) * (1-Params.tt_PctANVWhoAreInCareOrART_1) ...
             + store_ChronicTreatCostsByDistinctContStage(Params.distinctContStage_LTC,:) * (Params.tt_PctANVWhoAreInCareOrART_1) * (1-Params.tt_PctInCareANVWhoAreOnART_1) ...
             + store_ChronicTreatCostsByDistinctContStage(Params.distinctContStage_ART,:) * (Params.tt_PctANVWhoAreInCareOrART_1) * (Params.tt_PctInCareANVWhoAreOnART_1);
         
         costsByModelsContStageChronic_2to4(Params.chronicContStage_ANV,:) = ...
               store_ChronicTreatCostsByDistinctContStage(Params.distinctContStage_Aware,:) * (1-Params.tt_PctANVWhoAreInCareOrART_2to4) ...
             + store_ChronicTreatCostsByDistinctContStage(Params.distinctContStage_LTC,:) * (Params.tt_PctANVWhoAreInCareOrART_2to4) * (1-Params.tt_PctInCareANVWhoAreOnART_2to4) ...
             + store_ChronicTreatCostsByDistinctContStage(Params.distinctContStage_ART,:) * (Params.tt_PctANVWhoAreInCareOrART_2to4) * (Params.tt_PctInCareANVWhoAreOnART_2to4);
         
         costsByModelsContStageChronic_5(Params.chronicContStage_ANV,:) = ...
               store_ChronicTreatCostsByDistinctContStage(Params.distinctContStage_Aware,:) * (1-Params.tt_PctANVWhoAreInCareOrART_5) ...
             + store_ChronicTreatCostsByDistinctContStage(Params.distinctContStage_LTC,:) * (Params.tt_PctANVWhoAreInCareOrART_5) * (1-Params.tt_PctInCareANVWhoAreOnART_5) ...
             + store_ChronicTreatCostsByDistinctContStage(Params.distinctContStage_ART,:) * (Params.tt_PctANVWhoAreInCareOrART_5) * (Params.tt_PctInCareANVWhoAreOnART_5);
    % PrEP costs %JCPrEPUpdate: updated to read in PrEP cost for each state
        costPrEPComparts = ExcelValues_AllParameters(parameter_Index:parameter_Index+Params.numPrEPStates-1); 
        parameter_Index = parameter_Index + numel(costPrEPComparts);
         
        % Pre-allocate
        Params.hiv_annualCostPerCompart_1 = zeros(Params.numComparts,1);
        Params.hiv_annualCostPerCompart_2to4 = zeros(Params.numComparts,1);
        Params.hiv_annualCostPerCompart_5 = zeros(Params.numComparts,1);
         
        % Apply values to parameters    
        % Acute
            Params.hiv_annualCostPerCompart_1(Params.B1:Params.B3) = costsByDistinctContStageAcute;
            Params.hiv_annualCostPerCompart_2to4(Params.B1:Params.B3) = costsByDistinctContStageAcute;
            Params.hiv_annualCostPerCompart_5(Params.B1:Params.B3) = costsByDistinctContStageAcute;
                
        % Chronic
            Params.hiv_annualCostPerCompart_1(Params.C1:Params.F5) = costsByModelsContStageChronic_1;
            Params.hiv_annualCostPerCompart_2to4(Params.C1:Params.F5) = costsByModelsContStageChronic_2to4; 
            Params.hiv_annualCostPerCompart_5(Params.C1:Params.F5) = costsByModelsContStageChronic_5; 
       
        % PrEP
            % Note: doesn't change over time %JCPrEPUpdate: Updated to
            % include all PrEP comparts
            Params.hiv_annualCostPerCompart_1(Params.PrEPComparts) = costPrEPComparts;
            Params.hiv_annualCostPerCompart_2to4(Params.PrEPComparts) = costPrEPComparts;
            Params.hiv_annualCostPerCompart_5(Params.PrEPComparts) = costPrEPComparts; 
            
% HIV CONTINUUM OF CARE TRANSITION COSTS (added by JC on 02/07/2018)            
    
    % Read in from Excel
    Params.transcostPP_addtlPrEPTestCostIfInf = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
    Params.transcostPP_LTCFirst = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
    Params.transcostPP_LTCAfter = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
    Params.transcostPP_ARTInitiation = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
    Params.transcostPP_BecomeVLSfromANV = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
    Params.anntranscostPP_RemainVLS = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;

% NUMBER OF YEARS IN EACH AGE GROUP
    % 13-17, 18-24, 25-34, 35-44, 45-64, 65+
    
    % Read in unadjusted duration of age in age group 
    % (needed for calibration to multiply with calibrated multiplier)
    Params.durAge_normal = ExcelValues_AllParameters(parameter_Index:parameter_Index+Params.numAge-1);
    parameter_Index = parameter_Index + numel(Params.durAge_normal);    
    
    % Read in from Excel [7x1]
    hiv_durAge = ExcelValues_AllParameters(parameter_Index:parameter_Index+Params.numAge-1);
    parameter_Index = parameter_Index + numel(hiv_durAge);

    % Apply to subpopulation matrix [273x1]
    Params.hiv_durAge_a(Params.numStrats) = 0;
    Params.hiv_durAge_a = Params.ageIndicator * hiv_durAge;

% NUMBER OF YEARS IN EACH DISEASE STAGE (modified by JC on 08/02/2018)

    % Natural history [5x273] 
        % Applied to all PLWH who are not [ART-not-VLS] or [VLS]
        % Acute, CD4>500, CD4 350-500, CD4 200-350, CD4<200
        
        % Pre-allocate
        Params.hiv_durHIVStage_NaturalHistory(numHIVstages,Params.numStrats) = 0;        
        
        hiv_durHIVStage_NaturalHistory = ExcelValues_AllParameters(parameter_Index: ...
                parameter_Index + numHIVstages - 1);
        parameter_Index = parameter_Index + numHIVstages;
        
        hiv_durHIVStage_Adjust_Age(1,:) = ExcelValues_AllParameters(parameter_Index: ...
                parameter_Index + Params.numAge - 1);
        parameter_Index = parameter_Index + Params.numAge;
        
        hiv_durHIVStage_Adjust_Age = hiv_durHIVStage_Adjust_Age * Params.ageIndicator';
        
        for s=1:numHIVstages          
            Params.hiv_durHIVStage_NaturalHistory(s,:) = hiv_durHIVStage_NaturalHistory(s) * hiv_durHIVStage_Adjust_Age;            
        end
                        
        
    % ART-not-VLS
        % Applied to PLWH in [ART-not-VLS] compartments
        
        % Note: Acute individuals not eligible for ART
        % Note: AIDS/ANV receive ART mortality (from NA-ACCORD)
        % Note: Excel inputs break this into three age groups (13-44,
        % 45-64,65+)
        
        numAgeGroups = 3;
        
        hiv_durHIVStage_ANV_LatentA_Age(1,:) = ExcelValues_AllParameters(parameter_Index: ...
                parameter_Index + numAgeGroups - 1);
        parameter_Index = parameter_Index + numAgeGroups;
        
        hiv_durHIVStage_ANV_LatentA(1,1:4) = hiv_durHIVStage_ANV_LatentA_Age(1,1);
        hiv_durHIVStage_ANV_LatentA(1,5:Params.numAge-1) = hiv_durHIVStage_ANV_LatentA_Age(1,2);
        hiv_durHIVStage_ANV_LatentA(1,Params.numAge) = hiv_durHIVStage_ANV_LatentA_Age(1,3);
        
        hiv_durHIVStage_ANV_LatentB_Age(1,:) = ExcelValues_AllParameters(parameter_Index: ...
                parameter_Index + numAgeGroups - 1);
        parameter_Index = parameter_Index + numAgeGroups;
        
        hiv_durHIVStage_ANV_LatentB(1,1:4) = hiv_durHIVStage_ANV_LatentB_Age(1,1);
        hiv_durHIVStage_ANV_LatentB(1,5:Params.numAge-1) = hiv_durHIVStage_ANV_LatentB_Age(1,2);
        hiv_durHIVStage_ANV_LatentB(1,Params.numAge) = hiv_durHIVStage_ANV_LatentB_Age(1,3);
        
        hiv_durHIVStage_ANV_Late_Age(1,:) = ExcelValues_AllParameters(parameter_Index: ...
                parameter_Index + numAgeGroups - 1);
        parameter_Index = parameter_Index + numAgeGroups;
        
        hiv_durHIVStage_ANV_Late(1,1:4) = hiv_durHIVStage_ANV_Late_Age(1,1);
        hiv_durHIVStage_ANV_Late(1,5:Params.numAge-1) = hiv_durHIVStage_ANV_Late_Age(1,2);
        hiv_durHIVStage_ANV_Late(1,Params.numAge) = hiv_durHIVStage_ANV_Late_Age(1,3);
        
        Params.hiv_durHIVStage_ANV_LatentA = hiv_durHIVStage_ANV_LatentA * Params.ageIndicator';
        Params.hiv_durHIVStage_ANV_LatentB = hiv_durHIVStage_ANV_LatentB * Params.ageIndicator';
        Params.hiv_durHIVStage_ANV_Late = hiv_durHIVStage_ANV_Late * Params.ageIndicator';
        
% DISEASE STAGE TRANSITION RATES WHILE VLS (modified on 08/01/2018 - JC to modify again)

    % Decrease in CD4 count 
        % Rate of transitioning to a lower HIV stage while VLS       
        
        % Read in from Excel
        Params.hiv_rateVLSCD4decr_LatentA = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1;
        
        Params.hiv_rateVLSCD4decr_LatentB = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1;
        
        Params.hiv_rateVLSCD4decr_Late = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1;
               
    % Increase in CD4 count
        % Rate of transitioning to a higher HIV stage while VLS
        % Note: Excel inputs break this into three age groups (13-44,
        % 45-64,65+)
        
        % Read in from Excel            
        hiv_rateVLSCD4incr_LatentB_Age(1,:) = ExcelValues_AllParameters(parameter_Index: ...
                parameter_Index + numAgeGroups - 1);
        parameter_Index = parameter_Index + numAgeGroups;
        
        hiv_rateVLSCD4incr_LatentB(1,1:4) = hiv_rateVLSCD4incr_LatentB_Age(1,1);
        hiv_rateVLSCD4incr_LatentB(1,5:Params.numAge-1) = hiv_rateVLSCD4incr_LatentB_Age(1,2);
        hiv_rateVLSCD4incr_LatentB(1,Params.numAge) = hiv_rateVLSCD4incr_LatentB_Age(1,3);
        
        hiv_rateVLSCD4incr_Late_Age(1,:) = ExcelValues_AllParameters(parameter_Index: ...
                parameter_Index + numAgeGroups - 1);
        parameter_Index = parameter_Index + numAgeGroups;
        
        hiv_rateVLSCD4incr_Late(1,1:4) = hiv_rateVLSCD4incr_Late_Age(1,1);
        hiv_rateVLSCD4incr_Late(1,5:Params.numAge-1) = hiv_rateVLSCD4incr_Late_Age(1,2);
        hiv_rateVLSCD4incr_Late(1,Params.numAge) = hiv_rateVLSCD4incr_Late_Age(1,3);
        
        hiv_rateVLSCD4incr_AIDS_Age(1,:) = ExcelValues_AllParameters(parameter_Index: ...
                parameter_Index + numAgeGroups - 1);
        parameter_Index = parameter_Index + numAgeGroups;
        
        hiv_rateVLSCD4incr_AIDS(1,1:4) = hiv_rateVLSCD4incr_AIDS_Age(1,1);
        hiv_rateVLSCD4incr_AIDS(1,5:Params.numAge-1) = hiv_rateVLSCD4incr_AIDS_Age(1,2);
        hiv_rateVLSCD4incr_AIDS(1,Params.numAge) = hiv_rateVLSCD4incr_AIDS_Age(1,3);
        
        Params.hiv_rateVLSCD4incr_LatentB = hiv_rateVLSCD4incr_LatentB * Params.ageIndicator';
        Params.hiv_rateVLSCD4incr_Late = hiv_rateVLSCD4incr_Late * Params.ageIndicator';
        Params.hiv_rateVLSCD4incr_AIDS = hiv_rateVLSCD4incr_AIDS * Params.ageIndicator';    
                    
% PROBABILITY OF DEATH        

    % Time-period specific values:
        % _1 refers to period 1
        % _2to5 refers to periods 2 to 5
        
    % PLWH Acute are not eligible for ART
        
  %Preallocate

     %prob_deathVLSNotAIDS_12(Params.numStrats)=1;
     %prob_deathVLSAIDS_12(Params.numStrats)=1;
    
     %prob_deathVLSNotAIDS_3(Params.numStrats)=1;
     %prob_deathVLSAIDS_3(Params.numStrats)=1;
    
  
    Params.rate_deathVLSNotAIDS_1(Params.numStrats)=1;
    Params.rate_deathVLSAIDS_1(Params.numStrats)=1;
    
    Params.rate_deathVLSNotAIDS_2to5(Params.numStrats)=1;
    Params.rate_deathVLSAIDS_2to5(Params.numStrats)=1;
    
    Params.rate_deathANVNotAIDS_1(Params.numStrats)=1;
    Params.rate_deathANVAIDS_1(Params.numStrats)=1;
    
    Params.rate_deathANVNotAIDS_2to5(Params.numStrats)=1;
    Params.rate_deathANVAIDS_2to5(Params.numStrats)=1;
    
    Params.rate_deathnoARTEffectsNotAIDS_1(Params.numStrats)=1;
    
    Params.rate_deathnoARTEffectsNotAIDS_2to5(Params.numStrats)=1;

    % Disease and period specific values
    
    for s = 1:4 
            
        % Preallocate    
        storeAllAge = zeros(Params.numStrats,1);
        
        
        

            % Pick off the 5 values (age-stratified)
            % and disease stage
            ageDiseaseSpecificRateDeath = ExcelValues_AllParameters(parameter_Index: ...
                parameter_Index + Params.numAge - 1);
            parameter_Index = parameter_Index+numel(ageDiseaseSpecificRateDeath);

            ageSpecific_273 = Params.ageIndicator * ageDiseaseSpecificRateDeath; 

            storeAllAge = storeAllAge + ageSpecific_273;
        
          

        switch s
            case 1              
                Params.rate_deathVLSNotAIDS_1 = storeAllAge;                
            case 2
                Params.rate_deathVLSAIDS_1 = storeAllAge;                
            case 3
                Params.rate_deathVLSNotAIDS_2to5 = storeAllAge;                
            case 4
                Params.rate_deathVLSAIDS_2to5 = storeAllAge;              
        end
        
        clear storeAllRace

    end

    %Relative risk of death vs. VLS
        %No ART effects = rows 1,2,3
        %ANV = row 4
    
    relRiskDeath_noARTEffectsvsVLS_1 = ExcelValues_AllParameters(parameter_Index);
            parameter_Index = parameter_Index + numel(relRiskDeath_noARTEffectsvsVLS_1);
            
    relRiskDeath_ANVvsVLS_1 = ExcelValues_AllParameters(parameter_Index);
            parameter_Index = parameter_Index + numel(relRiskDeath_ANVvsVLS_1);
    
    relRiskDeath_noARTEffectsvsVLS_2to5 = ExcelValues_AllParameters(parameter_Index);
            parameter_Index = parameter_Index + numel(relRiskDeath_noARTEffectsvsVLS_2to5);
            
    relRiskDeath_ANVvsVLS_2to5 = ExcelValues_AllParameters(parameter_Index);
            parameter_Index = parameter_Index + numel(relRiskDeath_ANVvsVLS_2to5);
            
    %Probability of Death if in rows 1,2,3 (no ART Effects) or row 4 (ANV)       
    
        %rows 1,2,3   %Note: Not calculated for AIDS in rows 1,2,3 because
        %progression is based on natural history
        Params.rate_deathnoARTEffectsNotAIDS_1 = Params.rate_deathVLSNotAIDS_1 * ...
            relRiskDeath_noARTEffectsvsVLS_1;
        
        Params.rate_deathnoARTEffectsNotAIDS_2to5 = Params.rate_deathVLSNotAIDS_2to5 * ...
            relRiskDeath_noARTEffectsvsVLS_2to5;
        
        %row 4
        Params.rate_deathANVNotAIDS_1 = Params.rate_deathVLSNotAIDS_1 * ...
            relRiskDeath_ANVvsVLS_1;
        
        Params.rate_deathANVAIDS_1 = Params.rate_deathVLSAIDS_1 * ...
            relRiskDeath_ANVvsVLS_1;
        
        Params.rate_deathANVNotAIDS_2to5 = Params.rate_deathVLSNotAIDS_2to5 * ...
            relRiskDeath_ANVvsVLS_2to5;
        
        Params.rate_deathANVAIDS_2to5 = Params.rate_deathVLSAIDS_2to5 * ...
            relRiskDeath_ANVvsVLS_2to5;
    

% UTILITIES
    % By disease status
    
    % Read in from Excel
        % Note: numHIVstages isn't subtracted by 1 because it also includes uninfected
        store_utilities = ExcelValues_AllParameters(parameter_Index:parameter_Index + numHIVstages);          
        parameter_Index = parameter_Index + numel(store_utilities);
        
    % Apply to parameters
        Params.hiv_utility_A = store_utilities(1);
        Params.hiv_utility_B = store_utilities(2);
        Params.hiv_utility_C = store_utilities(3);
        Params.hiv_utility_D = store_utilities(4);
        Params.hiv_utility_E = store_utilities(5);
        Params.hiv_utility_F = store_utilities(6);

       
%% 12. Define Infectivity inputs
 
% OVERALL SEXUAL TRANSMISSION PROBABILITIES (BY CONTACT TYPE)
    % Per unprotected contact between HIV+ and (uncircumcised, if male) HIV- partners

    % Read in from Excel
    
        % Vaginal Insertive
            inf_overallProb_VaginalIns = ExcelValues_AllParameters(parameter_Index);
            parameter_Index = parameter_Index + numel(inf_overallProb_VaginalIns);

        % Vaginal Receptive
            inf_overallProb_VaginalRec = ExcelValues_AllParameters(parameter_Index);
            parameter_Index = parameter_Index + numel(inf_overallProb_VaginalRec);

        % Anal Insertive
            inf_overallProb_AnalIns = ExcelValues_AllParameters(parameter_Index);
            parameter_Index = parameter_Index + numel(inf_overallProb_AnalIns);

        % Anal Receptive
            inf_overallProb_AnalRec = ExcelValues_AllParameters(parameter_Index);
            parameter_Index = parameter_Index + numel(inf_overallProb_AnalRec);
            
        
% RELATIVE RISK OF TRANSMISSION BY HIV STAGE
    % Note: assumed applied to the probabilities
        
    % Pre-allocate
    Params.store_relRiskSexTransmission(numHIVstages,1)=0;
    
    % Read in from Excel
        % Note: this is included as a param because it will also be used in
        % UpdateParams.m duringa calibration
    Params.store_relRiskSexTransmission = ExcelValues_AllParameters(parameter_Index:...
        parameter_Index+numHIVstages-1);
    parameter_Index = parameter_Index + numel(Params.store_relRiskSexTransmission);
    
    
% BASE TRANSMISSION RATE
    % Calculated as the overall transmission rate by contact type *
    % relative risk of transmission by disease stage

    % Vaginal Insertive
    
        % Pre-allocate
        inf_basePerActProb_VaginalIns = zeros(Params.numStrats, Params.numStrats, Params.numComparts);

        % Apply relative risk by disease stage
        inf_basePerActProb_VaginalIns(:,:,AcuteComparts)=inf_overallProb_VaginalIns * Params.store_relRiskSexTransmission(Params.stage_Acute);
        inf_basePerActProb_VaginalIns(:,:,LatentAComparts)=inf_overallProb_VaginalIns * Params.store_relRiskSexTransmission(Params.stage_LatentA);
        inf_basePerActProb_VaginalIns(:,:,LatentBComparts)=inf_overallProb_VaginalIns * Params.store_relRiskSexTransmission(Params.stage_LatentB);
        inf_basePerActProb_VaginalIns(:,:,LateComparts)=inf_overallProb_VaginalIns * Params.store_relRiskSexTransmission(Params.stage_Late);
        inf_basePerActProb_VaginalIns(:,:,AIDSComparts)=inf_overallProb_VaginalIns * Params.store_relRiskSexTransmission(Params.stage_AIDS);

    % Anal Insertive

        % Pre-allocate
        inf_basePerActProb_AnalIns = zeros(Params.numStrats, Params.numStrats, Params.numComparts);

        % Apply relative risk by disease stage
        inf_basePerActProb_AnalIns(:,:,AcuteComparts)=inf_overallProb_AnalIns * Params.store_relRiskSexTransmission(Params.stage_Acute);
        inf_basePerActProb_AnalIns(:,:,LatentAComparts)=inf_overallProb_AnalIns * Params.store_relRiskSexTransmission(Params.stage_LatentA);
        inf_basePerActProb_AnalIns(:,:,LatentBComparts)=inf_overallProb_AnalIns * Params.store_relRiskSexTransmission(Params.stage_LatentB);
        inf_basePerActProb_AnalIns(:,:,LateComparts)=inf_overallProb_AnalIns * Params.store_relRiskSexTransmission(Params.stage_Late);
        inf_basePerActProb_AnalIns(:,:,AIDSComparts)=inf_overallProb_AnalIns * Params.store_relRiskSexTransmission(Params.stage_AIDS);
    
   % Vaginal receptive

        % Pre-allocate
        inf_basePerActProb_VaginalRec = zeros(Params.numStrats, Params.numStrats, Params.numComparts);

        % Apply relative risk by disease stage
        inf_basePerActProb_VaginalRec(:,:,AcuteComparts)=inf_overallProb_VaginalRec * Params.store_relRiskSexTransmission(Params.stage_Acute);
        inf_basePerActProb_VaginalRec(:,:,LatentAComparts)=inf_overallProb_VaginalRec * Params.store_relRiskSexTransmission(Params.stage_LatentA);
        inf_basePerActProb_VaginalRec(:,:,LatentBComparts)=inf_overallProb_VaginalRec * Params.store_relRiskSexTransmission(Params.stage_LatentB);
        inf_basePerActProb_VaginalRec(:,:,LateComparts)=inf_overallProb_VaginalRec * Params.store_relRiskSexTransmission(Params.stage_Late);
        inf_basePerActProb_VaginalRec(:,:,AIDSComparts)=inf_overallProb_VaginalRec * Params.store_relRiskSexTransmission(Params.stage_AIDS);
    
   % Anal receptive
        
        % Pre-allocate
        inf_basePerActProb_AnalRec = zeros(Params.numStrats, Params.numStrats, Params.numComparts);
 
        % Apply relative risk by disease stage
        inf_basePerActProb_AnalRec(:,:,AcuteComparts)=inf_overallProb_AnalRec * Params.store_relRiskSexTransmission(Params.stage_Acute);
        inf_basePerActProb_AnalRec(:,:,LatentAComparts)=inf_overallProb_AnalRec * Params.store_relRiskSexTransmission(Params.stage_LatentA);
        inf_basePerActProb_AnalRec(:,:,LatentBComparts)=inf_overallProb_AnalRec * Params.store_relRiskSexTransmission(Params.stage_LatentB);
        inf_basePerActProb_AnalRec(:,:,LateComparts)=inf_overallProb_AnalRec * Params.store_relRiskSexTransmission(Params.stage_Late);
        inf_basePerActProb_AnalRec(:,:,AIDSComparts)=inf_overallProb_AnalRec * Params.store_relRiskSexTransmission(Params.stage_AIDS);
 
% NEEDLE TRANSMISSION PROBABILITY

    % Read in from Excel
    store_needleTransmissionProb = ExcelValues_AllParameters(parameter_Index : parameter_Index + numHIVstages-1);
    parameter_Index = parameter_Index + numel(store_needleTransmissionProb);


    % Apply to IDU-IDU combinations, by disease stage

        % Preallocate
        inf_indicator_IDU_IDU = zeros(Params.numStrats,Params.numStrats,Params.numComparts);

        % Only allow IDU-IDU combinations
            %     for c = (1 + Params.numHIVNegStates):(Params.numComparts - Params.numAbsorbingStates)
            for c = Params.B1:Params.F5
                inf_indicator_IDU_IDU(:,:,c) = Params.idx_IDU_IDU;
            end

        % Adjust from 273x273x30 to 30x273x273
        inf_indicator_IDU_IDU = permute(inf_indicator_IDU_IDU, [3 1 2]);

        % Save indicator for use in UpdateParams.m
        Params.idx_3d_IDU_IDU = inf_indicator_IDU_IDU;

        % Use indicator to apply the tranmission probabilities to a full
        % matrix [30x273x273]
        
            % Pre-allocate
            inf_basePerNeedleTransProb = zeros(Params.numComparts,Params.numStrats, Params.numStrats);

        inf_basePerNeedleTransProb(AcuteComparts,:,:) = inf_indicator_IDU_IDU(AcuteComparts,:,:)*store_needleTransmissionProb(Params.stage_Acute);
        inf_basePerNeedleTransProb(LatentAComparts,:,:) = inf_indicator_IDU_IDU(LatentAComparts,:,:)*store_needleTransmissionProb(Params.stage_LatentA);
        inf_basePerNeedleTransProb(LatentBComparts,:,:) = inf_indicator_IDU_IDU(LatentBComparts,:,:)*store_needleTransmissionProb(Params.stage_LatentB);
        inf_basePerNeedleTransProb(LateComparts,:,:) = inf_indicator_IDU_IDU(LateComparts,:,:)*store_needleTransmissionProb(Params.stage_Late);
        inf_basePerNeedleTransProb(AIDSComparts,:,:) = inf_indicator_IDU_IDU(AIDSComparts,:,:)*store_needleTransmissionProb(Params.stage_AIDS);

    
% TRANSMISSION REDUCTION DUE TO CIRCUMCISION
    %Reduction in HIV transmission per contact if circumcised versus 
    %uncircumcised
    
    % VI   
        % Read in scalar from Excel
        inf_circReductVI_MF_1 = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(inf_circReductVI_MF_1);

       % Apply scalar to full 273x273x30 matrix where population 1 is
       % circumcised
       inf_circReductVI_MF_273x273(Params.numStrats,Params.numStrats) = 0;
        for s2 = 1:Params.numStrats
           inf_circReductVI_MF_273x273(:,s2) = inf_circReductVI_MF_1 * Params.circIndicator(:,2);
        end

        Params.inf_circReductVI_MF = repmat(inf_circReductVI_MF_273x273,[1,1,Params.numComparts]);


    % AI
        % Male-Male
        Params.inf_circReductAI_MM = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(Params.inf_circReductAI_MM);

        % Male-Female
        Params.inf_circReductAI_MF = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(Params.inf_circReductAI_MF);

    % Populate the AI circumcsion 
        % Preallocate
            circReductAI_MF_273_273=zeros(Params.numStrats,Params.numStrats);

        % Apply M-F AI circ reduct factor
            % This is initially applied to all circumcised men; this will be
            % replaced by the MSM-MSM factor in the next step
        for s2 = 1:Params.numStrats
           circReductAI_MF_273_273(:,s2) = Params.inf_circReductAI_MF * Params.circIndicator(:,2);
        end

    % Apply a zero to any MSM-MSM combinations
     circReductAI_MF_273_273(Params.idx_MSMandMSM)=0;
        
    % Calculate an index for combinations of partner 1: circumcised MSM and 
    % partner 2: any MSM 
        % We will apply the circ reduct for MSM-MSM contacts to this index
        
        % Preallocate
        idx_circMSM_273_273 = zeros(Params.numStrats,Params.numStrats);
        
        for s2 = 1:Params.numStrats
            idx_circMSM_273_273(:,s2) = Params.idx_MSM_MSM(:,s2).*Params.circIndicator(:,2);
        end

    % Determine the final circumcision reduction for all AI contacts
        % Calculated by adding the MSM-MSM reduction to the current MF reduction (circReduct_273_273)
    circReductAI_273_273 = idx_circMSM_273_273 * Params.inf_circReductAI_MM + circReductAI_MF_273_273;
    
    %<-- ET, update UpdateParams
    Params.circReductAI_273x273x30 = repmat(circReductAI_273_273,[1,1,Params.numComparts]);

% TRANSMISSION REDUCTION DUE TO SEP

    % Read in from Excel
    Params.SEPTransmissionReduct = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + numel(Params.SEPTransmissionReduct);    
    
% TRANSMISSION REDUCTION DUE TO PrEP % JCPrEPUpdate: modified to read in
% both oral and injectable and for high and low adherence
    
    % Read in from Excel - Oral - High Adherence [3x1]
    store_prepTransmissionReduct_Oral_HighAdherence(:,1) = ExcelValues_AllParameters(parameter_Index:parameter_Index+numPop-1);
    parameter_Index = parameter_Index + numel(store_prepTransmissionReduct_Oral_HighAdherence);
    
    % Read in from Excel - Oral - Low Adherence [3x1]
    store_prepTransmissionReduct_Oral_LowAdherence(:,1) = ExcelValues_AllParameters(parameter_Index:parameter_Index+numPop-1);
    parameter_Index = parameter_Index + numel(store_prepTransmissionReduct_Oral_LowAdherence);
    
    % Read in from Excel - Injectable - High Adherence [3x1]
    store_prepTransmissionReduct_Inject_HighAdherence(:,1) = ExcelValues_AllParameters(parameter_Index:parameter_Index+numPop-1);
    parameter_Index = parameter_Index + numel(store_prepTransmissionReduct_Inject_HighAdherence);
    
    % Read in from Excel - Injectable - Low Adherence [3x1]
    store_prepTransmissionReduct_Inject_LowAdherence(:,1) = ExcelValues_AllParameters(parameter_Index:parameter_Index+numPop-1);
    parameter_Index = parameter_Index + numel(store_prepTransmissionReduct_Inject_LowAdherence);

    % Apply to full Param [273x1]
    Params.inf_prepIncidenceReduct_Oral_HighAdherence = ...
        Params.popIndicator * store_prepTransmissionReduct_Oral_HighAdherence;
    
    Params.inf_prepIncidenceReduct_Oral_LowAdherence = ...
        Params.popIndicator * store_prepTransmissionReduct_Oral_LowAdherence;
    
    Params.inf_prepIncidenceReduct_Inject_HighAdherence = ...
        Params.popIndicator * store_prepTransmissionReduct_Inject_HighAdherence;
    
    Params.inf_prepIncidenceReduct_Inject_LowAdherence = ...
        Params.popIndicator * store_prepTransmissionReduct_Inject_LowAdherence;

% TRANSMISSION REDUCTION DUE TO VLS

    % Sexual Reduction
        % Reduction in HIV transmission per sexual contact if VLS (vs. not VLS)
        
        % Read in from Excel: stratified by transmission group [3x1]
            store_vlsReductS(:,1)=ExcelValues_AllParameters(parameter_Index:parameter_Index+numPop-1);
            parameter_Index = parameter_Index + numel(store_vlsReductS);
        
        % Apply to full subpopulation matrix [273x1]
            vlsReductS = Params.popIndicator * store_vlsReductS;
        
        % Apply to full partner matrix [273x273]
            Params.inf_vlsReductS_273x273 = repmat(vlsReductS,[1,Params.numStrats]);
        

    % Needle Reduction
        %Reduction in HIV transmission per shared needle if VLS (vs. not VLS)
        
        Params.inf_vlsReductN = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + numel(Params.inf_vlsReductN);

        

% PROPORTION RECEPTIVE
    % HETF, IDUF 100% receptive
    % MSM-MSM: stratified by age; MSM-Female: 0% receptive
    % [273x273x1]

    % Male-Female
    storeMF_rec(Params.numStrats,Params.numStrats,Params.numComparts)=0;
    for c = 1:Params.numComparts % compartments
        storeMF_rec(:,:,c) = Params.idx_Male_Female*0;
    end

    %Male-Male-defined above as [273 x 1]
    storeMM_rec(Params.numStrats,Params.numStrats,Params.numComparts)=0;
    %Converting to [273x273]
    behav_pctMSMRec = repmat(behav_pctMSMRec,[1,Params.numStrats]);
    %Converting to [273 x 273 x 30]
    for c = 1:Params.numComparts % compartments
        storeMM_rec(:,:,c) = Params.idx_MSM_MSM .* behav_pctMSMRec;
    end

    % Female-Male
    storeFM_rec(Params.numStrats,Params.numStrats,Params.numComparts)=0;
    for c = 1:Params.numComparts % compartments
        storeFM_rec(:,:,c) = Params.idx_Female_Male;
    end

    % Full matrix: all partnerships
  Params.factor_proportionRec = storeMF_rec + storeMM_rec + storeFM_rec;

%% 13. Beta Section

    % 13.i Calculate values that go into the betas
    for betaFactorAdjustmentSection = 1:1
       

% ADJUSTED TRANSMISSION RISKS
        
    % Vaginal Insertive Risk
    
        % Reductions due to
            % Circumcision
            % VLS
        
        % Apply circumcision reduction
        BetaFactor.adjPerActVaginalInsProb = ...
            inf_basePerActProb_VaginalIns .* ...
            (1-Params.inf_circReductVI_MF);
        
        % Apply VLS reduction
            for c = VLSComparts
                BetaFactor.adjPerActVaginalInsProb(:,:,c)= ...
                    BetaFactor.adjPerActVaginalInsProb(:,:,c) .*...
                    (1 - Params.inf_vlsReductS_273x273);
            end
            

    % Vaginal Receptive Risk
        % Adjusted per-unprotected-receptive-act sexual risk 

        % Reduction due to 
            % VLS
   
        % Preallocate: apply base risk to full 273x273x30 matrix
        BetaFactor.adjPerActVaginalRecProb = inf_basePerActProb_VaginalRec;
        
        % Apply VLS reduction
        for c = VLSComparts
            BetaFactor.adjPerActVaginalRecProb(:,:,c) = ...
                inf_basePerActProb_VaginalRec(:,:,c) .* ...
                (1 - Params.inf_vlsReductS_273x273);
        end


    % Anal Insertive Risk
    
        % Reductions due to 
            % Circumcision
            % VLS
    
        % Apply circumcision reduction
        BetaFactor.adjPerActAnalInsProb = ...
            inf_basePerActProb_AnalIns .* ...
            (1 -Params.circReductAI_273x273x30);
         
        % Apply VLS Reduction
        for c = VLSComparts
           BetaFactor.adjPerActAnalInsProb(:,:,c) = ...     
                BetaFactor.adjPerActAnalInsProb(:,:,c) .* ...
                (1 - Params.inf_vlsReductS_273x273);
        end

    
    % Anal Receptive Risk
    
        % Reduction due to 
            % VLS
            
        % Preallocate: apply base risk to full 273x273x30 matrix
        BetaFactor.adjPerActAnalRecProb = inf_basePerActProb_AnalRec;
        
         % Apply VLS Reduction
        for c = VLSComparts
            BetaFactor.adjPerActAnalRecProb(:,:,c) = ...
                inf_basePerActProb_AnalRec(:,:,c) .* ...
                (1 - Params.inf_vlsReductS_273x273);
        end

    end

    % 13.ii Calculate non-needle betas   
    
        % Note: the input argument "1" designates that the function is called
        % from InitializeParams (rather than during an Calib_updateParams
    Params = CalcBetas(Params,BetaFactor,1);
    
    % 13.iii Calculate needle betas
    for betaNSection = 1:1
    
        
% NUMBER OF NEEDLES SHARED (ADJUSTED)
    % Adjustment based on awareness of status
    
    % Pre-allocate
    inf_adjNeedlesShared_NoPrEP_1to4 = zeros(Params.numComparts,Params.numStrats,Params.numStrats);
    inf_adjNeedlesShared_NoPrEP_5 = zeros(Params.numComparts,Params.numStrats,Params.numStrats);
    
    % Adjust for awareness
        behav_needlesSharedUnaware = behav_baseNeedlesPerPartner;
        behav_needlesSharedAware_1to4 = behav_baseNeedlesPerPartner * (1-Params.behav_needleReductDueToAware_1to4);
        behav_needlesSharedAware_5 = behav_baseNeedlesPerPartner * (1-Params.behav_needleReductDueToAware_5);

        % Apply to parameter
                % Note: It's ok that they're all filled because the betas will remove the
                % compartments and people who won't be infected via needle.
            inf_adjNeedlesShared_NoPrEP_1to4(:,:,:) = behav_needlesSharedUnaware;
            inf_adjNeedlesShared_NoPrEP_5(:,:,:) = behav_needlesSharedUnaware;
                %Overwrite for the people who are aware of the their HIV+ status
            inf_adjNeedlesShared_NoPrEP_1to4(AwareComparts,:,:) = behav_needlesSharedAware_1to4;
            inf_adjNeedlesShared_NoPrEP_5(AwareComparts,:,:) = behav_needlesSharedAware_5;

    % Adjust for PrEP
        inf_adjNeedlesShared_PrEP_1to4 = ...
            inf_adjNeedlesShared_NoPrEP_1to4 .* (1 + Params.behav_incrNeedleShareDueToPrEP);
        inf_adjNeedlesShared_PrEP_5 = ...
            inf_adjNeedlesShared_NoPrEP_5 .* (1 + Params.behav_incrNeedleShareDueToPrEP);
    
        
% PER-NEEDLE PARTNERSHIP RISK OF INFECTION (BETA N)

    % Pre-allocate adjusted needle risk (set equal to base risk)
    inf_adjPerNeedleRisk = inf_basePerNeedleTransProb;

    % Reduce transmission risk for PLWH who are VLS
    inf_adjPerNeedleRisk(Params.VLSComparts,:,:) = inf_basePerNeedleTransProb(Params.VLSComparts,:,:)...
        *(1-Params.inf_vlsReductN);
    
    
    % Calculate betas
    
        % No PrEP
        inf_betaN_NoPrEP_1to4 = 1 - (1 - inf_adjPerNeedleRisk) .^ inf_adjNeedlesShared_NoPrEP_1to4;
        Params.oneMinusBetaN_NoPrEP_1to4 = 1 - inf_betaN_NoPrEP_1to4;
        
        inf_betaN_NoPrEP_5 = 1 - (1 - inf_adjPerNeedleRisk) .^ inf_adjNeedlesShared_NoPrEP_5;
        Params.oneMinusBetaN_NoPrEP_5 = 1 - inf_betaN_NoPrEP_5;

        % With PrEP
            % Only calculate when PrEP is run with behavior changes
        if  Params.CalcPrEPSpecificBetas == Params.ind_PrEP            
            inf_betaN_PrEP_1to4 = 1 - (1 - inf_adjPerNeedleRisk) .^ inf_adjNeedlesShared_PrEP_1to4;
            Params.oneMinusBetaN_PrEP_1to4 = 1 - inf_betaN_PrEP_1to4;
            inf_betaN_PrEP_5 = 1 - (1 - inf_adjPerNeedleRisk) .^ inf_adjNeedlesShared_PrEP_5;
            Params.oneMinusBetaN_PrEP_5 = 1 - inf_betaN_PrEP_5;
        end
    end

%% 14. Define Intervention inputs
   
% COST PER PERSON

    % Testing
    % Cost per initial test and confirmatory test in first testing period
        % Negative Rapid Tests
        costPP_TestNegRapid_TestPeriod1 = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1; 

        % Positive Rapid Tests
        costPP_TestPosRapid_TestPeriod1 = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1; 

        % Negative Conventional Tests
        costPP_TestNegConv_TestPeriod1 = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1;

        % Positive Conventional Tests
        costPP_TestPosConv_TestPeriod1 = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1;
    
    
    % Cost per initial test and confirmatory test in second testing period
        % Negative Rapid Tests
        costPP_TestNegRapid_TestPeriod2 = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1; 

        % Positive Rapid Tests
        costPP_TestPosRapid_TestPeriod2 = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1; 

        % Negative Conventional Tests
        costPP_TestNegConv_TestPeriod2 = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1;

        % Positive Conventional Tests
        costPP_TestPosConv_TestPeriod2 = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1;
    
    
    % Cost of NAT test
        % Applied when the confirmatory screen is negative
    Params.costPP_NATTest = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
    % Test outreach cost
        % Low-risk HETs
        costPP_TestOutreach_HETs_Low = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1;

        % High-risk HETs
        costPP_TestOutreach_HETs_High = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1;
    
    % Calculated: Apply test outreach costs to HET populations
    
        % Pre-allocate
        Params.costPP_TestOutreach(1,Params.numStrats)=0;
    
        % Apply it
        for i = 1:Params.numStrats

            % HETs
            if Params.popIndicator(i,Params.pop_HET)==1 %HET
                
                % Low-risk
                if Params.riskLevelIndicator(i,Params.risk_Main)==1 %Low risk
                    Params.costPP_TestOutreach(1,i) = costPP_TestOutreach_HETs_Low;
                
                % High-risk
                else 
                    Params.costPP_TestOutreach(1,i) = costPP_TestOutreach_HETs_High;
                end
            end
        end
   
    % Additional cost per test performed in non-clinical (vs. clinical) setting
    costAddl_nonclinical = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
    % Calculated: Weighted costs per test (based on percent that are
    % non-clinical)
    
    Params.costPP_TestNegRapid_TestPeriod1 = costPP_TestNegRapid_TestPeriod1 + (Params.tt_pctnonclinical * costAddl_nonclinical);
    Params.costPP_TestPosRapid_TestPeriod1 = costPP_TestPosRapid_TestPeriod1 + (Params.tt_pctnonclinical * costAddl_nonclinical);
    Params.costPP_TestNegConv_TestPeriod1 = costPP_TestNegConv_TestPeriod1 + (Params.tt_pctnonclinical * costAddl_nonclinical);
    Params.costPP_TestPosConv_TestPeriod1 = costPP_TestPosConv_TestPeriod1 + (Params.tt_pctnonclinical * costAddl_nonclinical);
    Params.costPP_TestNegRapid_TestPeriod2 = costPP_TestNegRapid_TestPeriod2 + (Params.tt_pctnonclinical * costAddl_nonclinical);
    Params.costPP_TestPosRapid_TestPeriod2 = costPP_TestPosRapid_TestPeriod2 + (Params.tt_pctnonclinical * costAddl_nonclinical);
    Params.costPP_TestNegConv_TestPeriod2 = costPP_TestNegConv_TestPeriod2 + (Params.tt_pctnonclinical * costAddl_nonclinical);
    Params.costPP_TestPosConv_TestPeriod2 = costPP_TestPosConv_TestPeriod2 + (Params.tt_pctnonclinical * costAddl_nonclinical);
        
    % Notification of negative rapid test
    Params.costPP_Notify_NegRapid = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
    % Notification of a positive rapid test
    Params.costPP_Notify_PosRapid = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
    % Notification of a negative conventional test
    Params.costPP_Notify_NegConv = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
    % Notification of a positive conventional test
    Params.costPP_Notify_PosConv = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;

    % PrEP costs 
    % ...are in the HIV section (10) with the costs by other compartments
    % Params.hiv_annualCostPerCompart_1(Params.A6), _2, _3 for 3 time periods
                
    % Linkage to care at diagnosis
    Params.costPP_LTCFirst = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;

    % Linkage to care after diagnosis
    Params.costPP_LTCAfter = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
    % Initiating ART
    Params.costPP_ARTInitiation = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;    
    
    % Treatment Adherence (to become VLS)
    Params.costPP_TxAdherence_BecomeVLS = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;

    % Treatment Compliance (to remain VLS)
    Params.costPP_TxAdherence_RemainVLS = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;
    
    % Syringe exchange program
    Params.costPP_SyringeExchange = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;

    % Years for applying each allocation
    Params.Intn_Allocation_StartYr(1) = Params.tt_periodFiveStartYear;
    
    for allocationPeriod = 2:Params.numMaxAllocationPeriods
        Params.Intn_Allocation_StartYr(allocationPeriod) = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1;    
    end
    
    % Determine number of allocation periods being used 
    TimeHorizonLastYear = Params.tt_modelStartYear + Params.tt_timeHorizon - 1;
    if TimeHorizonLastYear <= Params.Intn_Allocation_StartYr(2)
        Params.numAllocationPeriods = 1;
    elseif TimeHorizonLastYear <= Params.Intn_Allocation_StartYr(3)
        Params.numAllocationPeriods = 2;
    else
        Params.numAllocationPeriods = 3;
    end    
    
    % Determine outcome indices associated with each allocation period 
    Params.OutcomesIndex_AllocationStartYr(1) =  ...
        Params.Intn_Allocation_StartYr(1) - FirstOutcomeYr + 1;
    Params.OutcomesIndex_AllocationEndYr(1) =  ...
        Params.OutcomesIndex_AllocationStartYr(1) + ...
        (min(Params.Intn_Allocation_StartYr(2)-1,LastOutcomeYr) - ...
            Params.Intn_Allocation_StartYr(1));
    if Params.numAllocationPeriods >= 2
        Params.OutcomesIndex_AllocationStartYr(2) =  ...
            Params.OutcomesIndex_AllocationEndYr(1) + 1;
        Params.OutcomesIndex_AllocationEndYr(2) =  ...
            Params.OutcomesIndex_AllocationStartYr(2) + ...
            (min(Params.Intn_Allocation_StartYr(3)-1,LastOutcomeYr) - ...
                Params.Intn_Allocation_StartYr(2));
    end
    if Params.numAllocationPeriods == 3
        Params.OutcomesIndex_AllocationStartYr(3) =  ...
            Params.OutcomesIndex_AllocationEndYr(2) + 1;
        Params.OutcomesIndex_AllocationEndYr(3) =  ...
            Params.OutcomesIndex_AllocationStartYr(3) + ...
            (LastOutcomeYr - Params.Intn_Allocation_StartYr(3));
    end
    

% CDC INVESTMENT IN INTERVENTIONS
    for allocationPeriod = 1:Params.numMaxAllocationPeriods
        
        for investmentInIntn = 1:Params.numIntns

            importedVariable = ExcelValues_AllParameters(parameter_Index);
            parameter_Index = parameter_Index + 1;

            switch investmentInIntn
                case 1
                    Params.intn_Testing_Investment_HET_Low(allocationPeriod) = importedVariable;
                case 2
                    Params.intn_Testing_Investment_HET_High(allocationPeriod) = importedVariable;
                case 3
                    Params.intn_Testing_Investment_MSM_Low(allocationPeriod) = importedVariable;
                case 4
                    Params.intn_Testing_Investment_MSM_High(allocationPeriod) = importedVariable;
                case 5
                    Params.intn_Testing_Investment_IDU(allocationPeriod) = importedVariable;
                case 6
                    Params.intn_LTCatDiag_Investment(allocationPeriod) = importedVariable;
                case 7
                    Params.intn_LTCafterDiag_Investment(allocationPeriod) = importedVariable;
                case 8
                    Params.intn_ARTInitiation_Investment(allocationPeriod) = importedVariable;
                case 9
                    Params.intn_ARTAdher5to4_Investment(allocationPeriod) = importedVariable;
                case 10
                    Params.intn_ARTAdher4to5_Investment(allocationPeriod) = importedVariable;
                case 11
                    Params.intn_SEP_Investment_B(allocationPeriod) = importedVariable;
                case 12
                    Params.intn_SEP_Investment_H(allocationPeriod) = importedVariable;
                case 13
                    Params.intn_SEP_Investment_O(allocationPeriod) = importedVariable;    
                case 14
                    Params.intn_PrEP_Oral_Investment_HETM_B(allocationPeriod) = importedVariable;
                case 15
                    Params.intn_PrEP_Oral_Investment_HETM_H(allocationPeriod) = importedVariable;
                case 16
                    Params.intn_PrEP_Oral_Investment_HETM_O(allocationPeriod) = importedVariable;
                case 17
                    Params.intn_PrEP_Oral_Investment_HETF_B(allocationPeriod) = importedVariable;
                case 18
                    Params.intn_PrEP_Oral_Investment_HETF_H(allocationPeriod) = importedVariable; 
                case 19
                    Params.intn_PrEP_Oral_Investment_HETF_O(allocationPeriod) = importedVariable; 
                case 20
                    Params.intn_PrEP_Oral_Investment_MSM_B(allocationPeriod) = importedVariable;
                case 21
                    Params.intn_PrEP_Oral_Investment_MSM_H(allocationPeriod) = importedVariable;
                case 22
                    Params.intn_PrEP_Oral_Investment_MSM_O(allocationPeriod) = importedVariable;
                case 23
                    Params.intn_PrEP_Oral_Investment_IDU_B(allocationPeriod) = importedVariable;
                case 24
                    Params.intn_PrEP_Oral_Investment_IDU_H(allocationPeriod) = importedVariable;
                case 25
                    Params.intn_PrEP_Oral_Investment_IDU_O(allocationPeriod) = importedVariable;
                case 26
                    Params.intn_PrEP_Inject_Investment_HETM_B(allocationPeriod) = importedVariable;
                case 27
                    Params.intn_PrEP_Inject_Investment_HETM_H(allocationPeriod) = importedVariable;
                case 28
                    Params.intn_PrEP_Inject_Investment_HETM_O(allocationPeriod) = importedVariable;
                case 29
                    Params.intn_PrEP_Inject_Investment_HETF_B(allocationPeriod) = importedVariable;
                case 30
                    Params.intn_PrEP_Inject_Investment_HETF_H(allocationPeriod) = importedVariable;
                case 31
                    Params.intn_PrEP_Inject_Investment_HETF_O(allocationPeriod) = importedVariable;
                case 32
                    Params.intn_PrEP_Inject_Investment_MSM_B(allocationPeriod) = importedVariable;
                case 33
                    Params.intn_PrEP_Inject_Investment_MSM_H(allocationPeriod) = importedVariable;
                case 34
                    Params.intn_PrEP_Inject_Investment_MSM_O(allocationPeriod) = importedVariable;
                case 35
                    Params.intn_PrEP_Inject_Investment_IDU_B(allocationPeriod) = importedVariable;
                case 36
                    Params.intn_PrEP_Inject_Investment_IDU_H(allocationPeriod) = importedVariable;
                case 37
                    Params.intn_PrEP_Inject_Investment_IDU_O(allocationPeriod) = importedVariable;
            end
        end
    end

% REACH LEVELS BEYOND T2 AND INCREASE IN INTERVENTION COSTS WHEN REACH BEYOND T2


for r= 1:2
    for HighReachCostIncr = 1:Params.numIntnIncrCostCategories                            
        
        importedVariable = ExcelValues_AllParameters(parameter_Index:parameter_Index+Params.numReachLevels-1);
        parameter_Index = parameter_Index + Params.numReachLevels;
        
        switch r
            case 1
                switch HighReachCostIncr
                    case 1
                        Params.intn_Testing_ReachLvls = importedVariable;
                    case 2
                        Params.intn_LTCatDiag_ReachLvls = importedVariable;
                    case 3
                        Params.intn_LTCafterDiag_ReachLvls = importedVariable;
                    case 4
                        Params.intn_ARTInitiation_ReachLvls = importedVariable;
                    case 5
                        Params.intn_TxAdherence_ReachLvls = importedVariable;
                    case 6
                        Params.intn_SEP_ReachLvls = importedVariable;
                    case 7
                        Params.intn_PrEP_Oral_ReachLvls = importedVariable;
                    case 8
                        Params.intn_PrEP_Inject_ReachLvls = importedVariable;
                end
            case 2
                switch HighReachCostIncr
                    case 1
                        Params.intn_Testing_PctCostatReachLvl = importedVariable;
                    case 2
                        Params.intn_LTCatDiag_PctCostatReachLvl = importedVariable;
                    case 3
                        Params.intn_LTCafterDiag_PctCostatReachLvl = importedVariable;
                    case 4
                        Params.intn_ARTInitiation_PctCostatReachLvl = importedVariable;
                    case 5
                        Params.intn_TxAdherence_PctCostatReachLvl = importedVariable;
                    case 6
                        Params.intn_SEP_PctCostatReachLvl = importedVariable;
                    case 7
                        Params.intn_PrEP_Oral_PctCostatReachLvl = importedVariable;
                    case 8
                        Params.intn_PrEP_Inject_PctCostatReachLvl = importedVariable;
                end
        end
    end
    
end    

    
% RESOURCE ALLOCATION OPTIMIZATION OBJECTIVE
    Params.optObjNum = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;  
    
% SETTINGS FOR WHEN DOING COST MINIMIZATION
    Params.minCost_Tgt1NumYears = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;  

    Params.minCost_Tgt1PctReduce = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;  
    
    Params.minCost_Tgt2NumYears = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1;  

    Params.minCost_Tgt2PctReduce = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1; 
    
% TOTAL BUDGET FOR IMPLEMENTING INTERVENTIONS IN EACH ALLOCATION PERIOD
    for allocationPeriod = 1:Params.numMaxAllocationPeriods
        Params.alloc_TotalIntnBudget(allocationPeriod) = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1;    
    end
    
% INDICATOR TO REDUCE FUNDING AVAILABLE FOR INTERVENTIONS BY HIV TREATMENT
% AND CARE COSTS
    Params.alloc_ConsiderARTAndCareCosts = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1; 

% TOTAL ESTIMATED SPENDING ON HIV TREATMENT AND CARE
    Params.alloc_TotalARTAndCareCosts = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1; 
        
% PERCENTAGE OF TREATMENT AND CARE COSTS TO BE COVERED BY THIS BUDGET (VS OTHER FUNDING SOURCES
    Params.alloc_PctARTAndCareCoveredByBudget = ExcelValues_AllParameters(parameter_Index);
    parameter_Index = parameter_Index + 1; 
    
% INTN REACH ACROSS ELIGIBLE POPULATIONS

    % Minimum reach
    for minReach = 1:Params.numIntnReachCategories 
        
        importedVariable = ExcelValues_AllParameters(parameter_Index);
        parameter_Index = parameter_Index + 1;

        switch minReach
            case 1
                Params.Testing_MinReach = importedVariable;
            case 2
                Params.intn_LTCatDiag_MinReach = importedVariable;
            case 3
                Params.intn_LTCafterDiag_MinReach = importedVariable;
            case 4
                Params.intn_ARTInitiation_MinReach = importedVariable;
            case 5
                Params.intn_TxAdherence_MinReach = importedVariable; 
            case 6
                Params.intn_SEP_MinReach = importedVariable;     
            case 7
                intn_PrEP_Oral_HET_MinReach = importedVariable .* Params.popIndicator(:,Params.pop_HET) .* Params.riskLevelIndicator(:,Params.risk_Casual);
            case 8
                intn_PrEP_Oral_MSM_MinReach = importedVariable .* Params.popIndicator(:,Params.pop_MSM) .* Params.riskLevelIndicator(:,Params.risk_Casual); 
            case 9
                intn_PrEP_Oral_IDU_MinReach = importedVariable .* Params.popIndicator(:,Params.pop_IDU); 
            case 10
                intn_PrEP_Inject_HET_MinReach = importedVariable .* Params.popIndicator(:,Params.pop_HET) .* Params.riskLevelIndicator(:,Params.risk_Casual);
            case 11
                intn_PrEP_Inject_MSM_MinReach = importedVariable .* Params.popIndicator(:,Params.pop_MSM) .* Params.riskLevelIndicator(:,Params.risk_Casual);
            case 12
                intn_PrEP_Inject_IDU_MinReach = importedVariable .* Params.popIndicator(:,Params.pop_IDU);
        end
    end

    Params.intn_PrEP_Oral_MinReach = intn_PrEP_Oral_HET_MinReach + intn_PrEP_Oral_MSM_MinReach + intn_PrEP_Oral_IDU_MinReach;
    Params.intn_PrEP_Inject_MinReach = intn_PrEP_Inject_HET_MinReach + intn_PrEP_Inject_MSM_MinReach + intn_PrEP_Inject_IDU_MinReach; 
    
    % Maximum reach
        for maxReach = 1:Params.numIntnReachCategories 
        
            importedVariable = ExcelValues_AllParameters(parameter_Index);
            parameter_Index = parameter_Index + 1;

            switch maxReach
                case 1
                    Params.intn_Testing_MaxReach = importedVariable;
                case 2
                   Params.intn_LTCatDiag_MaxReach = importedVariable;
                case 3
                   Params.intn_LTCafterDiag_MaxReach = importedVariable;
                case 4
                    Params.intn_ARTInitiation_MaxReach = importedVariable;
                case 5
                    Params.intn_TxAdherence_MaxReach = importedVariable; 
                case 6
                    Params.intn_SEP_MaxReach = importedVariable;     
                case 7
                    Params.intn_PrEP_Oral_HET_MaxReach = importedVariable .* Params.popIndicator(:,Params.pop_HET) .* Params.riskLevelIndicator(:,Params.risk_Casual);
                case 8
                    Params.intn_PrEP_Oral_MSM_MaxReach = importedVariable .* Params.popIndicator(:,Params.pop_MSM) .* Params.riskLevelIndicator(:,Params.risk_Casual); 
                case 9
                    Params.intn_PrEP_Oral_IDU_MaxReach = importedVariable .* Params.popIndicator(:,Params.pop_IDU); 
                case 10
                    Params.intn_PrEP_Inject_HET_MaxReach = importedVariable .* Params.popIndicator(:,Params.pop_HET) .* Params.riskLevelIndicator(:,Params.risk_Casual);
                case 11
                    Params.intn_PrEP_Inject_MSM_MaxReach = importedVariable .* Params.popIndicator(:,Params.pop_MSM) .* Params.riskLevelIndicator(:,Params.risk_Casual);
                case 12
                    Params.intn_PrEP_Inject_IDU_MaxReach = importedVariable .* Params.popIndicator(:,Params.pop_IDU);
             end
        end

        Params.intn_PrEP_Oral_MaxReach = Params.intn_PrEP_Oral_HET_MaxReach + Params.intn_PrEP_Oral_MSM_MaxReach + Params.intn_PrEP_Oral_IDU_MaxReach;
        Params.intn_PrEP_Inject_MaxReach = Params.intn_PrEP_Inject_HET_MaxReach + Params.intn_PrEP_Inject_MSM_MaxReach + Params.intn_PrEP_Inject_IDU_MaxReach; 
    
    % Allocation to interventions
        % Minimum allocation to interventions
        for allocationPeriod = 1:Params.numMaxAllocationPeriods
            Params.minalloc(:,allocationPeriod) = ExcelValues_AllParameters(parameter_Index:parameter_Index+Params.numIntns-1,1);
            parameter_Index = parameter_Index + Params.numIntns;    
        end


        % Maximum allocation to interventions
        for allocationPeriod = 1:Params.numMaxAllocationPeriods
            Params.maxalloc(:,allocationPeriod) = ExcelValues_AllParameters(parameter_Index:parameter_Index+Params.numIntns-1,1);
            parameter_Index = parameter_Index + Params.numIntns;    
        end         
    
% MINIMUM ALLOCATION (BY SUBPOPULATION)
         
    % By Sex
    Params.alloc_Sex_Min(Params.sex_Male:Params.sex_Female) = ExcelValues_AllParameters(parameter_Index ...
        : parameter_Index + 1);
    parameter_Index = parameter_Index + numel(Params.alloc_Sex_Min);

    % By Transmission Group
    Params.alloc_Population_Min(Params.pop_HET:Params.pop_IDU) = ExcelValues_AllParameters(parameter_Index ...
        : parameter_Index + (Params.pop_IDU-Params.pop_HET)); % IDU(3) - HET(1) = 2 (the additional number of pops to add)
            % where 3 and 1 are the indices of IDU and HET respectively
    parameter_Index = parameter_Index + numel(Params.alloc_Population_Min);

    % By Race
    Params.alloc_Race_Min(Params.race_B:Params.race_O) = ExcelValues_AllParameters(parameter_Index ...
        : parameter_Index + (Params.race_O-Params.race_B)); % Other(3) - Black(1) = 2 (the additional number of pops to add)
    parameter_Index = parameter_Index + numel(Params.alloc_Race_Min);
    
    % By Age
    Params.alloc_Age_Min(Params.age_1:Params.numAge) = ExcelValues_AllParameters(parameter_Index ...
        : parameter_Index + Params.numAge - 1); 
    parameter_Index = parameter_Index + numel(Params.alloc_Age_Min);
   
    % By Risk Level
    Params.alloc_Risk_Min(Params.risk_Main:Params.risk_Casual) = ExcelValues_AllParameters(parameter_Index ...
        : parameter_Index + (Params.risk_Casual-Params.risk_Main));
    parameter_Index = parameter_Index + numel(Params.alloc_Risk_Min);
 
    
% MAXIMUM ALLOCATION (BY SUBPOPULATION)
    
    % By Sex
    Params.alloc_Sex_Max(Params.sex_Male:Params.sex_Female) = ExcelValues_AllParameters(parameter_Index ...
        : parameter_Index + 1);
    parameter_Index = parameter_Index + numel(Params.alloc_Sex_Max);

    % By Transmission Group
    Params.alloc_Population_Max(Params.pop_HET:Params.pop_IDU) = ExcelValues_AllParameters(parameter_Index ...
        : parameter_Index + (Params.pop_IDU-Params.pop_HET)); % IDU(3) - HET(1) = 2 (the additional number of pops to add)
            % where 3 and 1 are the indices of IDU and HET respectively
    parameter_Index = parameter_Index + numel(Params.alloc_Population_Max);

    % By Race
    Params.alloc_Race_Max(Params.race_B:Params.race_O) = ExcelValues_AllParameters(parameter_Index ...
        : parameter_Index + (Params.race_O-Params.race_B)); % Other(3) - Black(1) = 2 (the additional number of pops to add)
    parameter_Index = parameter_Index + numel(Params.alloc_Race_Max);
    
    % By Age
    Params.alloc_Age_Max(Params.age_1:Params.numAge) = ExcelValues_AllParameters(parameter_Index ...
        : parameter_Index + Params.numAge - 1); 
    parameter_Index = parameter_Index + numel(Params.alloc_Age_Max);
   
    % By Risk Level
    Params.alloc_Risk_Max(Params.risk_Main:Params.risk_Casual) = ExcelValues_AllParameters(parameter_Index ...
        : parameter_Index + (Params.risk_Casual-Params.risk_Main));
    parameter_Index = parameter_Index + numel(Params.alloc_Risk_Max);  
    
    
% INITIAL ALLOCATION
    % Used for the optimization
    for allocationPeriod = 1:Params.numMaxAllocationPeriods
        Params.initalloc(:,allocationPeriod) = ExcelValues_AllParameters(parameter_Index:parameter_Index+Params.numIntns-1,1);
        parameter_Index = parameter_Index + Params.numIntns;    
    end
       
%% 15. Define Initial and New Populations (from Import_InitPop and Import_NewPop sheets)
  
% INITIAL POPULATION
    % Initial distribution of the modeled population
    % No one in absorbing states

    % Initialize index
    InitAndNewPop_Index = 1;

    % Pre-allocate
    store_initPop(Params.F5,Params.numStrats) = 0;
    Params.initPop = zeros(Params.numComparts,Params.numStrats);

    % Import values
    for nCompart = Params.NonAbsorbingComparts
        store_initPop(nCompart, :) = ExcelValues_Populations(InitAndNewPop_Index + Params.numStrats*(nCompart-1): ...
            (InitAndNewPop_Index+Params.numStrats-1)+Params.numStrats*(nCompart-1));
    end
    
    % Apply initial population to main parameter
    Params.initPop(Params.NonAbsorbingComparts,:) = store_initPop;
 
% NEW POPULATION
    % Number of people aging into the youngest age group of each
    % subpopulation every year
    
    % Moved to top of section 6 (because variables related to this area read from param list now)
    % And code that crates newPop matrix is at the bottom of Aging.m
   
%% 16. Conversion Functions

function OutRate = ProbToRate(InAnnualProb)
    if (InAnnualProb >= 1) 
        then InAnnualProb = 0.99999;
    end
            
    OutRate = -log(1 - InAnnualProb);
end

function OutProb = RateToProb(InAnnualRate)
    OutProb = 1 - exp(-InAnnualRate);
end

end