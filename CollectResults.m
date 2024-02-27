function Results = CollectResults(Params, Compartments, Year,InfRate, InfRate_PrEP_Oral_High, InfRate_PrEP_Oral_Low, InfRate_PrEP_Inject_High, InfRate_PrEP_Inject_Low, ContCounterOrDiscStep,Transitions, Results, TTProg)
%% Purpose: A function which stores results from each time step and

% summarizes outcomes in the last run of the model

% Called from: StepStateDisc (Discrete) and StepStateCont (Continuous)

% Nomenclature
    %set_   setting
    %ts_    data is available for every time step 
    %ann_   data is available each year


    
%% 1. Determine Years for which Outcomes Are Collected

    FirstOutcomeYr = max(Params.tt_modelStartYear,Params.outcomeCollectionStartYr);
    LastOutcomeYr = min(Params.tt_modelStartYear+Params.tt_timeHorizon-1,Params.outcomeCollectionEndYr);

if and(Year >= FirstOutcomeYr, floor(Year) <= LastOutcomeYr)
%% 2. Define Runtime Variables

    numYears = LastOutcomeYr - FirstOutcomeYr +1;
    numSteps = numYears/Params.tt_tstep;
    numStrats = Params.numStrats;
    FirstContCounterOrDiscStep = (FirstOutcomeYr - Params.tt_modelStartYear)/Params.tt_tstep + 1;
    
    % 2.i. Version-specific parameters
    for sectionVersion = 1:1
    
     % Set the increments
     
        % Continuous
        if Params.SolveCont == 1 % if continuous run
            
            % Set the increment (index comes from the ode timestep)
            if isfield(Results,'First_t')==0
                Results.First_t = TTProg.k;
            end
            
            odeCounter = ContCounterOrDiscStep;
            t = TTProg.k - Results.First_t + 1; 
            YrsOmitted = FirstOutcomeYr - Params.tt_modelStartYear;
            stepsPerYear = TTProg.stepsPerYear((1 + YrsOmitted) : (1 + Year - FirstOutcomeYr + YrsOmitted));
            tstep = TTProg.timePerODEStep(TTProg.k);
            
        % Discrete
        else 
            % Set the increment (index comes from the model tstep set by
            % the user in the Excel file under 'Model Settings')
            
            t = ContCounterOrDiscStep - FirstContCounterOrDiscStep + 1; 
            
            odeCounter = 1; % this doesn't do anything in the disc version
                % but needs to be included in the code so it's compatible with the Cont version
            stepsPerYear = numSteps/numYears;
            tstep = Params.tt_tstep;
        end
    end        
    
    % 2.ii. Period-specific parameters    
    for sectionPeriodSpec = 1:1
    if Year < Params.tt_periodTwoStartYear 
        PctANVWhoAreInCareInclART = Params.tt_PctANVWhoAreInCareOrART_1;
        PctANVInCareWhoAreOnART = Params.tt_PctInCareANVWhoAreOnART_1;
        AnnualHIVCostsByCompartment = Params.hiv_annualCostPerCompart_1;
        
    elseif Year < Params.tt_periodFiveStartYear
        PctANVWhoAreInCareInclART = Params.tt_PctANVWhoAreInCareOrART_2to4;
        PctANVInCareWhoAreOnART = Params.tt_PctInCareANVWhoAreOnART_2to4;
        AnnualHIVCostsByCompartment = Params.hiv_annualCostPerCompart_2to4;
        
    else
        PctANVWhoAreInCareInclART = Params.tt_PctANVWhoAreInCareOrART_5;
        PctANVInCareWhoAreOnART = Params.tt_PctInCareANVWhoAreOnART_5;
        AnnualHIVCostsByCompartment = Params.hiv_annualCostPerCompart_5;
    end
    end

%% 3. Define reference indices (to be used to collect data from the appropriate compartments)

    %TT Continuum Comparts
    Comparts_CCStage1 = Params.UnawareComparts;
    Comparts_CCStage2 = [Params.B2, Params.C2, Params.D2,Params.E2, Params.F2];
    Comparts_CCStages2345 = Params.AwareComparts;
    Comparts_CCStages345 = Params.CCStages345;
    Comparts_CCStage3 = [Params.B3,Params.C3,Params.D3,Params.E3,Params.F3];
    Comparts_CCStage4 = Params.ANVComparts;
    Comparts_CCStages45 = Params.ANVandVLSComparts;
    Comparts_CCStage5 = Params.VLSComparts;
    HIVComparts = Params.HIVComparts;
    Comparts_SusAndPLWH = Params.NonAbsorbingComparts;

    %Disease Progression Comparts
    UninfectedComparts = Params.UninfectedComparts;
    AcuteComparts = Params.AcuteComparts;
    LatentAComparts = Params.LatentAComparts;
    LatentBComparts = Params.LatentBComparts;
    LateComparts = Params.LateComparts;
    AIDSComparts = Params.AIDSComparts;

    %Absorbing State Comparts
    NormDeathCompart = Params.Compart_NormDeath;
    AIDSDeathCompart = Params.Compart_AIDSDeath;
    
    % Create indicator to apply costs only to interventions
    % targeted by allocations (all subpops or YMSM only)
    if Params.TargetIntnsToYMSM == 1
                       
       AllocPopsTargetedIndicator = Params.YoungMSMIndicator;
                       
    else
                       
       AllocPopsTargetedIndicator = ones(Params.numStrats,1);
                       
    end

%% 4. Collect Model Settings and Preallocate Variables
    
    if t == 1 % first timestep or first year of outcomes collection
        
        % Collect model settings that don't change over time (%Added time % periods 3-5 - Clinkscales 5/13/2022)
        Results.set_ModelStartYear = Params.tt_modelStartYear;
        Results.set_PeriodTwoStartYear = Params.tt_periodTwoStartYear;
        Results.set_PeriodThreeStartYear = Params.tt_periodThreeStartYear;
        Results.set_PeriodFourStartYear = Params.tt_periodFourStartYear;
        Results.set_PeriodFiveStartYear = Params.tt_periodFiveStartYear;
        Results.set_TimeHorizon = Params.tt_timeHorizon;
        Results.set_ModelEndYear = Params.tt_modelStartYear + Params.tt_timeHorizon-1;
        Results.set_HETTestFreqIntervalsInMonths_byRisk = Params.tt_HETFreq_Intervals_l;
        Results.set_HETTestFreqAllHET = Params.tt_HETFreq_ApplyLRH;
        Results.set_FirstOutcomeYr = FirstOutcomeYr;
        Results.set_LastOutcomeYr = LastOutcomeYr;
        Results.set_DiagRateYr = Params.diagRateYr;
        Results.set_minCost_Tgt1NumYears = Params.minCost_Tgt1NumYears; 
        Results.set_minCost_Tgt1PctReduce = Params.minCost_Tgt1PctReduce; 
        Results.set_minCost_Tgt2NumYears = Params.minCost_Tgt2NumYears;  
        Results.set_minCost_Tgt2PctReduce = Params.minCost_Tgt2PctReduce; 
        
        switch Params.tt_progressionSource
            case 1
                Results.set_ProgressionSource = 'Direct entry';
            case 2
                Results.set_ProgressionSource = 'HET testing frequency';
            case 3
                Results.set_ProgressionSource = 'Allocation-based';
        end
        switch Params.ModelRunType
            case 1
                Results.set_ModelRunType = 'Epi model';
            case 2
                Results.set_ModelRunType = 'Calibration using LHS';
            case 3
                Results.set_ModelRunType = 'Calibration using optimization';
            case 4
                Results.set_ModelRunType = 'Resource allocation optimization';
        end
        
            Results.set_MinAgeIncluded = Params.MinAgeIncluded;
            Results.set_Circumcision = Params.pop_pctCirc_r;
        
        if LastOutcomeYr >= Results.set_PeriodFiveStartYear
            Results.annCostT1andReachLvlRate_Test_LowRiskHETs_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_Test_HighRiskHETs_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_Test_LowRiskMSM_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_Test_HighRiskMSM_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_Test_IDU_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlProb_LTCatDiag_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;               
            Results.annCostT1andReachLvlRate_LTCafterDiag_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;                                              
            Results.annCostT1andReachLvlRate_ARTInit_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;                
            Results.annCostT1andReachLvlRate_ARTAdhere_BecomeVLS_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;                
            Results.annCostT1andReachLvlRate_ARTAdhere_RemainVLS_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlPctServed_SEP_B_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls)  = 0;
            Results.annCostT1andReachLvlPctServed_SEP_H_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls)  = 0;
            Results.annCostT1andReachLvlPctServed_SEP_O_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls)  = 0;
            Results.annCostT1andReachLvlRate_PrEP_Oral_HETM_B_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Oral_HETM_H_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Oral_HETM_O_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Oral_HETF_B_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Oral_HETF_H_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Oral_HETF_O_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Oral_MSM_B_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Oral_MSM_H_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Oral_MSM_O_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Oral_IDU_B_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Oral_IDU_H_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Oral_IDU_O_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Inject_HETM_B_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Inject_HETM_H_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Inject_HETM_O_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Inject_HETF_B_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Inject_HETF_H_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Inject_HETF_O_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Inject_MSM_B_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Inject_MSM_H_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Inject_MSM_O_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Inject_IDU_B_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Inject_IDU_H_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;
            Results.annCostT1andReachLvlRate_PrEP_Inject_IDU_O_inM(Results.set_ModelEndYear - Results.set_PeriodThreeStartYear + 1,Params.numReachLvls) = 0;            
        end
       % Pre-allocate
        % Results.ann_numPosTestsofHIVPosPeopleCDCFunds=zeros(numYears,numStrats);
        % Results.ann_numPosTestsofHIVPosPeople=zeros(numYears,numStrats);
        % Results.ann_numNegTestsCDCFunds=zeros(numYears,numStrats);
        % Results.ann_numNegTestss=zeros(numYears,numStrats);
        % Results.ann_numNotifiedFromCDCFunds=zeros(numYears,numStrats);
        % Results.ann_numLinkedFirstFromCDCFunds=zeros(numYears,numStrats);
        % Results.ann_numberInCareEligForART=zeros(numYears,numStrats);
        % Results.ann_numAdherentFromCDCFunds=zeros(numYears,numStrats);
        % Results.ann_numARTInitFromCDCFunds=zeros(numYears,numStrats);
        % Results.ann_numPosRapidAcuteTestsCDCFunds=zeros(numYears,numStrats);
        % Results.ann_numPosConvAcuteTestsCDCFunds=zeros(numYears,numStrats);
        % Results.ann_numPosRapidChronicTestsCDCFunds=zeros(numYears,numStrats);
        % Results.ann_numPosConvChronicTestsCDCFunds=zeros(numYears,numStrats);
        % Results.ann_numNegRapidTestsCDCFunds=zeros(numYears,numStrats);
        % Results.ann_numNegConvTestsCDCFunds=zeros(numYears,numStrats);
        % Results.ann_numConvPosNotifiedFromCDCFunds=zeros(numYears,numStrats);
        % Results.ann_numRapidPosNotifiedFromCDCFunds=zeros(numYears,numStrats);

    end

%% 5. Collect Results for Each Time Step

% TRANSITIONS/COUNTS: Infection-related outcomes  
    for sectionInfects = 1:1

    % Infection rates, by force of infection
        
        Results.ts_infRateVonly_NoPrEP(t,:) = InfRate.lambdaVonly(:, odeCounter)';
        Results.ts_infRateAonly_NoPrEP(t,:) = InfRate.lambdaAonly(:, odeCounter)';
        Results.ts_infRateVsomeAI_NoPrEP(t,:) = InfRate.lambdaVsomeAI(:, odeCounter)';
        Results.ts_infRateAsomeVI_NoPrEP(t,:) = InfRate.lambdaAsomeVI(:, odeCounter)';
        Results.ts_infRateN_NoPrEP(t,:) = InfRate.lambdaN(:, odeCounter)';
        
        % JCPrEPUpdate: modified PrEP inf rate outcomes
        Results.ts_infRateVonly_OnPrEP_Oral_High(t,:) = InfRate_PrEP_Oral_High.lambdaVonly(:, odeCounter)';
        Results.ts_infRateAonly_OnPrEP_Oral_High(t,:) = InfRate_PrEP_Oral_High.lambdaAonly(:, odeCounter)';
        Results.ts_infRateVsomeAI_OnPrEP_Oral_High(t,:) = InfRate_PrEP_Oral_High.lambdaVsomeAI(:, odeCounter)';
        Results.ts_infRateAsomeVI_OnPrEP_Oral_High(t,:) = InfRate_PrEP_Oral_High.lambdaAsomeVI(:, odeCounter)';
        Results.ts_infRateN_OnPrEP_Oral_High(t,:) = InfRate_PrEP_Oral_High.lambdaN(:, odeCounter)';
        
        Results.ts_infRateVonly_OnPrEP_Oral_Low(t,:) = InfRate_PrEP_Oral_Low.lambdaVonly(:, odeCounter)';
        Results.ts_infRateAonly_OnPrEP_Oral_Low(t,:) = InfRate_PrEP_Oral_Low.lambdaAonly(:, odeCounter)';
        Results.ts_infRateVsomeAI_OnPrEP_Oral_Low(t,:) = InfRate_PrEP_Oral_Low.lambdaVsomeAI(:, odeCounter)';
        Results.ts_infRateAsomeVI_OnPrEP_Oral_Low(t,:) = InfRate_PrEP_Oral_Low.lambdaAsomeVI(:, odeCounter)';
        Results.ts_infRateN_OnPrEP_Oral_Low(t,:) = InfRate_PrEP_Oral_Low.lambdaN(:, odeCounter)';
        
        Results.ts_infRateVonly_OnPrEP_Inject_High(t,:) = InfRate_PrEP_Inject_High.lambdaVonly(:, odeCounter)';
        Results.ts_infRateAonly_OnPrEP_Inject_High(t,:) = InfRate_PrEP_Inject_High.lambdaAonly(:, odeCounter)';
        Results.ts_infRateVsomeAI_OnPrEP_Inject_High(t,:) = InfRate_PrEP_Inject_High.lambdaVsomeAI(:, odeCounter)';
        Results.ts_infRateAsomeVI_OnPrEP_Inject_High(t,:) = InfRate_PrEP_Inject_High.lambdaAsomeVI(:, odeCounter)';
        Results.ts_infRateN_OnPrEP_Inject_High(t,:) = InfRate_PrEP_Inject_High.lambdaN(:, odeCounter)';
        
        Results.ts_infRateVonly_OnPrEP_Inject_Low(t,:) = InfRate_PrEP_Inject_Low.lambdaVonly(:, odeCounter)';
        Results.ts_infRateAonly_OnPrEP_Inject_Low(t,:) = InfRate_PrEP_Inject_Low.lambdaAonly(:, odeCounter)';
        Results.ts_infRateVsomeAI_OnPrEP_Inject_Low(t,:) = InfRate_PrEP_Inject_Low.lambdaVsomeAI(:, odeCounter)';
        Results.ts_infRateAsomeVI_OnPrEP_Inject_Low(t,:) = InfRate_PrEP_Inject_Low.lambdaAsomeVI(:, odeCounter)';
        Results.ts_infRateN_OnPrEP_Inject_Low(t,:) = InfRate_PrEP_Inject_Low.lambdaN(:, odeCounter)';
        
    % Total new infections, by compartment % JCPrEPUpdate: added infection
    % outcomes from new PrEP states
        Results.ts_NewInfections_A1(t,:) = Transitions(Params.A1,Params.B1,:, odeCounter);
        Results.ts_NewInfections_A6(t,:) = Transitions(Params.A6,Params.B3,:, odeCounter);
        Results.ts_NewInfections_A7(t,:) = Transitions(Params.A7,Params.B3,:, odeCounter);
        Results.ts_NewInfections_A8(t,:) = Transitions(Params.A8,Params.B3,:, odeCounter);
        Results.ts_NewInfections_A9(t,:) = Transitions(Params.A9,Params.B3,:, odeCounter);
        
        Results.ts_NewInfections_OnPrEP(t,:) = Results.ts_NewInfections_A6(t,:) + Results.ts_NewInfections_A7(t,:) + Results.ts_NewInfections_A8(t,:) + Results.ts_NewInfections_A9(t,:);
        
        Results.ts_TotalNewInfections(t,:) = Results.ts_NewInfections_A1(t,:) + Results.ts_NewInfections_A6(t,:) + Results.ts_NewInfections_A7(t,:) + Results.ts_NewInfections_A8(t,:) + Results.ts_NewInfections_A9(t,:);
            
        
    % New infections, by force of infection
    
        % Calculated as:
            % [Total new infections] * [Proportion of infection rate that is due to that force of infection]
        
        % No PrEP
            Results.ts_NewInfectionsVonly_NoPrEP(t,:) = Results.ts_NewInfections_A1(t,:) .* InfRate.RelativeInfRates(:,Params.force_Vonly,odeCounter)';
            Results.ts_NewInfectionsAonly_NoPrEP(t,:) = Results.ts_NewInfections_A1(t,:) .* InfRate.RelativeInfRates(:,Params.force_Aonly,odeCounter)';
            Results.ts_NewInfectionsVsomeAI_NoPrEP(t,:) = Results.ts_NewInfections_A1(t,:) .* InfRate.RelativeInfRates(:,Params.force_VsomeA,odeCounter)';
            Results.ts_NewInfectionsAsomeVI_NoPrEP(t,:) = Results.ts_NewInfections_A1(t,:) .* InfRate.RelativeInfRates(:,Params.force_AsomeV,odeCounter)';
            Results.ts_NewInfectionsN_NoPrEP(t,:) = Results.ts_NewInfections_A1(t,:) .* InfRate.RelativeInfRates(:,Params.force_N,odeCounter)';

        % Oral PrEP, high adherence
            Results.ts_NewInfectionsVonly_OnPrEP_Oral_High(t,:) = Results.ts_NewInfections_A6(t,:) .* InfRate_PrEP_Oral_High.RelativeInfRates(:,Params.force_Vonly,odeCounter)';
            Results.ts_NewInfectionsAonly_OnPrEP_Oral_High(t,:) = Results.ts_NewInfections_A6(t,:) .* InfRate_PrEP_Oral_High.RelativeInfRates(:,Params.force_Aonly,odeCounter)';
            Results.ts_NewInfectionsVsomeAI_OnPrEP_Oral_High(t,:) = Results.ts_NewInfections_A6(t,:) .* InfRate_PrEP_Oral_High.RelativeInfRates(:,Params.force_VsomeA,odeCounter)';
            Results.ts_NewInfectionsAsomeVI_OnPrEP_Oral_High(t,:) = Results.ts_NewInfections_A6(t,:) .* InfRate_PrEP_Oral_High.RelativeInfRates(:,Params.force_AsomeV,odeCounter)';
            Results.ts_NewInfectionsN_OnPrEP_Oral_High(t,:) = Results.ts_NewInfections_A6(t,:) .* InfRate_PrEP_Oral_High.RelativeInfRates(:,Params.force_N,odeCounter)';
        
        % Oral PrEP, Low adherence
            Results.ts_NewInfectionsVonly_OnPrEP_Oral_Low(t,:) = Results.ts_NewInfections_A7(t,:) .* InfRate_PrEP_Oral_Low.RelativeInfRates(:,Params.force_Vonly,odeCounter)';
            Results.ts_NewInfectionsAonly_OnPrEP_Oral_Low(t,:) = Results.ts_NewInfections_A7(t,:) .* InfRate_PrEP_Oral_Low.RelativeInfRates(:,Params.force_Aonly,odeCounter)';
            Results.ts_NewInfectionsVsomeAI_OnPrEP_Oral_Low(t,:) = Results.ts_NewInfections_A7(t,:) .* InfRate_PrEP_Oral_Low.RelativeInfRates(:,Params.force_VsomeA,odeCounter)';
            Results.ts_NewInfectionsAsomeVI_OnPrEP_Oral_Low(t,:) = Results.ts_NewInfections_A7(t,:) .* InfRate_PrEP_Oral_Low.RelativeInfRates(:,Params.force_AsomeV,odeCounter)';
            Results.ts_NewInfectionsN_OnPrEP_Oral_Low(t,:) = Results.ts_NewInfections_A7(t,:) .* InfRate_PrEP_Oral_Low.RelativeInfRates(:,Params.force_N,odeCounter)';    
            
        % Injectable PrEP, high adherence
            Results.ts_NewInfectionsVonly_OnPrEP_Inject_High(t,:) = Results.ts_NewInfections_A8(t,:) .* InfRate_PrEP_Inject_High.RelativeInfRates(:,Params.force_Vonly,odeCounter)';
            Results.ts_NewInfectionsAonly_OnPrEP_Inject_High(t,:) = Results.ts_NewInfections_A8(t,:) .* InfRate_PrEP_Inject_High.RelativeInfRates(:,Params.force_Aonly,odeCounter)';
            Results.ts_NewInfectionsVsomeAI_OnPrEP_Inject_High(t,:) = Results.ts_NewInfections_A8(t,:) .* InfRate_PrEP_Inject_High.RelativeInfRates(:,Params.force_VsomeA,odeCounter)';
            Results.ts_NewInfectionsAsomeVI_OnPrEP_Inject_High(t,:) = Results.ts_NewInfections_A8(t,:) .* InfRate_PrEP_Inject_High.RelativeInfRates(:,Params.force_AsomeV,odeCounter)';
            Results.ts_NewInfectionsN_OnPrEP_Inject_High(t,:) = Results.ts_NewInfections_A8(t,:) .* InfRate_PrEP_Inject_High.RelativeInfRates(:,Params.force_N,odeCounter)';
        
        % Injectable PrEP, Low adherence
            Results.ts_NewInfectionsVonly_OnPrEP_Inject_Low(t,:) = Results.ts_NewInfections_A9(t,:) .* InfRate_PrEP_Inject_Low.RelativeInfRates(:,Params.force_Vonly,odeCounter)';
            Results.ts_NewInfectionsAonly_OnPrEP_Inject_Low(t,:) = Results.ts_NewInfections_A9(t,:) .* InfRate_PrEP_Inject_Low.RelativeInfRates(:,Params.force_Aonly,odeCounter)';
            Results.ts_NewInfectionsVsomeAI_OnPrEP_Inject_Low(t,:) = Results.ts_NewInfections_A9(t,:) .* InfRate_PrEP_Inject_Low.RelativeInfRates(:,Params.force_VsomeA,odeCounter)';
            Results.ts_NewInfectionsAsomeVI_OnPrEP_Inject_Low(t,:) = Results.ts_NewInfections_A9(t,:) .* InfRate_PrEP_Inject_Low.RelativeInfRates(:,Params.force_AsomeV,odeCounter)';
            Results.ts_NewInfectionsN_OnPrEP_Inject_Low(t,:) = Results.ts_NewInfections_A9(t,:) .* InfRate_PrEP_Inject_Low.RelativeInfRates(:,Params.force_N,odeCounter)';    
            
        % Overall
            Results.ts_NewInfectionsVonly(t,:) = Results.ts_NewInfectionsVonly_NoPrEP(t,:) + Results.ts_NewInfectionsVonly_OnPrEP_Oral(t,:) + Results.ts_NewInfectionsVonly_OnPrEP_Inject(t,:);
            Results.ts_NewInfectionsAonly(t,:) = Results.ts_NewInfectionsAonly_NoPrEP(t,:) + Results.ts_NewInfectionsAonly_OnPrEP_Oral(t,:) + Results.ts_NewInfectionsAonly_OnPrEP_Inject(t,:);
            Results.ts_NewInfectionsVsomeAI(t,:) = Results.ts_NewInfectionsVsomeAI_NoPrEP(t,:) + Results.ts_NewInfectionsVsomeAI_OnPrEP_Oral(t,:) + Results.ts_NewInfectionsVsomeAI_OnPrEP_Inject(t,:);
            Results.ts_NewInfectionsAsomeVI(t,:) = Results.ts_NewInfectionsAsomeVI_NoPrEP(t,:) + Results.ts_NewInfectionsAsomeVI_OnPrEP_Oral(t,:) + Results.ts_NewInfectionsAsomeVI_OnPrEP_Inject(t,:);
            Results.ts_NewInfectionsN(t,:) = Results.ts_NewInfectionsN_NoPrEP(t,:) + Results.ts_NewInfectionsN_OnPrEP_Oral(t,:) + Results.ts_NewInfectionsN_OnPrEP_Inject(t,:);
            
    % New infections, by source of infection (only collect if selected to collect infections by source on 'Model Settings' sheet)     
           
            if Params.CollectInfbySource == 1
                Results.ts_NewInfectionsbySource(:,:,t) = InfRate.pctInfSourceBySubpop(:, :, odeCounter) .* Results.ts_TotalNewInfections(t,:)';
            end            
    end
        
% TRANSITIONS/COUNTS: Test and treat progression outcomes
    for sectionTTprog = 1:1
        
    % New Diagnoses, not LTC
        Results.ts_NewDiagnosesNotLTC_B1(t,:)= Transitions(Params.B1,Params.B2,:, odeCounter);
        Results.ts_NewDiagnosesNotLTC_C1(t,:)= Transitions(Params.C1,Params.C2,:, odeCounter);
        Results.ts_NewDiagnosesNotLTC_D1(t,:)= Transitions(Params.D1,Params.D2,:, odeCounter);
        Results.ts_NewDiagnosesNotLTC_E1(t,:)= Transitions(Params.E1,Params.E2,:, odeCounter);
        Results.ts_NewDiagnosesNotLTC_F1(t,:)= Transitions(Params.F1,Params.F2,:, odeCounter);


    % Linked to care first
        Results.ts_LinkToCareFirst_B1(t,:) = Transitions(Params.B1,Params.B3,:, odeCounter);
        Results.ts_LinkToCareFirst_C1(t,:) = Transitions(Params.C1,Params.C3,:, odeCounter);
        Results.ts_LinkToCareFirst_D1(t,:) = Transitions(Params.D1,Params.D3,:, odeCounter);
        Results.ts_LinkToCareFirst_E1(t,:) = Transitions(Params.E1,Params.E3,:, odeCounter);
        Results.ts_LinkToCareFirst_F1(t,:) = Transitions(Params.F1,Params.F3,:, odeCounter);
        

    % Linked to care (no ART effects) after diagnosis (vs. at diagnosis)
        Results.ts_LinkToCare_B2(t,:) = Transitions(Params.B2,Params.B3,:, odeCounter);
        Results.ts_LinkToCare_C2(t,:) = Transitions(Params.C2,Params.C3,:, odeCounter);
        Results.ts_LinkToCare_D2(t,:) = Transitions(Params.D2,Params.D3,:, odeCounter);
        Results.ts_LinkToCare_E2(t,:) = Transitions(Params.E2,Params.E3,:, odeCounter);
        Results.ts_LinkToCare_F2(t,:) = Transitions(Params.F2,Params.F3,:, odeCounter);

        
    % Drop out
    
        % [Care] -> [Aware]
            Results.ts_DropOutOfInCare_B3(t,:) = Transitions(Params.B3,Params.B2,:, odeCounter);
            Results.ts_DropOutOfInCare_C3(t,:) = Transitions(Params.C3,Params.C2,:, odeCounter);
            Results.ts_DropOutOfInCare_D3(t,:) = Transitions(Params.D3,Params.D2,:, odeCounter);
            Results.ts_DropOutOfInCare_E3(t,:) = Transitions(Params.E3,Params.E2,:, odeCounter);
            Results.ts_DropOutOfInCare_F3(t,:) = Transitions(Params.F3,Params.F2,:, odeCounter);
        
        % [ANV] -> [Aware]
            Results.ts_DropOff_ANVToAware_C4(t,:) = Transitions(Params.C4,Params.C2,:, odeCounter);
            Results.ts_DropOff_ANVToAware_D4(t,:) = Transitions(Params.D4,Params.D2,:, odeCounter);
            Results.ts_DropOff_ANVToAware_E4(t,:) = Transitions(Params.E4,Params.E2,:, odeCounter);   
            Results.ts_DropOff_ANVToAware_F4(t,:) = Transitions(Params.F4,Params.F2,:, odeCounter); 

        % [ANV] -> [Care]
            Results.ts_DropOff_ANVToCare_C4(t,:) = Transitions(Params.C4,Params.C3,:, odeCounter);
            Results.ts_DropOff_ANVToCare_D4(t,:) = Transitions(Params.D4,Params.D3,:, odeCounter);
            Results.ts_DropOff_ANVToCare_E4(t,:) = Transitions(Params.E4,Params.E3,:, odeCounter);   
            Results.ts_DropOff_ANVToCare_F4(t,:) = Transitions(Params.F4,Params.F3,:, odeCounter);    
        
        % [VLS] -> [ANV]
            Results.ts_DropOff_VLSToANV_C5(t,:) = Transitions(Params.C5,Params.C4,:, odeCounter);
            Results.ts_DropOff_VLSToANV_D5(t,:) = Transitions(Params.D5,Params.D4,:, odeCounter);
            Results.ts_DropOff_VLSToANV_E5(t,:) = Transitions(Params.E5,Params.E4,:, odeCounter);   
            Results.ts_DropOff_VLSToANV_F5(t,:) = Transitions(Params.F5,Params.F4,:, odeCounter);    
        
        % [PrEP] -> [Sus-no-PrEP]
            Results.ts_DropOff_PrEP_A6(t,:) = Transitions(Params.A6,Params.A1,:, odeCounter);
            Results.ts_DropOff_PrEP_A7(t,:) = Transitions(Params.A7,Params.A1,:, odeCounter);
            Results.ts_DropOff_PrEP_A8(t,:) = Transitions(Params.A8,Params.A1,:, odeCounter);
            Results.ts_DropOff_PrEP_A9(t,:) = Transitions(Params.A9,Params.A1,:, odeCounter);
            
            Results.ts_DropOff_PrEP(t,:) = Results.ts_DropOff_PrEP_A6(t,:) + Results.ts_DropOff_PrEP_A7(t,:) + Results.ts_DropOff_PrEP_A8(t,:) + Results.ts_DropOff_PrEP_A9(t,:);
            Results.ts_DropOff_PrEP_Oral(t,:) = Results.ts_DropOff_PrEP_A6(t,:) + Results.ts_DropOff_PrEP_A7(t,:);
            Results.ts_DropOff_PrEP_Inject(t,:) = Results.ts_DropOff_PrEP_A8(t,:) + Results.ts_DropOff_PrEP_A9(t,:);
            
    % Initiate ART
    
        % Don't become VLS
    
            % Note: Acute PLWH who initiate ART move to LatentA/ART
        Results.ts_StartARTNotVLS_B3(t,:) = Transitions(Params.B3,Params.C4,:, odeCounter);
        Results.ts_StartARTNotVLS_C3(t,:) = Transitions(Params.C3,Params.C4,:, odeCounter);
        Results.ts_StartARTNotVLS_D3(t,:) = Transitions(Params.D3,Params.D4,:, odeCounter);
        Results.ts_StartARTNotVLS_E3(t,:) = Transitions(Params.E3,Params.E4,:, odeCounter);
        Results.ts_StartARTNotVLS_F3(t,:) = Transitions(Params.F3,Params.F4,:, odeCounter);
        
        
        % Become VLS
        
            % From LTC, no ART effects
            Results.ts_BecomeVLS_B3(t,:) = Transitions(Params.B3,Params.C5,:, odeCounter);
            Results.ts_BecomeVLS_C3(t,:) = Transitions(Params.C3,Params.C5,:, odeCounter);
            Results.ts_BecomeVLS_D3(t,:) = Transitions(Params.D3,Params.D5,:, odeCounter);
            Results.ts_BecomeVLS_E3(t,:) = Transitions(Params.E3,Params.E5,:, odeCounter);
            Results.ts_BecomeVLS_F3(t,:) = Transitions(Params.F3,Params.F5,:, odeCounter);
        
            % From ANV (not part of core model, only included for race/eth
            % analysis)
            Results.ts_BecomeVLS_C4(t,:) = Transitions(Params.C4,Params.C5,:, odeCounter);
            Results.ts_BecomeVLS_D4(t,:) = Transitions(Params.D4,Params.D5,:, odeCounter);
            Results.ts_BecomeVLS_E4(t,:) = Transitions(Params.E4,Params.E5,:, odeCounter);
            Results.ts_BecomeVLS_F4(t,:) = Transitions(Params.F4,Params.F5,:, odeCounter);
   
    % Initiate PrEP
        % Note: even though the initiate PrEP transition is defined under the
        % field TransRates.StartART, I'm not considering it a "StartART"
        % transition for the sake of outcomes
        % JCPrEPUpdate: Incorporate new transitions to PrEP into total
        % starting PrEP
            Results.ts_StartPrEP_A6(t,:) = Transitions(Params.A1,Params.A6,:, odeCounter);
            Results.ts_StartPrEP_A7(t,:) = Transitions(Params.A1,Params.A7,:, odeCounter);
            Results.ts_StartPrEP_A8(t,:) = Transitions(Params.A1,Params.A8,:, odeCounter);
            Results.ts_StartPrEP_A9(t,:) = Transitions(Params.A1,Params.A9,:, odeCounter);
            
            Results.ts_StartPrEP(t,:) = Results.ts_StartPrEP_A6(t,:) + Results.ts_StartPrEP_A7(t,:) + Results.ts_StartPrEP_A8(t,:) + Results.ts_StartPrEP_A9(t,:);
            Results.ts_StartPrEP_Oral(t,:) = Results.ts_StartPrEP_A6(t,:) + Results.ts_StartPrEP_A7(t,:);
            Results.ts_StartPrEP_Inject(t,:) = Results.ts_StartPrEP_A8(t,:) + Results.ts_StartPrEP_A9(t,:);
    end
        
% TRANSITIONS/COUNTS: Deaths and aged out
    for sectionDeaths = 1:1
        
    % Normal death
        % Eligible: All modeled individuals except AIDS and affected by ART
        
        
        Results.ts_numNormDeaths_A1(t,:) = Transitions(Params.A1,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_A6(t,:) = Transitions(Params.A6,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_A7(t,:) = Transitions(Params.A7,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_A8(t,:) = Transitions(Params.A8,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_A9(t,:) = Transitions(Params.A9,Params.Compart_NormDeath,:, odeCounter);

        Results.ts_numNormDeaths_B1(t,:) = Transitions(Params.B1,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_C1(t,:) = Transitions(Params.C1,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_D1(t,:) = Transitions(Params.D1,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_E1(t,:) = Transitions(Params.E1,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_F1(t,:) = Transitions(Params.F1,Params.Compart_NormDeath,:, odeCounter);

        Results.ts_numNormDeaths_B2(t,:) = Transitions(Params.B2,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_C2(t,:) = Transitions(Params.C2,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_D2(t,:) = Transitions(Params.D2,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_E2(t,:) = Transitions(Params.E2,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_F2(t,:) = Transitions(Params.F2,Params.Compart_NormDeath,:, odeCounter);

        Results.ts_numNormDeaths_B3(t,:) = Transitions(Params.B3,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_C3(t,:) = Transitions(Params.C3,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_D3(t,:) = Transitions(Params.D3,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_E3(t,:) = Transitions(Params.E3,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_F3(t,:) = Transitions(Params.F3,Params.Compart_NormDeath,:, odeCounter);

        % Note: no Acute ANV state
        Results.ts_numNormDeaths_C4(t,:) = Transitions(Params.C4,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_D4(t,:) = Transitions(Params.D4,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_E4(t,:) = Transitions(Params.E4,Params.Compart_NormDeath,:, odeCounter);
         % Note: AIDS/ANV only die of AIDS
          
        % Note: no Acute VLS state
        Results.ts_numNormDeaths_C5(t,:) = Transitions(Params.C5,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_D5(t,:) = Transitions(Params.D5,Params.Compart_NormDeath,:, odeCounter);
        Results.ts_numNormDeaths_E5(t,:) = Transitions(Params.E5,Params.Compart_NormDeath,:, odeCounter);
        % Note: AIDS/VLS only die of AIDS
            
        
    % AIDS death
        % Eligible: all PLWH with AIDS

        Results.ts_numAIDSDeaths_F1(t,:) = Transitions(Params.F1,Params.Compart_AIDSDeath,:, odeCounter);
        Results.ts_numAIDSDeaths_F2(t,:) = Transitions(Params.F2,Params.Compart_AIDSDeath,:, odeCounter);
        Results.ts_numAIDSDeaths_F3(t,:) = Transitions(Params.F3,Params.Compart_AIDSDeath,:, odeCounter);
        Results.ts_numAIDSDeaths_F4(t,:) = Transitions(Params.F4,Params.Compart_AIDSDeath,:, odeCounter);
        Results.ts_numAIDSDeaths_F5(t,:) = Transitions(Params.F5,Params.Compart_AIDSDeath,:, odeCounter);
     
      
    end
        
% STATE VARIABLES: Population sizes
    for sectionPop = 1:1
        
    % Total Population
    Results.ts_PopulationSize(t,:) = sum(Compartments(Comparts_SusAndPLWH,:,odeCounter))'; 
    
    % All Compartments
    Results.store_Compartments(:,:,t) = Compartments(:,:,odeCounter);
    
    % HIV Prevalence
    Results.ts_HIVPrevalence(t,:)=sum(Compartments(HIVComparts,:,odeCounter));
    
    end
    
% STATE VARIABLES: Distibution along distinct continuum stages
    % E.g., Unaware/Aware/LTC/ART/VLS
    for sectionContDist = 1:1
        
    % On and eligible for PrEP %JCPrEPUpdate: updated num on PrEP calc to
    % include all PrEP states
        Results.ts_NumberOnPrEP(t,:) = sum(Compartments(Params.PrEPComparts,:,odeCounter),1);
        Results.ts_NumberOnPrEP_Oral_HighAdh(t,:) = Compartments(Params.A6,:,odeCounter);
        Results.ts_NumberOnPrEP_Oral_LowAdh(t,:) = Compartments(Params.A7,:,odeCounter);
        Results.ts_NumberOnPrEP_Oral(t,:) = sum(Compartments([Params.A6,Params.A7],:,odeCounter),1);
        Results.ts_NumberOnPrEP_Inject_HighAdh(t,:) = Compartments(Params.A8,:,odeCounter);
        Results.ts_NumberOnPrEP_Inject_LowAdh(t,:) = Compartments(Params.A9,:,odeCounter);
        Results.ts_NumberOnPrEP_Inject(t,:) = sum(Compartments([Params.A8,Params.A9],:,odeCounter),1);

        Results.ts_NumberEligForPrEP(t,:) = sum(Compartments(Params.UninfectedComparts,:,odeCounter),1)'.*Params.tt_prepPctEligible;
        
         for counter_strats = 1:Params.numStrats 
            Results.ts_NumberEligForPrEPNotOnPrEP(t,counter_strats) = ...
                max((Results.ts_NumberEligForPrEP(t,counter_strats) - sum(Compartments(Params.PrEPComparts,counter_strats,odeCounter))'),0);
         end
         
         % Time on PrEP (for calculating annual/cumulative person-years on PrEP [added by JC on 2022/03/01])
        Results.ts_TimeOnPrEP(t,:) = sum(Compartments(Params.PrEPComparts,:,odeCounter),1) * tstep;
        Results.ts_TimeOnPrEP_Oral(t,:) = sum(Compartments([Params.A6,Params.A7],:,odeCounter),1) * tstep;
        Results.ts_TimeOnPrEP_Inject(t,:) = sum(Compartments([Params.A8,Params.A9],:,odeCounter),1) * tstep;

    % PLWH Unaware
        Results.ts_NumberUnaware(t,:) = sum(Compartments(Comparts_CCStage1,:,odeCounter)); 

    % Aware
        % Everyone in model diagram rows 2-5
        Results.ts_NumberAware(t,:) = sum(Compartments(Comparts_CCStages2345,:,odeCounter));

    % In Care
        % All CC stage 3 (In Care)
        % In-care portion of CC stage 4 (ANV)
        % All CC stage 5 (VLS)   

        Results.ts_NumberInCare(t,:) = ...  
              sum(Compartments(Comparts_CCStage3,:,odeCounter)) ...
            + sum(Compartments(Comparts_CCStage4,:,odeCounter)) * PctANVWhoAreInCareInclART ... 
            + sum(Compartments(Comparts_CCStage5,:,odeCounter));

    % ART
        % ART portion of CC stage 4 (ANV)
        % All CC stage 5 (VLS)

        % Doesn't include PrEP
        
        Results.ts_NumberOnART(t,:) = ...
              sum(Compartments(Comparts_CCStage4,:,odeCounter)) * PctANVWhoAreInCareInclART * PctANVInCareWhoAreOnART ... 
            + sum(Compartments(Comparts_CCStage5,:,odeCounter));


    % VLS
        Results.ts_NumberVLS(t,:)=sum(Compartments(Comparts_CCStage5,:,odeCounter));
    end
                    
% STATE VARIABLES: Distribution along model's continuum (r)
    % E.g., Stage 1, stage 2
   for sectionComp = 1:1         
    % Diagnosed but not in care, or affected by ART
        Results.ts_NumberCCStage2(t,:)=sum(Compartments(Comparts_CCStage2,:,odeCounter));
        Results.ts_NumberAIDSUndiag(t,:)=(Compartments(Params.F1,:,odeCounter));
        Results.ts_NumberAIDSCCStage2(t,:)=(Compartments(Params.F2,:,odeCounter));
        Results.ts_NumberAIDSCCStage3(t,:)=(Compartments(Params.F3,:,odeCounter));
        Results.ts_NumberAIDSCCStage4(t,:)=(Compartments(Params.F4,:,odeCounter));
        Results.ts_NumberAIDSVLS(t,:)=(Compartments(Params.F5,:,odeCounter));
        
        Results.ts_NumberAIDSAware(t,:)=(Compartments(Params.F2,:,odeCounter))...
                                       +(Compartments(Params.F3,:,odeCounter))...
                                       +(Compartments(Params.F4,:,odeCounter))...
                                       +(Compartments(Params.F5,:,odeCounter));

    % LTC, not affected  by ART          
        Results.ts_NumberCCStage3(t,:)=sum(Compartments(Comparts_CCStage3,:,odeCounter));
            
    % ANV
        Results.ts_NumberCCStage4(t,:)=sum(Compartments(Comparts_CCStage4,:,odeCounter));
            
    % ANV, on ART
        % Calculated as:
            % [Number in ANV stage]*[Percent ANV in care]*[Percent in-care ANV who are on ART]
        Results.ts_NumberCCStage4_OnART(t,:) = ...
            sum(Compartments(Comparts_CCStage4,:,odeCounter)) * PctANVWhoAreInCareInclART * PctANVInCareWhoAreOnART; 
            
    % ANV, in Care (incl ART)
         % Calculated as:
            % [Number in ANV stage]*[Percent ANV in care]
        Results.ts_NumberCCStage4_InCareInclART(t,:) = ...
            sum(Compartments(Comparts_CCStage4,:,odeCounter)) * PctANVWhoAreInCareInclART; 
                               
    % VLS
        % NumberCCStage5 = NumberVLS
        % Note: this has been defined earlier
        
    % Stages 3-4-5
        % In Care or affected by ART
             
         Results.ts_NumberCCStages345(t,:)=sum(Compartments(Comparts_CCStages345,:,odeCounter)); 
            
    % ANV + VLS 
         % Affected by ART or VLS
            Results.ts_NumberCCStages45(t,:)=sum(Compartments(Comparts_CCStages45,:,odeCounter));
   end
   
% STATE VARIABLES: Distribution along disease stages (h)
   for sectionHIV = 1:1
        % Uninfected
        Results.ts_NumberUninfected(t,:) = sum(Compartments(UninfectedComparts,:,odeCounter));

        % Acute
        Results.ts_NumberAcute(t,:) = sum(Compartments(AcuteComparts,:,odeCounter));

        % LatentA (CD4>500)
        Results.ts_NumberLatentA(t,:) = sum(Compartments(LatentAComparts,:,odeCounter));

        % LatentB (500>CD4>350)
        Results.ts_NumberLatentB(t,:) = sum(Compartments(LatentBComparts,:,odeCounter));

        % Late (350>CD4>200)
        Results.ts_NumberLate(t,:) = sum(Compartments(LateComparts,:,odeCounter));

        % AIDS (CD4<200)
        Results.ts_NumberAIDS(t,:) = sum(Compartments(AIDSComparts,:,odeCounter));
   end
 
% STATE VARIABLE: PrEP
    
    % Coverage
        Results.ts_PctSusOnPrEP = Results.ts_NumberOnPrEP ./ Results.ts_NumberUninfected;
        
        
% STATE VARIABLES: Absorbing states
    for sectionAbsorb = 1:1
        
        % Cumulative Normal Deaths
        Results.ts_CumulativeNumNormDeaths(t,:)=Compartments(NormDeathCompart,:, odeCounter);

        % Cumulative AIDS Deaths
        Results.ts_CumulativeNumAIDSDeaths(t,:)=Compartments(AIDSDeathCompart,:, odeCounter);

    end

% TRANSITIONS/COUNTS: Progression Counts
    for sectionTransitions = 1:1
        
    % ET note: These are all applied in CollectResults as non-CDC funded tests.
    % Be careful once you get into the allocation code.
        
    % Acute
        Results.ts_numPosRapidAcuteTests(t,:) = TTProg.numTSPosRapidAcuteTests';
        Results.ts_numPosConvAcuteTests(t,:) = TTProg.numTSPosConvAcuteTests';
    % Chronic
        Results.ts_numPosRapidChronicTests(t,:)= TTProg.numTSPosRapidChronicTests';
        Results.ts_numPosConvChronicTests(t,:) = TTProg.numTSPosConvChronicTests';
    % Negative
        Results.ts_numNegRapidTests(t,:) = TTProg.numTSNegRapidTests';
        Results.ts_numNegConvTests(t,:)= TTProg.numTSNegConvTests';
    % Number of Notifications
        Results.ts_numConvPosNotify(t,:) = TTProg.numTSConvPosNotify';
        Results.ts_numConvNegNotify(t,:) = TTProg.numTSConvNegNotify';
        Results.ts_numRapidPosNotify(t,:) = TTProg.numTSRapidPosNotify';
        Results.ts_numRapidNegNotify(t,:) = TTProg.numTSRapidNegNotify';
        
        
%     % Calculate number of non-CDC transiitons
%         Results.ts_numLinkedFirstFromFunds(t,:) = TTProg.numTSLinkedFirstFromFunds';
%         Results.ts_numAdherentFromFunds(t,:) = TTProg.numTSAdherentFromFunds;
%         Results.ts_numARTInitFromFunds(t,:) = TTProg.numTSARTInitFromFunds;
         
    end
        
% TRANSITIONS/COUNTS: Health state costs
    for sectionHealthStates = 1:1
        
    % PrEP
        Results.ts_healthStateCost_A6(t,:) = Compartments(Params.A6,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.A6);
        Results.ts_healthStateCost_A7(t,:) = Compartments(Params.A7,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.A7);
        Results.ts_healthStateCost_A8(t,:) = Compartments(Params.A8,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.A8);
        Results.ts_healthStateCost_A9(t,:) = Compartments(Params.A9,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.A9);
        
        Results.ts_healthStateCost_PrEP_Oral_Undisc = ...
            Results.ts_healthStateCost_A6 + ...
            Results.ts_healthStateCost_A7;
        
        Results.ts_healthStateCost_PrEP_Inject_Undisc = ...            
            Results.ts_healthStateCost_A8 + ...
            Results.ts_healthStateCost_A9;
        
        Results.ts_healthStateCost_PrEP_Undisc = ...
            Results.ts_healthStateCost_A6 + ...
            Results.ts_healthStateCost_A7 + ...
            Results.ts_healthStateCost_A8 + ...
            Results.ts_healthStateCost_A9;
        
    % Acute    
        Results.ts_healthStateCost_B1(t,:) = Compartments(Params.B1,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.B1);
        Results.ts_healthStateCost_B2(t,:) = Compartments(Params.B2,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.B2);
        Results.ts_healthStateCost_B3(t,:) = Compartments(Params.B3,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.B3);
        
    % LatentA
        Results.ts_healthStateCost_C1(t,:) = Compartments(Params.C1,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.C1);
        Results.ts_healthStateCost_C2(t,:) = Compartments(Params.C2,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.C2);
        Results.ts_healthStateCost_C3(t,:) = Compartments(Params.C3,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.C3);
        Results.ts_healthStateCost_C4(t,:) = Compartments(Params.C4,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.C4);
        Results.ts_healthStateCost_C5(t,:) = Compartments(Params.C5,:,odeCounter)*tstep * ... 
            AnnualHIVCostsByCompartment(Params.C5);

    % LatentB
        Results.ts_healthStateCost_D1(t,:) = Compartments(Params.D1,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.D1);
        Results.ts_healthStateCost_D2(t,:) = Compartments(Params.D2,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.D2);
        Results.ts_healthStateCost_D3(t,:) = Compartments(Params.D3,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.D3);
        Results.ts_healthStateCost_D4(t,:) = Compartments(Params.D4,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.D4);
        Results.ts_healthStateCost_D5(t,:) = Compartments(Params.D5,:,odeCounter)*tstep * ... 
            AnnualHIVCostsByCompartment(Params.D5);

    % Late
        Results.ts_healthStateCost_E1(t,:) = Compartments(Params.E1,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.E1);
        Results.ts_healthStateCost_E2(t,:) = Compartments(Params.E2,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.E2);
        Results.ts_healthStateCost_E3(t,:) = Compartments(Params.E3,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.E3);
        Results.ts_healthStateCost_E4(t,:) = Compartments(Params.E4,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.E4);
        Results.ts_healthStateCost_E5(t,:) = Compartments(Params.E5,:,odeCounter)*tstep * ... 
            AnnualHIVCostsByCompartment(Params.E5);

    % AIDS
        Results.ts_healthStateCost_F1(t,:) = Compartments(Params.F1,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.F1);
        Results.ts_healthStateCost_F2(t,:) = Compartments(Params.F2,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.F2);
        Results.ts_healthStateCost_F3(t,:) = Compartments(Params.F3,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.F3);
        Results.ts_healthStateCost_F4(t,:) = Compartments(Params.F4,:,odeCounter)*tstep * ...
            AnnualHIVCostsByCompartment(Params.F4);
        Results.ts_healthStateCost_F5(t,:) = Compartments(Params.F5,:,odeCounter)*tstep * ...  
            AnnualHIVCostsByCompartment(Params.F5);

        Results.ts_TransCost_RemainVLS(t,:) = ...
            max(Results.ts_NumberVLS(t,:) - ...
                (Results.ts_DropOff_VLSToANV_C5(t,:) ...
                + Results.ts_DropOff_VLSToANV_D5(t,:) ...
                + Results.ts_DropOff_VLSToANV_E5(t,:) ...
                + Results.ts_DropOff_VLSToANV_F5(t,:)),0) * ...
                tstep * Params.anntranscostPP_RemainVLS;
        
    % Total undiscounted
        Results.ts_TotalARTAndCareCost_Undisc = ...
            Results.ts_healthStateCost_B1 + ...
            Results.ts_healthStateCost_B2 + ...
            Results.ts_healthStateCost_B3 + ...
            Results.ts_healthStateCost_C1 + ...
            Results.ts_healthStateCost_C2 + ...
            Results.ts_healthStateCost_C3 + ...
            Results.ts_healthStateCost_C4 + ...
            Results.ts_healthStateCost_C5 + ...
            Results.ts_healthStateCost_D1 + ...
            Results.ts_healthStateCost_D2 + ...
            Results.ts_healthStateCost_D3 + ...
            Results.ts_healthStateCost_D4 + ...
            Results.ts_healthStateCost_D5 + ...
            Results.ts_healthStateCost_E1 + ...
            Results.ts_healthStateCost_E2 + ...
            Results.ts_healthStateCost_E3 + ...
            Results.ts_healthStateCost_E4 + ...
            Results.ts_healthStateCost_E5 + ...
            Results.ts_healthStateCost_F1 + ...
            Results.ts_healthStateCost_F2 + ...
            Results.ts_healthStateCost_F3 + ...
            Results.ts_healthStateCost_F4 + ...
            Results.ts_healthStateCost_F5;

        Results.ts_ARTCarePrEPCost_Undisc = ...
            Results.ts_healthStateCost_PrEP_Undisc + ...
            Results.ts_TotalARTAndCareCost_Undisc;
                                    
    % Note: total health state cost is calculated during the last time step 
    % rather than every time step
    
    end
   
% TRANSITIONS/COUNTS: Life years
    for sectionLY = 1:1
    % Calculated as:
        % [Number in each HIV stage] * [length of time step]

        % HIV-
        Results.ts_LifeYears_A(t,:) = sum(Compartments(UninfectedComparts,:,odeCounter))*tstep;
        % Acute
        Results.ts_LifeYears_B(t,:) = sum(Compartments(AcuteComparts,:,odeCounter))*tstep;
        % LatentA
        Results.ts_LifeYears_C(t,:) = sum(Compartments(LatentAComparts,:,odeCounter))*tstep; 
        % LatentB
        Results.ts_LifeYears_D(t,:) = sum(Compartments(LatentBComparts,:,odeCounter))*tstep; 
        % Late
        Results.ts_LifeYears_E(t,:) = sum(Compartments(LateComparts,:,odeCounter))*tstep; 
        % AIDS
        Results.ts_LifeYears_F(t,:) = sum(Compartments(AIDSComparts,:,odeCounter))*tstep; 

        % Note: QALYs not calculated every tstep - only calc'd at the end
    end

%% 6. Collect Results at the Beginning of Each Year
   
    if (Year-floor(Year)) == 0 % beginning of each year
    
        columnIndex = Year - FirstOutcomeYr + 1;
 

           % Other Test and Treat Progression Results
            Results.ann_pctRapid(columnIndex,:)=TTProg.PctTestsRapid';

        
            % Test sensitivites over time (set up to have different values
            % in periods 1/2 and period 3)
            Results.ann_TestSens_Rapid_Acute(columnIndex,1) = TTProg.TestSens_Rapid_Acute;
            Results.ann_TestSens_Conv_Acute(columnIndex,1) = TTProg.TestSens_Conv_Acute;
            Results.ann_TestSens_Rapid_Chronic(columnIndex,1) = TTProg.TestSens_Rapid_Chronic;
            Results.ann_TestSens_Conv_Chronic(columnIndex,1) = TTProg.TestSens_Conv_Chronic;
            Results.ann_TestSens_Confirm_Acute(columnIndex,1) = TTProg.TestSens_Confirm_Acute;
            Results.ann_TestSens_Confirm_Chronic(columnIndex,1) = TTProg.TestSens_Confirm_Chronic;
        
            % Costs per test over time 
                % by type of test and test result
                % Note: NAT test unit cost doesn't vary over time.
                    % However the percent of time it is used varies over
                    % time because the sensitivity of the initial screens
                    % and confirmatory tests vary over time.
                    % This is applied later on in the code. 
                        % (search for costPP_NATTest)
            
            Results.ann_costPP_TestNegRapid(columnIndex,1) = TTProg.costPP_TestNegRapid; 
            Results.ann_costPP_TestPosRapid(columnIndex,1) = TTProg.costPP_TestPosRapid; 
            Results.ann_costPP_TestNegConv(columnIndex,1) = TTProg.costPP_TestNegConv; 
            Results.ann_costPP_TestPosConv(columnIndex,1) = TTProg.costPP_TestPosConv;
            
            % Testing rates
                    TTProg.WeightedTestRateAllEligByCohort;
                    TTProg.WeightedTestRateHIVPosByCohort; 
                    TTProg.RaceSpecificAllEligTestRate;
                    TTProg.RaceSpecificHIVPosTestRate;% added 9Oct
                    TTProg.RiskGroupAllEligTestRate;
                    TTProg.RiskGroupHIVPosTestRate;
                    TTProg.HRHAllEligTestRate;
                    TTProg.HRHHIVPosTestRate;
                    TTProg.LRHAllEligTestRate;
                    TTProg.LRHHIVPosTestRate;
                    TTProg.HRMSMAllEligTestRate;
                    TTProg.LRMSMAllEligTestRate;
                    TTProg.IDUAllEligTestRate;
                    TTProg.AllEligTestRateByRiskGpLevel; % = [LRHTestAllEligRate; HRHTestAllEligRate; LRMSMTestAllEligRate; HRMSMTestAllEligRate; RiskGroupAllEligTestRate(Params.pop_IDU)];
                    TTProg.DiseaseStageTestRate;
                    TTProg.OverallUninfectedTestRate;
             
                    
             % Num served by syringe exchange program (SEP)
             Results.ann_NumServedbySEP_B(columnIndex,1) = TTProg.NumPWIDServedbySEP(1,1);
             Results.ann_NumServedbySEP_H(columnIndex,1) = TTProg.NumPWIDServedbySEP(2,1);
             Results.ann_NumServedbySEP_O(columnIndex,1) = TTProg.NumPWIDServedbySEP(3,1);
             Results.ann_NumServedbySEP(columnIndex,1) = sum(TTProg.NumPWIDServedbySEP);


            
        if Params.tt_progressionSource == 3  % Allocation-based progression

            % Only collect these in the 5th period
            if Year >= Params.tt_periodFiveStartYear

                Results.annCostT1andReachLvlRate_Test_LowRiskHETs_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_Test_LowRiskHETs ./ 1000000;
                Results.annCostT1andReachLvlRate_Test_HighRiskHETs_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_Test_HighRiskHETs ./ 1000000;
                Results.annCostT1andReachLvlRate_Test_LowRiskMSM_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_Test_LowRiskMSM ./ 1000000;
                Results.annCostT1andReachLvlRate_Test_HighRiskMSM_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_Test_HighRiskMSM ./ 1000000;
                Results.annCostT1andReachLvlRate_Test_IDU_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_Test_IDU ./ 1000000;
                Results.annCostT1andReachLvlProb_LTCatDiag_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlProb_LTCatDiag ./ 1000000;               
                Results.annCostT1andReachLvlRate_LTCafterDiag_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_LTCafterDiag ./ 1000000;                                              
                Results.annCostT1andReachLvlRate_ARTInit_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_ARTInit ./ 1000000;                
                Results.annCostT1andReachLvlRate_ARTAdhere_BecomeVLS_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_ARTAdhere_BecomeVLS ./ 1000000;                
                Results.annCostT1andReachLvlRate_ARTAdhere_RemainVLS_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_ARTAdhere_RemainVLS ./ 1000000;
                Results.annCostT1andReachLvlPctServed_SEP_B_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlPctServed_SEP(1,:) ./ 1000000;
                Results.annCostT1andReachLvlPctServed_SEP_H_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlPctServed_SEP(2,:) ./ 1000000;
                Results.annCostT1andReachLvlPctServed_SEP_O_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlPctServed_SEP(3,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Oral_HETM_B_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(1,1,1,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Oral_HETM_H_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(1,1,2,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Oral_HETM_O_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(1,1,3,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Oral_HETF_B_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(1,2,1,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Oral_HETF_H_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(1,2,2,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Oral_HETF_O_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(1,2,3,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Oral_MSM_B_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(1,3,1,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Oral_MSM_H_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(1,3,2,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Oral_MSM_O_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(1,3,3,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Oral_IDU_B_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(1,4,1,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Oral_IDU_H_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(1,4,2,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Oral_IDU_O_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(1,4,3,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Inject_HETM_B_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(2,1,1,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Inject_HETM_H_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(2,1,2,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Inject_HETM_O_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(2,1,3,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Inject_HETF_B_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(2,2,1,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Inject_HETF_H_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(2,2,2,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Inject_HETF_O_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(2,2,3,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Inject_MSM_B_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(2,3,1,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Inject_MSM_H_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(2,3,2,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Inject_MSM_O_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(2,3,3,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Inject_IDU_B_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(2,4,1,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Inject_IDU_H_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(2,4,2,:) ./ 1000000;
                Results.annCostT1andReachLvlRate_PrEP_Inject_IDU_O_inM(Year - Params.tt_periodFiveStartYear + 1, :)  = TTProg.CostT1andReachLvlRate_PrEP(2,4,3,:) ./ 1000000;
                
                % Allocation-funded progression
                    
                    %<-- update when allocation outcomes are confirmed
%                 Results.ann_numPosRapidAcuteTestsCDCFunds(columnIndex,:) = TTProg.numPosRapidAcuteTestsCDCFunds';
%                 Results.ann_numPosConvAcuteTestsCDCFunds(columnIndex,:) = TTProg.numPosConvAcuteTestsCDCFunds';
%                 Results.ann_numPosRapidChronicTestsCDCFunds(columnIndex,:) = TTProg.numPosRapidChronicTestsCDCFunds';
%                 Results.ann_numPosConvChronicTestsCDCFunds(columnIndex,:) = TTProg.numPosConvChronicTestsCDCFunds';
%                 Results.ann_numNegRapidTestsCDCFunds(columnIndex,:) = TTProg.numNegRapidTestsCDCFunds';
%                 Results.ann_numNegConvTestsCDCFunds(columnIndex,:) = TTProg.numNegConvTestsCDCFunds';
                
                %Results.ann_numPosTestsofHIVPosPeopleCDCFunds(columnIndex,:)=TTProg.numPosTestsofHIVPosPeopleCDCFunds';
                %Results.ann_numNegTestsCDCFunds(columnIndex,:)=TTProg.numNegTestsCDCFunds';
%                 %Results.ann_numNotifiedFromCDCFunds(columnIndex,:)=TTProg.numNotifiedFromCDCFunds';
%                 Results.ann_numConvPosNotifiedFromCDCFunds(columnIndex,:)=TTProg.numConvPosNotifyFromCDC'; % added 12Nov2014 ET 
%                 Results.ann_numRapidPosNotifiedFromCDCFunds(columnIndex,:)=TTProg.numRapidPosNotifyFromCDC'; % added 12Nov2014 ET
%                 Results.ann_numLinkedFirstFromCDCFunds(columnIndex,:)=TTProg.numLinkedFirstFromCDCFunds';
%               
%                 
%                 
%                 Results.ann_numberInCareEligForART(columnIndex,:)=TTProg.numberInCareEligForART';
%                 Results.ann_numAdherentFromCDCFunds(columnIndex,:)=TTProg.numAdherentFromCDCFunds';
%                 Results.ann_numARTInitFromCDCFunds(columnIndex,:)=TTProg.numARTInitFromCDCFunds';
            end
        end
    end

%% 7. Reformat Results: Last Run of the Model

    %Adjust collected outcomes to be in the desired results format
    %Move from local variables to Results struct

    if (Params.SolveCont == 1 && numYears == numel(stepsPerYear) && odeCounter == stepsPerYear(end)) ...
            || (Params.SolveCont ~= 1 && t == numSteps)
     % If continuous and the total years is the number of elements in
            % steps per year and it's the last timestep of that year
            % OR if not continuous and t is equal to the total time steps in the model 

    % 7.i. Calculate totals from time step outcomes
        for sectionTSTotals = 1:1
            clear Results.First_t

            
       % DIAGNOSES
            % Not immediately linked to care + immediately linked to care
            
            % Total diagnoses by disease stage
            Results.ts_NewDiagnoses_B1 = Results.ts_NewDiagnosesNotLTC_B1+Results.ts_LinkToCareFirst_B1;
            Results.ts_NewDiagnoses_C1 = Results.ts_NewDiagnosesNotLTC_C1+Results.ts_LinkToCareFirst_C1;
            Results.ts_NewDiagnoses_D1 = Results.ts_NewDiagnosesNotLTC_D1+Results.ts_LinkToCareFirst_D1;
            Results.ts_NewDiagnoses_E1 = Results.ts_NewDiagnosesNotLTC_E1+Results.ts_LinkToCareFirst_E1;
            Results.ts_NewDiagnoses_F1 = Results.ts_NewDiagnosesNotLTC_F1+Results.ts_LinkToCareFirst_F1;
         
            % Diagnoses not immediately linked to care
            Results.ts_TotalNewDiagnosesNotLTC = ...
                  Results.ts_NewDiagnosesNotLTC_B1 ...
                + Results.ts_NewDiagnosesNotLTC_C1 ...
                + Results.ts_NewDiagnosesNotLTC_D1 ...
                + Results.ts_NewDiagnosesNotLTC_E1 ...
                + Results.ts_NewDiagnosesNotLTC_F1;

            % Diagnoses who are immediately linked to care
            Results.ts_TotalLinkToCareFirst = ...
                  Results.ts_LinkToCareFirst_B1 ...
                + Results.ts_LinkToCareFirst_C1 ...
                + Results.ts_LinkToCareFirst_D1 ...
                + Results.ts_LinkToCareFirst_E1 ...
                + Results.ts_LinkToCareFirst_F1;

            % Overall total new diagnoses
            Results.ts_TotalNewDiagnoses = ...
                  Results.ts_TotalNewDiagnosesNotLTC ...
                + Results.ts_TotalLinkToCareFirst;
            
            
       % LINK TO CARE
    
            % Linked to care (not immediately)
            Results.ts_TotalLinkToCareAfterFirst = ...
                  Results.ts_LinkToCare_B2 ...
                + Results.ts_LinkToCare_C2 ...
                + Results.ts_LinkToCare_D2 ...
                + Results.ts_LinkToCare_E2 ...
                + Results.ts_LinkToCare_F2;

            % Overall linked to care
            Results.ts_TotalLinkToCare = ...
                  Results.ts_TotalLinkToCareFirst ...
                + Results.ts_TotalLinkToCareAfterFirst;

            
      % DROP OUT
        
            % [In care] -> [Aware]
            Results.ts_TotalDropOutofInCare = ...
                  Results.ts_DropOutOfInCare_B3 ...
                + Results.ts_DropOutOfInCare_C3 ...
                + Results.ts_DropOutOfInCare_D3 ...
                + Results.ts_DropOutOfInCare_E3 ...
                + Results.ts_DropOutOfInCare_F3;

            % [ANV] -> [Aware]
            Results.ts_TotalDropOff_ANVToAware = ...
                + Results.ts_DropOff_ANVToAware_C4 ...
                + Results.ts_DropOff_ANVToAware_D4 ...
                + Results.ts_DropOff_ANVToAware_E4 ...
                + Results.ts_DropOff_ANVToAware_F4;
            
            % [ANV] -> [Care]
            Results.ts_TotalDropOff_ANVToCare = ...
                + Results.ts_DropOff_ANVToCare_C4 ...
                + Results.ts_DropOff_ANVToCare_D4 ...
                + Results.ts_DropOff_ANVToCare_E4 ...
                + Results.ts_DropOff_ANVToCare_F4;
            
            % [VLS] -> [ANV]
            Results.ts_TotalDropOff_VLSToANV = ...
                + Results.ts_DropOff_VLSToANV_C5 ...
                + Results.ts_DropOff_VLSToANV_D5 ...
                + Results.ts_DropOff_VLSToANV_E5 ...
                + Results.ts_DropOff_VLSToANV_F5;  
            
            
     % START ART  
        
            Results.ts_TotalStartARTNotVLS = ...
                + Results.ts_StartARTNotVLS_B3 ...
                + Results.ts_StartARTNotVLS_C3 ...
                + Results.ts_StartARTNotVLS_D3 ...
                + Results.ts_StartARTNotVLS_E3 ...
                + Results.ts_StartARTNotVLS_F3;
     
          
     % BECOME VLS
        % From LTC, no ART effects
            Results.ts_BecomeVLSFromLTC = ...
                + Results.ts_BecomeVLS_B3 ...
                + Results.ts_BecomeVLS_C3 ...
                + Results.ts_BecomeVLS_D3 ...
                + Results.ts_BecomeVLS_E3 ...
                + Results.ts_BecomeVLS_F3;
            
        % From ANV
            Results.ts_BecomeVLSFromANV = ...
                + Results.ts_BecomeVLS_C4 ...
                + Results.ts_BecomeVLS_D4 ...
                + Results.ts_BecomeVLS_E4 ...
                + Results.ts_BecomeVLS_F4;
     
     % DEATHS
        
        % Normal deaths for PLWH not on ART, by HIV stage
            Results.ts_numNormDeathsAcuteNoART = ...
                Results.ts_numNormDeaths_B1 + ...
                Results.ts_numNormDeaths_B2 + Results.ts_numNormDeaths_B3;
            
            Results.ts_numNormDeathsLatentANoART = ...
                Results.ts_numNormDeaths_C1 + ...
                Results.ts_numNormDeaths_C2 + Results.ts_numNormDeaths_C3;
             
            Results.ts_numNormDeathsLatentBNoART = ...
                Results.ts_numNormDeaths_D1 + ...
                Results.ts_numNormDeaths_D2 + Results.ts_numNormDeaths_D3;
            
            Results.ts_numNormDeathsLateNoART = ...
                Results.ts_numNormDeaths_E1 + ...
                Results.ts_numNormDeaths_E2 + Results.ts_numNormDeaths_E3;
            
            Results.ts_numNormDeathsAIDSNoART = ...
                Results.ts_numNormDeaths_F1 + ...
                Results.ts_numNormDeaths_F2 + Results.ts_numNormDeaths_F3;
            
            
        % Total normal deaths for PLWH, no ART

            Results.ts_numNormDeathsHIVPosNoART =  ...
                  Results.ts_numNormDeathsAcuteNoART ...
                + Results.ts_numNormDeathsLatentANoART ...
                + Results.ts_numNormDeathsLatentBNoART ...
                + Results.ts_numNormDeathsLateNoART ...
                + Results.ts_numNormDeathsAIDSNoART;


        % Total normal deaths for PLWH on ART or VLS

            Results.ts_numNormDeathsANVandVLS = ...
               + Results.ts_numNormDeaths_C4 ...
               + Results.ts_numNormDeaths_D4 ...
               + Results.ts_numNormDeaths_E4 ...
               + Results.ts_numNormDeaths_C5 ...
               + Results.ts_numNormDeaths_D5 ...
               + Results.ts_numNormDeaths_E5;

            Results.ts_numDeathsVLS = ...
               + Results.ts_numNormDeaths_C5 ...
               + Results.ts_numNormDeaths_D5 ...
               + Results.ts_numNormDeaths_E5...
               + Results.ts_numAIDSDeaths_F5...
               ;

            Results.ts_numAIDSDeathsVLS = ...
               + Results.ts_numAIDSDeaths_F5;

           
            Results.ts_numDeathsCCStage4 = ...
               + Results.ts_numNormDeaths_C4 ...
               + Results.ts_numNormDeaths_D4 ...
               + Results.ts_numNormDeaths_E4...
               + Results.ts_numAIDSDeaths_F4;

            Results.ts_numAIDSDeathsCCStage4 = ...
               + Results.ts_numAIDSDeaths_F4;
           
            Results.ts_numDeathsCCStage2 = ...
               + Results.ts_numNormDeaths_B2 ...
               + Results.ts_numNormDeaths_C2 ...
               + Results.ts_numNormDeaths_D2 ...
               + Results.ts_numNormDeaths_E2...
               + Results.ts_numNormDeaths_F2...
               + Results.ts_numAIDSDeaths_F2;

             Results.ts_numAIDSDeathsCCStage2 = ...
               + Results.ts_numAIDSDeaths_F2;
           
           Results.ts_numDeathsCCStage3 = ...
               + Results.ts_numNormDeaths_B3 ...
               + Results.ts_numNormDeaths_C3 ...
               + Results.ts_numNormDeaths_D3 ...
               + Results.ts_numNormDeaths_E3...
               + Results.ts_numNormDeaths_F3...
               + Results.ts_numAIDSDeaths_F3;

            Results.ts_numAIDSDeathsCCStage3 = ...
               + Results.ts_numAIDSDeaths_F3;

        % Total deaths PLWH aware
            Results.ts_numDeathsPLWHAware = ...
                + Results.ts_numNormDeaths_B2 ...
                + Results.ts_numNormDeaths_C2 ...
                + Results.ts_numNormDeaths_D2 ...
                + Results.ts_numNormDeaths_E2 ...
                + Results.ts_numNormDeaths_F2 ...
                + Results.ts_numNormDeaths_B3 ...
                + Results.ts_numNormDeaths_C3 ...
                + Results.ts_numNormDeaths_D3 ...
                + Results.ts_numNormDeaths_E3 ...
                + Results.ts_numNormDeaths_F3 ...
                + Results.ts_numNormDeaths_C4 ...
                + Results.ts_numNormDeaths_D4 ...
                + Results.ts_numNormDeaths_E4 ...
                + Results.ts_numNormDeaths_C5 ...
                + Results.ts_numNormDeaths_D5 ...
                + Results.ts_numNormDeaths_E5 ...
                + Results.ts_numAIDSDeaths_F2 ...
                + Results.ts_numAIDSDeaths_F3 ...
                + Results.ts_numAIDSDeaths_F4 ...
                + Results.ts_numAIDSDeaths_F5;
           
        % Normal deaths for PLWH Unaware
        
            Results.ts_numNormDeathsHIVPosUndiag = ...  
                  Results.ts_numNormDeaths_B1 ...
                + Results.ts_numNormDeaths_C1 ...
                + Results.ts_numNormDeaths_D1 ...
                + Results.ts_numNormDeaths_E1 ...
                + Results.ts_numNormDeaths_F1;
  
       % All deaths for PLWH Unaware
        
            Results.ts_numDeathsHIVPosUndiag = ...  
                  Results.ts_numNormDeathsHIVPosUndiag ...
                + Results.ts_numAIDSDeaths_F1;
            
            Results.ts_numAIDSDeathsHIVPosUndiag = ...
               + Results.ts_numAIDSDeaths_F1;

        % AIDS Deaths no ART
        
            Results.ts_numAIDSDeathsNoART = ...
                  Results.ts_numAIDSDeaths_F1 ...
                + Results.ts_numAIDSDeaths_F2 ...
                + Results.ts_numAIDSDeaths_F3;
            
            
        % AIDS Deaths on ART or VLS
            Results.ts_numAIDSDeathsANVandVLS = ...
                  Results.ts_numAIDSDeaths_F4 ...
                + Results.ts_numAIDSDeaths_F5;
            

        % AIDS Deaths Diagnosed
            Results.ts_numAIDSDeathsAware = ...
                  Results.ts_numAIDSDeaths_F2 ...
                + Results.ts_numAIDSDeaths_F3 ...
                + Results.ts_numAIDSDeaths_F4 ...
                + Results.ts_numAIDSDeaths_F5;            

        end

    % 7.ii. Calculate QALYs, Life Years, Discounted Infections and
    % discounted treatment costs
        for sectionQALYs = 1:1
            
    % QALYs, LYs, and discounted infections    
    %        DiscRates_QALYs(numSteps,1) = 0;
    %        DiscRates_Costs(numSteps,1) = 0;

            % Set up vector of discount factors
             if Params.SolveCont == 1 % continuous

                 for i = 1:numYears
                    YrsFromDollarYear = (FirstOutcomeYr + i - 1) - Params.DollarYear;
%                     DiscRates_QALYs(i,1)= 1/((1+Params.Disc_QALYs)^(YrsFromDollarYear));
%                     DiscRates_Costs(i,1)= 1/((1+Params.Disc_Costs)^(YrsFromDollarYear));
%                     DiscRates_Infs(i,1)= 1/((1+Params.alloc_discRateNewInfections)^(YrsFromDollarYear));
                    DiscRates_QALYs(stepsPerYear(i),1)=0;
                    DiscRates_LifeYears(stepsPerYear(i),1)=0;
                    DiscRates_Costs(stepsPerYear(i),1)=0;
                    DiscRates_Infs(stepsPerYear(i),1)=0;
                    if i == 1 % first year
                        for j = 1:stepsPerYear(i)
                            DiscRates_QALYs(j,1)= 1/((1+Params.Disc_QALYs)^(YrsFromDollarYear));
                            DiscRates_LifeYears(j,1)=1/((1+Params.Disc_LifeYears)^(YrsFromDollarYear));
                            DiscRates_Costs(j,1)= 1/((1+Params.Disc_Costs)^(YrsFromDollarYear));
                            DiscRates_Infs(j,1)= 1/((1+Params.alloc_discRateNewInfections)^(YrsFromDollarYear));                       
                        end
                    else
                        for j = 1:stepsPerYear(i)
                            DiscRates_QALYs(sum(stepsPerYear(1:i-1))+j,1)= 1/((1+Params.Disc_QALYs)^(YrsFromDollarYear));
                            DiscRates_LifeYears(sum(stepsPerYear(1:i-1))+j,1)= 1/((1+Params.Disc_LifeYears)^(YrsFromDollarYear));
                            DiscRates_Costs(sum(stepsPerYear(1:i-1))+j,1)= 1/((1+Params.Disc_Costs)^(YrsFromDollarYear));
                            DiscRates_Infs(sum(stepsPerYear(1:i-1))+j,1)= 1/((1+Params.alloc_discRateNewInfections)^(YrsFromDollarYear));
                        end
                    end
                 end    

             else % discrete

                for i = 1:numYears
                    YrsFromDollarYear = (FirstOutcomeYr + i - 1) - Params.DollarYear;
                    for j = 1:stepsPerYear
                        DiscRates_QALYs((i-1)*stepsPerYear+j,1)= 1/((1+Params.Disc_QALYs)^(YrsFromDollarYear));
                        DiscRates_LifeYears((i-1)*stepsPerYear+j,1)= 1/((1+Params.Disc_LifeYears)^(YrsFromDollarYear));
                        DiscRates_Costs((i-1)*stepsPerYear+j,1)= 1/((1+Params.Disc_Costs)^(YrsFromDollarYear));
                        DiscRates_Infs((i-1)*stepsPerYear+j,1)= 1/((1+Params.alloc_discRateNewInfections)^(YrsFromDollarYear));
                    end
                end    

             end

                      
                % QALYS by tstep
                    % Calculated using the undiscounted life years because
                    % QALYs may have a different discount factor applied to
                    % them

                % HIV-
                Results.ts_QALYs_A = Results.ts_LifeYears_A.*(DiscRates_QALYs*ones(1,numStrats)) * Params.hiv_utility_A;
                % Acute
                Results.ts_QALYs_B = Results.ts_LifeYears_B.*(DiscRates_QALYs*ones(1,numStrats)) * Params.hiv_utility_B;
                % LatentA
                Results.ts_QALYs_C = Results.ts_LifeYears_C.*(DiscRates_QALYs*ones(1,numStrats)) * Params.hiv_utility_C;
                % LatentB
                Results.ts_QALYs_D = Results.ts_LifeYears_D.*(DiscRates_QALYs*ones(1,numStrats)) * Params.hiv_utility_D;
                % Late
                Results.ts_QALYs_E = Results.ts_LifeYears_E.*(DiscRates_QALYs*ones(1,numStrats)) * Params.hiv_utility_E;
                % AIDS
                Results.ts_QALYs_F = Results.ts_LifeYears_F.*(DiscRates_QALYs*ones(1,numStrats)) * Params.hiv_utility_F;

                % Total QALYs by TStep
                Results.ts_QALYs = Results.ts_QALYs_A + ...
                    Results.ts_QALYs_B + Results.ts_QALYs_C + ...
                    Results.ts_QALYs_D + Results.ts_QALYs_E + ...
                    Results.ts_QALYs_F;
                
                % Life years by tstep
                    % Calculated using a different discount factor than the
                    % QALYs
                
                % HIV-
                Results.ts_LifeYears_A = Results.ts_LifeYears_A.*(DiscRates_LifeYears*ones(1,numStrats));
                % Acute
                Results.ts_LifeYears_B = Results.ts_LifeYears_B.*(DiscRates_LifeYears*ones(1,numStrats));
                % LatentA
                Results.ts_LifeYears_C = Results.ts_LifeYears_C.*(DiscRates_LifeYears*ones(1,numStrats));
                % LatentB
                Results.ts_LifeYears_D = Results.ts_LifeYears_D.*(DiscRates_LifeYears*ones(1,numStrats));
                % Late
                Results.ts_LifeYears_E = Results.ts_LifeYears_E.*(DiscRates_LifeYears*ones(1,numStrats));
                % AIDS
                Results.ts_LifeYears_F = Results.ts_LifeYears_F.*(DiscRates_LifeYears*ones(1,numStrats));
                

                % Total Life Years
                    % 1x273
                Results.total_LifeYears_A = sum(Results.ts_LifeYears_A);
                Results.total_LifeYears_B = sum(Results.ts_LifeYears_B);
                Results.total_LifeYears_C = sum(Results.ts_LifeYears_C);
                Results.total_LifeYears_D = sum(Results.ts_LifeYears_D);
                Results.total_LifeYears_E = sum(Results.ts_LifeYears_E);
                Results.total_LifeYears_F = sum(Results.ts_LifeYears_F);

                Results.total_LifeYears_byCohort = Results.total_LifeYears_A + ...
                    Results.total_LifeYears_B + Results.total_LifeYears_C + ...
                    Results.total_LifeYears_D + Results.total_LifeYears_E + ...
                    Results.total_LifeYears_F;
                
                Results.total_LifeYears = sum(Results.total_LifeYears_byCohort);
                
                
                % Total LifeYears by tstep
                Results.ts_LifeYears = Results.ts_LifeYears_A + ...
                    Results.ts_LifeYears_B + Results.ts_LifeYears_C + ...
                    Results.ts_LifeYears_D + Results.ts_LifeYears_E + ...
                    Results.ts_LifeYears_F;
                
                
                % Total QALYs
                    % [1x273]
                Results.total_QALYs_A = sum(Results.ts_QALYs_A);
                Results.total_QALYs_B = sum(Results.ts_QALYs_B);
                Results.total_QALYs_C = sum(Results.ts_QALYs_C);
                Results.total_QALYs_D = sum(Results.ts_QALYs_D);
                Results.total_QALYs_E = sum(Results.ts_QALYs_E);
                Results.total_QALYs_F = sum(Results.ts_QALYs_F);

                
                % Overall QALYs 
                    % [1x273]
                Results.total_QALYs_byCohort = Results.total_QALYs_A + ...
                    Results.total_QALYs_B + Results.total_QALYs_C + ...
                    Results.total_QALYs_D + Results.total_QALYs_E + ...
                    Results.total_QALYs_F;
                
                
                % Total QALYs [1x1]
                Results.total_QALYs = sum(Results.total_QALYs_byCohort);

                
                % Discounted infections by tstep
                Results.TotalNewInfections_Disc = sum(Results.ts_TotalNewInfections,2)'*DiscRates_Infs;            
                
                % Discounted infections by tstep for YMSM
                Results.TotalNewInfections_YMSM_Disc = sum(Results.ts_TotalNewInfections .* Params.YoungMSMIndicator',2)'*DiscRates_Infs;
                
                % Length of each timestep

                if Params.SolveCont == 1 % if continuous run
                    Results.ts_timeStepLength = TTProg.timePerODEStep;
                    Results.ann_stepsPerYear = TTProg.stepsPerYear;
                else % discrete
                    Results.ts_timeStepLength = Params.tt_tstep;
                end
                
                
            % Health state costs (discounted)
                    Results.ts_healthStateCost_A6 = Results.ts_healthStateCost_A6 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_A7 = Results.ts_healthStateCost_A7 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_A8 = Results.ts_healthStateCost_A8 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_A9 = Results.ts_healthStateCost_A9 .* (DiscRates_Costs * ones(1,numStrats));
                    
                    Results.ts_healthStateCost_PrEP_Oral = Results.ts_healthStateCost_A6 + Results.ts_healthStateCost_A7;
                    Results.ts_healthStateCost_PrEP_Inject = Results.ts_healthStateCost_A8 + Results.ts_healthStateCost_A9;
                    
                    Results.ts_healthStateCost_PrEP = Results.ts_healthStateCost_A6 + Results.ts_healthStateCost_A7 + Results.ts_healthStateCost_A8 + Results.ts_healthStateCost_A9;
            
                    Results.ts_healthStateCost_B1 = Results.ts_healthStateCost_B1 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_B2 = Results.ts_healthStateCost_B2 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_B3 = Results.ts_healthStateCost_B3 .* (DiscRates_Costs * ones(1,numStrats));

                    Results.ts_healthStateCost_C1 = Results.ts_healthStateCost_C1 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_C2 = Results.ts_healthStateCost_C2 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_C3 = Results.ts_healthStateCost_C3 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_C4 = Results.ts_healthStateCost_C4 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_C5 = Results.ts_healthStateCost_C5 .* (DiscRates_Costs * ones(1,numStrats));

                    Results.ts_healthStateCost_D1 = Results.ts_healthStateCost_D1 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_D2 = Results.ts_healthStateCost_D2 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_D3 = Results.ts_healthStateCost_D3 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_D4 = Results.ts_healthStateCost_D4 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_D5 = Results.ts_healthStateCost_D5 .* (DiscRates_Costs * ones(1,numStrats));

                    Results.ts_healthStateCost_E1 = Results.ts_healthStateCost_E1 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_E2 = Results.ts_healthStateCost_E2 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_E3 = Results.ts_healthStateCost_E3 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_E4 = Results.ts_healthStateCost_E4 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_E5 = Results.ts_healthStateCost_E5 .* (DiscRates_Costs * ones(1,numStrats));

                    Results.ts_healthStateCost_F1 = Results.ts_healthStateCost_F1 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_F2 = Results.ts_healthStateCost_F2 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_F3 = Results.ts_healthStateCost_F3 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_F4 = Results.ts_healthStateCost_F4 .* (DiscRates_Costs * ones(1,numStrats));
                    Results.ts_healthStateCost_F5 = Results.ts_healthStateCost_F5 .* (DiscRates_Costs * ones(1,numStrats));  
                    
                    % Collect cost of remaining VLS by timestep
                    Results.ts_TransCost_RemainVLS_Undisc = Results.ts_TransCost_RemainVLS;
                    
                    Results.ts_TransCost_RemainVLS = Results.ts_TransCost_RemainVLS .* (DiscRates_Costs * ones(1,numStrats));
 
           % Overall care and treatment costs
                    Results.ts_TotalARTAndCareCost = ...
                        Results.ts_healthStateCost_B1 + ...
                        Results.ts_healthStateCost_B2 + ...
                        Results.ts_healthStateCost_B3 + ...
                        Results.ts_healthStateCost_C1 + ...
                        Results.ts_healthStateCost_C2 + ...
                        Results.ts_healthStateCost_C3 + ...
                        Results.ts_healthStateCost_C4 + ...
                        Results.ts_healthStateCost_C5 + ...
                        Results.ts_healthStateCost_D1 + ...
                        Results.ts_healthStateCost_D2 + ...
                        Results.ts_healthStateCost_D3 + ...
                        Results.ts_healthStateCost_D4 + ...
                        Results.ts_healthStateCost_D5 + ...
                        Results.ts_healthStateCost_E1 + ...
                        Results.ts_healthStateCost_E2 + ...
                        Results.ts_healthStateCost_E3 + ...
                        Results.ts_healthStateCost_E4 + ...
                        Results.ts_healthStateCost_E5 + ...
                        Results.ts_healthStateCost_F1 + ...
                        Results.ts_healthStateCost_F2 + ...
                        Results.ts_healthStateCost_F3 + ...
                        Results.ts_healthStateCost_F4 + ...
                        Results.ts_healthStateCost_F5;
                    
                    Results.ts_ARTCarePrEPCost = ...
                        Results.ts_healthStateCost_PrEP + ...
                        Results.ts_TotalARTAndCareCost;

        end
         
    % 7.iii. Calculate state variable outcomes annually
    
        %Pre-allocate
        finalVariable(numYears,Params.numStrats) = 0;
        ann_DiscRates_Costs(numYears,1)=0;
        
        for sectionTS = 1:1
            
            % This code takes applies the value of the variable in the last time
            % step of the year as the annual outcome.
                                    
            for yearCount = 1:numYears
            % Calculate discount factor annually
               if Params.SolveCont == 1 % if continuous
                            ann_DiscRates_Costs(yearCount,1) = ...
                                DiscRates_Costs(sum(stepsPerYear(1:yearCount)),1);
                        else % discrete
                            % the discrete version has the same number of steps
                            % each year while the continuous version varies
                            ann_DiscRates_Costs(yearCount,1) = DiscRates_Costs(stepsPerYear*yearCount,1);
               end
        
            end
            
            for u = 2:39
                for yearCount = 1:numYears


        % Set up the initial outputs
                   switch u
                        case 2 
                            initialVariable = Results.ts_NumberCCStages345;
                        case 3 
                            initialVariable = Results.ts_NumberAware;
                        case 4
                            initialVariable = Results.ts_NumberUnaware;
                        case 5
                            initialVariable = Results.ts_HIVPrevalence;
                        case 6
                            initialVariable = Results.ts_NumberCCStage2;
                        case 7
                            initialVariable = Results.ts_NumberInCare;
                        case 8
                            initialVariable = Results.ts_NumberCCStage4;
                        case 9
                            initialVariable = Results.ts_NumberUninfected;
                        case 10
                            initialVariable = Results.ts_NumberAcute;
                        case 11
                            initialVariable = Results.ts_NumberLatentA;
                        case 12
                            initialVariable = Results.ts_NumberLatentB;
                        case 13
                            initialVariable = Results.ts_NumberLate;
                        case 14
                            initialVariable = Results.ts_NumberAIDS;
                        case 15
                            initialVariable = Results.ts_CumulativeNumNormDeaths;
                        case 16
                            initialVariable = Results.ts_CumulativeNumAIDSDeaths;                               
                        case 17
                            initialVariable = Results.ts_PopulationSize;
                        case 18
                            initialVariable = Results.ts_NumberCCStages45;
                        case 19
                            initialVariable = Results.ts_NumberVLS;
                        case 20
                            initialVariable = Results.ts_NumberOnART;
                        case 21
                            initialVariable = Results.ts_NumberCCStage3;
                        case 22
                            initialVariable = Results.ts_NumberCCStage4_OnART;
                        case 23
                            initialVariable = Results.ts_NumberCCStage4_InCareInclART;
                       case 24
                            initialVariable = Results.ts_NumberOnPrEP;
                       case 25
                            initialVariable = Results.ts_PctSusOnPrEP;
                       case 26
                            initialVariable = Results.ts_NumberEligForPrEP;
                       case 27
                            initialVariable = Results.ts_NumberEligForPrEPNotOnPrEP;
                       case 28
                            initialVariable = Results.ts_NumberOnPrEP_Oral;
                       case 29
                            initialVariable = Results.ts_NumberOnPrEP_Inject;
                       case 30
                            initialVariable = Results.ts_NumberOnPrEP_Oral_HighAdh;
                       case 31
                            initialVariable = Results.ts_NumberOnPrEP_Oral_LowAdh;
                       case 32
                            initialVariable = Results.ts_NumberOnPrEP_Inject_HighAdh;
                       case 33
                            initialVariable = Results.ts_NumberOnPrEP_Inject_LowAdh;
                        case 34
                            initialVariable = Results.ts_NumberAIDSUndiag;
                        case 35 
                            initialVariable = Results.ts_NumberAIDSCCStage2;
                        case 36
                            initialVariable = Results.ts_NumberAIDSCCStage3;
                        case 37 
                            initialVariable = Results.ts_NumberAIDSCCStage4;
                        case 38 
                            initialVariable = Results.ts_NumberAIDSVLS;
                        case 39
                            initialVariable = Results.ts_NumberAIDSAware;
                    end 

                    % takes the value from the last timestep of each year 
                    % and applies it as the value each year

                    if Params.SolveCont == 1 % if continuous
                        finalVariable(yearCount,:) = ...
                            initialVariable(sum(stepsPerYear(1:yearCount)),:);
                    else % discrete
                        % the discrete version has the same number of steps
                        % each year while the continuous version varies
                        finalVariable(yearCount,:) = initialVariable(stepsPerYear*yearCount,:);
                    end

                    switch u
                        case 2
                            Results.ann_NumberCCStages345 = finalVariable;
                        case 3
                            Results.ann_NumberAware = finalVariable;
                        case 4
                            Results.ann_NumberUnaware = finalVariable;
                        case 5
                            Results.ann_HIVPrevalence = finalVariable;
                        case 6   
                            Results.ann_NumberCCStage2 = finalVariable;
                        case 7
                            Results.ann_NumberInCare = finalVariable;
                        case 8    
                            Results.ann_NumberCCStage4 = finalVariable;
                        case 9
                            Results.ann_NumberUninfected = finalVariable;
                        case 10    
                            Results.ann_NumberAcute = finalVariable;
                        case 11
                            Results.ann_NumberLatentA = finalVariable;
                        case 12
                            Results.ann_NumberLatentB = finalVariable;
                        case 13    
                            Results.ann_NumberLate = finalVariable;
                        case 14
                            Results.ann_NumberAIDS = finalVariable;
                        case 15
                            Results.ann_CumulativeNumNormDeath = finalVariable;
                        case 16
                            Results.ann_CumulativeNumAIDSDeath = finalVariable;                       
                        case 17
                            Results.ann_PopulationSize = finalVariable;
                        case 18
                            Results.ann_NumberCCStages45 = finalVariable;
                        case 19
                            Results.ann_NumberVLS = finalVariable;
                        case 20
                            Results.ann_NumberOnART = finalVariable;
                        case 21
                            Results.ann_NumberCCStage3 = finalVariable;
                        case 22
                            Results.ann_NumberCCStage4_OnART = finalVariable;
                        case 23
                            Results.ann_NumberCCStage4_InCareInclART = finalVariable;
                        case 24
                            Results.ann_NumberOnPrEP = finalVariable;
                        case 25
                            Results.ann_PctSusOnPrEP = finalVariable;
                       case 26
                            Results.ann_NumberEligForPrEP = finalVariable; 
                       case 27
                            Results.ann_NumberEligForPrEPNotOnPrEP = finalVariable;    
                        case 28
                            Results.ann_NumberOnPrEP_Oral = finalVariable;
                        case 29
                            Results.ann_NumberOnPrEP_Inject = finalVariable;
                        case 30
                            Results.ann_NumberOnPrEP_Oral_HighAdh = finalVariable;
                        case 31
                            Results.ann_NumberOnPrEP_Oral_LowAdh = finalVariable;
                        case 32
                            Results.ann_NumberOnPrEP_Inject_HighAdh = finalVariable;
                        case 33
                            Results.ann_NumberOnPrEP_Inject_LowAdh = finalVariable;
                        case 34
                            Results.ann_NumberAIDSUndiag = finalVariable;
                        case 35 
                            Results.ann_NumberAIDSCCStage2 = finalVariable;
                        case 36
                            Results.ann_NumberAIDSCCStage3 = finalVariable;
                        case 37 
                            Results.ann_NumberAIDSCCStage4 = finalVariable;
                        case 38 
                            Results.ann_NumberAIDSVLS = finalVariable;
                        case 39
                            Results.ann_NumberAIDSAware = finalVariable;                            
                    end                                   
                    
                end 
            end
        end

    % 7.iv. Calculate transition/count outcomes annually
        for sectionTransition = 1:1

            % This code takes sums the values of the variable across the
            % entire year and applies it as the annual value.
            
            % Re-arrange 3d outcome for infections by subpop source (did
            % not collect correctly with original arrangement of dims)
            if Params.CollectInfbySource == 1
                Results.ts_NewInfectionsbySource = permute(Results.ts_NewInfectionsbySource,[3 1 2]);
                Results.ann_NewInfectionsbySource_Age(Params.numAge,Params.numAge,numYears) = 0;
                ann_NewInfectionsbySource_Age(Params.numStrats,Params.numStrats,numYears) = 0;
            end
            
            for u = 22:110 
                
                % ET Note: the counter starts at 22 because I originally
                % thought it couldn't overlap with the previous counter.
                % That's not the case so feel free to begin it at any number.
                % Added NewDiagnoses in first 2 disease stages as cases 99 and 100: Clinkscales 04/21/2022
                
                for yearCount = 1:numYears
                switch u

                    % Set up the initial outputs
                    % For variables that are cumulative
                    case 22
                        initialVariable = Results.ts_TotalNewInfections;
                    case 23
                        initialVariable = Results.ts_NewInfectionsVonly;
                    case 24
                        initialVariable = Results.ts_NewInfectionsAonly;
                    case 25
                        initialVariable = Results.ts_NewInfectionsVsomeAI;
                    case 26
                        initialVariable = Results.ts_NewInfectionsAsomeVI;
                    case 27
                        initialVariable = Results.ts_NewInfectionsN;
                    case 28
                        initialVariable = Results.ts_TotalNewDiagnoses;
                    case 29
                        initialVariable = Results.ts_TotalLinkToCareFirst;
                    case 30
                        initialVariable = Results.ts_TotalLinkToCareAfterFirst;
                    case 31
                        initialVariable = Results.ts_TotalLinkToCare;
                    case 32
                        initialVariable = Results.ts_TotalDropOutofInCare;
                    case 33
                        initialVariable = Results.ts_TotalDropOff_ANVToAware;
                    case 34
                        initialVariable = Results.ts_TotalStartARTNotVLS;  
                    case 35
                        initialVariable = Results.ts_ARTCarePrEPCost;
                    case 36
                        initialVariable = Results.ts_TotalNewDiagnosesNotLTC;
                    case 37
                        initialVariable = Results.ts_QALYs;
                    case 38
                        initialVariable = Results.ts_LifeYears;
                    case 39
                        initialVariable = Results.ts_QALYs_A;
                    case 40
                        initialVariable = Results.ts_QALYs_B;
                    case 41
                        initialVariable = Results.ts_QALYs_C;
                    case 42
                        initialVariable = Results.ts_QALYs_D;
                    case 43
                        initialVariable = Results.ts_QALYs_E;
                    case 44
                        initialVariable = Results.ts_QALYs_F;
                    case 45
                        initialVariable = Results.ts_LifeYears_A;
                    case 46
                        initialVariable = Results.ts_LifeYears_B;
                    case 47
                        initialVariable = Results.ts_LifeYears_C;
                    case 48
                        initialVariable = Results.ts_LifeYears_D;
                    case 49
                        initialVariable = Results.ts_LifeYears_E;
                    case 50
                        initialVariable = Results.ts_LifeYears_F;
                    case 51
                        initialVariable = Results.ts_numPosRapidAcuteTests;
                    case 52
                        initialVariable = Results.ts_numPosConvAcuteTests;
                    case 53
                        initialVariable = Results.ts_numPosRapidChronicTests;
                    case 54
                        initialVariable = Results.ts_numPosConvChronicTests;
                    case 55
                        initialVariable = Results.ts_numNegRapidTests;
                    case 56
                        initialVariable = Results.ts_numNegConvTests;
                    case 57
                        initialVariable = Results.ts_numConvPosNotify;
                    case 58
                        initialVariable = Results.ts_numConvNegNotify;
                    case 59
                        initialVariable = Results.ts_numRapidPosNotify;
                    case 60
                        initialVariable = Results.ts_numRapidNegNotify;
                    case 61
                        initialVariable = Results.ts_numNormDeathsHIVPosNoART;
                    case 62
                        initialVariable = Results.ts_numNormDeathsANVandVLS;
                    case 63
                        initialVariable = Results.ts_numNormDeathsHIVPosUndiag;
                    case 64
                        initialVariable = Results.ts_numAIDSDeathsNoART;
                    case 65
                        initialVariable = Results.ts_numAIDSDeathsANVandVLS;
                    case 66
                        initialVariable = Results.ts_TotalDropOff_ANVToCare;
                    case 67
                        initialVariable = Results.ts_TotalDropOff_VLSToANV;
                    case 68
                        initialVariable = Results.ts_BecomeVLSFromLTC;
                    case 69
                        initialVariable = Results.ts_numDeathsPLWHAware;
                    case 70
                        initialVariable = Results.ts_BecomeVLSFromANV;
                    case 71
                        initialVariable = Results.ts_DropOff_PrEP;
                    case 72
                        initialVariable = Results.ts_StartPrEP;
                    case 73
                        initialVariable = Results.ts_NewInfections_A1;
                    case 74
                        initialVariable = Results.ts_NewInfections_OnPrEP;
                    case 75
                        initialVariable = Results.ts_healthStateCost_PrEP;
                    case 76
                        initialVariable = Results.ts_numDeathsHIVPosUndiag;
                    case 77
                        initialVariable = Results.ts_NumberCCStage2;
                    case 78
                        initialVariable = Results.ts_NumberCCStage3;
                    case 79
                        initialVariable = Results.ts_NumberCCStage4;                        
                    case 80
                        initialVariable = Results.ts_NumberVLS; 
                    case 81
                        initialVariable = Results.ts_NumberEligForPrEPNotOnPrEP;
                    case 82
                        initialVariable = Results.ts_TotalARTAndCareCost;
                    case 83
                        initialVariable = Results.ts_ARTCarePrEPCost_Undisc;
                    case 84
                        initialVariable = Results.ts_TotalARTAndCareCost_Undisc;
                    case 85
%                         initialVariable = Results.ts_TransCost_RemainVLS;
                    case 86
%                         initialVariable = Results.ts_TransCost_RemainVLS_Undisc;
                    case 87    
                        initialVariable = Results.ts_healthStateCost_PrEP_Oral_Undisc;                       
                    case 88
                        initialVariable = Results.ts_healthStateCost_PrEP_Inject_Undisc;    
                    case 89    
                        initialVariable = Results.ts_healthStateCost_PrEP_Undisc;                       
                    case 90
                        initialVariable = Results.ts_DropOff_PrEP_Oral;
                    case 91
                        initialVariable = Results.ts_DropOff_PrEP_Inject;
                    case 92
                        initialVariable = Results.ts_StartPrEP_Oral;
                    case 93
                        initialVariable = Results.ts_StartPrEP_Inject;
                    case 94
                        initialVariable = Results.ts_healthStateCost_PrEP_Oral;    
                    case 95
                        initialVariable = Results.ts_healthStateCost_PrEP_Inject;
                    case 96
                        initialVariable = Results.ts_TimeOnPrEP;
                    case 97
                        initialVariable = Results.ts_TimeOnPrEP_Oral;
                    case 98
                        initialVariable = Results.ts_TimeOnPrEP_Inject;
                    case 99
                        initialVariable = Results.ts_NewDiagnoses_B1;
                    case 100
                        initialVariable = Results.ts_NewDiagnoses_C1;
                    case 101
                        initialVariable = Results.ts_numDeathsVLS;
                    case 102
                        initialVariable = Results.ts_numAIDSDeathsVLS;
                    case 103
                        initialVariable = Results.ts_numDeathsCCStage2;
                    case 104
                        initialVariable = Results.ts_numAIDSDeathsCCStage2;
                    case 105
                        initialVariable = Results.ts_numDeathsCCStage3;
                    case 106
                        initialVariable = Results.ts_numAIDSDeathsCCStage3;
                    case 107
                        initialVariable = Results.ts_numDeathsCCStage4;
                    case 108
                        initialVariable = Results.ts_numAIDSDeathsCCStage4;
                    case 109
                        initialVariable = Results.ts_numAIDSDeathsHIVPosUndiag;
                    case 110
                        initialVariable = Results.ts_numAIDSDeathsAware;
                end


                if Params.SolveCont == 1 % continuous

                    if yearCount == 1 % first year in the loop
                        finalVariable(yearCount,:)=sum(initialVariable(1:stepsPerYear(1),:));
                    else
                       finalVariable(yearCount,:)= sum(initialVariable(1 + ...
                           sum(stepsPerYear(1:yearCount-1)): sum(stepsPerYear(1:yearCount)),:));
                    end

                else % discrete

                    finalVariable(yearCount,:) = sum(initialVariable(1+stepsPerYear...
                        *(yearCount-1) : stepsPerYear*yearCount, :));       
                end

                switch u

                    case 22
                        Results.ann_TotalNewInfections = finalVariable;
                    case 23
                        Results.ann_NewInfectionsVonly = finalVariable;
                    case 24
                        Results.ann_NewInfectionsAonly = finalVariable;
                    case 25
                        Results.ann_NewInfectionsVsomeAI = finalVariable;
                    case 26
                        Results.ann_NewInfectionsAsomeVI = finalVariable;
                    case 27
                        Results.ann_NewInfectionsN = finalVariable;
                    case 28
                        Results.ann_TotalNewDiagnoses= finalVariable;
                    case 29
                        Results.ann_TotalLinkToCareFirst = finalVariable;
                    case 30
                        Results.ann_TotalLinkToCareAfterFirst = finalVariable;
                    case 31
                        Results.ann_TotalLinkToCare = finalVariable;
                    case 32
                        Results.ann_TotalDropOutofInCare = finalVariable;
                    case 33
                        Results.ann_TotalDropOff_ANVToAware =  finalVariable;
                    case 34
                        Results.ann_TotalStartARTNotVLS = finalVariable;
                    case 35
                        Results.ann_ARTCarePrEPCost_Disc = finalVariable;
                    case 36
                        Results.ann_TotalNewDiagnosesNotLTC = finalVariable;
                    case 37
                        Results.ann_QALYs = finalVariable;
                    case 38
                        Results.ann_LifeYears = finalVariable;
                    case 39
                        Results.ann_QALYs_A = finalVariable;
                    case 40
                        Results.ann_QALYs_B = finalVariable;
                    case 41
                        Results.ann_QALYs_C = finalVariable;
                    case 42
                        Results.ann_QALYs_D = finalVariable;
                    case 43
                        Results.ann_QALYs_E = finalVariable;
                    case 44
                        Results.ann_QALYs_F = finalVariable;
                    case 45
                        Results.ann_LifeYears_A = finalVariable;
                    case 46
                        Results.ann_LifeYears_B = finalVariable;
                    case 47
                        Results.ann_LifeYears_C = finalVariable;
                    case 48
                        Results.ann_LifeYears_D = finalVariable;
                    case 49
                        Results.ann_LifeYears_E = finalVariable;
                    case 50
                        Results.ann_LifeYears_F = finalVariable;
                    case 51
                        Results.ann_numPosRapidAcuteTests = finalVariable;
                    case 52
                        Results.ann_numPosConvAcuteTests = finalVariable;
                    case 53
                        Results.ann_numPosRapidChronicTests = finalVariable;
                    case 54
                        Results.ann_numPosConvChronicTests = finalVariable;
                    case 55
                         Results.ann_numNegTests_Rapid = finalVariable;
                    case 56
                        Results.ann_numNegTests_Conv = finalVariable;
                    case 57
                        Results.ann_numNotifiedConvPos = finalVariable;
                    case 58
                        Results.ann_numNotifiedConvNeg = finalVariable;
                    case 59
                        Results.ann_numNotifiedRapidPos = finalVariable;
                    case 60
                        Results.ann_numNotifiedRapidNeg = finalVariable;
                    case 61 
                        Results.ann_numNormDeathsHIVPosNoART = finalVariable;
                    case 62
                        Results.ann_numNormDeathsANVandVLS = finalVariable;
                    case 63 
                        Results.ann_numNormDeathsHIVPosUndiag = finalVariable;
                    case 64
                        Results.ann_numAIDSDeathsNoART = finalVariable;
                    case 65
                        Results.ann_numAIDSDeathsANVandVLS = finalVariable;
                    case 66
                        Results.ann_TotalDropOff_ANVToCare = finalVariable;
                    case 67
                        Results.ann_TotalDropOff_VLSToANV = finalVariable;
                    case 68
                        Results.ann_BecomeVLSFromLTC = finalVariable;
                    case 69
                        Results.ann_numDeathsPLWHAware = finalVariable;
                    case 70
                        Results.ann_BecomeVLSFromANV = finalVariable;
                    case 71                            
                        Results.ann_DropOff_PrEP = finalVariable;
                    case 72
                        Results.ann_StartPrEP = finalVariable;
                    case 73
                        Results.ann_NewInfections_NoPrEP = finalVariable;
                    case 74
                        Results.ann_NewInfections_OnPrEP = finalVariable;
                    case 75
                        Results.ann_HealthStateCost_PrEP_Disc = finalVariable .* AllocPopsTargetedIndicator';
                    case 76
                        Results.ann_numDeathsHIVPosUndiag = finalVariable;
                    case 77
                        Results.ann_sumNumberCCStage2 = finalVariable;
                    case 78
                        Results.ann_sumNumberCCStage3 = finalVariable;
                    case 79
                        Results.ann_sumNumberCCStage4 = finalVariable;
                    case 80
                        Results.ann_sumNumberVLS = finalVariable;
                    case 81
                        Results.ann_sumTotalNumberEligForPrEPNotOnPrEP = finalVariable;
                    case 82
                        Results.ann_TotalARTAndCareCost_Disc = finalVariable .* AllocPopsTargetedIndicator'; 
                    case 83
                        Results.ann_ARTCarePrEPCost_Undisc = finalVariable .* AllocPopsTargetedIndicator';
                    case 84
                        Results.ann_TotalARTAndCareCost_Undisc = finalVariable .* AllocPopsTargetedIndicator';
                    case 85
%                         Results.ann_TransCost_RemainVLS_Disc = finalVariable .* AllocPopsTargetedIndicator';
                    case 86
%                         Results.ann_TransCost_RemainVLS_Undisc = finalVariable .* AllocPopsTargetedIndicator';
                    case 87
                        Results.ann_healthStateCost_PrEP_Oral_Undisc = finalVariable .* AllocPopsTargetedIndicator';
                    case 88
                        Results.ann_healthStateCost_PrEP_Inject_Undisc = finalVariable .* AllocPopsTargetedIndicator';    
                    case 89
                        Results.ann_healthStateCost_PrEP_Undisc = finalVariable .* AllocPopsTargetedIndicator';
                    case 90
                        Results.ann_DropOff_PrEP_Oral = finalVariable;
                    case 91
                        Results.ann_DropOff_PrEP_Inject = finalVariable;
                    case 92
                        Results.ann_StartPrEP_Oral = finalVariable;
                    case 93
                        Results.ann_StartPrEP_Inject = finalVariable;
                    case 94
                        Results.ann_HealthStateCost_PrEP_Oral_Disc = finalVariable;    
                    case 95
                        Results.ann_HealthStateCost_PrEP_Inject_Disc = finalVariable;  
                    case 96
                        Results.ann_PersonYears_onPrEP = finalVariable;
                    case 97
                        Results.ann_PersonYears_onPrEP_Oral = finalVariable;
                    case 98
                        Results.ann_PersonYears_onPrEP_Inject = finalVariable;
                    case 99
                        Results.ann_NewDiagnoses_B1 = finalVariable;
                    case 100
                        Results.ann_NewDiagnoses_C1 = finalVariable; 
                    case 101
                        Results.ann_numDeathsVLS = finalVariable;
                    case 102
                        Results.ann_numAIDSDeathsVLS = finalVariable;
                    case 103
                        Results.ann_numDeathsCCStage2 = finalVariable;
                    case 104
                        Results.ann_numAIDSDeathsCCStage2 = finalVariable;
                    case 105
                        Results.ann_numDeathsCCStage3 = finalVariable;
                    case 106
                        Results.ann_numAIDSDeathsCCStage3 = finalVariable;
                    case 107
                        Results.ann_numDeathsCCStage4 = finalVariable;
                    case 108
                        Results.ann_numAIDSDeathsCCStage4 = finalVariable;
                    case 109
                        Results.ann_numAIDSDeathsHIVPosUndiag = finalVariable;
                    case 110
                        Results.ann_numAIDSDeathsAware = finalVariable;
                end
                
                    % Calculate annual outcome for new infections, by source of infection
                    if Params.CollectInfbySource == 1
                        
                        ann_NewInfectionsbySource(:,:) = sum(Results.ts_NewInfectionsbySource(1+stepsPerYear...
                                *(yearCount-1) : stepsPerYear*yearCount,:,:));   

                        ann_NewInfectionsbySource_Age = ann_NewInfectionsbySource * Params.ageIndicator;
                        Results.ann_NewInfectionsbySource_Age(:,:,yearCount) = ann_NewInfectionsbySource_Age' * Params.ageIndicator;                   
                    end
                    
                end
            end
            
            % Re-arrange 3d outcomes for infections by subpop source (for
            % some reason, if you do not re-arrange back to
            % 273x273xnumtsteps, it gives you trouble viewing outcomes
            % after the model run)
            
            if Params.CollectInfbySource == 1
                Results.ts_NewInfectionsbySource = permute(Results.ts_NewInfectionsbySource,[2 3 1]);
            end
            
            % Calculate compartments annually
            permutedCompartments = permute(Results.store_Compartments,[3 1 2]);
            
                % takes the value from the last timestep of each year 
                % and applies it as the value each year

                if Params.SolveCont == 1 % if continuous
                    for yearCount = 1:numYears
                        Results.ann_Compartments(yearCount,:,:) = ...
                            permutedCompartments(sum(stepsPerYear(1:yearCount)),:,:);
                    end
                else % discrete
                    % the discrete version has the same number of steps
                    % each year while the continuous version varies
                    for yearCount = 1:numYears
                        Results.ann_Compartments(yearCount,:,:) = permutedCompartments(stepsPerYear*yearCount,:,:);
                    end
                end
                
                Results.ann_NumInCompartments = sum(Results.ann_Compartments,3);
                               
        end
           
    % 7.v. Calculate Costs
    
        % Pre-allocate
        costPosRapidAcuteTests(numYears,Params.numStrats) = 0;
        costPosRapidChronicTests(numYears,Params.numStrats) = 0;
        costPosConvAcuteTests(numYears,Params.numStrats) = 0;
        costPosConvChronicTests(numYears,Params.numStrats) = 0;
        
        for sectionCosts = 1:1

            % 7v.i) Cost of transitions each year (except testing from
            % unaware not on PrEP)
            for sectionTransitionCosts = 1:1
                   
                   
                
                   % TRANSITIONS/COUNTS: Continuum of care transition costs (all discounted) 
           
                   % Undiscounted
                    Results.ann_TransCost_DiagOnPrEP_Undisc = ...
                        Results.ann_NewInfections_OnPrEP * Params.transcostPP_addtlPrEPTestCostIfInf;
                    Results.ann_TransCost_LTCFirst_Undisc = ...
                        Results.ann_TotalLinkToCareFirst * Params.transcostPP_LTCFirst;
                    Results.ann_TransCost_LTCAfter_Undisc = ...
                        Results.ann_TotalLinkToCareAfterFirst * Params.transcostPP_LTCAfter;
                    Results.ann_TransCost_ARTInitiation_Undisc = ...
                        (Results.ann_TotalStartARTNotVLS + Results.ann_BecomeVLSFromLTC) * Params.transcostPP_ARTInitiation;
                    Results.ann_TransCost_BecomeVLSfromANV_Undisc = ...
                        Results.ann_BecomeVLSFromANV * Params.transcostPP_BecomeVLSfromANV;                   
                    %Results.ann_TransCost_RemainVLS_Undisc = calculated
                    %earlier in code
                    Results.ann_TransCost_RemainVLS_Undisc = ...
                        Results.ann_NumberVLS .* ...
                        (1- RateToProb(Results.ann_TotalDropOff_VLSToANV ./ max(Results.ann_NumberVLS,1) )) ...
                        * Params.anntranscostPP_RemainVLS;
                   
                   % Undiscounted - YMSM only      
                    Results.ann_TransCost_DiagOnPrEP_YMSM_Undisc = ...
                        Results.ann_NewInfections_OnPrEP * Params.transcostPP_addtlPrEPTestCostIfInf .* Params.YoungMSMIndicator';
                    Results.ann_TransCost_LTCFirst_YMSM_Undisc = ...
                        Results.ann_TotalLinkToCareFirst * Params.transcostPP_LTCFirst .* Params.YoungMSMIndicator';
                    Results.ann_TransCost_LTCAfter_YMSM_Undisc = ...
                        Results.ann_TotalLinkToCareAfterFirst * Params.transcostPP_LTCAfter .* Params.YoungMSMIndicator';
                    Results.ann_TransCost_ARTInitiation_YMSM_Undisc = ...
                        (Results.ann_TotalStartARTNotVLS + Results.ann_BecomeVLSFromLTC) * Params.transcostPP_ARTInitiation .* Params.YoungMSMIndicator';
                    Results.ann_TransCost_BecomeVLSfromANV_YMSM_Undisc = ...
                        Results.ann_BecomeVLSFromANV * Params.transcostPP_BecomeVLSfromANV .* Params.YoungMSMIndicator'; 
                    Results.ann_TransCost_RemainVLS_YMSM_Undisc = ...
                        Results.ann_TransCost_RemainVLS_Undisc .* Params.YoungMSMIndicator';       

                   % Discounted
                    Results.ann_TransCost_DiagOnPrEP_Disc = ...
                        Results.ann_NewInfections_OnPrEP * Params.transcostPP_addtlPrEPTestCostIfInf .* (ann_DiscRates_Costs * ones(1,numStrats));
                    Results.ann_TransCost_LTCFirst_Disc = ...
                        Results.ann_TotalLinkToCareFirst * Params.transcostPP_LTCFirst .* (ann_DiscRates_Costs * ones(1,numStrats));
                    Results.ann_TransCost_LTCAfter_Disc = ...
                        Results.ann_TotalLinkToCareAfterFirst * Params.transcostPP_LTCAfter .* (ann_DiscRates_Costs * ones(1,numStrats));
                    Results.ann_TransCost_ARTInitiation_Disc = ...
                        (Results.ann_TotalStartARTNotVLS + Results.ann_BecomeVLSFromLTC) * Params.transcostPP_ARTInitiation .* (ann_DiscRates_Costs * ones(1,numStrats));
                    Results.ann_TransCost_BecomeVLSfromANV_Disc = ...
                        Results.ann_BecomeVLSFromANV * Params.transcostPP_BecomeVLSfromANV .* (ann_DiscRates_Costs * ones(1,numStrats));
% KH code below added on 16Feb2022 - need to check
                    Results.ann_TransCost_RemainVLS_Disc = ...
                        Results.ann_NumberVLS .* ...
                        ((1- RateToProb(Results.ann_TotalDropOff_VLSToANV ./ max(Results.ann_NumberVLS,1) )) ...
                        * Params.anntranscostPP_RemainVLS) .* (ann_DiscRates_Costs * ones(1,numStrats));                    
                    
                    % Total including testing of unaware not on PrEP
                    % calculated below testing costs calcs (Results.ann_TotalTransCost_Disc)
                
            end
            
      
            % 7v.ii) Intervention and Testing Costs
                        
            for sectionIntnCosts = 1:1
                

                    % Costs for positive tests + notify + NAT if
                    % applicable
                        % NAT is only applied if confirmatory test was
                        % negative
  
                        % [Total cost of pos rapid/conv tests for Acute/Chronic HIV] = 
                        % [Num pos rapid/conv acute/chronic tests ] 
                        % * ([Cost per rapid/conv test and confirm test]
                        % + [Outreach cost]
                        % + (1-sens of confirm for acute/chronic)
                        % * [NAT Cost] )
                
                     
                        % Total number of non-CDC tests
                            % Positive tests
                         Results.ann_numPosTests = ...
                             Results.ann_numPosRapidAcuteTests + ...
                             Results.ann_numPosRapidChronicTests + ...
                             Results.ann_numPosConvAcuteTests + ...
                             Results.ann_numPosConvChronicTests;
                            % Negative tests
                         Results.ann_numNegTests = ...
                             Results.ann_numNegTests_Rapid + ...
                             Results.ann_numNegTests_Conv;

                         Results.ann_numTotalTests = ...
                             Results.ann_numPosTests + ...
                             Results.ann_numNegTests;
                         
                        for yearCount = 1:numYears
                       % non-allocation
                            costPosRapidAcuteTests(yearCount,:) = ...
                                Results.ann_numPosRapidAcuteTests(yearCount,:) .* ...
                                (Results.ann_costPP_TestPosRapid(yearCount) + ...
                                Params.costPP_TestOutreach + ...
                                (1 - Results.ann_TestSens_Confirm_Acute(yearCount)) ...
                                * Params.costPP_NATTest);
                            
                             costPosRapidChronicTests(yearCount,:) = ...
                                Results.ann_numPosRapidChronicTests(yearCount,:) .* ...
                                (Results.ann_costPP_TestPosRapid(yearCount) + ...
                                Params.costPP_TestOutreach + ...
                                (1 - Results.ann_TestSens_Confirm_Chronic(yearCount)) ...
                                * Params.costPP_NATTest);
                            
                            costPosConvAcuteTests(yearCount,:) = ...
                                Results.ann_numPosConvAcuteTests(yearCount,:) .* ...
                                (Results.ann_costPP_TestPosConv(yearCount) + ...
                                Params.costPP_TestOutreach + ...
                                (1 - Results.ann_TestSens_Confirm_Acute(yearCount)) ...
                                * Params.costPP_NATTest);
                            
                             costPosConvChronicTests(yearCount,:) = ...
                                Results.ann_numPosConvChronicTests(yearCount,:) .* ...
                                (Results.ann_costPP_TestPosConv(yearCount) + ...
                                Params.costPP_TestOutreach + ...
                                (1 - Results.ann_TestSens_Confirm_Chronic(yearCount)) ...
                                * Params.costPP_NATTest);
                            
                            
                            % allocation
             
                            %costPosRapidAcuteTestsFromCDC(yearCount,:) = ...
                            %    Results.ann_numPosRapidAcuteTestsCDCFunds(yearCount,:) .* ...
                            %    (Results.ann_costPP_TestPosRapid(yearCount) + ...
                            %    Params.costPP_TestOutreach + ...
                            %    (1 - Results.ann_TestSens_Confirm_Acute(yearCount)) ...
                            %    * Params.costPP_NATTest);
                            
                             %costPosRapidChronicTestsFromCDC(yearCount,:) = ...
                             %  Results.ann_numPosRapidChronicTestsCDCFunds(yearCount,:) .* ...
                             %   (Results.ann_costPP_TestPosRapid(yearCount) + ...
                             %   Params.costPP_TestOutreach + ...
                             %   (1 - Results.ann_TestSens_Confirm_Chronic(yearCount)) ...
                             %   * Params.costPP_NATTest);
                            
                            %costPosConvAcuteTestsFromCDC(yearCount,:) = ...
                            %    Results.ann_numPosConvAcuteTestsCDCFunds(yearCount,:) .* ...
                            %    (Results.ann_costPP_TestPosConv(yearCount) + ...
                            %    Params.costPP_TestOutreach + ...
                            %    (1 - Results.ann_TestSens_Confirm_Acute(yearCount)) ...
                            %    * Params.costPP_NATTest);
                            
                             %costPosConvChronicTestsFromCDC(yearCount,:) = ...
                             %   Results.ann_numPosConvChronicTestsCDCFunds(yearCount,:) .* ...
                             %   (Results.ann_costPP_TestPosConv(yearCount) + ...
                             %   Params.costPP_TestOutreach + ...
                             %   (1 - Results.ann_TestSens_Confirm_Chronic(yearCount)) ...
                             %   * Params.costPP_NATTest);
                            
                        end
                        
                    % Costs for negative rapid tests (Undiscounted)
                    ann_costNegRapidTests_Undisc = ...
                         Results.ann_numNegTests_Rapid ...
                         .* (Results.ann_costPP_TestNegRapid * ones(1,numStrats) ...
                         + ones(numYears,1)*Params.costPP_TestOutreach);
                    
                     % Costs for negative rapid tests
                     ann_costNegRapidTests = ...
                         Results.ann_numNegTests_Rapid ...
                         .* (Results.ann_costPP_TestNegRapid * ones(1,numStrats) ...
                         + ones(numYears,1)*Params.costPP_TestOutreach)...
                         .* (ann_DiscRates_Costs*ones(1,numStrats));
                     
                     % ann_costNegRapidTests_FromCDC = ...
                     %     Results.ann_numNegRapidTestsCDCFunds ...
                     %    .* (Results.ann_costPP_TestNegRapid * ones(1,numStrats) ...
                     %    + ones(numYears,1)*Params.costPP_TestOutreach)...
                     %    .* (ann_DiscRates_Costs*ones(1,numStrats));

                     
                    % Cost negative conv tests (Undiscounted)
                       ann_costNegConvTests_Undisc = ...
                         Results.ann_numNegTests_Conv ...
                         .* (Results.ann_costPP_TestNegConv * ones(1,numStrats) ...
                         + ones(numYears,1)*Params.costPP_TestOutreach);
                     
                     % Cost negative conv tests
                       ann_costNegConvTests = ...
                         Results.ann_numNegTests_Conv ...
                         .* (Results.ann_costPP_TestNegConv * ones(1,numStrats) ...
                         + ones(numYears,1)*Params.costPP_TestOutreach)...
                         .* (ann_DiscRates_Costs*ones(1,numStrats));
                     
                       %ann_costNegConvTests_FromCDC = ...
                       %  Results.ann_numNegConvTestsCDCFunds ...
                       %  .* (Results.ann_costPP_TestNegConv * ones(1,numStrats) ...
                       %  + ones(numYears,1)*Params.costPP_TestOutreach)...
                       %  .* (ann_DiscRates_Costs*ones(1,numStrats));
                     ann_costNegTests_Undisc = ...
                            ann_costNegRapidTests_Undisc + ann_costNegConvTests_Undisc;  
                       
                     Results.ann_costNegTests_Disc = ...
                            ann_costNegRapidTests + ann_costNegConvTests;
                     %Results.ann_costNegTests_FromCDC = ...
                     %       ann_costNegRapidTests_FromCDC + ann_costNegConvTests_FromCDC;
                            
                            
            % Costs for tests and notification for rapid tests
                    % Costs for positive rapid tests + notify
                    % (Undiscounted)
                    ann_costAllPosTests_Undisc = ...
                         (costPosRapidAcuteTests + costPosRapidChronicTests ...
                         + costPosConvAcuteTests + costPosConvChronicTests);
                     
                     % Costs for positive rapid tests + notify
                    Results.ann_costAllPosTests_Disc = ...
                         (costPosRapidAcuteTests + costPosRapidChronicTests ...
                         + costPosConvAcuteTests + costPosConvChronicTests) ...
                         .* (ann_DiscRates_Costs * ones(1,numStrats));
                     
                    %Results.ann_costAllPosTests_FromCDC = ...
                    %     (costPosRapidAcuteTestsFromCDC + costPosRapidChronicTestsFromCDC ...
                    %     + costPosConvAcuteTestsFromCDC + costPosConvChronicTestsFromCDC) ...
                    %     .* (ann_DiscRates_Costs * ones(1,numStrats));
                    
                    
                   % Pull summary results (Undiscounted)
                    ann_costAllTesting_Undisc = ...
                        ann_costAllPosTests_Undisc + ann_costNegTests_Undisc;
                    
                    % Pull summary results
                    Results.ann_costAllTesting_Disc = ...
                        Results.ann_costAllPosTests_Disc + Results.ann_costNegTests_Disc;
                       
                    %Results.ann_costAllTesting_FromCDC = ...
                    %    Results.ann_costAllPosTests_FromCDC + Results.ann_costNegTests_FromCDC;
                    
                    Results.ann_costTotalTesting_Disc = ...
                        Results.ann_costAllTesting_Disc; % + Results.ann_costAllTesting_FromCDC;
                    

                  
                        %Conv Pos

                         Results.ann_costNotifyConvPos_Disc = ...
                            Results.ann_numNotifiedConvPos ...
                            * Params.costPP_Notify_PosConv .* (ann_DiscRates_Costs*ones(1,numStrats));
                        %Conv Neg
                        Results.ann_costNotifyConvNeg_Disc = ...
                            Results.ann_numNotifiedConvNeg ...
                            * Params.costPP_Notify_NegConv .* (ann_DiscRates_Costs*ones(1,numStrats));               
                        % Rapid pos
                        Results.ann_costNotifyRapidPos_Disc = ...
                            Results.ann_numNotifiedRapidPos ...
                            * Params.costPP_Notify_PosRapid .* (ann_DiscRates_Costs*ones(1,numStrats));          
                        %Rapid Neg
                        Results.ann_costNotifyRapidNeg_Disc = ...
                            Results.ann_numNotifiedRapidNeg ...
                            * Params.costPP_Notify_NegRapid .* (ann_DiscRates_Costs*ones(1,numStrats));
                        
                        % Undisc Notification costs
                        %Conv Pos

                        ann_costNotifyConvPos_Undisc = ...
                            Results.ann_numNotifiedConvPos ...
                            * Params.costPP_Notify_PosConv;
                        %Conv Neg
                        ann_costNotifyConvNeg_Undisc = ...
                            Results.ann_numNotifiedConvNeg ...
                            * Params.costPP_Notify_NegConv;               
                        % Rapid pos
                        ann_costNotifyRapidPos_Undisc = ...
                            Results.ann_numNotifiedRapidPos ...
                            * Params.costPP_Notify_PosRapid;          
                        %Rapid Neg
                        ann_costNotifyRapidNeg_Undisc = ...
                            Results.ann_numNotifiedRapidNeg ...
                            * Params.costPP_Notify_NegRapid;
                        
                        %CDC Rapid positive
                        %Results.ann_costNotifyRapidPosFromCDC = ...
                        %    Results.ann_numRapidPosNotifiedFromCDCFunds ...
                        %    * Params.costPP_Notify_PosRapid .* (ann_DiscRates_Costs*ones(1,numStrats));
                        
                        %Results.ann_costNotifyConvPosFromCDC = ...
                        %    Results.ann_numConvPosNotifiedFromCDCFunds ...
                        %    * Params.costPP_Notify_PosConv .* (ann_DiscRates_Costs*ones(1,numStrats));
                        
                        % Total cost to notify (undiscounted)
                        ann_totalCostNotify_Undisc = ...
                            ann_costNotifyConvPos_Undisc + ann_costNotifyConvNeg_Undisc + ...
                            ann_costNotifyRapidPos_Undisc + ann_costNotifyRapidNeg_Undisc;
                        
                        % Total cost to notify 
                        Results.ann_totalCostNotify_Disc = ...
                            Results.ann_costNotifyConvPos_Disc + Results.ann_costNotifyConvNeg_Disc + ...
                            Results.ann_costNotifyRapidPos_Disc + Results.ann_costNotifyRapidNeg_Disc;
                        
                        % Cost to notify CDC
                        %Results.ann_totalCostNotifyPos_CDCFunds = ...
                        %    Results.ann_costNotifyRapidPosFromCDC + Results.ann_costNotifyConvPosFromCDC;

                        
                        % added 4Nov2014 by ET
                % Total test and notification costs (Undisc)
                Results.ann_TransCost_CostTestAndNotify_Undisc = ...
                     (ann_totalCostNotify_Undisc + ... 
                     ann_costAllTesting_Undisc);
                
                % Total test and notification costs (Undisc) - YMSM only
                Results.ann_TransCost_CostTestAndNotify_YMSM_Undisc = ...
                     (ann_totalCostNotify_Undisc + ... 
                     ann_costAllTesting_Undisc) .* Params.YoungMSMIndicator';
                 
                % Total test and notification costs
                Results.ann_TransCost_CostTestAndNotify_Disc = ...
                     (Results.ann_totalCostNotify_Disc + ... 
                     Results.ann_costTotalTesting_Disc);
                
                % Total test and notification costs by transmission and
                % risk group (discounted)
                for yearCount = 1:numYears
                
                    % Total test and notification costs by transmission and
                    % risk group (undiscounted)
                    Results.ann_TransCost_CostTestAndNotify_HighRiskHETs_Undisc(yearCount,:) = ...
                        Results.ann_TransCost_CostTestAndNotify_Undisc(yearCount,:) .* Params.HRHIndicator';
                    Results.ann_TransCost_CostTestAndNotify_LowRiskHETs_Undisc(yearCount,:) = ...
                        Results.ann_TransCost_CostTestAndNotify_Undisc(yearCount,:) .* Params.LowRiskHETIndicator';
                    Results.ann_TransCost_CostTestAndNotify_HighRiskMSM_Undisc(yearCount,:) = ...
                        Results.ann_TransCost_CostTestAndNotify_Undisc(yearCount,:) .* Params.HighRiskMSMIndicator';
                    Results.ann_TransCost_CostTestAndNotify_LowRiskMSM_Undisc(yearCount,:) = ...
                        Results.ann_TransCost_CostTestAndNotify_Undisc(yearCount,:) .* Params.riskLevelIndicator(:,Params.risk_Main)' .* Params.popIndicator(:,Params.pop_MSM)';
                    Results.ann_TransCost_CostTestAndNotify_PWID_Undisc(yearCount,:) = ...
                        Results.ann_TransCost_CostTestAndNotify_Undisc(yearCount,:) .* Params.popIndicator(:,Params.pop_IDU)';    
                    
                    % Total test and notification costs by transmission and
                    % risk group (undiscounted) - YMSM
                    Results.ann_TransCost_CostTestAndNotify_HighRiskYMSM_Undisc(yearCount,:) = ...
                        Results.ann_TransCost_CostTestAndNotify_Undisc(yearCount,:) .* Params.HighRiskMSMIndicator' .* Params.YoungMSMIndicator';
                    Results.ann_TransCost_CostTestAndNotify_LowRiskYMSM_Undisc(yearCount,:) = ...
                        Results.ann_TransCost_CostTestAndNotify_Undisc(yearCount,:) .* Params.riskLevelIndicator(:,Params.risk_Main)' .* Params.popIndicator(:,Params.pop_MSM)' .* Params.YoungMSMIndicator';

                    % Total test and notification costs by transmission and
                    % risk group (discounted)
                    Results.ann_CostTestAndNotify_HighRiskMSM_Disc(yearCount,:) = ...
                        Results.ann_TransCost_CostTestAndNotify_Disc(yearCount,:) .* Params.HighRiskMSMIndicator';
                    Results.ann_CostTestAndNotify_LowRiskMSM_Disc(yearCount,:) = ...
                        Results.ann_TransCost_CostTestAndNotify_Disc(yearCount,:) .* Params.riskLevelIndicator(:,Params.risk_Main)' .* Params.popIndicator(:,Params.pop_MSM)';
                    Results.ann_CostTestAndNotify_PWID_Disc(yearCount,:) = ...
                        Results.ann_TransCost_CostTestAndNotify_Disc(yearCount,:) .* Params.popIndicator(:,Params.pop_IDU)';
                
                end
                 
                Results.ann_TotalTransCost_Disc = Results.ann_TransCost_DiagOnPrEP_Disc + Results.ann_TransCost_LTCFirst_Disc + ...
                    Results.ann_TransCost_LTCAfter_Disc + Results.ann_TransCost_ARTInitiation_Disc + ...
                    Results.ann_TransCost_BecomeVLSfromANV_Disc + Results.ann_TransCost_RemainVLS_Disc + Results.ann_TransCost_CostTestAndNotify_Disc;                             
                
                Results.ann_SEPCost_B_Undisc = Results.ann_NumServedbySEP_B * Params.costPP_SyringeExchange;
                Results.ann_SEPCost_H_Undisc = Results.ann_NumServedbySEP_H * Params.costPP_SyringeExchange;
                Results.ann_SEPCost_O_Undisc = Results.ann_NumServedbySEP_O * Params.costPP_SyringeExchange;
                Results.ann_SEPCost_Undisc = Results.ann_NumServedbySEP * Params.costPP_SyringeExchange;
                
                Results.ann_SEPCost_B_Disc = Results.ann_SEPCost_B_Undisc .* ann_DiscRates_Costs;
                Results.ann_SEPCost_H_Disc = Results.ann_SEPCost_H_Undisc .* ann_DiscRates_Costs;
                Results.ann_SEPCost_O_Disc = Results.ann_SEPCost_O_Undisc .* ann_DiscRates_Costs;
                Results.ann_SEPCost_Disc = Results.ann_SEPCost_Undisc .* ann_DiscRates_Costs;
                
                % Total annual undiscounted intervention spending (for YMSM
                % budget optimization analysis)
                Results.ann_TotalTransCost_Undisc = sum(Results.ann_TransCost_DiagOnPrEP_Undisc + Results.ann_TransCost_LTCFirst_Undisc + ...
                    Results.ann_TransCost_LTCAfter_Undisc + Results.ann_TransCost_ARTInitiation_Undisc + ...
                    Results.ann_TransCost_BecomeVLSfromANV_Undisc + Results.ann_TransCost_RemainVLS_Undisc + Results.ann_TransCost_CostTestAndNotify_Undisc,2);
                    
                 %<-- update when allocation outcomes finalized
    %                  % Linked First from CDC Intervention
%                     Results.ann_costLinkedFirstFromCDCFunds = ...
%                         Results.ann_numLinkedFirstFromCDCFunds * Params.costPP_LTCFirst .* (ann_DiscRates_Costs*ones(1,numStrats));
% %     
% %                 % Linked First from Non-CDC funds
% %                     Results.ann_costLinkedFirst = ...
% %                         Results.ann_numLinkedFirst * Params.costPP_LTCFirst;
% %     

%
%                 % Number adherent from CDC intervention
%                     Results.ann_costAdherentFromCDCFunds = ...
%                         Results.ann_numAdherentFromCDCFunds * Params.costPP_TxAdherence .* (ann_DiscRates_Costs*ones(1,numStrats));
%     
% %                 % Cost of adherence 
% %                     Results.ann_costAdherent = ...
% %                         Results.ann_numAdherence * Params.costPP_TxAdherence;
%     
%                  % Cost of Initiating ART from CDC inv
%                     Results.ann_costARTInitFromCDCFunds = ...
%                         Results.ann_numARTInitFromCDCFunds * Params.costPP_ARTInitiation.* (ann_DiscRates_Costs*ones(1,numStrats));
%     
%                  % Cost of Initiating ART not from CDC inv
%                     Results.ann_costARTInitFrom = ...
%                         Results.ann_numARTInit * Params.costPP_ARTInitiation;
%             
            end     
        end
        
    % 7.vi. Calculate Subpopulation-specific Results
    
        for sectionSubpop = 1:1
            
        if Params.SolveCont == 1
            totalSteps = numYears;

        else
            totalSteps = numSteps;
        end

        
        % New infections by transmission group, by time step
        for stepCount = 1:totalSteps
            Results.ts_TotalNewInfectionsHET(stepCount,:) = Results.ts_TotalNewInfections(stepCount,:) .* ...
                Params.popIndicator(:,Params.pop_HET)';
            Results.ts_TotalNewInfectionsMSM(stepCount,:) = Results.ts_TotalNewInfections(stepCount,:) .* ...
                Params.popIndicator(:,Params.pop_MSM)';
            Results.ts_TotalNewInfectionsIDU(stepCount,:) = Results.ts_TotalNewInfections(stepCount,:) .* ...
                Params.popIndicator(:,Params.pop_IDU)';

            % NewInfections by transmission group and sex
            Results.ts_TotalNewInfectionsHETM(stepCount,:) = Results.ts_TotalNewInfections(stepCount,:) .* ...
                Params.popSexIndicator(:,Params.popSex_HETM)';
            Results.ts_TotalNewInfectionsHETF(stepCount,:) = Results.ts_TotalNewInfections(stepCount,:) .* ...
                Params.popSexIndicator(:,Params.popSex_HETF)';
            Results.ts_TotalNewInfectionsMSM(stepCount,:) = Results.ts_TotalNewInfections(stepCount,:) .* ...
                Params.popSexIndicator(:,Params.popSex_MSM)';
            Results.ts_TotalNewInfectionsIDUM(stepCount,:) = Results.ts_TotalNewInfections(stepCount,:) .* ...
                Params.popSexIndicator(:,Params.popSex_IDUM)';
            Results.ts_TotalNewInfectionsIDUF(stepCount,:) = Results.ts_TotalNewInfections(stepCount,:) .* ...
                Params.popSexIndicator(:,Params.popSex_IDUF)';
        end
        
        %Chris added 2017/07/18 (IDK if this should go here)
        %Annual population in model by race/sex
        if 1==1
        end
        Results.ann_populationMaleBlk = squeeze(sum(Results.ann_PopulationSize .* ...
            repmat(Params.sexIndicator(:,Params.sex_Male)',[numYears,1]) .* ...
            repmat(Params.raceIndicator(:,Params.race_B)',[numYears,1]),2));

        Results.ann_populationMaleHisp = squeeze(sum(Results.ann_PopulationSize .* ...
            repmat(Params.sexIndicator(:,Params.sex_Male)',[numYears,1]) .* ...
            repmat(Params.raceIndicator(:,Params.race_H)',[numYears,1]),2));
        
        Results.ann_populationMaleOth = squeeze(sum(Results.ann_PopulationSize .* ...
            repmat(Params.sexIndicator(:,Params.sex_Male)',[numYears,1]) .* ...
            repmat(Params.raceIndicator(:,Params.race_O)',[numYears,1]),2));
        
        Results.ann_populationFemaleBlk = squeeze(sum(Results.ann_PopulationSize .* ...
            repmat(Params.sexIndicator(:,Params.sex_Female)',[numYears,1]) .* ...
            repmat(Params.raceIndicator(:,Params.race_B)',[numYears,1]),2));
        
        Results.ann_populationFemaleHisp = squeeze(sum(Results.ann_PopulationSize .* ...
            repmat(Params.sexIndicator(:,Params.sex_Female)',[numYears,1]) .* ...
            repmat(Params.raceIndicator(:,Params.race_H)',[numYears,1]),2));
        
        Results.ann_populationFemaleOth = squeeze(sum(Results.ann_PopulationSize .* ...
            repmat(Params.sexIndicator(:,Params.sex_Female)',[numYears,1]) .* ...
            repmat(Params.raceIndicator(:,Params.race_O)',[numYears,1]),2));
        %check equality
        %sum(Results.ann_PopulationSize,2) - ...
        %(Results.ann_populationMaleBlk + Results.ann_populationMaleHisp + Results.ann_populationMaleOth + ...
        %Results.ann_populationFemaleBlk + Results.ann_populationFemaleHisp + Results.ann_populationFemaleOth)
                
    % Total undiagnosed
    Results.ann_UndiagnosedTotal = Results.ann_HIVPrevalence - Results.ann_NumberAware;
    
    % Initialize
    Results.ann_HIVPrevalence_HET(numYears,Params.numStrats)=0;
    Results.ann_HIVPrevalence_MSM(numYears,Params.numStrats)=0;
    Results.ann_HIVPrevalence_IDU(numYears,Params.numStrats)=0;

    % New infections and continuum of care by various subpopulations
    
        % QALYs by Transmission Group, annually
        Results.ann_QALYs_pop = Results.ann_QALYs * Params.popIndicator;
        
        % Number Aware by risk group (IDUs calculated in for loop)
        Results.ann_NumberAware_LowRiskHETs = ...
                Results.ann_NumberAware * Params.LowRiskHETIndicator;
        Results.ann_NumberAware_HighRiskHETs = ...
                Results.ann_NumberAware * Params.HRHIndicator;
        Results.ann_NumberAware_LowRiskMSMs = ...
                Results.ann_NumberAware * (Params.riskLevelIndicator(:,Params.risk_Main).*Params.popIndicator(:,Params.pop_MSM));
        Results.ann_NumberAware_HighRiskMSMs = ...
                Results.ann_NumberAware * (Params.riskLevelIndicator(:,Params.risk_Casual).*Params.popIndicator(:,Params.pop_MSM));
        
        %HIV Prevalence by risk group (HETs calculated in section 7.viii  and IDUs calculated in for loop)
        Results.ann_HIVPrevalence_LowRiskMSMs = Results.ann_HIVPrevalence * ...
                (Params.riskLevelIndicator(:,Params.risk_Main).*Params.popIndicator(:,Params.pop_MSM));
        Results.ann_HIVPrevalence_HighRiskMSMs = Results.ann_HIVPrevalence * ...
                (Params.riskLevelIndicator(:,Params.risk_Casual).*Params.popIndicator(:,Params.pop_MSM));
        
        %Number VLS by risk group(IDUs calculated in for loop)
        Results.ann_NumberVLS_LowRiskHETs = ...
                Results.ann_NumberVLS * Params.LowRiskHETIndicator;   
        Results.ann_NumberVLS_HighRiskHETs = ...
                Results.ann_NumberVLS * Params.HRHIndicator;
        Results.ann_NumberVLS_LowRiskMSMs = ...
                Results.ann_NumberVLS * (Params.riskLevelIndicator(:,Params.risk_Main).*Params.popIndicator(:,Params.pop_MSM));
        Results.ann_NumberVLS_HighRiskMSMs = ...
                Results.ann_NumberVLS * (Params.riskLevelIndicator(:,Params.risk_Casual).*Params.popIndicator(:,Params.pop_MSM));
        
        
        for yearCount = 1:numYears

            % NewInfections by Transmission group, annually
        Results.ann_NewInfectionsHET(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popIndicator(:,Params.pop_HET)';
        Results.ann_NewInfectionsMSM(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popIndicator(:,Params.pop_MSM)';
        Results.ann_NewInfectionsIDU(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popIndicator(:,Params.pop_IDU)';

            % NewInfections by transmission group and sex
        Results.ann_NewInfectionsHETM(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_HETM)';
        Results.ann_NewInfectionsHETF(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_HETF)';
        Results.ann_NewInfectionsMSM(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_MSM)';
        Results.ann_NewInfectionsIDUM(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_IDUM)';
        Results.ann_NewInfectionsIDUF(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_IDUF)';
        
            % New Infections by race
        Results.ann_NewInfections_Blk(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .*...
            Params.raceIndicator(:,Params.race_B)';
        Results.ann_NewInfections_Hisp(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .*...
            Params.raceIndicator(:,Params.race_H)';
        Results.ann_NewInfections_Oth(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .*...
            Params.raceIndicator(:,Params.race_O)';
        
            % New Infections by transmission group, sex and race
            
        %Black
        Results.ann_NewInfections_Blk_HET_M(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_HETM)' .* Params.raceIndicator(:,Params.race_B)';
        Results.ann_NewInfections_Blk_HET_F(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_HETF)' .* Params.raceIndicator(:,Params.race_B)';
        Results.ann_NewInfections_Blk_MSM_M(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_MSM)' .* Params.raceIndicator(:,Params.race_B)';
        Results.ann_NewInfections_Blk_IDU_M(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_IDUM)' .* Params.raceIndicator(:,Params.race_B)';
        Results.ann_NewInfections_Blk_IDU_F(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_IDUF)' .* Params.raceIndicator(:,Params.race_B)';  
            
        %Hispanic/Latino
        Results.ann_NewInfections_Hisp_HET_M(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_HETM)' .* Params.raceIndicator(:,Params.race_H)';
        Results.ann_NewInfections_Hisp_HET_F(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_HETF)' .* Params.raceIndicator(:,Params.race_H)';
        Results.ann_NewInfections_Hisp_MSM_M(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_MSM)' .* Params.raceIndicator(:,Params.race_H)';
        Results.ann_NewInfections_Hisp_IDU_M(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_IDUM)' .* Params.raceIndicator(:,Params.race_H)';
        Results.ann_NewInfections_Hisp_IDU_F(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_IDUF)' .* Params.raceIndicator(:,Params.race_H)';  
                    
        %Other
        Results.ann_NewInfections_Oth_HET_M(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_HETM)' .* Params.raceIndicator(:,Params.race_O)';
        Results.ann_NewInfections_Oth_HET_F(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_HETF)' .* Params.raceIndicator(:,Params.race_O)';
        Results.ann_NewInfections_Oth_MSM_M(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_MSM)' .* Params.raceIndicator(:,Params.race_O)';
        Results.ann_NewInfections_Oth_IDU_M(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_IDUM)' .* Params.raceIndicator(:,Params.race_O)';
        Results.ann_NewInfections_Oth_IDU_F(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
            Params.popSexIndicator(:,Params.popSex_IDUF)' .* Params.raceIndicator(:,Params.race_O)';  


       % ANV
        Results.ann_NumberCCStage4_HET(yearCount,:) = ...
            Results.ann_NumberCCStage4(yearCount,:)' .* Params.popIndicator(:,Params.pop_HET);
        Results.ann_NumberCCStage4_MSM(yearCount,:)  = ...
            Results.ann_NumberCCStage4(yearCount,:)' .* Params.popIndicator(:,Params.pop_MSM);
        Results.ann_NumberCCStage4_IDU(yearCount,:) = ...
            Results.ann_NumberCCStage4(yearCount,:)' .* Params.popIndicator(:,Params.pop_IDU);

        % Just on VLS
        Results.ann_NumberVLS_HET(yearCount,:) = ...
            Results.ann_NumberVLS(yearCount,:)' .* Params.popIndicator(:,Params.pop_HET);
        Results.ann_NumberVLS_MSM(yearCount,:)  = ...
            Results.ann_NumberVLS(yearCount,:)' .* Params.popIndicator(:,Params.pop_MSM);
        Results.ann_NumberVLS_IDU(yearCount,:) = ...
            Results.ann_NumberVLS(yearCount,:)' .* Params.popIndicator(:,Params.pop_IDU);

        
        % Aware, by transmission group
            Results.ann_NumberAware_HET(yearCount,:) = ...
                Results.ann_NumberAware(yearCount,:)'.*Params.popIndicator(:,Params.pop_HET);
            Results.ann_NumberAware_MSM(yearCount,:) = ...
                Results.ann_NumberAware(yearCount,:)'.*Params.popIndicator(:,Params.pop_MSM);
            Results.ann_NumberAware_IDU(yearCount,:) = ...
                Results.ann_NumberAware(yearCount,:)'.*Params.popIndicator(:,Params.pop_IDU);
            
        
        % HIV Prevalence, by transmission group
            Results.ann_HIVPrevalence_HET(yearCount,:) = ...
                Results.ann_HIVPrevalence(yearCount,:)'.*Params.popIndicator(:,Params.pop_HET);
            Results.ann_HIVPrevalence_MSM(yearCount,:) = ...
                Results.ann_HIVPrevalence(yearCount,:)'.*Params.popIndicator(:,Params.pop_MSM);
            Results.ann_HIVPrevalence_IDU(yearCount,:) = ...
                Results.ann_HIVPrevalence(yearCount,:)'.*Params.popIndicator(:,Params.pop_IDU);
        
       % HIV Prevalence, by race/eth
            Results.ann_HIVPrevalence_Blk(yearCount,:) = ...
                Results.ann_HIVPrevalence(yearCount,:)' .* Params.raceIndicator(:,Params.race_B);
            Results.ann_HIVPrevalence_Hisp(yearCount,:) = ...
                Results.ann_HIVPrevalence(yearCount,:)' .* Params.raceIndicator(:,Params.race_H);
            Results.ann_HIVPrevalence_Oth(yearCount,:) = ...
                Results.ann_HIVPrevalence(yearCount,:)' .* Params.raceIndicator(:,Params.race_O);
            
     % CONTINUUM OF CARE BY RACE
            
            % INCLUSIVE CONTINUUM, BY RACE/ETH
            
                % Undiagnosed
                Results.ann_Undiagnosed_Blk(yearCount,:) = ...
                    Results.ann_UndiagnosedTotal(yearCount,:)'.* Params.raceIndicator(:,Params.race_B);
                Results.ann_Undiagnosed_Hisp(yearCount,:) = ...
                    Results.ann_UndiagnosedTotal(yearCount,:)'.* Params.raceIndicator(:,Params.race_H);
                Results.ann_Undiagnosed_Oth(yearCount,:) = ...
                    Results.ann_UndiagnosedTotal(yearCount,:)'.* Params.raceIndicator(:,Params.race_O);


                % Aware
                Results.ann_NumberAware_Blk(yearCount,:) = ...
                    Results.ann_NumberAware(yearCount,:)' .* Params.raceIndicator(:,Params.race_B);
                Results.ann_NumberAware_Hisp(yearCount,:) = ...
                    Results.ann_NumberAware(yearCount,:)' .* Params.raceIndicator(:,Params.race_H);
                Results.ann_NumberAware_Oth(yearCount,:) = ...
                    Results.ann_NumberAware(yearCount,:)' .* Params.raceIndicator(:,Params.race_O);

                % In care
                Results.ann_NumberInCare_Blk(yearCount,:) = ...
                    Results.ann_NumberInCare(yearCount,:)' .* Params.raceIndicator(:,Params.race_B);
                Results.ann_NumberInCare_Hisp(yearCount,:) = ...
                    Results.ann_NumberInCare(yearCount,:)' .* Params.raceIndicator(:,Params.race_H);
                Results.ann_NumberInCare_Oth(yearCount,:) = ...
                    Results.ann_NumberInCare(yearCount,:)' .* Params.raceIndicator(:,Params.race_O);

                % On ART
                    % Includes VLS
                 Results.ann_NumberOnART_Blk(yearCount,:) = ...
                    Results.ann_NumberOnART(yearCount,:)' .* Params.raceIndicator(:,Params.race_B);
                 Results.ann_NumberOnART_Hisp(yearCount,:) = ...
                    Results.ann_NumberOnART(yearCount,:)' .* Params.raceIndicator(:,Params.race_H);
                 Results.ann_NumberOnART_Oth(yearCount,:) = ...
                    Results.ann_NumberOnART(yearCount,:)' .* Params.raceIndicator(:,Params.race_O);

                % VLS
                Results.ann_NumberVLS_Blk(yearCount,:) = ...
                    Results.ann_NumberVLS(yearCount,:)' .*Params.raceIndicator(:,Params.race_B);
                Results.ann_NumberVLS_Hisp(yearCount,:) = ...
                    Results.ann_NumberVLS(yearCount,:)' .*Params.raceIndicator(:,Params.race_H);
                Results.ann_NumberVLS_Oth(yearCount,:) = ...
                    Results.ann_NumberVLS(yearCount,:)' .*Params.raceIndicator(:,Params.race_O);
                
                
          % MUTUALLY EXCLUSIVE CONTINNUM, BY RACE/ETH
        
                % Aware, not in care
                    % Includes both ART and non-ART affected
                Results.ann_NumberAwareNotInCare_Blk(yearCount,:) = ...
                    (Results.ann_NumberAware(yearCount,:) - Results.ann_NumberInCare(yearCount,:))' ...
                    .* Params.raceIndicator(:,Params.race_B);
                Results.ann_NumberAwareNotInCare_Hisp(yearCount,:) = ...
                    (Results.ann_NumberAware(yearCount,:) - Results.ann_NumberInCare(yearCount,:))' ...
                    .* Params.raceIndicator(:,Params.race_H);
               Results.ann_NumberAwareNotInCare_Oth(yearCount,:) = ...
                    (Results.ann_NumberAware(yearCount,:) - Results.ann_NumberInCare(yearCount,:))' ...
                    .* Params.raceIndicator(:,Params.race_O);

                % In Care, not on treatment 
                    % Includes both ART and non-ART affected
                Results.ann_NumberInCareNotART_Blk(yearCount,:) = ...
                    (Results.ann_NumberInCare(yearCount,:) - Results.ann_NumberOnART(yearCount,:))' ...
                    .* Params.raceIndicator(:,Params.race_B);
                Results.ann_NumberInCareNotART_Hisp(yearCount,:) = ...
                    (Results.ann_NumberInCare(yearCount,:) - Results.ann_NumberOnART(yearCount,:))' ...
                    .* Params.raceIndicator(:,Params.race_H);
               Results.ann_NumberInCareNotART_Oth(yearCount,:) = ...
                    (Results.ann_NumberInCare(yearCount,:) - Results.ann_NumberOnART(yearCount,:))' ...
                    .* Params.raceIndicator(:,Params.race_O);

                % On ART, not VLS
                Results.ann_NumberOnARTNotVLS_Blk(yearCount,:) = ...
                    (Results.ann_NumberOnART(yearCount,:) - Results.ann_NumberVLS(yearCount,:))' ...
                    .* Params.raceIndicator(:,Params.race_B);
                Results.ann_NumberOnARTNotVLS_Hisp(yearCount,:) = ...
                    (Results.ann_NumberOnART(yearCount,:) - Results.ann_NumberVLS(yearCount,:))' ...
                    .* Params.raceIndicator(:,Params.race_H);
               Results.ann_NumberOnARTNotVLS_Oth(yearCount,:) = ...
                    (Results.ann_NumberOnART(yearCount,:) - Results.ann_NumberVLS(yearCount,:))' ...
                    .* Params.raceIndicator(:,Params.race_O);

        % CONTINUUM BY MODEL'S CC STAGES, BY RACE/ETH
            

            % Continuum stage 2 (aware not care, not affected by ART)
            Results.ann_NumberCCStage2_Blk(yearCount,:) = ...
                Results.ann_NumberCCStage2(yearCount,:)' .* Params.raceIndicator(:,Params.race_B);
            Results.ann_NumberCCStage2_Hisp(yearCount,:) = ...
                Results.ann_NumberCCStage2(yearCount,:)' .* Params.raceIndicator(:,Params.race_H);
            Results.ann_NumberCCStage2_Oth(yearCount,:) = ...
                Results.ann_NumberCCStage2(yearCount,:)' .* Params.raceIndicator(:,Params.race_O);

            % In Care (not affected by ART)
            Results.ann_NumberCCStage3_Blk(yearCount,:) = ...
                Results.ann_NumberCCStage3(yearCount,:)' .* Params.raceIndicator(:,Params.race_B);
            Results.ann_NumberCCStage3_Hisp(yearCount,:) = ...
                Results.ann_NumberCCStage3(yearCount,:)' .* Params.raceIndicator(:,Params.race_H);
            Results.ann_NumberCCStage3_Oth(yearCount,:) = ...
                Results.ann_NumberCCStage3(yearCount,:)' .* Params.raceIndicator(:,Params.race_O);     

            % Affected by ART - not VLS
            Results.ann_NumberCCStage4_Blk(yearCount,:) = ...
                Results.ann_NumberCCStage4(yearCount,:)' .* Params.raceIndicator(:,Params.race_B);
            Results.ann_NumberCCStage4_Hisp(yearCount,:) = ...
                Results.ann_NumberCCStage4(yearCount,:)' .* Params.raceIndicator(:,Params.race_H);
            Results.ann_NumberCCStage4_Oth(yearCount,:) = ...
                Results.ann_NumberCCStage4(yearCount,:)' .* Params.raceIndicator(:,Params.race_O);

            % Affected by ART    
            Results.ann_NumberCCStages45_Blk(yearCount,:) = ...
                Results.ann_NumberCCStages45(yearCount,:)' .* Params.raceIndicator(:,Params.race_B);
            Results.ann_NumberCCStages45_Hisp(yearCount,:) = ...
                Results.ann_NumberCCStages45(yearCount,:)' .* Params.raceIndicator(:,Params.race_H);
            Results.ann_NumberCCStages45_Oth(yearCount,:) = ...
                Results.ann_NumberCCStages45(yearCount,:)' .* Params.raceIndicator(:,Params.race_O);
            
            % Note: VLS is already recorded by race/eth above.

        end
       
        end
        
    % 7.vii. Calculate Calibration Outcomes   
        for sectionCalib = 1:1
            
        % Gather new diagnoses per year for calibration results
         %added by age - 6/20/17. Laurel    
            for yearCount = 1:numYears
                
                % by Transmission group
                Results.ann_TotalNewDiagnosesHET(yearCount,:) = Results.ann_TotalNewDiagnoses(yearCount,:) .* ...
                    Params.popIndicator(:,Params.pop_HET)';
                Results.ann_TotalNewDiagnosesMSM(yearCount,:) = Results.ann_TotalNewDiagnoses(yearCount,:) .* ...
                    Params.popIndicator(:,Params.pop_MSM)';
                Results.ann_TotalNewDiagnosesIDU(yearCount,:) = Results.ann_TotalNewDiagnoses(yearCount,:) .* ...
                    Params.popIndicator(:,Params.pop_IDU)';
                
                Results.ann_TotalNewDiagnoses_BlackMSM(yearCount,:) = Results.ann_TotalNewDiagnoses(yearCount,:) .* ...
                    Params.popIndicator(:,Params.pop_MSM)' .* ...
                    Params.raceIndicator(:,Params.race_B)';

                % by Transmission group and sex
                Results.ann_TotalNewDiagnosesHETM(yearCount,:) = Results.ann_TotalNewDiagnoses(yearCount,:) .* ...
                    Params.popSexIndicator(:,Params.popSex_HETM)';
                Results.ann_TotalNewDiagnosesHETF(yearCount,:) = Results.ann_TotalNewDiagnoses(yearCount,:) .* ...
                    Params.popSexIndicator(:,Params.popSex_HETF)';
                Results.ann_TotalNewDiagnosesMSM(yearCount,:) = Results.ann_TotalNewDiagnoses(yearCount,:) .* ...
                    Params.popSexIndicator(:,Params.popSex_MSM)';
                Results.ann_TotalNewDiagnosesIDUM(yearCount,:) = Results.ann_TotalNewDiagnoses(yearCount,:) .* ...
                    Params.popSexIndicator(:,Params.popSex_IDUM)';
                Results.ann_TotalNewDiagnosesIDUF(yearCount,:) = Results.ann_TotalNewDiagnoses(yearCount,:) .* ...
                    Params.popSexIndicator(:,Params.popSex_IDUF)';
            
                % by Race/Ethnicity
                Results.ann_TotalNewDiagnoses_Blk(yearCount,:) = ...
                    Results.ann_TotalNewDiagnoses(yearCount,:) .* Params.raceIndicator(:,Params.race_B)';
                Results.ann_TotalNewDiagnoses_Hisp(yearCount,:) = ...
                   Results.ann_TotalNewDiagnoses(yearCount,:) .* Params.raceIndicator(:,Params.race_H)';
                Results.ann_TotalNewDiagnoses_Oth(yearCount,:) = ...
                   Results.ann_TotalNewDiagnoses(yearCount,:) .* Params.raceIndicator(:,Params.race_O)';
               
                % by Age (Need to account for 13-24, 18-24, and 25-34)
                Results.ann_TotalNewDiagnoses_13_17(yearCount,:) = ...
                    Results.ann_TotalNewDiagnoses(yearCount,:) .* Params.ageIndicator(:,Params.age_1)';
                Results.ann_TotalNewDiagnoses_18_24(yearCount,:) = ...
                   Results.ann_TotalNewDiagnoses(yearCount,:) .* Params.ageIndicator(:,Params.age_2)';
                Results.ann_TotalNewDiagnoses_25_34(yearCount,:) = ...
                   Results.ann_TotalNewDiagnoses(yearCount,:) .* Params.ageIndicator(:,Params.age_3)';

            end    
                    
            %NEW CALIB RESULTS START HERE. 
            %Changed start year to 2010 instead of 2006. Bates. 6/20/19
        % 2010 outcomes - used for calculating other targets only
        if FirstOutcomeYr <= 2010 && (LastOutcomeYr >= 2010)    
             
            % HIV Prevalence in 2010
            calib_HIVPrevalence2010 = sum(Results.ann_HIVPrevalence(Params.year_2010,:));
           
            %HIV Prevalence for LR and HR HETS in 2010
            calib_LRHETPrevalence2010 = Results.ann_HIVPrevalence_HET(Params.year_2010,:)*...
                Params.riskLevelIndicator(:,Params.risk_Main);

            calib_HRHETPrevalence2010 = Results.ann_HIVPrevalence_HET(Params.year_2010,:)*...
                Params.riskLevelIndicator(:,Params.risk_Casual);         
            
        end
        %Removed 24 2010 targets. Bates 6/11/19.
        %Removed 2012 targets and values. Bates 6/19/16

       % 2014 outcomes Moved 2014 calibration targets to 2019. Bates
       % 2/10/23
        
       % 2015 outcomes
        if FirstOutcomeYr <= 2015 && (LastOutcomeYr >= 2015)    
            %Many 2012 values moved to 2015 values. 4/7/16. Laurel Bates
            
        %Number aware is needed for both 2015 and 2016, as diagnosed is
        %2016 and %VLS among diagnosed is 2015. 
            numaware2015_13_17 = Results.ann_NumberAware(Params.year_2015,:)*...
                Params.ageIndicator(:,Params.age_1);
            numaware2015_18_24 = Results.ann_NumberAware(Params.year_2015,:)*...
                Params.ageIndicator(:,Params.age_2);     
            
            
            numaware2015_25_34 = Results.ann_NumberAware(Params.year_2015,:)*...
                Params.ageIndicator(:,Params.age_3);
            numaware2015_35_44 = Results.ann_NumberAware(Params.year_2015,:)*...
                Params.ageIndicator(:,Params.age_4);
            numaware2015_45_54 = Results.ann_NumberAware(Params.year_2015,:)*...
                Params.ageIndicator(:,Params.age_5);
            numaware2015_55_64 = Results.ann_NumberAware(Params.year_2015,:)*...
                Params.ageIndicator(:,Params.age_6);
            numaware2015_65 = Results.ann_NumberAware(Params.year_2015,:)*...
                Params.ageIndicator(:,Params.age_7);
         
     
        % Continuum of Care in 2015
            % By race
            
            %Unaware-removed 6/13/19. Bates

            %VLS among Diagnosed- by race moved to 2017 4/20/20
           
            Results.calib_TTdist2015_VLSamongdiag_13_17 = ...
                Results.ann_NumberVLS(Params.year_2015,:)*Params.ageIndicator(:,Params.age_1)...
                /numaware2015_13_17;
            Results.calib_TTdist2015_VLSamongdiag_18_24 = ...
                Results.ann_NumberVLS(Params.year_2015,:)*Params.ageIndicator(:,Params.age_2)...
                /numaware2015_18_24;
            Results.calib_TTdist2015_VLSamongdiag_25_34 = ...
                Results.ann_NumberVLS(Params.year_2015,:)*Params.ageIndicator(:,Params.age_3)...
                /numaware2015_25_34;
            Results.calib_TTdist2015_VLSamongdiag_35_44 = ...
                Results.ann_NumberVLS(Params.year_2015,:)*Params.ageIndicator(:,Params.age_4)...
                /numaware2015_35_44;
            Results.calib_TTdist2015_VLSamongdiag_45_54 = ...
                Results.ann_NumberVLS(Params.year_2015,:)*Params.ageIndicator(:,Params.age_5)...
                /numaware2015_45_54;
            Results.calib_TTdist2015_VLSamongdiag_55_64 = ...
                Results.ann_NumberVLS(Params.year_2015,:)*Params.ageIndicator(:,Params.age_6)...
                /numaware2015_55_64;
            Results.calib_TTdist2015_VLSamongdiag_65 = ...
                Results.ann_NumberVLS(Params.year_2015,:)*Params.ageIndicator(:,Params.age_7)...
                /numaware2015_65;
            
            Results.calib_TTdist2015_VLSamongdiag_13_24 = ...
                (Results.ann_NumberVLS(Params.year_2015,:)*Params.ageIndicator(:,Params.age_1)...
                + Results.ann_NumberVLS(Params.year_2015,:)*Params.ageIndicator(:,Params.age_2))...
                /(numaware2015_13_17 + numaware2015_18_24);
            
            %VLS among Diagnosed for young MSM - added by JC on 8/15/2017.
            %Moved to 2016 by Laurel 6/19/19
            %VLS total among diagnosed total- Moved to 2017 4/20/20
      
                   % PLWH Aware Deaths - moved to 2016
            
        end  
        
        % 2016 outcomes
        %New infections among YMSM groups moved to 2016 6/20/19. Laurel

        %New infections total and among older age groups moved to 2016.
        %6/20/19 Laurel

               % New Diagnoses in 2015. Removed outcomes by race
               % 6/17/19. Moved outcomes by transmission group and
               % total to 2016

               % New Infections in 2015- moved to 2016 6/19/19 Bates

        %HIV Prevalence for LR and HR HETS in 2015 moved to 2016. Bates 6/20/19  
        
        if FirstOutcomeYr <= 2016 && (LastOutcomeYr >= 2016)
        % Age / Race population values

        
        %Prevalence overall and by race needed for 2016 to calculate other
        %calibration values
          Results.calib_HIVPrevalence2016 = sum(Results.ann_HIVPrevalence(Params.year_2016,:));
 
          Results.calib_HIVPrevalence2016_B = Results.ann_HIVPrevalence(Params.year_2016,:)*...
                Params.raceIndicator(:,Params.race_B);
          Results.calib_HIVPrevalence2016_H = Results.ann_HIVPrevalence(Params.year_2016,:)*...
                Params.raceIndicator(:,Params.race_H);
          Results.calib_HIVPrevalence2016_O = Results.ann_HIVPrevalence(Params.year_2016,:)*...
                Params.raceIndicator(:,Params.race_O); 
            
                            
          %Prevalence by age- Updated to 2016 6/13/19, Bates
          Results.calib_HIVPrevalence2016_13_17 = Results.ann_HIVPrevalence(Params.year_2016,:)*...
                Params.ageIndicator(:,Params.age_1);
          Results.calib_HIVPrevalence2016_18_24 = Results.ann_HIVPrevalence(Params.year_2016,:)*...
                Params.ageIndicator(:,Params.age_2);
          Results.calib_HIVPrevalence2016_25_34 = Results.ann_HIVPrevalence(Params.year_2016,:)*...
                Params.ageIndicator(:,Params.age_3);
          Results.calib_HIVPrevalence2016_35_44 = Results.ann_HIVPrevalence(Params.year_2016,:)*...
                Params.ageIndicator(:,Params.age_4);  
          Results.calib_HIVPrevalence2016_45_54 = Results.ann_HIVPrevalence(Params.year_2016,:)*...
                Params.ageIndicator(:,Params.age_5);  
          Results.calib_HIVPrevalence2016_55_64 = Results.ann_HIVPrevalence(Params.year_2016,:)*...
                Params.ageIndicator(:,Params.age_6);  
          Results.calib_HIVPrevalence2016_65 = Results.ann_HIVPrevalence(Params.year_2016,:)*...
                Params.ageIndicator(:,Params.age_7);  
            
          Results.calib_HIVPrevalence2016_13_24 = Results.calib_HIVPrevalence2016_13_17 + Results.calib_HIVPrevalence2016_18_24; 
          
          %Prevalence by race and age. Added 6/18/19. Laurel Bates
   
          Results.calib_HIVPrevalence2016_B_13_24 = sum(Results.ann_HIVPrevalence(Params.year_2016,:)'.* ...
                 Params.raceIndicator(:,Params.race_B).*...
                 sum(Params.ageIndicator(:,[Params.age_1,Params.age_2]),2));
       
          Results.calib_HIVPrevalence2016_H_13_24 = sum(Results.ann_HIVPrevalence(Params.year_2016,:)'.* ...
                 Params.raceIndicator(:,Params.race_H).*...
                 sum(Params.ageIndicator(:,[Params.age_1,Params.age_2]),2));
             
          Results.calib_HIVPrevalence2016_O_13_24 = sum(Results.ann_HIVPrevalence(Params.year_2016,:)'.* ...
                 Params.raceIndicator(:,Params.race_O).*...
                 sum(Params.ageIndicator(:,[Params.age_1,Params.age_2]),2));    
             
          Results.calib_HIVPrevalence2016_B_25_34 = sum(Results.ann_HIVPrevalence(Params.year_2016,:)'.* ...
                 Params.raceIndicator(:,Params.race_B).*...
                 Params.ageIndicator(:,Params.age_3));
            
          Results.calib_HIVPrevalence2016_H_25_34 = sum(Results.ann_HIVPrevalence(Params.year_2016,:)'.* ...
                 Params.raceIndicator(:,Params.race_H).*...
                 Params.ageIndicator(:,Params.age_3));
             
          Results.calib_HIVPrevalence2016_O_25_34 = sum(Results.ann_HIVPrevalence(Params.year_2016,:)'.* ...
                 Params.raceIndicator(:,Params.race_O).*...
                 Params.ageIndicator(:,Params.age_3)); 
       
            %Aware by age -  Added 2016. Bates. 6/17/19
            numaware2016_13_17 = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_1));
            numaware2016_18_24 = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_2));     
            
            numaware2016_13_24 = numaware2016_13_17 + numaware2016_18_24;
            
            numaware2016_25_34 = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_3));
            numaware2016_35_44 = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_4));
            numaware2016_45_54 = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_5));
            numaware2016_55_64 = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_6));
            numaware2016_65 = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_7));
            
            %Aware YMSM- Added 2016. Bates. 6/17/19
            
            numaware2016_13_17_MSM_B = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_1) ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator (:, Params.race_B));
            numaware2016_13_17_MSM_H = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_1) ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator (:, Params.race_H));
             numaware2016_13_17_MSM_O = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_1) ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator (:, Params.race_O));
            
           numaware2016_18_24_MSM_B = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_2) ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator (:, Params.race_B));
            numaware2016_18_24_MSM_H = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_2) ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator (:, Params.race_H));
             numaware2016_18_24_MSM_O = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_2) ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator (:, Params.race_O));
            
            numaware2016_13_24_MSM_B = numaware2016_13_17_MSM_B + numaware2016_18_24_MSM_B;
            numaware2016_13_24_MSM_H = numaware2016_13_17_MSM_H + numaware2016_18_24_MSM_H;
            numaware2016_13_24_MSM_O = numaware2016_13_17_MSM_O + numaware2016_18_24_MSM_O;
            
            
            numaware2016_25_34_MSM_B = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_3) ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator (:, Params.race_B));
            numaware2016_25_34_MSM_H = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_3) ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator (:, Params.race_H));
             numaware2016_25_34_MSM_O = sum(Results.ann_NumberAware(Params.year_2016,:)'.*...
                Params.ageIndicator(:,Params.age_3) ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator (:, Params.race_O));
            
           % Diagnosed targets for percent aware for MSM by race and age.
           % [This code used to work, but had to re-write when we changed
           % age groups from 18-24 to 13-24.]
            % (cats 2 and 3)
            %Results.calib_TTdist2015_Diagnosed_MSM_B_18_24 = percentAware_pop_race_age('MSM','B','2');
            %Results.calib_TTdist2015_Diagnosed_MSM_H_18_24 = percentAware_pop_race_age('MSM','H','2');
           % Results.calib_TTdist2015_Diagnosed_MSM_O_18_24 = percentAware_pop_race_age('MSM','O','2');
           % Results.calib_TTdist2015_Diagnosed_MSM_B_25_34 = percentAware_pop_race_age('MSM','B','3');
           % Results.calib_TTdist2015_Diagnosed_MSM_H_25_34 = percentAware_pop_race_age('MSM','H','3');
           % Results.calib_TTdist2015_Diagnosed_MSM_O_25_34 = percentAware_pop_race_age('MSM','O','3');
        
           %Diagnosed for MSM by race and age. Updated to 2016 6/17/19.
           %Bates
          % Changed 18-24 age group to 13-24 
            Results.calib_TTdist2016_Diagnosed_MSM_B_13_24 = ...
              numaware2016_13_24_MSM_B...
               / ...
              sum(Results.ann_HIVPrevalence(Params.year_2016,:)' ...
                .* Params.popIndicator(:,Params.pop_MSM)...
                .* Params.raceIndicator(:,Params.race_B)...
                .* sum(Params.ageIndicator(:,[Params.age_1,Params.age_2]),2));
         
           
            Results.calib_TTdist2016_Diagnosed_MSM_H_13_24 = ...
              numaware2016_13_24_MSM_H...
               / ...
              sum(Results.ann_HIVPrevalence(Params.year_2016,:)' ...
                .* Params.popIndicator(:,Params.pop_MSM)...
                .* Params.raceIndicator(:,Params.race_H)...
                .* sum(Params.ageIndicator(:,[Params.age_1,Params.age_2]),2));
           
           Results.calib_TTdist2016_Diagnosed_MSM_O_13_24 = ...
              numaware2016_13_24_MSM_O...
               / ...
              sum(Results.ann_HIVPrevalence(Params.year_2016,:)' ...
                .* Params.popIndicator(:,Params.pop_MSM)...
                .* Params.raceIndicator(:,Params.race_O)...
                .* sum(Params.ageIndicator(:,[Params.age_1,Params.age_2]),2));
           
           Results.calib_TTdist2016_Diagnosed_MSM_B_25_34 = ...
              numaware2016_25_34_MSM_B...
              /...
              sum(Results.ann_HIVPrevalence(Params.year_2016,:)' ...
                .* Params.popIndicator(:,Params.pop_MSM)...
                .* Params.raceIndicator(:,Params.race_B)...
                .* Params.ageIndicator(:,Params.age_3));
           
            Results.calib_TTdist2016_Diagnosed_MSM_H_25_34 = ...
              numaware2016_25_34_MSM_H...
              /...
              sum(Results.ann_HIVPrevalence(Params.year_2016,:)' ...
                .* Params.popIndicator(:,Params.pop_MSM)...
                .* Params.raceIndicator(:,Params.race_H)...
                .* Params.ageIndicator(:,Params.age_3));
           
           Results.calib_TTdist2016_Diagnosed_MSM_O_25_34 = ...
              numaware2016_25_34_MSM_O...
              /...
              sum(Results.ann_HIVPrevalence(Params.year_2016,:)' ...
                .* Params.popIndicator(:,Params.pop_MSM)...
                .* Params.raceIndicator(:,Params.race_O)...
                .* Params.ageIndicator(:,Params.age_3));
           
        %Diagnosed
            %By race-updated to 2016. Bates. 6/14/19. Moved to 2019, Bates,
            %11/22/22

          %Diagnosed by transmission group. Updated 6/17/19- Moved to 2018
          %4/21/20
          
            %Diagnosed by age          
             Results.calib_TTdist2016_Diagnosed_13_24 = ...
              numaware2016_13_24...
               /Results.calib_HIVPrevalence2016_13_24;   
           
           %Results.calib_TTdist2016_Diagnosed_13_24 = ...
              % sum(Results.ann_HIVPrevalence(Params.year_2016,:)'.*...
               % Params.ageIndicator(:,Params.age_1));
            
            Results.calib_TTdist2016_Diagnosed_25_34 = ...
                numaware2016_25_34...
                /Results.calib_HIVPrevalence2016_25_34;
            
            Results.calib_TTdist2016_Diagnosed_35_44 = ...
                numaware2016_35_44...
                /Results.calib_HIVPrevalence2016_35_44;
            
            Results.calib_TTdist2016_Diagnosed_45_54 = ...
                numaware2016_45_54...
                /Results.calib_HIVPrevalence2016_45_54;
            
            Results.calib_TTdist2016_Diagnosed_55_64 = ...
                numaware2016_55_64...
                /Results.calib_HIVPrevalence2016_55_64;
            
            Results.calib_TTdist2016_Diagnosed_65 = ...
                numaware2016_65...
                /Results.calib_HIVPrevalence2016_65;
            
            %Diagnosed overall- moved to 2018

            %In Care - Updated to 2016 (Allaire 6-14-2019)
            %Removed. Bates 6/17/19. Moved to 2019 target. Bates, 11/4/22
            
            %No Care - Updated prevalence to 2016 (Allaire 6-14-2019)
            %Removed. Bates 6/17/19

            %On ART, Not VLS- removed 4/7/17 Laurel Bates
            
            %VLS among PLWH-  removed 6/19/19 Laurel Bates 
            
            %VLS among Diagnosed for young MSM - added by JC on 8/15/2017.
            %Moved to 2016 by Laurel 6/19/19. Changed 18-24 age group to 13-24. 
               
            Results.calib_TTdist2016_VLSamongdiag_13_24_MSM_B = ...
                sum(Results.ann_NumberVLS(Params.year_2016,:)'.*...
                Params.popSexIndicator(:,Params.popSex_MSM).*...
                Params.raceIndicator(:, Params.race_B).*...
            sum(Params.ageIndicator(:,[Params.age_1, Params.age_2]),2))...
                /numaware2016_13_24_MSM_B;
            
           Results.calib_TTdist2016_VLSamongdiag_13_24_MSM_H = ...
                sum(Results.ann_NumberVLS(Params.year_2016,:)'.*...
                Params.popSexIndicator(:,Params.popSex_MSM).*...
                Params.raceIndicator(:, Params.race_H).*...
            sum(Params.ageIndicator(:,[Params.age_1, Params.age_2]),2))...
                /numaware2016_13_24_MSM_H;  
            
           Results.calib_TTdist2016_VLSamongdiag_13_24_MSM_O = ...
                sum(Results.ann_NumberVLS(Params.year_2016,:)'.*...
                Params.popSexIndicator(:,Params.popSex_MSM).*...
                Params.raceIndicator(:, Params.race_O).*...
            sum(Params.ageIndicator(:,[Params.age_1, Params.age_2]),2))...
                /numaware2016_13_24_MSM_O;             
           
            Results.calib_TTdist2016_VLSamongdiag_25_34_MSM_B = ...
                sum(Results.ann_NumberVLS(Params.year_2016,:)'.*Params.ageIndicator(:,Params.age_3).*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator(:, Params.race_B))...
                /numaware2016_25_34_MSM_B;
            Results.calib_TTdist2016_VLSamongdiag_25_34_MSM_H = ...
                sum(Results.ann_NumberVLS(Params.year_2016,:)'.*Params.ageIndicator(:,Params.age_3).*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator(:, Params.race_H))...
                /numaware2016_25_34_MSM_H;
            Results.calib_TTdist2016_VLSamongdiag_25_34_MSM_O = ...
                sum(Results.ann_NumberVLS(Params.year_2016,:)'.*Params.ageIndicator(:,Params.age_3).*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator(:, Params.race_O))...
                /numaware2016_25_34_MSM_O;
            
            
            Results.calib_AwarePLWHDeaths2016_45_54 = sum(Results.ann_numDeathsPLWHAware(Params.year_2016,:)' ...                
                .*Params.ageIndicator (:,Params.age_5));
            Results.calib_AwarePLWHDeaths2016_55_64 = sum(Results.ann_numDeathsPLWHAware(Params.year_2016,:)' ...                
                .*Params.ageIndicator (:,Params.age_6));
            Results.calib_AwarePLWHDeaths2016_65 = sum(Results.ann_numDeathsPLWHAware(Params.year_2016,:)' ...                
                .*Params.ageIndicator (:,Params.age_7));
          
            %Aware by age: Note that the 13-24 group is used here, a
            %combination of the 13-17 and 18-24. Since they are both
            %defined above, this should be OK
            

            %New diagnoses, 2016. Moved 6/19/19 Laurel
            
            Results.calib_NewDiagnoses2016_HETM = sum(Results.ann_TotalNewDiagnosesHETM(Params.year_2016,:));
                       
            Results.calib_NewDiagnoses2016_HETF = sum(Results.ann_TotalNewDiagnosesHETF(Params.year_2016,:));
            
            Results.calib_NewDiagnoses2016_MSM = sum(Results.ann_TotalNewDiagnosesMSM(Params.year_2016,:));
            
            Results.calib_NewDiagnoses2016_IDUM = sum(Results.ann_TotalNewDiagnosesIDUM(Params.year_2016,:));
            
            Results.calib_NewDiagnoses2016_IDUF = sum(Results.ann_TotalNewDiagnosesIDUF(Params.year_2016,:));
            
            Results.calib_NewDiagnoses2016_Total = sum(Results.ann_TotalNewDiagnoses(Params.year_2016,:));  
               
            % New Infections in 2015- moved to 2016 6/19/19 Bates. Moved to
            % 2018 4/21/19 Bates

            %New infections among YMSM- changed to 2016. 6/20/19 Bates. Changed 18-24 age group to 13-24   
             Results.calib_NewInfections_2016_MSM_13_24 = ...
                sum(Results.ann_TotalNewInfections(Params.year_2016,:)' ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*sum(Params.ageIndicator(:,[Params.age_1, Params.age_2]),2));
            
            Results.calib_NewInfections_2016_MSM_13_24_B = ...
                sum(Results.ann_TotalNewInfections(Params.year_2016,:)' ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator (:, Params.race_B)...
                .*sum(Params.ageIndicator(:,[Params.age_1, Params.age_2]),2));
            
            Results.calib_NewInfections_2016_MSM_13_24_H = ...
                sum(Results.ann_TotalNewInfections(Params.year_2016,:)' ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.ageIndicator (:,Params.age_2) ...
                .*Params.raceIndicator (:, Params.race_H)...
                .*sum(Params.ageIndicator(:,[Params.age_1, Params.age_2]),2));
            
            Results.calib_NewInfections_2016_MSM_13_24_O = ...
                sum(Results.ann_TotalNewInfections(Params.year_2016,:)' ...
                .*Params.popSexIndicator(:,Params.popSex_MSM) ...
                .*Params.ageIndicator (:,Params.age_2) ...
                .*Params.raceIndicator (:, Params.race_O)...
                .*sum(Params.ageIndicator(:,[Params.age_1, Params.age_2]),2));
            
            Results.calib_NewInfections_2016_MSM_25_34 = ...
                sum(Results.ann_TotalNewInfections(Params.year_2016,:)' ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.ageIndicator (:,Params.age_3));
            
            Results.calib_NewInfections_2016_MSM_25_34_B = ...
                sum(Results.ann_TotalNewInfections(Params.year_2016,:)' ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.ageIndicator (:,Params.age_3)... 
                .*Params.raceIndicator (:,Params.race_B));
            
            Results.calib_NewInfections_2016_MSM_25_34_H = ...
                sum(Results.ann_TotalNewInfections(Params.year_2016,:)' ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.ageIndicator (:,Params.age_3)... 
                .*Params.raceIndicator (:,Params.race_H));
            
            Results.calib_NewInfections_2016_MSM_25_34_O = ...
                sum(Results.ann_TotalNewInfections(Params.year_2016,:)' ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.ageIndicator (:,Params.age_3)... 
                .*Params.raceIndicator (:,Params.race_O));
            
            %New infections among older age groups moved 6/20/19
            
            Results.calib_NewInfections2016_45_54 = ...
                sum(Results.ann_TotalNewInfections(Params.year_2016,:)' ...                
                .*Params.ageIndicator (:,Params.age_5));
            
            Results.calib_NewInfections2016_55_64 = ...
                sum(Results.ann_TotalNewInfections(Params.year_2016,:)' ...                
                .*Params.ageIndicator (:,Params.age_6));
            
            Results.calib_NewInfections2016_65 = ...
                sum(Results.ann_TotalNewInfections(Params.year_2016,:)' ...                
                .*Params.ageIndicator (:,Params.age_7));
            
            %Population values 2016- Updated to 2019 11/22/22, Bates
              
            % populations by age values- moved to 2018 

        end    
        
            %2017 targets- Added 4/20/20 then moved to 2018 5/29/20
                  
        %2018 targets- some 2018 targets were being collected in 2019
        %section - JC moved here 0n 2021/08/02
        if FirstOutcomeYr <= 2018 && (LastOutcomeYr >= 2018)
            
            % Calculate number aware for the denominators            
             numaware2018_HET = sum(Results.ann_NumberAware(Params.year_2018,:)'.*...
                Params.popIndicator(:,Params.pop_HET));                         
            numaware2018_MSM = sum(Results.ann_NumberAware(Params.year_2018,:)'.*...
                Params.popIndicator(:,Params.pop_MSM));                         
            numaware2018_IDU = sum(Results.ann_NumberAware(Params.year_2018,:)'.*...
                Params.popIndicator(:,Params.pop_IDU));           
              
            %VLS by race (moved from 2015), transmission group (added
            %4/20/20), and total (moved from 2015). Moved by race and
            %overall to 2019 11/22/22
            
            Results.calib_TTdist2018_VLSamongdiag_HET = ...
                sum(Results.ann_NumberVLS(Params.year_2018,:)'.*Params.popIndicator(:,Params.pop_HET))...
                /numaware2018_HET;
            Results.calib_TTdist2018_VLSamongdiag_MSM = ...
                sum(Results.ann_NumberVLS(Params.year_2018,:)'.*Params.popIndicator(:,Params.pop_MSM))...
                /numaware2018_MSM;
            Results.calib_TTdist2018_VLSamongdiag_PWID = ...
                sum(Results.ann_NumberVLS(Params.year_2018,:)'.*Params.popIndicator(:,Params.pop_IDU))...
                /numaware2018_IDU;
            %Moved number on PrEP to 2021. Bates 02/10/2023
        end
            
        %2019 targets- Added 4/20/20 Bates
        if FirstOutcomeYr <= 2019 && (LastOutcomeYr >= 2019)

        % Calculate number aware for the denominators: LB added 11/17/22
            numaware2019_B = sum(Results.ann_NumberAware(Params.year_2019,:)'.*...
                Params.raceIndicator(:,Params.race_B));                         
            numaware2019_H = sum(Results.ann_NumberAware(Params.year_2019,:)'.*...
                Params.raceIndicator(:,Params.race_H));                         
            numaware2019_O = sum(Results.ann_NumberAware(Params.year_2019,:)'.*...
                Params.raceIndicator(:,Params.race_O));
            
            % HIV Prevalence: updated to 2018, 4/21/20 Bates
            % Prevalence updated to 2019. Clinkscales 6/29/2021
            Results.calib_HIVPrevalence2019 = sum(Results.ann_HIVPrevalence(Params.year_2019,:));
            
         % By race updated to 2018, 4/21/20. Values updated to 2019,
         % 6/29/2021
          Results.calib_HIVPrevalence2019_B = Results.ann_HIVPrevalence(Params.year_2019,:)*...
                Params.raceIndicator(:,Params.race_B);
          Results.calib_HIVPrevalence2019_H = Results.ann_HIVPrevalence(Params.year_2019,:)*...
                Params.raceIndicator(:,Params.race_H);
          Results.calib_HIVPrevalence2019_O = Results.ann_HIVPrevalence(Params.year_2019,:)*...
                Params.raceIndicator(:,Params.race_O);   
            
            
            %Prevalence by transmission group- Added 4/20/20. Values
            %updated to 2019. Clinkscales 6/29/2021
          Results.calib_HIVPrevalence2019_HETM = Results.ann_HIVPrevalence(Params.year_2019,:)*...
              Params.popSexIndicator(:,Params.popSex_HETM);
          Results.calib_HIVPrevalence2019_HETF = Results.ann_HIVPrevalence(Params.year_2019,:)*...
              Params.popSexIndicator(:,Params.popSex_HETF);
           Results.calib_HIVPrevalence2019_MSM = Results.ann_HIVPrevalence(Params.year_2019,:)*...
              Params.popSexIndicator(:,Params.popSex_MSM);
          Results.calib_HIVPrevalence2019_IDUM = Results.ann_HIVPrevalence(Params.year_2019,:)*...
              Params.popSexIndicator(:,Params.popSex_IDUM);
          Results.calib_HIVPrevalence2019_IDUF = Results.ann_HIVPrevalence(Params.year_2019,:)*...
              Params.popSexIndicator(:,Params.popSex_IDUF);     
          
          %Diagnosed overall- Moved from 2016. Updated to 2019 - Clinkcales
          %6/30/2021
          Results.calib_TTdist2019_Diagnosed_Total = ...
               sum(Results.ann_NumberAware(Params.year_2019,:)) ...
                /sum(Results.ann_HIVPrevalence(Params.year_2019,:));
            
          %By race-Moved to 2019, Bates, 11/22/22
            Results.calib_TTdist2019_Diagnosed_B = ...
                sum(Results.ann_NumberAware(Params.year_2019,:)'.*Params.raceIndicator(:,Params.race_B))...
                /Results.calib_HIVPrevalence2019_B;
            Results.calib_TTdist2019_Diagnosed_H = ...
                sum(Results.ann_NumberAware(Params.year_2019,:)'.*Params.raceIndicator(:,Params.race_H))...
                /Results.calib_HIVPrevalence2019_H;
            Results.calib_TTdist2019_Diagnosed_O = ...
                sum(Results.ann_NumberAware(Params.year_2019,:)'.*Params.raceIndicator(:,Params.race_O))...
                /Results.calib_HIVPrevalence2019_O;   
          %Diagnosed by transmission group- Moved from 2016 4/21/20.
          %Updated to 2019 - Clinkscales 6/30/2021
          
             Results.calib_TTdist2019_Diagnosed_HET = ...
                (Results.ann_NumberAware(Params.year_2019,:) ...
                    * Params.popIndicator(:,Params.pop_HET)) ...
                / ...
                (Results.ann_HIVPrevalence(Params.year_2019,:) ...
                    * Params.popIndicator(:,Params.pop_HET));
            
            Results.calib_TTdist2019_Diagnosed_MSM = ...
                (Results.ann_NumberAware(Params.year_2019,:) ...
                    * Params.popIndicator(:,Params.pop_MSM)) ...
                / ...
                (Results.ann_HIVPrevalence(Params.year_2019,:) ...
                    * Params.popIndicator(:,Params.pop_MSM));
            
            Results.calib_TTdist2019_Diagnosed_PWID = ...
                (Results.ann_NumberAware(Params.year_2019,:) ...
                    * Params.popIndicator(:,Params.pop_IDU)) ...
                / ...
                (Results.ann_HIVPrevalence(Params.year_2019,:) ...
                    * Params.popIndicator(:,Params.pop_IDU));
            
          %In care by race/ethnicity and overall. Added to model 11/4/22 after previously removed. Bates   

            Results.calib_TTdist2019_InCareAmongDiag_B = ...
                sum(Results.ann_NumberInCare(Params.year_2019,:)'.*Params.raceIndicator(:,Params.race_B))...
                /numaware2019_B;
            Results.calib_TTdist2019_InCareAmongDiag_H = ...
                sum(Results.ann_NumberInCare(Params.year_2019,:)'.*Params.raceIndicator(:,Params.race_H))...
                /numaware2019_H;
            Results.calib_TTdist2019_InCareAmongDiag_O = ...
                sum(Results.ann_NumberInCare(Params.year_2019,:)'.*Params.raceIndicator(:,Params.race_O))...
                /numaware2019_O;
            Results.calib_TTdist2019_InCareAmongDiag_Total = ...
                sum(Results.ann_NumberInCare(Params.year_2019,:))...
                /sum(Results.ann_NumberAware(Params.year_2019,:));

           %VLS by race and overall (moved from 2018) 11/22/22
            Results.calib_TTdist2019_VLSamongdiag_B = ...
                sum(Results.ann_NumberVLS(Params.year_2019,:)'.*Params.raceIndicator(:,Params.race_B))...
                /numaware2019_B;
            Results.calib_TTdist2019_VLSamongdiag_H = ...
                sum(Results.ann_NumberVLS(Params.year_2019,:)'.*Params.raceIndicator(:,Params.race_H))...
                /numaware2019_H;
            Results.calib_TTdist2019_VLSamongdiag_O = ...
                sum(Results.ann_NumberVLS(Params.year_2019,:)'.*Params.raceIndicator(:,Params.race_O))...
                /numaware2019_O;

             Results.calib_TTdist2019_VLSamongdiag_Total = ...
                sum(Results.ann_NumberVLS(Params.year_2019,:)) ...
                /sum(Results.ann_NumberAware(Params.year_2019,:));
            
          %Incidence by race- Moved from 2016 4/21/20 (renamed to B, H, O).
          %Updated to 2019 - Clinkscales 6/30/2021
                      
            Results.calib_NewInfections2019_B = ...
                Results.ann_TotalNewInfections(Params.year_2019,:) ...
                *Params.raceIndicator(:,Params.race_B);
            
            Results.calib_NewInfections2019_H = ...
                Results.ann_TotalNewInfections(Params.year_2019,:) ...
                *Params.raceIndicator(:,Params.race_H);
            
            Results.calib_NewInfections2019_O = ...
                Results.ann_TotalNewInfections(Params.year_2019,:) ...
                *Params.raceIndicator(:,Params.race_O);
            
          
          %Incidence by transmission group- Added 4/20/20
          Results.calib_NewInfections2019_HETM = Results.ann_TotalNewInfections(Params.year_2019,:)*...
              Params.popSexIndicator(:,Params.popSex_HETM);
          Results.calib_NewInfections2019_HETF = Results.ann_TotalNewInfections(Params.year_2019,:)*...
              Params.popSexIndicator(:,Params.popSex_HETF);
           Results.calib_NewInfections2019_MSM = Results.ann_TotalNewInfections(Params.year_2019,:)*...
              Params.popSexIndicator(:,Params.popSex_MSM);
          Results.calib_NewInfections2019_IDUM = Results.ann_TotalNewInfections(Params.year_2019,:)*...
              Params.popSexIndicator(:,Params.popSex_IDUM);
          Results.calib_NewInfections2019_IDUF = Results.ann_TotalNewInfections(Params.year_2019,:)*...
              Params.popSexIndicator(:,Params.popSex_IDUF); 
          
          %Incidence overall- Moved from 2016 4/21/20. Moved to 2019 -
          %Clinkscales 6/30/2021
            Results.calib_NewInfections2019_Total = ...
                sum(Results.ann_TotalNewInfections(Params.year_2019,:));  
            
            %HIV Prevalence for LR and HR HETS in 2015 moved to 2016- moved
            %to 2018 4/21/20. Moved t0 2019 - Clinkcales 6/30/2021
            calib_LRHETPrevalence2019 = Results.ann_HIVPrevalence_HET(Params.year_2019,:)*...
                Params.riskLevelIndicator(:,Params.risk_Main);

            calib_HRHETPrevalence2019 = Results.ann_HIVPrevalence_HET(Params.year_2019,:)*...
                Params.riskLevelIndicator(:,Params.risk_Casual);

            
            % PLWH Aware Deaths - moved to 2016. 6/14/19. Bates. Updated to
            % 2019 - Clinkscales 6/30/2021
            Results.calib_AwarePLWHDeaths2019 = sum(Results.ann_numDeathsPLWHAware(Params.year_2019,:));

            % PLWH AIDS deaths in 2019, Viguerie. Aware from surveillance, and unaware is a guess (it should be low). Added 12/9/2022.
            Results.calib_AwarePLWHAIDSDeaths2019 = sum(Results.ann_numAIDSDeathsAware(Params.year_2019,:));
            Results.calib_UnawarePLWHAIDSDeaths2019 = sum(Results.ann_numAIDSDeathsHIVPosUndiag(Params.year_2019,:));


            %Added New Diagnoses in 2019 - Clinkscales 04/21/2022
            Results.calib_NewDiagnoses2019_Total = sum(Results.ann_TotalNewDiagnoses(Params.year_2019,:));  
            Results.calib_NewDiagnoses2019_PctB1C1 = ...
                sum(Results.ann_NewDiagnoses_B1(Params.year_2019,:) +...
                Results.ann_NewDiagnoses_C1(Params.year_2019,:)) /...
                Results.calib_NewDiagnoses2019_Total;

            %Number on PrEP in 2019- Bates 02/10/23
            Results.calib_NumOnPrEP2019_Male = Results.ann_NumberOnPrEP(Params.year_2019,:) * ...
                    Params.sexIndicator(:,Params.sex_Male);
            Results.calib_NumOnPrEP2019_Female = Results.ann_NumberOnPrEP(Params.year_2019,:) * ...
                    Params.sexIndicator(:,Params.sex_Female);
            Results.calib_NumOnPrEP2019_Black = Results.ann_NumberOnPrEP(Params.year_2019,:) * ...
                    Params.raceIndicator(:,Params.race_B);
            Results.calib_NumOnPrEP2019_Hispanic = Results.ann_NumberOnPrEP(Params.year_2019,:) * ...
                    Params.raceIndicator(:,Params.race_H);
            Results.calib_NumOnPrEP2019_Other = Results.ann_NumberOnPrEP(Params.year_2019,:) * ...
                    Params.raceIndicator(:,Params.race_O);
                
            Results.calib_NumOnPrEP2019_Total = sum(Results.ann_NumberOnPrEP(Params.year_2019,:));
            
            % Targets that use 2010 outcomes
            if FirstOutcomeYr <= 2010 && (LastOutcomeYr >= 2010)

                % Ratio of HIV prevalence 2018/2010. Updated to 2019.
                % Clinkscales 6/29/2021
                Results.calib_OverallPrev2019v2010 = ...
                    Results.calib_HIVPrevalence2019 / calib_HIVPrevalence2010;     

                Results.calib_HETPrev2019v2010_LR = ...
                    calib_LRHETPrevalence2019 / calib_LRHETPrevalence2010;

                Results.calib_HETPrev2019v2010_HR = ...
                    calib_HRHETPrevalence2019 / calib_HRHETPrevalence2010; 

            end  
            
          %Population values 2019 (Updated from 2016) 11/22/22, Bates
            Results.calib_populationMaleBlk2019 = sum(Results.ann_PopulationSize(Params.year_2019,:)' ...
                    .* Params.sexIndicator(:,Params.sex_Male) ...
                    .* Params.raceIndicator(:,Params.race_B));
                
            Results.calib_populationMaleHisp2019 = sum(Results.ann_PopulationSize(Params.year_2019,:)' ...
                    .* Params.sexIndicator(:,Params.sex_Male) ...
                    .* Params.raceIndicator(:,Params.race_H));

            Results.calib_populationMaleOth2019 = sum(Results.ann_PopulationSize(Params.year_2019,:)' ...
                    .* Params.sexIndicator(:,Params.sex_Male) ...
                    .* Params.raceIndicator(:,Params.race_O));
                
            Results.calib_populationFemaleBlk2019 = sum(Results.ann_PopulationSize(Params.year_2019,:)' ...
                    .* Params.sexIndicator(:,Params.sex_Female) ...
                    .* Params.raceIndicator(:,Params.race_B));
                
            Results.calib_populationFemaleHisp2019 = sum(Results.ann_PopulationSize(Params.year_2019,:)' ...
                    .* Params.sexIndicator(:,Params.sex_Female) ...
                    .* Params.raceIndicator(:,Params.race_H));

            Results.calib_populationFemaleOth2019 = sum(Results.ann_PopulationSize(Params.year_2019,:)' ...
                    .* Params.sexIndicator(:,Params.sex_Female) ...
                    .* Params.raceIndicator(:,Params.race_O));

            %Population targets by age. Moved from 2016 4/21/20. Updated to
            %2019 - Clinkscales 6/30/2021
            Results.calib_population_1334_2019 = Results.ann_PopulationSize(Params.year_2019,:) ...
                    * sum(Params.ageIndicator(:,[Params.age_1,Params.age_2,Params.age_3]),2);
            Results.calib_population_3564_2019 = Results.ann_PopulationSize(Params.year_2019,:) ...
                    * sum(Params.ageIndicator(:,[Params.age_4,Params.age_5,Params.age_6]),2);
            Results.calib_population_65_2019 = Results.ann_PopulationSize(Params.year_2019,:) ...
                    * Params.ageIndicator(:,[Params.age_7]);             
        end

        %2021 outcomes- added 02/10/23. Bates
        if FirstOutcomeYr <= 2021 && (LastOutcomeYr >= 2021)
            
            %Number on PrEP in 2021- Bates 02/10/23
            Results.calib_NumOnPrEP2021_Male = Results.ann_NumberOnPrEP(Params.year_2021,:) * ...
                    Params.sexIndicator(:,Params.sex_Male);
            Results.calib_NumOnPrEP2021_Female = Results.ann_NumberOnPrEP(Params.year_2021,:) * ...
                    Params.sexIndicator(:,Params.sex_Female);
            Results.calib_NumOnPrEP2021_Black = Results.ann_NumberOnPrEP(Params.year_2021,:) * ...
                    Params.raceIndicator(:,Params.race_B);
            Results.calib_NumOnPrEP2021_Hispanic = Results.ann_NumberOnPrEP(Params.year_2021,:) * ...
                    Params.raceIndicator(:,Params.race_H);
            Results.calib_NumOnPrEP2021_Other = Results.ann_NumberOnPrEP(Params.year_2021,:) * ...
                    Params.raceIndicator(:,Params.race_O);
                
            Results.calib_NumOnPrEP2021_Total = sum(Results.ann_NumberOnPrEP(Params.year_2021,:));
        end

        end
        


    % 7.viii. Analysis-specific Outcomes
    
        % Pre-allocate
        HighRiskOutput(numYears,Params.numStrats)=0;
        LowRiskOutput(numYears,Params.numStrats)=0;
        store_age(numYears,Params.numStrats)=0;
        store_HET(numYears,Params.numStrats)=0;
        store_MSM(numYears,Params.numStrats)=0;
        store_IDU(numYears,Params.numStrats)=0;
        
        
        for sectionAnalysis = 1:1
            
        % Population Sizes by transmission group, risk level, sex
        for sectionPopSize = 1:1
            Results.ann_popSize_Blk(numYears,1)=0;
            Results.ann_popSize_Hisp(numYears,1)=0;
            Results.ann_popSize_Oth(numYears,1)=0;
            Results.ann_popSize_HETLowMale(numYears,1)=0;
            Results.ann_popSize_HETLowFemale(numYears,1)=0;
            Results.ann_popSize_HETHighMale(numYears,1)=0;
            Results.ann_popSize_HETHighFemale(numYears,1)=0;
            Results.ann_popSize_MSMLow(numYears,1)=0;
            Results.ann_popSize_MSMHigh(numYears,1)=0;
            Results.ann_popSize_IDUHighMale(numYears,1)=0;

        for yearCount = 1:numYears
            
            Results.ann_popSize_Blk(yearCount,1)=sum(Results.ann_PopulationSize(yearCount,:)' ...
                .* Params.raceIndicator(:,Params.race_B));
            
            Results.ann_popSize_Hisp(yearCount,1)=sum(Results.ann_PopulationSize(yearCount,:)' ...
                .* Params.raceIndicator(:,Params.race_H));
            
            Results.ann_popSize_Oth(yearCount,1)=sum(Results.ann_PopulationSize(yearCount,:)' ...
                .* Params.raceIndicator(:,Params.race_O));
            
            Results.ann_popSize_HETLowMale(yearCount,1)=sum(Results.ann_PopulationSize(yearCount,:)' ...
                .* Params.popSexIndicator(:,Params.popSex_HETM) .* ...
                Params.riskLevelIndicator(:,Params.risk_Main));
            
            Results.ann_popSize_HETLowFemale(yearCount,1)=sum(Results.ann_PopulationSize(yearCount,:)' ...
                .* Params.popSexIndicator(:,Params.popSex_HETF) .* ...
                Params.riskLevelIndicator(:,Params.risk_Main));
            
            Results.ann_popSize_HETHighMale(yearCount,1)=sum(Results.ann_PopulationSize(yearCount,:)' ...
                .* Params.popSexIndicator(:,Params.popSex_HETM) .* ...
                Params.riskLevelIndicator(:,Params.risk_Casual));
            
            Results.ann_popSize_HETHighFemale(yearCount,1)=sum(Results.ann_PopulationSize(yearCount,:)' ...
                .* Params.popSexIndicator(:,Params.popSex_HETF) .* ...
                Params.riskLevelIndicator(:,Params.risk_Casual));
            
            Results.ann_popSize_MSMLow(yearCount,1)=sum(Results.ann_PopulationSize(yearCount,:)' ...
                .* Params.popSexIndicator(:,Params.popSex_MSM) .* ...
                Params.riskLevelIndicator(:,Params.risk_Main));
            
            Results.ann_popSize_MSMHigh(yearCount,1)=sum(Results.ann_PopulationSize(yearCount,:)' ...
                .* Params.popSexIndicator(:,Params.popSex_MSM) .* ...
                Params.riskLevelIndicator(:,Params.risk_Casual));
            
            Results.ann_popSize_IDUHighMale(yearCount,1)=sum(Results.ann_PopulationSize(yearCount,:)' ...
                .* Params.popSexIndicator(:,Params.popSex_IDUM) .* ...
                Params.riskLevelIndicator(:,Params.risk_Casual));
            
            Results.ann_popSize_IDUHighFemale(yearCount,1)=sum(Results.ann_PopulationSize(yearCount,:)' ...
                .* Params.popSexIndicator(:,Params.popSex_IDUF) .* ...
                Params.riskLevelIndicator(:,Params.risk_Casual)); 
            
        end
        end

        % Outcomes for testing frequency analysis
        
        for testFreqAnalysis = 1:1
             for u = 1:5
                 
                 % Calculate the high and low risk, HET-specific outcomes
                 % for each of the following variables:

                    % Main input
                 switch u                 
                     case 1
                         Input = Results.ann_TotalNewInfections;
                     case 2 
                         Input = Results.ann_LifeYears;
                     case 3 
                         Input = Results.ann_QALYs;     
                     case 4
                         Input = Results.ann_TransCost_CostTestAndNotify_Disc;
                     case 5 
                         Input = Results.ann_ARTCarePrEPCost_Disc; 
                 end
                    % Preallocate
                 
                 
                 % Generic calculation using main input
                 for yearCount = 1:numYears
                     HighRiskOutput(yearCount,:) = Input(yearCount,:) .* Params.HRHIndicator';
                     LowRiskOutput(yearCount,:) = Input(yearCount,:) .* Params.LowRiskHETIndicator';
                 end
                 
                 % Applying the calculated output (from the generic
                 % formula) to the relevant outcome
                switch u
                    case 1
                        Results.ann_NewInfections_HighRiskHETs = HighRiskOutput;
                        Results.ann_NewInfections_LowRiskHETs = LowRiskOutput;
                    case 2
                        Results.ann_LifeYears_HighRiskHETs = HighRiskOutput;
                        Results.ann_LifeYears_LowRiskHETs = LowRiskOutput;
                    case 3
                        Results.ann_QALYs_HighRiskHETs = HighRiskOutput;
                        Results.ann_QALYs_LowRiskHETs = LowRiskOutput;
                    case 4
                        Results.ann_CostTestAndNotify_HighRiskHETs_Disc = HighRiskOutput;
                        Results.ann_CostTestAndNotify_LowRiskHETs_Disc = LowRiskOutput;
                    case 5
                        Results.ann_ARTCarePrEPCost_Disc_HighRiskHETs = HighRiskOutput;
                        Results.ann_ARTCarePrEPCost_Disc_LowRiskHETs = LowRiskOutput;
                end
                
                clear HighRiskOutput
                clear LowRiskOutput
             end
             
            % Calculate values for overall HET population
            Results.ann_LifeYears_HET = ...
                  Results.ann_LifeYears_HighRiskHETs ...
                + Results.ann_LifeYears_LowRiskHETs;
            
            Results.ann_QALYs_HET = ...
                  Results.ann_QALYs_HighRiskHETs ...
                + Results.ann_QALYs_LowRiskHETs;
            
            Results.ann_CostTestAndNotify_HET_Disc = ...
                  Results.ann_CostTestAndNotify_HighRiskHETs_Disc ...
                + Results.ann_CostTestAndNotify_LowRiskHETs_Disc;
            
            Results.ann_ARTCarePrEPCost_Disc_HET = ...
                  Results.ann_ARTCarePrEPCost_Disc_HighRiskHETs ...
                + Results.ann_ARTCarePrEPCost_Disc_LowRiskHETs;
           
            % Calculate total cost (without SEP [SEP not stratified by subpopulation])
            Results.ann_ARTCarePrEPTransitionAndSEPCost_Disc_noSEP = ...
                  Results.ann_ARTCarePrEPCost_Disc ...
                + Results.ann_TotalTransCost_Disc;
            
            % Calculate total cost overall and by race/ethnicity 
            Results.ann_ARTCarePrEPTransitionAndSEPCost_Disc = ...
                  sum(Results.ann_ARTCarePrEPCost_Disc,2) ...
                + sum(Results.ann_TotalTransCost_Disc,2) ...
                + Results.ann_SEPCost_Disc; 
                        
            Results.ann_ARTCarePrEPTransitionAndSEPCost_Disc_Blk = ...
                  sum(Results.ann_ARTCarePrEPCost_Disc .* Params.raceIndicator(:,Params.race_B)',2) ...
                + sum(Results.ann_TotalTransCost_Disc .* Params.raceIndicator(:,Params.race_B)',2) ...
                + Results.ann_SEPCost_B_Disc;
            
            Results.ann_ARTCarePrEPTransitionAndSEPCost_Disc_Hisp = ...
                  sum(Results.ann_ARTCarePrEPCost_Disc .* Params.raceIndicator(:,Params.race_H)',2) ...
                + sum(Results.ann_TotalTransCost_Disc .* Params.raceIndicator(:,Params.race_H)',2) ...
                + Results.ann_SEPCost_H_Disc;
            
            Results.ann_ARTCarePrEPTransitionAndSEPCost_Disc_Oth = ...
                  sum(Results.ann_ARTCarePrEPCost_Disc .* Params.raceIndicator(:,Params.race_O)',2) ...
                + sum(Results.ann_TotalTransCost_Disc .* Params.raceIndicator(:,Params.race_O)',2) ...
                + Results.ann_SEPCost_O_Disc;
                
                % Cumulative total cost
                Results.total_ARTCarePrEPTransitionAndSEPCost_Disc = sum(Results.ann_ARTCarePrEPTransitionAndSEPCost_Disc);
            
            if and(FirstOutcomeYr < Params.tt_periodTwoStartYear,LastOutcomeYr >= Params.tt_periodTwoStartYear - 1)
                
                % Distribution of people between compartments at end of time period 1 
                DistCompartsLastYrPeriod1_store1 = squeeze(Results.ann_Compartments(Params.tt_periodTwoStartYear - FirstOutcomeYr,HIVComparts,:));
                DistCompartsLastYrPeriod1_store2 = sum(DistCompartsLastYrPeriod1_store1,2) ./ sum(sum(DistCompartsLastYrPeriod1_store1,2));
                Results.DistCompartsLastYrPeriod1 = zeros(Params.numModelsContStagesIfChronic, Params.numHIVstages);
                Results.DistCompartsLastYrPeriod1(1:3,1) = DistCompartsLastYrPeriod1_store2(1:3,1);
                Results.DistCompartsLastYrPeriod1(1:5,2) = DistCompartsLastYrPeriod1_store2(4:8,1);
                Results.DistCompartsLastYrPeriod1(1:5,3) = DistCompartsLastYrPeriod1_store2(9:13,1);
                Results.DistCompartsLastYrPeriod1(1:5,4) = DistCompartsLastYrPeriod1_store2(14:18,1);
                Results.DistCompartsLastYrPeriod1(1:5,5) = DistCompartsLastYrPeriod1_store2(19:23,1);
                
            end
            
            if and(and(FirstOutcomeYr < Params.tt_periodThreeStartYear, ...
                    Params.tt_periodTwoStartYear < Params.tt_periodThreeStartYear), ...
                    LastOutcomeYr >= Params.tt_periodThreeStartYear - 1)
                
                % Distribution of people between compartments at end of time period 2 
                DistCompartsLastYrPeriod2_store1 = squeeze(Results.ann_Compartments(Params.tt_periodThreeStartYear - FirstOutcomeYr,HIVComparts,:));
                DistCompartsLastYrPeriod2_store2 = sum(DistCompartsLastYrPeriod2_store1,2) ./ sum(sum(DistCompartsLastYrPeriod2_store1,2));
                Results.DistCompartsLastYrPeriod2 = zeros(Params.numModelsContStagesIfChronic, Params.numHIVstages);
                Results.DistCompartsLastYrPeriod2(1:3,1) = DistCompartsLastYrPeriod2_store2(1:3,1);
                Results.DistCompartsLastYrPeriod2(1:5,2) = DistCompartsLastYrPeriod2_store2(4:8,1);
                Results.DistCompartsLastYrPeriod2(1:5,3) = DistCompartsLastYrPeriod2_store2(9:13,1);
                Results.DistCompartsLastYrPeriod2(1:5,4) = DistCompartsLastYrPeriod2_store2(14:18,1);
                Results.DistCompartsLastYrPeriod2(1:5,5) = DistCompartsLastYrPeriod2_store2(19:23,1);
                
            end
            
            if and(and(FirstOutcomeYr < Params.tt_periodFourStartYear, ...
                    Params.tt_periodThreeStartYear < Params.tt_periodFourStartYear), ...
                    LastOutcomeYr >= Params.tt_periodFourStartYear - 1)
                
                % Distribution of people between compartments at end of time period 3 
                DistCompartsLastYrPeriod3_store1 = squeeze(Results.ann_Compartments(Params.tt_periodFourStartYear - FirstOutcomeYr,HIVComparts,:));
                DistCompartsLastYrPeriod3_store2 = sum(DistCompartsLastYrPeriod3_store1,2) ./ sum(sum(DistCompartsLastYrPeriod3_store1,2));
                Results.DistCompartsLastYrPeriod3 = zeros(Params.numModelsContStagesIfChronic, Params.numHIVstages);
                Results.DistCompartsLastYrPeriod3(1:3,1) = DistCompartsLastYrPeriod3_store2(1:3,1);
                Results.DistCompartsLastYrPeriod3(1:5,2) = DistCompartsLastYrPeriod3_store2(4:8,1);
                Results.DistCompartsLastYrPeriod3(1:5,3) = DistCompartsLastYrPeriod3_store2(9:13,1);
                Results.DistCompartsLastYrPeriod3(1:5,4) = DistCompartsLastYrPeriod3_store2(14:18,1);
                Results.DistCompartsLastYrPeriod3(1:5,5) = DistCompartsLastYrPeriod3_store2(19:23,1);
                
            end
            
            if and(and(FirstOutcomeYr < Params.tt_periodFiveStartYear, ...
                    Params.tt_periodFourStartYear < Params.tt_periodFiveStartYear), ...
                    LastOutcomeYr >= Params.tt_periodFiveStartYear - 1)
                
                % Distribution of people between compartments at end of time period 4 
                DistCompartsLastYrPeriod4_store1 = squeeze(Results.ann_Compartments(Params.tt_periodFiveStartYear - FirstOutcomeYr,HIVComparts,:));
                DistCompartsLastYrPeriod4_store2 = sum(DistCompartsLastYrPeriod4_store1,2) ./ sum(sum(DistCompartsLastYrPeriod4_store1,2));
                Results.DistCompartsLastYrPeriod4 = zeros(Params.numModelsContStagesIfChronic, Params.numHIVstages);
                Results.DistCompartsLastYrPeriod4(1:3,1) = DistCompartsLastYrPeriod4_store2(1:3,1);
                Results.DistCompartsLastYrPeriod4(1:5,2) = DistCompartsLastYrPeriod4_store2(4:8,1);
                Results.DistCompartsLastYrPeriod4(1:5,3) = DistCompartsLastYrPeriod4_store2(9:13,1);
                Results.DistCompartsLastYrPeriod4(1:5,4) = DistCompartsLastYrPeriod4_store2(14:18,1);
                Results.DistCompartsLastYrPeriod4(1:5,5) = DistCompartsLastYrPeriod4_store2(19:23,1);
                
            end
            
            % Distribution of people between compartments at end of last year 
                DistCompartsLastYr_store1 = squeeze(Results.ann_Compartments(numYears,HIVComparts,:));
                DistCompartsLastYr_store2 = sum(DistCompartsLastYr_store1,2) ./ sum(sum(DistCompartsLastYr_store1,2));
                Results.DistCompartsLastYr = zeros(Params.numModelsContStagesIfChronic, Params.numHIVstages);
                Results.DistCompartsLastYr(1:3,1) = DistCompartsLastYr_store2(1:3,1);
                Results.DistCompartsLastYr(1:5,2) = DistCompartsLastYr_store2(4:8,1);
                Results.DistCompartsLastYr(1:5,3) = DistCompartsLastYr_store2(9:13,1);
                Results.DistCompartsLastYr(1:5,4) = DistCompartsLastYr_store2(14:18,1);
                Results.DistCompartsLastYr(1:5,5) = DistCompartsLastYr_store2(19:23,1);
            
            % Distribution of HR-HET PLWH between compartments
            DistCompartsLastYr = squeeze(Results.ann_Compartments(numYears,:,:));
            Results.ann_DistHRHByComparts_LastYr = DistCompartsLastYr*Params.HRHIndicator;
                
     
            % Only pull test/diag rates for just one year in discrete version
            % if Params.SolveCont~= 1
            rateYear = Params.diagRateYr;
            
            if FirstOutcomeYr <= rateYear && (LastOutcomeYr > rateYear)
            
                YrIndex = max(floor(rateYear - FirstOutcomeYr + 1),0);


                % Testing Rates - only for selected year
                %All Eligible (HIV- and HIV+)
                Results.tr_WeightedAllEligTestRateByCohort = TTProg.WeightedTestRateAllEligByCohort(YrIndex,:);
                Results.tr_raceAllEligTestRate = TTProg.RaceSpecificAllEligTestRate(YrIndex,:);
                Results.tr_riskGroupAllEligTestRate = TTProg.RiskGroupAllEligTestRate(YrIndex,:);
                Results.tr_HRHAllEligTestRate = TTProg.HRHAllEligTestRate(YrIndex,:);
                Results.tr_LRHAllEligTestRate = TTProg.LRHAllEligTestRate(YrIndex,:);
                Results.tr_HRMSMAllEligTestRate = TTProg.HRMSMAllEligTestRate(YrIndex,:);
                Results.tr_LRMSMAllEligTestRate = TTProg.LRMSMAllEligTestRate(YrIndex,:);
                Results.tr_IDUAllEligTestRate = TTProg.IDUAllEligTestRate(YrIndex,:);
                Results.tr_AllEligTestRateByRiskGpLevel = TTProg.AllEligTestRateByRiskGpLevel(YrIndex,:);

                % HIV+ only       
                Results.tr_WeightedHIVPosTestRateByCohort = TTProg.WeightedTestRateHIVPosByCohort(YrIndex,:);
                Results.tr_raceHIVPosTestRate = TTProg.RaceSpecificHIVPosTestRate(YrIndex,:);
                Results.tr_riskGroupHIVPosTestRate = TTProg.RiskGroupHIVPosTestRate(YrIndex,:);
                Results.tr_HRHHIVPosTestRate = TTProg.HRHHIVPosTestRate(YrIndex,:);
                Results.tr_LRHHIVPosTestRate = TTProg.LRHHIVPosTestRate(YrIndex,:);
                Results.tr_uninfectedTestRate = TTProg.OverallUninfectedTestRate(YrIndex,:);
                Results.tr_diseaseStageTestRate = TTProg.DiseaseStageTestRate(YrIndex,:);

                % Diagnosis rates
                Results.tr_DiagRateAll = TTProg.DiagRateAll(YrIndex,:);
                Results.tr_raceDiagRate = TTProg.RaceSpecificDiagRate(YrIndex,:);
                Results.tr_riskGroupDiagRate = TTProg.RiskGroupDiagRate(YrIndex,:);
                Results.tr_diseaseStageDiagRate = TTProg.DiseaseStageDiagRate(YrIndex,:);
                Results.tr_LRHDiagRate = TTProg.LowRiskHETDiagRate(YrIndex,:);
                Results.tr_HRHDiagRate = TTProg.HRHDiagRate(YrIndex,:);
            % end
            end         
        end
        end
             
        % Outcomes for HET AI Analysis
        for HetAIAnalysis = 1:1
             
             % Total Infections from Vaginal Sex and total Infections from
             % Anal Sex
            Results.ann_NewInfectionsV = ...
                  Results.ann_NewInfectionsVonly ...
                + Results.ann_NewInfectionsVsomeAI;
            
            Results.ann_NewInfectionsA = ...
                  Results.ann_NewInfectionsAonly ...
                + Results.ann_NewInfectionsAsomeVI;
             
            % Stratify new infections by force of infection and transission group
            
                % IDU Needle Infections
                Results.ann_NewInfectionsN_IDU = Results.ann_NewInfectionsN;
                
                % All others
            for yearCount = 1:numYears
                
                % HET
                Results.ann_NewInfectionsV_HET(yearCount,:) = ...
                    Results.ann_NewInfectionsV(yearCount,:) .* ...
                    Params.popIndicator(:,Params.pop_HET)';
                
                Results.ann_NewInfectionsA_HET(yearCount,:) = ...
                    Results.ann_NewInfectionsA(yearCount,:) .* ...
                    Params.popIndicator(:,Params.pop_HET)';        
                
                % MSM   
                Results.ann_NewInfectionsV_MSM(yearCount,:) = ...
                    Results.ann_NewInfectionsV(yearCount,:) .* ...
                    Params.popIndicator(:,Params.pop_MSM)';
                
                Results.ann_NewInfectionsA_MSM(yearCount,:) = ...
                    Results.ann_NewInfectionsA(yearCount,:) .* ...
                    Params.popIndicator(:,Params.pop_MSM)';    
                
                % IDU
                Results.ann_NewInfectionsV_IDU(yearCount,:) = ...
                    Results.ann_NewInfectionsV(yearCount,:) .* ...
                    Params.popIndicator(:,Params.pop_IDU)';
                
                Results.ann_NewInfectionsA_IDU(yearCount,:) = ...
                    Results.ann_NewInfectionsA(yearCount,:) .* ...
                    Params.popIndicator(:,Params.pop_IDU)';    
            end
            
            % HET AI infections by age
            for yearCount = 1:numYears
                
                % 13-17
               Results.ann_NewInfectionsA_HET_13_17(yearCount,:)= ...
                   Results.ann_NewInfectionsA_HET(yearCount,:).* ...
                   Params.ageIndicator(:,Params.age_1)';
               %18-24
               Results.ann_NewInfectionsA_HET_18_24(yearCount,:)= ...
                   Results.ann_NewInfectionsA_HET(yearCount,:).* ...
                   Params.ageIndicator(:,Params.age_2)';
               %25-34
               Results.ann_NewInfectionsA_HET_25_34(yearCount,:)= ...
                   Results.ann_NewInfectionsA_HET(yearCount,:).* ...
                   Params.ageIndicator(:,Params.age_3)';
               %35-44
               Results.ann_NewInfectionsA_HET_35_44(yearCount,:)= ...
                   Results.ann_NewInfectionsA_HET(yearCount,:).* ...
                   Params.ageIndicator(:,Params.age_4)';
               %45-64
               Results.ann_NewInfectionsA_HET_45_54(yearCount,:)= ...
                   Results.ann_NewInfectionsA_HET(yearCount,:).* ...
                   Params.ageIndicator(:,Params.age_5)';
               %55-64
               Results.ann_NewInfectionsA_HET_55_64(yearCount,:)= ...
                   Results.ann_NewInfectionsA_HET(yearCount,:).* ...
                   Params.ageIndicator(:,Params.age_6)';
               %65+
               Results.ann_NewInfectionsA_HET_65(yearCount,:)= ...
                   Results.ann_NewInfectionsA_HET(yearCount,:).* ...
                   Params.ageIndicator(:,Params.age_7)';
            end
            
            % Overall HET new infections
            Results.ann_NewInfections_HET = ...
                  Results.ann_NewInfectionsV_HET ...
                + Results.ann_NewInfectionsA_HET;
            
            % HET Infections by Age
            for yearCount = 1:numYears
                %13-17
                Results.ann_NewInfections_HET_13_17(yearCount,:) = ...
                    Results.ann_NewInfections_HET(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_1)';
                %18-24
                Results.ann_NewInfections_HET_18_24(yearCount,:) = ...
                    Results.ann_NewInfections_HET(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_2)';
                %25-34
                Results.ann_NewInfections_HET_25_34(yearCount,:) = ...
                    Results.ann_NewInfections_HET(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_3)';
                %35-44
                Results.ann_NewInfections_HET_35_44(yearCount,:) = ...
                    Results.ann_NewInfections_HET(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_4)';
                %45-54
                Results.ann_NewInfections_HET_45_54(yearCount,:) = ...
                    Results.ann_NewInfections_HET(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_5)';
                %55-64
                Results.ann_NewInfections_HET_55_64(yearCount,:) = ...
                    Results.ann_NewInfections_HET(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_6)';
                %65+
                Results.ann_NewInfections_HET_65(yearCount,:) = ...
                    Results.ann_NewInfections_HET(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_7)';
            end
            
            % Overall Infections by Age
            for yearCount = 1:numYears
                %13-17
                Results.ann_NewInfections_13_17(yearCount,:) = ...
                    Results.ann_TotalNewInfections(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_1)';
                %18-24
                Results.ann_NewInfections_18_24(yearCount,:) = ...
                    Results.ann_TotalNewInfections(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_2)';
                %25-34
                Results.ann_NewInfections_25_34(yearCount,:) = ...
                    Results.ann_TotalNewInfections(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_3)';
                %35-44
                Results.ann_NewInfections_35_44(yearCount,:) = ...
                    Results.ann_TotalNewInfections(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_4)';
                %45-54
                Results.ann_NewInfections_45_54(yearCount,:) = ...
                    Results.ann_TotalNewInfections(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_5)';
                %55-64
                Results.ann_NewInfections_55_64(yearCount,:) = ...
                    Results.ann_TotalNewInfections(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_6)';
                %65+
                Results.ann_NewInfections_65(yearCount,:) = ...
                    Results.ann_TotalNewInfections(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_7)';
                
            end
            
            % PLWH Deaths by Age            
            for yearCount = 1:numYears
                %13-17
                Results.ann_numDeathsPLWHAware_13(yearCount,:) = ...
                    Results.ann_numDeathsPLWHAware(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_1)';
                %18-24
                Results.ann_numDeathsPLWHAware_18(yearCount,:) = ...
                    Results.ann_numDeathsPLWHAware(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_2)';
                %25-34
                Results.ann_numDeathsPLWHAware_25(yearCount,:) = ...
                    Results.ann_numDeathsPLWHAware(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_3)';
                %35-44
                Results.ann_numDeathsPLWHAware_35(yearCount,:) = ...
                    Results.ann_numDeathsPLWHAware(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_4)';
                %45-54
                Results.ann_numDeathsPLWHAware_45(yearCount,:) = ...
                    Results.ann_numDeathsPLWHAware(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_5)';
                %55-64
                Results.ann_numDeathsPLWHAware_55(yearCount,:) = ...
                    Results.ann_numDeathsPLWHAware(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_6)';
                %65+
                Results.ann_numDeathsPLWHAware_65(yearCount,:) = ...
                    Results.ann_numDeathsPLWHAware(yearCount,:) .* ...
                    Params.ageIndicator(:,Params.age_7)';
                
            end

            % HIV Prevalence
            for age = 1:Params.numAge
                
                for yearCount = 1:numYears
                    store_age(yearCount,:) = Results.ann_HIVPrevalence(yearCount,:) .* ...
                        Params.ageIndicator(:,age)';
                    
                    store_HET(yearCount,:) = Results.ann_HIVPrevalence(yearCount,:) .* ...
                        Params.ageIndicator(:,age)' .* Params.popIndicator(:,Params.pop_HET)';
                    
                    store_MSM(yearCount,:) = Results.ann_HIVPrevalence(yearCount,:) .* ...
                        Params.ageIndicator(:,age)' .* Params.popIndicator(:,Params.pop_MSM)';
                    
                    store_IDU(yearCount,:) = Results.ann_HIVPrevalence(yearCount,:) .* ...
                        Params.ageIndicator(:,age)' .* Params.popIndicator(:,Params.pop_IDU)';
                end
                
                switch age
                    case Params.age_1
                        Results.ann_HIVPrevalence_13_17 = store_age;
                        Results.ann_HIVPrevalenceHET_13_17 = store_HET;
                        Results.ann_HIVPrevalenceMSM_13_17 = store_MSM;
                        Results.ann_HIVPrevalenceIDU_13_17 = store_IDU;
                        
                    case Params.age_2
                        Results.ann_HIVPrevalence_18_24 = store_age;
                        Results.ann_HIVPrevalenceHET_18_24 = store_HET;
                        Results.ann_HIVPrevalenceMSM_18_24 = store_MSM;
                        Results.ann_HIVPrevalenceIDU_18_24 = store_IDU;
                        
                    case Params.age_3
                        Results.ann_HIVPrevalence_25_34 = store_age;
                        Results.ann_HIVPrevalenceHET_25_34 = store_HET;
                        Results.ann_HIVPrevalenceMSM_25_34 = store_MSM;
                        Results.ann_HIVPrevalenceIDU_25_34 = store_IDU;
                        
                    case Params.age_4
                        Results.ann_HIVPrevalence_35_44 = store_age;
                        Results.ann_HIVPrevalenceHET_35_44 = store_HET;
                        Results.ann_HIVPrevalenceMSM_35_44 = store_MSM;
                        Results.ann_HIVPrevalenceIDU_35_44 = store_IDU;
                        
                    case Params.age_5
                        Results.ann_HIVPrevalence_45_54 = store_age;
                        Results.ann_HIVPrevalenceHET_45_54 = store_HET;
                        Results.ann_HIVPrevalenceMSM_45_54 = store_MSM;
                        Results.ann_HIVPrevalenceIDU_45_54 = store_IDU;
                        
                    case Params.age_6
                        Results.ann_HIVPrevalence_55_64 = store_age;
                        Results.ann_HIVPrevalenceHET_55_64 = store_HET;
                        Results.ann_HIVPrevalenceMSM_55_64 = store_MSM;
                        Results.ann_HIVPrevalenceIDU_55_64 = store_IDU;
                        
                    case Params.age_7
                        Results.ann_HIVPrevalence_65 = store_age;
                        Results.ann_HIVPrevalenceHET_65 = store_HET;
                        Results.ann_HIVPrevalenceMSM_65 = store_MSM;
                        Results.ann_HIVPrevalenceIDU_65 = store_IDU;
                end
            end
            
         end
      
        % High Risk HETs versus Low Risk HETs
        for yearCount = 1:numYears    
             
             % HRH
            Results.ann_HIVPrevalence_HighRiskHETs(yearCount,:) = Results.ann_HIVPrevalence(yearCount,:) .* ...
                Params.HRHIndicator';
            
            Results.ann_NewInfections_HighRiskHETs(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
                Params.HRHIndicator';
            
            Results.ann_NewInfectionsA_HighRiskHETs(yearCount,:) = Results.ann_NewInfectionsA(yearCount,:) .* ...
                Params.HRHIndicator';
            
            Results.ann_NewInfectionsV_HighRiskHETs(yearCount,:) = Results.ann_NewInfectionsV(yearCount,:) .* ...
                Params.HRHIndicator';
            
            % Low Risk HETs
            Results.ann_HIVPrevalence_LowRiskHETs(yearCount,:) = Results.ann_HIVPrevalence(yearCount,:) .* ...
                Params.LowRiskHETIndicator';
            
            Results.ann_NewInfections_LowRiskHETs(yearCount,:) = Results.ann_TotalNewInfections(yearCount,:) .* ...
                Params.LowRiskHETIndicator';
            
            Results.ann_NewInfectionsA_LowRiskHETs(yearCount,:) = Results.ann_NewInfectionsA(yearCount,:) .* ...
                Params.LowRiskHETIndicator';
            
            Results.ann_NewInfectionsV_LowRiskHETs(yearCount,:) = Results.ann_NewInfectionsV(yearCount,:) .* ...
                Params.LowRiskHETIndicator';
            
         end  
     
        % NHAS analysis
        for sectionNHAS = 1:1
            
            % To collect
            
                % New Infections annually, by r/e and total
                % Prevalance annually, by r/e and total
                % Continuum-of-care in 2020
                    % Percentage: Among PLWH:
                        % Percent Aware among HIV+, by r/e and total
                        % Percent LTC among HIV+, by r/e and total
                        % Percent On ART among HIV+, by r/e and total
                        % Percent VLS among HIV+, by r/e and total
                    % Percentage: Among diagnosed:
                        % Percent LTC among diagnosed, by r/e and total
                        % Percent On ART among diagnosed, by r/e and total
                        % Percent VLS among diagnosed, by r/e and total
                % Counts: Continuum-of-care in 2015 and 2020
                    % Number Aware, by r/e and total
                    % Number LTC, by r/e and total
                    % Number on ART, by r/e and total
                    % Number VLS, by r/e and total
          
       %  Check if 2015 is included in the time horizon
         if and(Params.year_2020 > 0, numYears >= Params.year_2020)
       
            Results.year_2010 = Params.year_2010;
            Results.year_2015 = Params.year_2015;
            Results.year_2020 = Params.year_2020;

           % Calculate denominators
            
                % Race-specific prevalence in 2020
                Prev2020_Blk = sum(Results.ann_HIVPrevalence_Blk(Params.year_2020,:));
                Prev2020_Hisp = sum(Results.ann_HIVPrevalence_Hisp(Params.year_2020,:));
                Prev2020_Oth = sum(Results.ann_HIVPrevalence_Oth(Params.year_2020,:));
                Prev2020_Total = sum(Results.ann_HIVPrevalence(Params.year_2020,:));
                
                % Race-specific number aware in 2020
                NumAware2020_Blk = sum(Results.ann_NumberAware_Blk(Params.year_2020,:));
                NumAware2020_Hisp = sum(Results.ann_NumberAware_Hisp(Params.year_2020,:));
                NumAware2020_Oth = sum(Results.ann_NumberAware_Oth(Params.year_2020,:));
                NumAware2020_Total = sum(Results.ann_NumberAware(Params.year_2020,:));
                        
                
                      
       % CONTINUUM-OF-CARE, COUNTS
             for sectionCC_Counts = 1:1
                    % Number Aware, by r/e and total
                    % Number LTC, by r/e and total
                    % Number on ART, by r/e and total
                    % Number VLS, by r/e and total
          
                
            % Number undiagnosed          
                Results.nhas_NumberUndiagnosed_2020_Blk = Prev2020_Blk - NumAware2020_Blk;
                Results.nhas_NumberUndiagnosed_2020_Hisp = Prev2020_Hisp - NumAware2020_Hisp;
                Results.nhas_NumberUndiagnosed_2020_Oth = Prev2020_Oth - NumAware2020_Oth;
                Results.nhas_NumberUndiagnosed_2020_Total = Prev2020_Total - NumAware2020_Total;
                
            % Number Aware, by r/e and total
                Results.nhas_NumberAware_2020_Blk = sum(Results.ann_NumberAware_Blk(Params.year_2020,:));
                Results.nhas_NumberAware_2020_Hisp = sum(Results.ann_NumberAware_Hisp(Params.year_2020,:));
                Results.nhas_NumberAware_2020_Oth = sum(Results.ann_NumberAware_Oth(Params.year_2020,:));
                Results.nhas_NumberAware_2020_Total = sum(Results.ann_NumberAware(Params.year_2020,:));

            % Number In Care, by r/e and total
                Results.nhas_NumberInCare_2020_Blk = sum(Results.ann_NumberInCare_Blk(Params.year_2020,:));
                Results.nhas_NumberInCare_2020_Hisp = sum(Results.ann_NumberInCare_Hisp(Params.year_2020,:));
                Results.nhas_NumberInCare_2020_Oth = sum(Results.ann_NumberInCare_Oth(Params.year_2020,:));
                Results.nhas_NumberInCare_2020_Total = sum(Results.ann_NumberInCare(Params.year_2020,:));

            % Number On ART, by r/e and total
                Results.nhas_NumberOnART_2020_Blk = sum(Results.ann_NumberOnART_Blk(Params.year_2020,:));
                Results.nhas_NumberOnART_2020_Hisp = sum(Results.ann_NumberOnART_Hisp(Params.year_2020,:));
                Results.nhas_NumberOnART_2020_Oth = sum(Results.ann_NumberOnART_Oth(Params.year_2020,:));
                Results.nhas_NumberOnART_2020_Total = sum(Results.ann_NumberOnART(Params.year_2020,:));
       
           % Number ANV+VLS, by r/e and total     
                Results.nhas_NumberAffByART_2020_Blk = sum(Results.ann_NumberCCStages45_Blk(Params.year_2020,:));
                Results.nhas_NumberAffByART_2020_Hisp = sum(Results.ann_NumberCCStages45_Hisp(Params.year_2020,:));
                Results.nhas_NumberAffByART_2020_Oth = sum(Results.ann_NumberCCStages45_Oth(Params.year_2020,:));
                Results.nhas_NumberAffByART_2020_Total = sum(Results.ann_NumberCCStages45(Params.year_2020,:));
                
           % Number VLS, by r/e and total
                Results.nhas_NumberVLS_2020_Blk = sum(Results.ann_NumberVLS_Blk(Params.year_2020,:));
                Results.nhas_NumberVLS_2020_Hisp = sum(Results.ann_NumberVLS_Hisp(Params.year_2020,:));
                Results.nhas_NumberVLS_2020_Oth = sum(Results.ann_NumberVLS_Oth(Params.year_2020,:));
                Results.nhas_NumberVLS_2020_Total = sum(Results.ann_NumberVLS(Params.year_2020,:));

             end
              
        % CONTINUUM-OF-CARE, AMONG PLWH     
          for sectionCCAmongPLWH = 1:1
            
             % Percent Unaware among HIV+, by r/e and total
                Results.nhas_pctOfPLWH_Undiagnosed_2020_Blk = Results.nhas_NumberUndiagnosed_2020_Blk / Prev2020_Blk;
                Results.nhas_pctOfPLWH_Undiagnosed_2020_Hisp = Results.nhas_NumberUndiagnosed_2020_Hisp / Prev2020_Hisp;
                Results.nhas_pctOfPLWH_Undiagnosed_2020_Oth = Results.nhas_NumberUndiagnosed_2020_Oth / Prev2020_Oth;
                Results.nhas_pctOfPLWH_Undiagnosed_2020_Total = Results.nhas_NumberUndiagnosed_2020_Total / Prev2020_Total;
                
            % Percent Aware among HIV+, by r/e and total
                Results.nhas_pctOfPLWH_Aware_2020_Blk = Results.nhas_NumberAware_2020_Blk / Prev2020_Blk;
                Results.nhas_pctOfPLWH_Aware_2020_Hisp = Results.nhas_NumberAware_2020_Hisp / Prev2020_Hisp;
                Results.nhas_pctOfPLWH_Aware_2020_Oth = Results.nhas_NumberAware_2020_Oth / Prev2020_Oth;
                Results.nhas_pctOfPLWH_Aware_2020_Total = Results.nhas_NumberAware_2020_Total / Prev2020_Total;

             % Percent In Care among HIV+, by r/e and total
                Results.nhas_pctOfPLWH_InCare_2020_Blk = Results.nhas_NumberInCare_2020_Blk / Prev2020_Blk;
                Results.nhas_pctOfPLWH_InCare_2020_Hisp = Results.nhas_NumberInCare_2020_Hisp / Prev2020_Hisp;
                Results.nhas_pctOfPLWH_InCare_2020_Oth = Results.nhas_NumberInCare_2020_Oth / Prev2020_Oth;
                Results.nhas_pctOfPLWH_InCare_2020_Total = Results.nhas_NumberInCare_2020_Total / Prev2020_Total;
                
            % Percent on ART among HIV+, by r/e and total
                Results.nhas_pctOfPLWH_OnART_2020_Blk = Results.nhas_NumberOnART_2020_Blk / Prev2020_Blk;
                Results.nhas_pctOfPLWH_OnART_2020_Hisp = Results.nhas_NumberOnART_2020_Hisp / Prev2020_Hisp;
                Results.nhas_pctOfPLWH_OnART_2020_Oth = Results.nhas_NumberOnART_2020_Oth / Prev2020_Oth;
                Results.nhas_pctOfPLWH_OnART_2020_Total = Results.nhas_NumberOnART_2020_Total / Prev2020_Total;       

            % Percent VLS among HIV+, by r/e and total
                Results.nhas_pctOfPLWH_VLS_2020_Blk = Results.nhas_NumberVLS_2020_Blk / Prev2020_Blk;
                Results.nhas_pctOfPLWH_VLS_2020_Hisp = Results.nhas_NumberVLS_2020_Hisp / Prev2020_Hisp;
                Results.nhas_pctOfPLWH_VLS_2020_Oth = Results.nhas_NumberVLS_2020_Oth / Prev2020_Oth;
                Results.nhas_pctOfPLWH_VLS_2020_Total = Results.nhas_NumberVLS_2020_Total / Prev2020_Total;        
                 
             % Percent VLS among ART-affected
                Results.nhas_pctOfAffByART_VLS_2020_Blk = Results.nhas_NumberVLS_2020_Blk / Results.nhas_NumberAffByART_2020_Blk;
                Results.nhas_pctOfAffByART_VLS_2020_Hisp = Results.nhas_NumberVLS_2020_Hisp / Results.nhas_NumberAffByART_2020_Hisp;
                Results.nhas_pctOfAffByART_VLS_2020_Oth = Results.nhas_NumberVLS_2020_Oth / Results.nhas_NumberAffByART_2020_Oth;
                Results.nhas_pctOfAffByART_VLS_2020_Total = Results.nhas_NumberVLS_2020_Total / Results.nhas_NumberAffByART_2020_Total;
          end
          
      % CONTINUUM-OF-CARE AMONG DIAGNOSED
          for contCareAmongDiag = 1:1
              
             % Percent In Care among diagnosed, by r/e and total
                Results.nhas_pctOfAware_InCare_2020_Blk = Results.nhas_NumberInCare_2020_Blk / NumAware2020_Blk;
                Results.nhas_pctOfAware_InCare_2020_Hisp = Results.nhas_NumberInCare_2020_Hisp / NumAware2020_Hisp;
                Results.nhas_pctOfAware_InCare_2020_Oth = Results.nhas_NumberInCare_2020_Oth / NumAware2020_Oth;
                Results.nhas_pctOfAware_InCare_2020_Total = Results.nhas_NumberInCare_2020_Total / NumAware2020_Total;
                
            % Percent on ART among diagnosed, by r/e and total
                Results.nhas_pctOfAware_OnART_2020_Blk = Results.nhas_NumberOnART_2020_Blk / NumAware2020_Blk;
                Results.nhas_pctOfAware_OnART_2020_Hisp = Results.nhas_NumberOnART_2020_Hisp / NumAware2020_Hisp;
                Results.nhas_pctOfAware_OnART_2020_Oth = Results.nhas_NumberOnART_2020_Oth / NumAware2020_Oth;
                Results.nhas_pctOfAware_OnART_2020_Total = Results.nhas_NumberOnART_2020_Total / NumAware2020_Total;       

            % Percent VLS among diagnosed, by r/e and total
                Results.nhas_pctOfAware_VLS_2020_Blk = Results.nhas_NumberVLS_2020_Blk / NumAware2020_Blk;
                Results.nhas_pctOfAware_VLS_2020_Hisp = Results.nhas_NumberVLS_2020_Hisp / NumAware2020_Hisp;
                Results.nhas_pctOfAware_VLS_2020_Oth = Results.nhas_NumberVLS_2020_Oth / NumAware2020_Oth;
                Results.nhas_pctOfAware_VLS_2020_Total = Results.nhas_NumberVLS_2020_Total / NumAware2020_Total;      
              
          end

            
            
        end
        end 
        
        % PrEP Analysis
        for sectionPrEP = 1:1
            
            
            for yearCount = 1:numYears

                % Number on PrEP by sex, by transmission group, and by race
                % (added sex and race 2/10/23)


                
                    Results.ann_NumberOnPrEP_HighRiskHETs(yearCount,:) = ...
                           Results.ann_NumberOnPrEP(yearCount,:) ...
                        .* Params.popIndicator(:,Params.pop_HET)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_NumberOnPrEP_HighRiskHETs_M(yearCount,:) = ...
                           Results.ann_NumberOnPrEP(yearCount,:) ...
                        .* Params.popSexIndicator(:,Params.popSex_HETM)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_NumberOnPrEP_HighRiskHETs_F(yearCount,:) = ...
                           Results.ann_NumberOnPrEP(yearCount,:) ...
                        .* Params.popSexIndicator(:,Params.popSex_HETF)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';

                    Results.ann_NumberOnPrEP_HighRiskMSM(yearCount,:) = ...
                           Results.ann_NumberOnPrEP(yearCount,:) ...
                        .* Params.popIndicator(:,Params.pop_MSM)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';

                    Results.ann_NumberOnPrEP_HighRiskIDUs(yearCount,:) = ...
                           Results.ann_NumberOnPrEP(yearCount,:) ...
                        .* Params.popIndicator(:,Params.pop_IDU)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_NumberOnPrEP_Blk(yearCount,:) = ...
                           Results.ann_NumberOnPrEP(yearCount,:) ...
                        .* Params.raceIndicator(:,Params.race_B)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_NumberOnPrEP_Hisp(yearCount,:) = ...
                           Results.ann_NumberOnPrEP(yearCount,:) ...
                        .* Params.raceIndicator(:,Params.race_H)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_NumberOnPrEP_Oth(yearCount,:) = ...
                           Results.ann_NumberOnPrEP(yearCount,:) ...
                        .* Params.raceIndicator(:,Params.race_O)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                 % Person-years on PrEP
                    Results.ann_PersonYears_onPrEP_HETM(yearCount,:) = Results.ann_PersonYears_onPrEP(yearCount,:) .* ...
                    Params.popSexIndicator(:,Params.popSex_HETM)';
                
                    Results.ann_PersonYears_onPrEP_HETF(yearCount,:) = Results.ann_PersonYears_onPrEP(yearCount,:) .* ...
                    Params.popSexIndicator(:,Params.popSex_HETF)';
                
                    Results.ann_PersonYears_onPrEP_MSM(yearCount,:) = Results.ann_PersonYears_onPrEP(yearCount,:) .* ...
                    Params.popSexIndicator(:,Params.popSex_MSM)';
                
                    Results.ann_PersonYears_onPrEP_IDU(yearCount,:) = Results.ann_PersonYears_onPrEP(yearCount,:) .* ...
                    Params.popIndicator(:,Params.pop_IDU)';
                    
                % PrEP costs (undiscounted) by risk group
                    Results.ann_healthStateCost_PrEP_Undisc_HighRiskHETs(yearCount,:) = ...
                        Results.ann_healthStateCost_PrEP_Undisc(yearCount,:) ...
                        .* Params.popIndicator(:,Params.pop_HET)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_healthStateCost_PrEP_Undisc_HighRiskMSM(yearCount,:) = ...
                        Results.ann_healthStateCost_PrEP_Undisc(yearCount,:) ...
                        .* Params.popIndicator(:,Params.pop_MSM)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_healthStateCost_PrEP_Undisc_HighRiskIDUs(yearCount,:) = ...
                        Results.ann_healthStateCost_PrEP_Undisc(yearCount,:) ...
                        .* Params.popIndicator(:,Params.pop_IDU)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
              
                % Number susceptible high-risk by risk group and race/ethnicity
                    Results.ann_NumberUninfected_HighRiskHETs(yearCount,:) = ...
                           Results.ann_NumberUninfected(yearCount,:) ...
                        .* Params.popIndicator(:,Params.pop_HET)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_NumberUninfected_HighRiskMSM(yearCount,:) = ...
                           Results.ann_NumberUninfected(yearCount,:) ...
                        .* Params.popIndicator(:,Params.pop_MSM)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_NumberUninfected_HighRiskIDUs(yearCount,:) = ...
                           Results.ann_NumberUninfected(yearCount,:) ...
                        .* Params.popIndicator(:,Params.pop_IDU)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_NumberUninfected_HighRiskBlk(yearCount,:) = ...
                           Results.ann_NumberUninfected(yearCount,:) ...
                        .* Params.raceIndicator(:,Params.race_B)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_NumberUninfected_HighRiskHisp(yearCount,:) = ...
                           Results.ann_NumberUninfected(yearCount,:) ...
                        .* Params.raceIndicator(:,Params.race_H)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_NumberUninfected_HighRiskOth(yearCount,:) = ...
                           Results.ann_NumberUninfected(yearCount,:) ...
                        .* Params.raceIndicator(:,Params.race_O)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                 % Start PrEP
                    Results.ann_StartPrEP_HighRiskHETs(yearCount,:) = ...
                           Results.ann_StartPrEP(yearCount,:) ...
                        .* Params.popIndicator(:,Params.pop_HET)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_StartPrEP_HighRiskMSM(yearCount,:) = ...
                           Results.ann_StartPrEP(yearCount,:) ...
                        .* Params.popIndicator(:,Params.pop_MSM)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_StartPrEP_HighRiskIDUs(yearCount,:) = ...
                           Results.ann_StartPrEP(yearCount,:) ...
                        .* Params.popIndicator(:,Params.pop_IDU)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                 % Drop-off PrEP
                    Results.ann_DropOff_PrEP_HighRiskHETs(yearCount,:) = ...
                           Results.ann_DropOff_PrEP(yearCount,:) ...
                        .* Params.popIndicator(:,Params.pop_HET)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_DropOff_PrEP_HighRiskMSM(yearCount,:) = ...
                           Results.ann_DropOff_PrEP(yearCount,:) ...
                        .* Params.popIndicator(:,Params.pop_MSM)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_DropOff_PrEP_HighRiskIDUs(yearCount,:) = ...
                           Results.ann_DropOff_PrEP(yearCount,:) ...
                        .* Params.popIndicator(:,Params.pop_IDU)' ...
                        .* Params.riskLevelIndicator(:,Params.risk_Casual)';
                
            end
               
            
            % Pct On-PrEP Among Sus
%                 Results.ann_PctSusOnPrEP_Total = ...
%                     sum(Results.ann_NumberOnPrEP,2) ./ max(sum(Results.ann_NumberUninfected,2),0.00000000001);
                
                Results.ann_PctSusOnPrEP_HighRiskHETs = ...
                    sum(Results.ann_NumberOnPrEP_HighRiskHETs,2) ./ max(sum(Results.ann_NumberUninfected_HighRiskHETs,2),0.00000000001);
                
                Results.ann_PctSusOnPrEP_HighRiskMSM = ...
                    sum(Results.ann_NumberOnPrEP_HighRiskMSM,2) ./ max(sum(Results.ann_NumberUninfected_HighRiskMSM,2),0.00000000001);
                
                Results.ann_PctSusOnPrEP_HighRiskIDUs = ...
                    sum(Results.ann_NumberOnPrEP_HighRiskIDUs,2) ./ max(sum(Results.ann_NumberUninfected_HighRiskIDUs,2),0.00000000001);
                
                Results.ann_PctSusOnPrEP_HighRiskBlk = ...
                    sum(Results.ann_NumberOnPrEP_Blk,2) ./ max(sum(Results.ann_NumberUninfected_HighRiskBlk,2),0.00000000001);
                
                Results.ann_PctSusOnPrEP_HighRiskHisp = ...
                    sum(Results.ann_NumberOnPrEP_Hisp,2) ./ max(sum(Results.ann_NumberUninfected_HighRiskHisp,2),0.00000000001);
                
                Results.ann_PctSusOnPrEP_HighRiskOth = ...
                    sum(Results.ann_NumberOnPrEP_Oth,2) ./ max(sum(Results.ann_NumberUninfected_HighRiskOth,2),0.00000000001);

                Results.ann_PctSusOnPrEP_HighRiskHETs_ByRaceEth = ...
                    Results.ann_NumberOnPrEP_HighRiskHETs * Params.raceIndicator(:,:)...
                    ./ max(Results.ann_NumberUninfected_HighRiskHETs * Params.raceIndicator(:,:),0.00000000001);

                Results.ann_PctSusOnPrEP_HighRiskMSM_ByRaceEth = ...
                    Results.ann_NumberOnPrEP_HighRiskMSM * Params.raceIndicator(:,:)...
                    ./ max(Results.ann_NumberUninfected_HighRiskMSM * Params.raceIndicator(:,:),0.00000000001);

                Results.ann_PctSusOnPrEP_HighRiskIDUs_ByRaceEth = ...
                    Results.ann_NumberOnPrEP_HighRiskIDUs * Params.raceIndicator(:,:)...
                    ./ max(Results.ann_NumberUninfected_HighRiskIDUs * Params.raceIndicator(:,:),0.00000000001);

                % Rate of PrEP initiation 
            
                ann_NumberEligForPrEPNotOnPrEP_HighRiskHETs = ...
                    Results.ann_sumTotalNumberEligForPrEPNotOnPrEP * Params.HRHIndicator;
                ann_NumberEligForPrEPNotOnPrEP_HighRiskMSM = ...
                    Results.ann_sumTotalNumberEligForPrEPNotOnPrEP * (Params.HighRiskMSMIndicator.* AllocPopsTargetedIndicator);
                ann_NumberEligForPrEPNotOnPrEP_IDUs = ...
                    Results.ann_sumTotalNumberEligForPrEPNotOnPrEP * Params.popIndicator(:,Params.pop_IDU);  
                
                Results.ann_PrEPInitProb_HighRiskHETs = ...
                    RateToProb(sum(Results.ann_StartPrEP_HighRiskHETs,2) ./ ...
                    max(ann_NumberEligForPrEPNotOnPrEP_HighRiskHETs,1) ./ ...
                    (1./stepsPerYear.'));                                    
                Results.ann_PrEPInitProb_HighRiskMSM = ...
                    RateToProb(sum(Results.ann_StartPrEP_HighRiskMSM .* AllocPopsTargetedIndicator',2) ./ ...
                    max(ann_NumberEligForPrEPNotOnPrEP_HighRiskMSM,1) ./ ...
                    (1./stepsPerYear.'));                   
                Results.ann_PrEPInitProb_IDUs = ...
                    RateToProb(sum(Results.ann_StartPrEP_HighRiskIDUs,2) ./ ...
                    max(ann_NumberEligForPrEPNotOnPrEP_IDUs,1) ./ ...
                    (1./stepsPerYear.'));  
        end
        
        % Outcomes for optimization for RA: Rates of progression targeted by interventions
        for sectionRAOptimization = 1:1
            
            % Allocation to each intervention
            
                Results.alloc_Testing_HET_Low = Params.intn_Testing_Investment_HET_Low;
                Results.alloc_Testing_HET_High = Params.intn_Testing_Investment_HET_High;
                Results.alloc_Testing_MSM_Low = Params.intn_Testing_Investment_MSM_Low;
                Results.alloc_Testing_MSM_High = Params.intn_Testing_Investment_MSM_High;
                Results.alloc_Testing_IDU = Params.intn_Testing_Investment_IDU;
                
                Results.alloc_LTCatDiag = Params.intn_LTCatDiag_Investment;
                Results.alloc_LTCafterDiag = Params.intn_LTCafterDiag_Investment;
                
                Results.alloc_ARTInitiation = Params.intn_ARTInitiation_Investment;
                Results.alloc_ARTAdher5to4 = Params.intn_ARTAdher5to4_Investment;
                Results.alloc_ARTAdher4to5 = Params.intn_ARTAdher4to5_Investment;
                
                Results.alloc_SEP_B = Params.intn_SEP_Investment_B;
                Results.alloc_SEP_H = Params.intn_SEP_Investment_H;
                Results.alloc_SEP_O = Params.intn_SEP_Investment_O;

                Results.alloc_PrEP_Oral_HETM_B = Params.intn_PrEP_Oral_Investment_HETM_B;
                Results.alloc_PrEP_Oral_HETM_H = Params.intn_PrEP_Oral_Investment_HETM_H;
                Results.alloc_PrEP_Oral_HETM_O = Params.intn_PrEP_Oral_Investment_HETM_O;
                Results.alloc_PrEP_Oral_HETF_B = Params.intn_PrEP_Oral_Investment_HETF_B;
                Results.alloc_PrEP_Oral_HETF_H = Params.intn_PrEP_Oral_Investment_HETF_H;
                Results.alloc_PrEP_Oral_HETF_O = Params.intn_PrEP_Oral_Investment_HETF_O;
                Results.alloc_PrEP_Oral_MSM_B = Params.intn_PrEP_Oral_Investment_MSM_B;
                Results.alloc_PrEP_Oral_MSM_H = Params.intn_PrEP_Oral_Investment_MSM_H;
                Results.alloc_PrEP_Oral_MSM_O = Params.intn_PrEP_Oral_Investment_MSM_O;
                Results.alloc_PrEP_Oral_IDU_B = Params.intn_PrEP_Oral_Investment_IDU_B;
                Results.alloc_PrEP_Oral_IDU_H = Params.intn_PrEP_Oral_Investment_IDU_H;
                Results.alloc_PrEP_Oral_IDU_O = Params.intn_PrEP_Oral_Investment_IDU_O;
                Results.alloc_PrEP_Inject_HETM_B = Params.intn_PrEP_Inject_Investment_HETM_B;
                Results.alloc_PrEP_Inject_HETM_H = Params.intn_PrEP_Inject_Investment_HETM_H;
                Results.alloc_PrEP_Inject_HETM_O = Params.intn_PrEP_Inject_Investment_HETM_O;
                Results.alloc_PrEP_Inject_HETF_B = Params.intn_PrEP_Inject_Investment_HETF_B;
                Results.alloc_PrEP_Inject_HETF_H = Params.intn_PrEP_Inject_Investment_HETF_H;
                Results.alloc_PrEP_Inject_HETF_O = Params.intn_PrEP_Inject_Investment_HETF_O;
                Results.alloc_PrEP_Inject_MSM_B = Params.intn_PrEP_Inject_Investment_MSM_B;
                Results.alloc_PrEP_Inject_MSM_H = Params.intn_PrEP_Inject_Investment_MSM_H;
                Results.alloc_PrEP_Inject_MSM_O = Params.intn_PrEP_Inject_Investment_MSM_O;
                Results.alloc_PrEP_Inject_IDU_B = Params.intn_PrEP_Inject_Investment_IDU_B;
                Results.alloc_PrEP_Inject_IDU_H = Params.intn_PrEP_Inject_Investment_IDU_H;
                Results.alloc_PrEP_Inject_IDU_O = Params.intn_PrEP_Inject_Investment_IDU_O;                
                
            % Incidence outcomes
            
                Results.ann_undiscNewInfectionsHET = sum(Results.ann_NewInfectionsHET,2);
                Results.ann_undiscNewInfectionsMSM = sum(Results.ann_NewInfectionsMSM,2);
                Results.ann_undiscNewInfectionsYMSM = sum(Results.ann_NewInfectionsMSM .* AllocPopsTargetedIndicator',2);
                Results.ann_undiscNewInfectionsIDU = sum(Results.ann_NewInfectionsIDU,2);
                Results.ann_undiscNewInfections = sum(Results.ann_TotalNewInfections,2);
                Results.total_undiscNewInfections = sum(sum(Results.ann_TotalNewInfections,2),2);
                
            % Incidence outcomes for cost minimization optimization for
            % reaching EHE goals (added 08/20/2020)
            OutcomesIndexTarget0=max(Params.OutcomesIndex_AllocationStartYr(1)-1,1);
            OutcomesIndexTarget1=Params.OutcomesIndex_AllocationStartYr(1)+Params.minCost_Tgt1NumYears-1;
            OutcomesIndexTarget2=Params.OutcomesIndex_AllocationStartYr(1)+Params.minCost_Tgt2NumYears-1;            
            
            if FirstOutcomeYr <= Params.tt_periodFiveStartYear - OutcomesIndexTarget0 && (LastOutcomeYr >= Params.tt_periodFiveStartYear + OutcomesIndexTarget1 - OutcomesIndexTarget0 - 1)
                
                IncBase = sum(Results.ann_TotalNewInfections(OutcomesIndexTarget0,:));
                IncTG1 = sum(Results.ann_TotalNewInfections(OutcomesIndexTarget1,:));
                
                Results.IncidenceReductionTarget1 = -(IncTG1 / IncBase - 1);                
                
            end            
            
            if FirstOutcomeYr <= Params.tt_periodFiveStartYear - OutcomesIndexTarget0 && (LastOutcomeYr >= Params.tt_periodFiveStartYear + OutcomesIndexTarget2 - OutcomesIndexTarget0 - 1)
                               
                IncTG2 = sum(Results.ann_TotalNewInfections(OutcomesIndexTarget2,:));
                
                Results.IncidenceReductionTarget2 = -(IncTG2 / IncBase - 1);
                
            end
            
                %Total discounted infections calculated earlier
                
            % QALYs    
                
                Results.total_discQALYs_pop = sum(Results.ann_QALYs_pop,1);
                Results.ann_discQALYs = sum(Results.ann_QALYs,2);
                Results.ann_discQALYs_YMSM = sum(Results.ann_discQALYs .* AllocPopsTargetedIndicator',2);
            
            % HIV Treatment and care costs
                
                Results.ann_ARTCarePrEPCost_Disc_agg = sum(Results.ann_ARTCarePrEPCost_Disc,2);
            
            % Testing rates (tr_LRHAllEligTestRate, etc.) calc'ed above in testing freq section
           % Testing Rates - only for selected year
            %All Eligible (HIV- and HIV+)
            Results.tr_WeightedAllEligTestRateByCohort_AllYrs = TTProg.WeightedTestRateAllEligByCohort;
            Results.tr_WeightedAllEligTestRate_AllYrs = TTProg.WeightedTestRateAllElig;
            Results.tr_raceAllEligTestRate_AllYrs = TTProg.RaceSpecificAllEligTestRate;
            Results.tr_riskGroupAllEligTestRate_AllYrs = TTProg.RiskGroupAllEligTestRate;
            Results.tr_HRHAllEligTestProb_AllYrs = RateToProb(TTProg.HRHAllEligTestRate);
            Results.tr_LRHAllEligTestProb_AllYrs = RateToProb(TTProg.LRHAllEligTestRate);
            Results.tr_HRMSMAllEligTestProb_AllYrs = RateToProb(TTProg.HRMSMAllEligTestRate);
            Results.tr_HRYMSMAllEligTestProb_AllYrs = RateToProb(TTProg.HRYMSMAllEligTestRate);
            Results.tr_LRMSMAllEligTestProb_AllYrs = RateToProb(TTProg.LRMSMAllEligTestRate);
            Results.tr_LRYMSMAllEligTestProb_AllYrs = RateToProb(TTProg.LRYMSMAllEligTestRate);
            Results.tr_IDUAllEligTestProb_AllYrs = RateToProb(TTProg.IDUAllEligTestRate);
            Results.tr_AllEligTestRateByRiskGpLevel_AllYrs = TTProg.AllEligTestRateByRiskGpLevel;

            % HIV+ only       
            Results.tr_WeightedHIVPosTestRateByCohort_AllYrs = TTProg.WeightedTestRateHIVPosByCohort;
            Results.tr_raceHIVPosTestRate_AllYrs = TTProg.RaceSpecificHIVPosTestRate;
            Results.tr_riskGroupHIVPosTestRate_AllYrs = TTProg.RiskGroupHIVPosTestRate;
            Results.tr_HRHHIVPosTestRate_AllYrs = TTProg.HRHHIVPosTestRate;
            Results.tr_LRHHIVPosTestRate_AllYrs = TTProg.LRHHIVPosTestRate;
            Results.tr_uninfectedTestRate_AllYrs = TTProg.OverallUninfectedTestRate;
            Results.tr_diseaseStageTestRate_AllYrs = TTProg.DiseaseStageTestRate;

            % Diagnosis rates
            Results.tr_DiagRateAll_AllYrs = TTProg.DiagRateAll;
            Results.tr_raceDiagRate_AllYrs = TTProg.RaceSpecificDiagRate;
            Results.tr_riskGroupDiagRate_AllYrs = TTProg.RiskGroupDiagRate;
            Results.tr_diseaseStageDiagRate_AllYrs = TTProg.DiseaseStageDiagRate;
            Results.tr_LRHDiagRate_AllYrs = TTProg.LowRiskHETDiagRate;
            Results.tr_HRHDiagRate_AllYrs = TTProg.HRHDiagRate;

            % Linkage rates at diagnosis
                Results.ann_PctLinkToCareFirst = ...
                    sum(Results.ann_TotalLinkToCareFirst .* AllocPopsTargetedIndicator',2) ./ ...
                    max(sum((Results.ann_numNotifiedConvPos + Results.ann_numNotifiedRapidPos) .* AllocPopsTargetedIndicator',2),0.000000000000000001);

            % Linkage Probs after diagnosis
                Results.ann_LinkageAfterProb = ...
                    RateToProb(sum(Results.ann_TotalLinkToCareAfterFirst .* AllocPopsTargetedIndicator',2) ./ ...
                    max(sum(Results.ann_sumNumberCCStage2 .* AllocPopsTargetedIndicator',2),1) ./ ...
                    (1./stepsPerYear.'));
                
            % Collect outcomes needed to observe ART prescription
            % probability by race
            
            for yearCount = 1:numYears
                
                Results.ann_TotalStartARTNotVLS_Blk(yearCount,:) = ...
                    Results.ann_TotalStartARTNotVLS(yearCount,:)' .* Params.raceIndicator(:,Params.race_B);
                Results.ann_TotalStartARTNotVLS_Hisp(yearCount,:) = ...
                    Results.ann_TotalStartARTNotVLS(yearCount,:)' .* Params.raceIndicator(:,Params.race_H);
                Results.ann_TotalStartARTNotVLS_Oth(yearCount,:) = ...
                    Results.ann_TotalStartARTNotVLS(yearCount,:)' .* Params.raceIndicator(:,Params.race_O);
                
                Results.ann_BecomeVLSFromLTC_Blk(yearCount,:) = ...
                    Results.ann_BecomeVLSFromLTC(yearCount,:)' .* Params.raceIndicator(:,Params.race_B);
                Results.ann_BecomeVLSFromLTC_Hisp(yearCount,:) = ...
                    Results.ann_BecomeVLSFromLTC(yearCount,:)' .* Params.raceIndicator(:,Params.race_H);
                Results.ann_BecomeVLSFromLTC_Oth(yearCount,:) = ...
                    Results.ann_BecomeVLSFromLTC(yearCount,:)' .* Params.raceIndicator(:,Params.race_O);
                
                Results.ann_sumNumberCCStage3_Blk(yearCount,:) = ...
                    Results.ann_sumNumberCCStage3(yearCount,:)' .* Params.raceIndicator(:,Params.race_B);
                Results.ann_sumNumberCCStage3_Hisp(yearCount,:) = ...
                    Results.ann_sumNumberCCStage3(yearCount,:)' .* Params.raceIndicator(:,Params.race_H);
                Results.ann_sumNumberCCStage3_Oth(yearCount,:) = ...
                    Results.ann_sumNumberCCStage3(yearCount,:)' .* Params.raceIndicator(:,Params.race_O);
                
            end
            
                
            % ART prescription Probs from LTC
                Results.ann_ARTPrescrProb = ...
                    RateToProb(sum((Results.ann_TotalStartARTNotVLS + Results.ann_BecomeVLSFromLTC) .* AllocPopsTargetedIndicator',2) ./ ...
                    max(sum(Results.ann_sumNumberCCStage3 .* AllocPopsTargetedIndicator',2),1) ./ ...
                    (1./stepsPerYear.'));
                
            % ART prescription Probs from LTC (by race)
                Results.ann_ARTPrescrProb_Blk = ...
                    RateToProb(sum(Results.ann_TotalStartARTNotVLS_Blk + Results.ann_BecomeVLSFromLTC_Blk,2) ./ ...
                    max(sum(Results.ann_sumNumberCCStage3_Blk,2),1) ./ ...
                    (1./stepsPerYear.'));
                
                Results.ann_ARTPrescrProb_Hisp = ...
                    RateToProb(sum(Results.ann_TotalStartARTNotVLS_Hisp + Results.ann_BecomeVLSFromLTC_Hisp,2) ./ ...
                    max(sum(Results.ann_sumNumberCCStage3_Hisp,2),1) ./ ...
                    (1./stepsPerYear.'));
                
                Results.ann_ARTPrescrProb_Oth = ...
                    RateToProb(sum(Results.ann_TotalStartARTNotVLS_Oth + Results.ann_BecomeVLSFromLTC_Oth,2) ./ ...
                    max(sum(Results.ann_sumNumberCCStage3_Oth,2),1) ./ ...
                    (1./stepsPerYear.'));
                                
            % Prob of departure from VLS
                Results.ann_VLSToANVProb = ...
                    RateToProb(sum(Results.ann_TotalDropOff_VLSToANV .* AllocPopsTargetedIndicator',2) ./ ...
                    max(sum(Results.ann_sumNumberVLS .* AllocPopsTargetedIndicator',2),1) ./ ...
                    (1./stepsPerYear.'));    

            % Prob of movement from ANV to VLS
                Results.ann_ANVToVLSProb = ...
                    RateToProb(sum(Results.ann_BecomeVLSFromANV .* AllocPopsTargetedIndicator',2) ./ ...
                    max(sum(Results.ann_sumNumberCCStage4 .* AllocPopsTargetedIndicator',2),1) ./ ...
                    (1./stepsPerYear.'));   
                
            % Percent of PWID being served by SEP. Updated 3/9/23 by Bates to only
            % include uninfected for overall. Names changed to reflect
                Results.ann_PctUninfPWIDServedbySEP = ...
                    min(Results.ann_NumServedbySEP ./ ...
                    max(sum(Results.ann_NumberUninfected .* repmat(Params.popIndicator(:,Params.pop_IDU)',numYears,1) * Params.behav_pctActivePWID,2),1),1);                
                
                Results.ann_PctUninfPWIDServedbySEP_Blk = ...
                    min(Results.ann_NumServedbySEP_B ./ ...
                    max(sum(Results.ann_NumberUninfected .* repmat((Params.popIndicator(:,Params.pop_IDU) .* Params.raceIndicator(:,Params.race_B))',numYears,1) * Params.behav_pctActivePWID,2),1),1);
                
                Results.ann_PctUninfPWIDServedbySEP_Hisp = ...
                    min(Results.ann_NumServedbySEP_H ./ ...
                    max(sum(Results.ann_NumberUninfected .* repmat((Params.popIndicator(:,Params.pop_IDU) .* Params.raceIndicator(:,Params.race_H))',numYears,1) * Params.behav_pctActivePWID,2),1),1);
                
                Results.ann_PctUninfPWIDServedbySEP_Oth = ...
                    min(Results.ann_NumServedbySEP_O ./ ...
                    max(sum(Results.ann_NumberUninfected .* repmat((Params.popIndicator(:,Params.pop_IDU) .* Params.raceIndicator(:,Params.race_O))',numYears,1) * Params.behav_pctActivePWID,2),1),1);
            

            % Rate of PrEP initiation - calc'ed above in PrEP section
            
            
            % Distribution of PLWH among continuum of HIV Care
            
                Results.ann_PctAware = ...
                    sum(Results.ann_NumberAware .* AllocPopsTargetedIndicator',2) ./ ...
                    sum(Results.ann_HIVPrevalence .* AllocPopsTargetedIndicator',2);
                
                for yearCount = 1:numYears
                    
                    %13-17
                    ann_NumberAware_13_17(yearCount,:) = ...
                        Results.ann_NumberAware(yearCount,:) .* ...
                        Params.ageIndicator(:,Params.age_1)';
                    
                    %18-24
                    ann_NumberAware_18_24(yearCount,:) = ...
                        Results.ann_NumberAware(yearCount,:) .* ...
                        Params.ageIndicator(:,Params.age_2)';
                    
                    %25-34
                    ann_NumberAware_25_34(yearCount,:) = ...
                        Results.ann_NumberAware(yearCount,:) .* ...
                        Params.ageIndicator(:,Params.age_3)';
                    
                    %35-44
                    ann_NumberAware_35_44(yearCount,:) = ...
                        Results.ann_NumberAware(yearCount,:) .* ...
                        Params.ageIndicator(:,Params.age_4)';
                    
                    %45-54
                    ann_NumberAware_45_54(yearCount,:) = ...
                        Results.ann_NumberAware(yearCount,:) .* ...
                        Params.ageIndicator(:,Params.age_5)';
                 
                    %55-64
                    ann_NumberAware_55_64(yearCount,:) = ...
                        Results.ann_NumberAware(yearCount,:) .* ...
                        Params.ageIndicator(:,Params.age_6)';
                    
                    %65+
                    ann_NumberAware_65(yearCount,:) = ...
                        Results.ann_NumberAware(yearCount,:) .* ...
                        Params.ageIndicator(:,Params.age_7)';
                end
                Results.ann_PctAware_13_17 = ...
                    sum(ann_NumberAware_13_17,2) ./ ...
                    sum(Results.ann_HIVPrevalence_13_17,2);
                Results.ann_PctAware_18_24 = ...
                    sum(ann_NumberAware_18_24,2) ./ ...
                    sum(Results.ann_HIVPrevalence_18_24,2);
                Results.ann_PctAware_25_34 = ...
                    sum(ann_NumberAware_25_34,2) ./ ...
                    sum(Results.ann_HIVPrevalence_25_34,2);
                Results.ann_PctAware_35_44 = ...
                    sum(ann_NumberAware_35_44,2) ./ ...
                    sum(Results.ann_HIVPrevalence_35_44,2);
                Results.ann_PctAware_45_54 = ...
                    sum(ann_NumberAware_45_54,2) ./ ...
                    sum(Results.ann_HIVPrevalence_45_54,2);
                Results.ann_PctAware_55_64 = ...
                    sum(ann_NumberAware_55_64,2) ./ ...
                    sum(Results.ann_HIVPrevalence_55_64,2);
                Results.ann_PctAware_65 = ...
                    sum(ann_NumberAware_65,2) ./ ...
                    sum(Results.ann_HIVPrevalence_65,2);
                
                % by race/ethnicity
                
                Results.ann_PctAware_Blk = ...
                    sum(Results.ann_NumberAware_Blk,2) ./ ...
                    sum(Results.ann_HIVPrevalence_Blk,2);
                
                Results.ann_PctAware_Hisp = ...
                    sum(Results.ann_NumberAware_Hisp,2) ./ ...
                    sum(Results.ann_HIVPrevalence_Hisp,2);
                
                Results.ann_PctAware_Oth = ...
                    sum(Results.ann_NumberAware_Oth,2) ./ ...
                    sum(Results.ann_HIVPrevalence_Oth,2);
                
                % by transmission/risk group
                
                Results.ann_PctAware_HRHET = ...
                    sum(Results.ann_NumberAware_HighRiskHETs,2) ./ ...
                    sum(Results.ann_HIVPrevalence_HighRiskHETs,2);
                
                Results.ann_PctAware_OHET = ...
                    sum(Results.ann_NumberAware_LowRiskHETs,2) ./ ...
                    sum(Results.ann_HIVPrevalence_LowRiskHETs,2);
                
                Results.ann_PctAware_HRMSM = ...
                    sum(Results.ann_NumberAware_HighRiskMSMs .* AllocPopsTargetedIndicator',2) ./ ...
                    sum(Results.ann_HIVPrevalence_HighRiskMSMs .* AllocPopsTargetedIndicator',2);
                
                Results.ann_PctAware_OMSM = ...
                    sum(Results.ann_NumberAware_LowRiskMSMs .* AllocPopsTargetedIndicator',2) ./ ...
                    sum(Results.ann_HIVPrevalence_LowRiskMSMs .* AllocPopsTargetedIndicator',2);
                
                Results.ann_PctAware_PWID = ...
                    sum(Results.ann_NumberAware_IDU,2) ./ ...
                    sum(Results.ann_HIVPrevalence_IDU,2);
                
                Results.ann_PctInCare = ...
                    sum(Results.ann_NumberInCare .* AllocPopsTargetedIndicator',2) ./ ...
                    sum(Results.ann_HIVPrevalence .* AllocPopsTargetedIndicator',2);
                
                Results.ann_PctOnART = ...
                    sum(Results.ann_NumberOnART .* AllocPopsTargetedIndicator',2) ./ ...
                    sum(Results.ann_HIVPrevalence .* AllocPopsTargetedIndicator',2);
                
                Results.ann_PctVLS = ...
                    sum(Results.ann_NumberVLS .* AllocPopsTargetedIndicator',2) ./ ...
                    sum(Results.ann_HIVPrevalence .* AllocPopsTargetedIndicator',2);
                
                Results.ann_PctVLSamongdiag = ...
                sum(Results.ann_NumberVLS .* AllocPopsTargetedIndicator',2) ./ ...
                sum(Results.ann_NumberAware .* AllocPopsTargetedIndicator',2);
            
                for yearCount = 1:numYears
                    
                    %13-17
                    ann_NumberVLS_13_17(yearCount,:) = ...
                        Results.ann_NumberVLS(yearCount,:) .* ...
                        Params.ageIndicator(:,Params.age_1)';
                    
                    %18-24
                    ann_NumberVLS_18_24(yearCount,:) = ...
                        Results.ann_NumberVLS(yearCount,:) .* ...
                        Params.ageIndicator(:,Params.age_2)';
                    
                    %25-34
                    ann_NumberVLS_25_34(yearCount,:) = ...
                        Results.ann_NumberVLS(yearCount,:) .* ...
                        Params.ageIndicator(:,Params.age_3)';
                    
                    %35-44
                    ann_NumberVLS_35_44(yearCount,:) = ...
                        Results.ann_NumberVLS(yearCount,:) .* ...
                        Params.ageIndicator(:,Params.age_4)';
                    
                    %45-54
                    ann_NumberVLS_45_54(yearCount,:) = ...
                        Results.ann_NumberVLS(yearCount,:) .* ...
                        Params.ageIndicator(:,Params.age_5)';
                 
                    %55-64
                    ann_NumberVLS_55_64(yearCount,:) = ...
                        Results.ann_NumberVLS(yearCount,:) .* ...
                        Params.ageIndicator(:,Params.age_6)';
                    
                    %65+
                    ann_NumberVLS_65(yearCount,:) = ...
                        Results.ann_NumberVLS(yearCount,:) .* ...
                        Params.ageIndicator(:,Params.age_7)';
                end
                Results.ann_PctVLSamongdiag_13_17 = ...
                    sum(ann_NumberVLS_13_17,2) ./ ...
                    sum(ann_NumberAware_13_17,2);
                Results.ann_PctVLSamongdiag_18_24 = ...
                    sum(ann_NumberVLS_18_24,2) ./ ...
                    sum(ann_NumberAware_18_24,2);
                Results.ann_PctVLSamongdiag_25_34 = ...
                    sum(ann_NumberVLS_25_34,2) ./ ...
                    sum(ann_NumberAware_25_34,2);
                Results.ann_PctVLSamongdiag_35_44 = ...
                    sum(ann_NumberVLS_35_44,2) ./ ...
                    sum(ann_NumberAware_35_44,2);
                Results.ann_PctVLSamongdiag_45_54 = ...
                    sum(ann_NumberVLS_45_54,2) ./ ...
                    sum(ann_NumberAware_45_54,2);
                Results.ann_PctVLSamongdiag_55_64 = ...
                    sum(ann_NumberVLS_55_64,2) ./ ...
                    sum(ann_NumberAware_55_64,2);
                Results.ann_PctVLSamongdiag_65 = ...
                    sum(ann_NumberVLS_65,2) ./ ...
                    sum(ann_NumberAware_65,2);
                
                % Percent VLS among diagnosed by race/ethnicity
                
                Results.ann_PctVLSamongdiag_Blk = ...
                    sum(Results.ann_NumberVLS_Blk,2) ./ ...
                    sum(Results.ann_NumberAware_Blk,2);
                Results.ann_PctVLSamongdiag_Hisp = ...
                    sum(Results.ann_NumberVLS_Hisp,2) ./ ...
                    sum(Results.ann_NumberAware_Hisp,2);
                Results.ann_PctVLSamongdiag_Oth = ...
                    sum(Results.ann_NumberVLS_Oth,2) ./ ...
                    sum(Results.ann_NumberAware_Oth,2);
                
                % Percent VLS among diagnosed by race/ethnicity
                
                Results.ann_PctVLSamongdiag_HET = ...
                    sum(Results.ann_NumberVLS_HET,2) ./ ...
                    sum(Results.ann_NumberAware_HET,2);
                Results.ann_PctVLSamongdiag_MSM = ...
                    sum(Results.ann_NumberVLS_MSM,2) ./ ...
                    sum(Results.ann_NumberAware_MSM,2);
                Results.ann_PctVLSamongdiag_IDU = ...
                    sum(Results.ann_NumberVLS_IDU,2) ./ ...
                    sum(Results.ann_NumberAware_IDU,2);
                
            % Percent Eligible on PrEP
            
                Results.ann_PctEligOnPrEP_HRHET = ...
                    sum(Results.ann_NumberOnPrEP_HighRiskHETs,2) ./ ...
                    sum(Results.ann_NumberEligForPrEP * Params.HRHIndicator,2);
                
                Results.ann_PctEligOnPrEP_HRMSM = ...
                    sum(Results.ann_NumberOnPrEP_HighRiskMSM .* AllocPopsTargetedIndicator',2) ./ ...
                    sum(Results.ann_NumberEligForPrEP .* AllocPopsTargetedIndicator' * (Params.riskLevelIndicator(:,Params.risk_Casual).*Params.popIndicator(:,Params.pop_MSM)),2);
                
                Results.ann_PctEligOnPrEP_PWID = ...
                    sum(Results.ann_NumberOnPrEP_HighRiskIDUs,2) ./ ...
                    sum(Results.ann_NumberEligForPrEP * Params.popIndicator(:,Params.pop_IDU),2);
                
            
        end
        
        % Outcomes for young MSM analysis
        for sectionYMSM = 1:1     
            
            % YMSM incidence by age group and race/eth
            Results.YMSM_NewInfections_MSM_18_24_B = ...
                sum(sum(Results.ann_TotalNewInfections,1)' ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.ageIndicator (:,Params.age_2) ...
                .*Params.raceIndicator (:, Params.race_B));
            
            Results.YMSM_NewInfections_MSM_18_24_H = ...
                sum(sum(Results.ann_TotalNewInfections,1)' ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.ageIndicator (:,Params.age_2) ...
                .*Params.raceIndicator (:, Params.race_H));
            
            Results.YMSM_NewInfections_MSM_18_24_O = ...
                sum(sum(Results.ann_TotalNewInfections,1)' ...
                .*Params.popSexIndicator(:,Params.popSex_MSM) ...
                .*Params.ageIndicator (:,Params.age_2) ...
                .*Params.raceIndicator (:, Params.race_O));
            

            Results.YMSM_NewInfections_MSM_25_34_B = ...
                sum(sum(Results.ann_TotalNewInfections,1)' ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.ageIndicator (:,Params.age_3)... 
                .*Params.raceIndicator (:,Params.race_B));
            
            Results.YMSM_NewInfections_MSM_25_34_H = ...
                sum(sum(Results.ann_TotalNewInfections,1)' ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.ageIndicator (:,Params.age_3)... 
                .*Params.raceIndicator (:,Params.race_H));
            
            Results.YMSM_NewInfections_MSM_25_34_O = ...
                sum(sum(Results.ann_TotalNewInfections,1)' ...
                .*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.ageIndicator (:,Params.age_3)... 
                .*Params.raceIndicator (:,Params.race_O));

            % Number of HIV-infected YMSM by age and race/eth in last year
            % of outcomes collected
            YMSM_NumAwareLastYr_18_24_MSM_B = ...
                sum(Results.ann_NumberAware(numYears,:)'.*...
                    Params.ageIndicator(:,Params.age_2) ...
                    .*Params.popSexIndicator(:,Params.popSex_MSM)...
                    .*Params.raceIndicator (:, Params.race_B));
            YMSM_NumAwareLastYr_18_24_MSM_H = ...
                sum(Results.ann_NumberAware(numYears,:)'.*...
                    Params.ageIndicator(:,Params.age_2) ...
                    .*Params.popSexIndicator(:,Params.popSex_MSM)...
                    .*Params.raceIndicator (:, Params.race_H));
            YMSM_NumAwareLastYr_18_24_MSM_O = ...
                sum(Results.ann_NumberAware(numYears,:)'.*...
                    Params.ageIndicator(:,Params.age_2) ...
                    .*Params.popSexIndicator(:,Params.popSex_MSM)...
                    .*Params.raceIndicator (:, Params.race_O));                

            YMSM_NumAwareLastYr_25_34_MSM_B = ...
                sum(Results.ann_NumberAware(numYears,:)'.*...
                    Params.ageIndicator(:,Params.age_3) ...
                    .*Params.popSexIndicator(:,Params.popSex_MSM)...
                    .*Params.raceIndicator (:, Params.race_B));                
            YMSM_NumAwareLastYr_25_34_MSM_H = ...
                sum(Results.ann_NumberAware(numYears,:)'.*...
                    Params.ageIndicator(:,Params.age_3) ...
                    .*Params.popSexIndicator(:,Params.popSex_MSM)...
                    .*Params.raceIndicator (:, Params.race_H));
            YMSM_NumAwareLastYr_25_34_MSM_O = ...
                sum(Results.ann_NumberAware(numYears,:)'.*...
                    Params.ageIndicator(:,Params.age_3) ...
                    .*Params.popSexIndicator(:,Params.popSex_MSM)...
                    .*Params.raceIndicator (:, Params.race_O)); 
            
            YMSM_NumAwareLastYr_HET = ...
                sum(Results.ann_NumberAware(numYears,:)'.*...
                    Params.popIndicator(:,Params.pop_HET));
            YMSM_NumAwareLastYr_MSM = ...
                sum(Results.ann_NumberAware(numYears,:)'.*...
                    Params.popIndicator(:,Params.pop_MSM));
            YMSM_NumAwareLastYr_IDU = ...
                sum(Results.ann_NumberAware(numYears,:)'.*...
                    Params.popIndicator(:,Params.pop_IDU));
                
            % Percentage of HIV-infected YMSM aware of status by age
            % and race/eth in last year of outcomes collected                
            Results.YMSM_PctAwareLastYr_18_24_MSM_B = ...
                YMSM_NumAwareLastYr_18_24_MSM_B / ...
                sum(Results.ann_HIVPrevalence(numYears,:)'.*...
                    Params.ageIndicator(:,Params.age_2)...
                    .*Params.popSexIndicator(:,Params.popSex_MSM)...
                    .*Params.raceIndicator (:, Params.race_B));
            Results.YMSM_PctAwareLastYr_18_24_MSM_H = ...
                YMSM_NumAwareLastYr_18_24_MSM_H / ...
                sum(Results.ann_HIVPrevalence(numYears,:)'.*...
                    Params.ageIndicator(:,Params.age_2)...
                    .*Params.popSexIndicator(:,Params.popSex_MSM)...
                    .*Params.raceIndicator (:, Params.race_H));
            Results.YMSM_PctAwareLastYr_18_24_MSM_O = ...
                YMSM_NumAwareLastYr_18_24_MSM_O / ...
                sum(Results.ann_HIVPrevalence(numYears,:)'.*...
                    Params.ageIndicator(:,Params.age_2)...
                    .*Params.popSexIndicator(:,Params.popSex_MSM)...
                    .*Params.raceIndicator (:, Params.race_O)); 
            Results.YMSM_PctAwareLastYr_25_34_MSM_B = ...
                YMSM_NumAwareLastYr_25_34_MSM_B / ...
                sum(Results.ann_HIVPrevalence(numYears,:)'.*...
                    Params.ageIndicator(:,Params.age_3)...
                    .*Params.popSexIndicator(:,Params.popSex_MSM)...
                    .*Params.raceIndicator (:, Params.race_B));
            Results.YMSM_PctAwareLastYr_25_34_MSM_H = ...
                YMSM_NumAwareLastYr_25_34_MSM_H / ...
                sum(Results.ann_HIVPrevalence(numYears,:)'.*...
                    Params.ageIndicator(:,Params.age_3)...
                    .*Params.popSexIndicator(:,Params.popSex_MSM)...
                    .*Params.raceIndicator (:, Params.race_H));
            Results.YMSM_PctAwareLastYr_25_34_MSM_O = ...
                YMSM_NumAwareLastYr_25_34_MSM_O / ...
                sum(Results.ann_HIVPrevalence(numYears,:)'.*...
                    Params.ageIndicator(:,Params.age_3)...
                    .*Params.popSexIndicator(:,Params.popSex_MSM)...
                    .*Params.raceIndicator (:, Params.race_O));

            Results.YMSM_PctAwareLastYr_HET = ...
                YMSM_NumAwareLastYr_HET / ...
                sum(Results.ann_HIVPrevalence(numYears,:)'.*...
                    Params.popIndicator(:,Params.pop_HET));
            Results.YMSM_PctAwareLastYr_MSM = ...
                YMSM_NumAwareLastYr_MSM / ...
                sum(Results.ann_HIVPrevalence(numYears,:)'.*...
                    Params.popIndicator(:,Params.pop_MSM));
            Results.YMSM_PctAwareLastYr_IDU = ...
                YMSM_NumAwareLastYr_IDU / ...
                sum(Results.ann_HIVPrevalence(numYears,:)'.*...
                    Params.popIndicator(:,Params.pop_IDU));

            % Percentage of HIV-infected YMSM aware of status by age
            % and race/eth in last year of outcomes collected                   
            Results.YMSM_VLSAmongDiagLastYr_18_24_MSM_B = ...
                sum(Results.ann_NumberVLS(numYears,:)'.*...
                Params.ageIndicator(:,Params.age_2).*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator(:, Params.race_B))...
                / YMSM_NumAwareLastYr_18_24_MSM_B;
            Results.YMSM_VLSAmongDiagLastYr_18_24_MSM_H = ...
                sum(Results.ann_NumberVLS(numYears,:)'.*...
                Params.ageIndicator(:,Params.age_2).*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator(:, Params.race_H))...
                / YMSM_NumAwareLastYr_18_24_MSM_H;
            Results.YMSM_VLSAmongDiagLastYr_18_24_MSM_O = ...
                sum(Results.ann_NumberVLS(numYears,:)'.*...
                Params.ageIndicator(:,Params.age_2).*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator(:, Params.race_O))...
                / YMSM_NumAwareLastYr_18_24_MSM_O;
            
            Results.YMSM_VLSAmongDiagLastYr_25_34_MSM_B = ...
                sum(Results.ann_NumberVLS(numYears,:)'.*...
                Params.ageIndicator(:,Params.age_3).*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator(:, Params.race_B))...
                / YMSM_NumAwareLastYr_25_34_MSM_B;
            Results.YMSM_VLSAmongDiagLastYr_25_34_MSM_H = ...
                sum(Results.ann_NumberVLS(numYears,:)'.*...
                Params.ageIndicator(:,Params.age_3).*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator(:, Params.race_H))...
                / YMSM_NumAwareLastYr_25_34_MSM_H;
            Results.YMSM_VLSAmongDiagLastYr_25_34_MSM_O = ...
                sum(Results.ann_NumberVLS(numYears,:)'.*...
                Params.ageIndicator(:,Params.age_3).*Params.popSexIndicator(:,Params.popSex_MSM)...
                .*Params.raceIndicator(:, Params.race_O))...
                / YMSM_NumAwareLastYr_25_34_MSM_O;
            
             Results.YMSM_VLSAmongDiagLastYr_HET = ...
                sum(Results.ann_NumberVLS(numYears,:)'.*...
                Params.popIndicator(:,Params.pop_HET))...
                / YMSM_NumAwareLastYr_HET;
            Results.YMSM_VLSAmongDiagLastYr_MSM = ...
                sum(Results.ann_NumberVLS(numYears,:)'.*...
                Params.popIndicator(:,Params.pop_MSM))...
                / YMSM_NumAwareLastYr_MSM;
            Results.YMSM_VLSAmongDiagLastYr_IDU = ...
                sum(Results.ann_NumberVLS(numYears,:)'.*...
                Params.popIndicator(:,Params.pop_IDU))...
                / YMSM_NumAwareLastYr_IDU;
            
            for yearCount = 1:numYears
                     
                    Results.ann_CostTestAndNotify_LRMSM_Disc(yearCount,:) = Results.ann_TransCost_CostTestAndNotify_Disc(yearCount,:) .* ...
                        Params.popSexIndicator(:,Params.popSex_MSM)'...
                        .*Params.riskLevelIndicator(:,Params.risk_Main)';
                    Results.ann_CostTestAndNotify_HRMSM_Disc(yearCount,:) = Results.ann_TransCost_CostTestAndNotify_Disc(yearCount,:) .* ...
                        Params.popSexIndicator(:,Params.popSex_MSM)'...
                        .*Params.riskLevelIndicator(:,Params.risk_Casual)';
                    
                    Results.ann_TotalNewDiagnoses_LRMSM(yearCount,:) = Results.ann_TotalNewDiagnoses(yearCount,:) .* ...
                        Params.popIndicator(:,Params.pop_MSM)' .* ...
                        Params.riskLevelIndicator(:,Params.risk_Main)';
                    Results.ann_TotalNewDiagnoses_HRMSM(yearCount,:) = Results.ann_TotalNewDiagnoses(yearCount,:) .* ...
                        Params.popIndicator(:,Params.pop_MSM)' .* ...
                        Params.riskLevelIndicator(:,Params.risk_Casual)';
                 
            end
        
        end
    end
end

%% Functions

% Conversion function
function OutProb = RateToProb(InAnnualRate)
        OutProb = 1 - exp(-InAnnualRate);
end
end

% define functions for use in section 7 (Chris G did this, I didn't
% know where else to put, the proble is our for ... 1:1 loops)

% % percent aware for MSM by age and race
% % create function then call 6 times to create vars
% function out = percentAware_pop_race_age(pop,race,age)
% 
%     % make these into pop_MSM
%     pop = ['pop_' pop];
%     race = ['race_' race];
%     age = ['age_' age];
% 
%     out = ...
%     sum(Results.ann_NumberAware(Params.year_2016,:)' ...
%         .* Params.popIndicator(:,Params.(pop)) ...
%         .* Params.raceIndicator(:,Params.(race)) ...
%         .* Params.ageIndicator(:,Params.(age))) ...
%     / ...
%     sum(Results.ann_HIVPrevalence(Params.year_2016,:)' ...
%         .* Params.popIndicator(:,Params.(pop)) ...
%         .* Params.raceIndicator(:,Params.(race)) ...
%         .* Params.ageIndicator(:,Params.(age)));
% 
% end
