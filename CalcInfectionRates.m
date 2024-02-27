function InfRate = CalcInfectionRates(Params, Compartments, TTProg, lambdaType, PrEPType)
%% Purpose: Calculate the annual per person infection rate
    
% Output:
    % InfRate.lambdaVonly [273x1]
    % InfRate.lambdaAonly [273x1]
    % InfRate.lambdaVsomeAI [273x1]
    % InfRate.lambdaAsomeVI [273x1]
    % InfRate.lambdaN [273x1]
    % InfRate.totalInfRate [273x1]
    % InfRate.RelativeRates [273x5]
    
%% 1. Preallocate local variables

    % Pre-allocate key variables
         NumInfectedPeopleInStage = zeros(Params.numComparts,Params.numStrats);
         PctOfInfectedPopulationInStage(Params.numComparts,Params.numStrats)=0;
         NumSexPartnersByCohort(Params.numStrats,Params.numStrats)=0;
         pctInfSourceBySubpop(Params.numStrats,Params.numStrats)=0;
    
    % Assign key variables to local variables
   
        NumPartnersSex = Params.partnersSex;

        % Assign betas
        if lambdaType == Params.ind_PrEP
            
            % PrEP
            oneMinusBetaVonly = TTProg.oneMinusBetaVonly_PrEP; % Defined in CalcTransRates.m
            oneMinusBetaAonly = TTProg.oneMinusBetaAonly_PrEP; % Defined in CalcTransRates.m
            oneMinusBetaN = TTProg.oneMinusBetaN_PrEP; % Defined in CalcTransRates.m
            oneMinusBetaAsomeVI = TTProg.oneMinusBetaAsomeVI_PrEP; % Defined in CalcTransRates.m 
            oneMinusBetaVsomeAI = TTProg.oneMinusBetaVsomeAI_PrEP; % Defined in CalcTransRates.m
           %disp('Calc PrEP lambdas')
        else
           
            % Non-PrEP
            oneMinusBetaN = TTProg.oneMinusBetaN_NoPrEP; % Defined in CalcTransRates.m
            oneMinusBetaVonly = TTProg.oneMinusBetaVonly_NoPrEP; % Defined in CalcTransRates.m
            oneMinusBetaAonly = TTProg.oneMinusBetaAonly_NoPrEP; % Defined in CalcTransRates.m
            oneMinusBetaAsomeVI = TTProg.oneMinusBetaAsomeVI_NoPrEP; % Defined in CalcTransRates.m
            oneMinusBetaVsomeAI = TTProg.oneMinusBetaVsomeAI_NoPrEP; % Defined in CalcTransRates.m
            
        end
        
        
        % Among people who have AI in M-F partnerships, percent of
        % partnerships that include AI
        pctMFPartnersA_IfSomeAInMF = TTProg.pctMFPartnersA_IfSomeAInMF; 
       
        % Pct of people who have A in M-F partnerships
            pctPeopleWithA_InMF = TTProg.pctPeopleWithA_InMF;
            pctPeopleVonly_InMF = 1 - pctPeopleWithA_InMF; % This can be calculated as a complement b/c it's 273x1 not 273x273

%% 2. Calculate percent of population infected per stage (compartment)

% Number of people infected in each compartment/subpop
    % [30x273]
        % Note: Params.A1 is 0 because they're uninfected
        % Absorbing states aren't included
    
    NumInfectedPeopleInStage(Params.HIVComparts,:)=Compartments(Params.HIVComparts,:);

% Total population by cohort
    % [1x273]
        % For each cohort, sum across all active compartments
    TotalPopulationByCohort = sum(Compartments(Params.NonAbsorbingComparts,:));

% Percent of total HIV+ population that is a particular cohort in a
% particular compartment
    %[30x273]
        % E.g., if 10 total people are HIV+, and 3 are Black and in compart B2,
        % 2 are Hispanic and in compart B2, and 5 are Other in B3 and the values for 
            % B2 and Black = 0.3; B2 and Hispanic = 0.2; B3 and Other = 0.5
       % Note: Sum across cohort at model start should equal HIV prevalence input
            
    % Formula: NumInfectedPeopleInStage/TotalPop
     
    for n = Params.NonAbsorbingComparts
        PctOfInfectedPopulationInStage(n,:) = NumInfectedPeopleInStage(n,:) ./ max(TotalPopulationByCohort,1);
    end
    
    InfRate.PctInfPerStage = PctOfInfectedPopulationInStage;
    
    %Removed lines here 4/3/2020 to speed up model- Bates
%% 3. Calculate the Number of Infected Partners by Partnership Type

    % 3.i. Sexual partners
    for sectionInfectedPartners = 1:1

    % Total number of sexual partners
        % [PartnersSex] * [Sexual Mixing]
            for a = 1:Params.numStrats
                NumSexPartnersByCohort(a,:) = NumPartnersSex(a) * Params.behav_SexualMixing(a,:);
            end

    % Diagnosed compartments vector
    diagnosed_vector = ones(Params.numComparts,1);
    for i = union(Params.UnawareComparts, Params.UninfectedComparts)
        diagnosed_vector(i,1) = 0;
    end
    
    % Undiagnosed compartments vector
    undiag_vector = zeros(Params.numComparts,1);
    for i = union(Params.UnawareComparts, Params.UninfectedComparts)
        undiag_vector(i,1) = 1;
    end
    
    % make compart x strats x strats (30x273x273) version of diag and undiag vector
    repm_diagnosed_vector = repmat(diagnosed_vector,[1,Params.numStrats,Params.numStrats]);
    repm_undiag_vector = repmat(undiag_vector,[1,Params.numStrats,Params.numStrats]);

    % Diagnosed reduction values
    HETDiagReduc_value = 1- Params.pctReductDiagPartners(Params.pop_HET);
    MSMDiagReduc_value = 1 - Params.pctReductDiagPartners(Params.pop_MSM);
    IDUDiagReduc_value = 1 - Params.pctReductDiagPartners(Params.pop_IDU);
       
    % Create HET reduction arrary
    HET_reduc_array = diagnosed_vector * Params.popIndicator(:,Params.pop_HET)';
    HET_reduc_array = HET_reduc_array * HETDiagReduc_value;
    
    % Create MSM reduction arrary
    MSM_reduc_array = diagnosed_vector * Params.popIndicator(:,Params.pop_MSM)';
    MSM_reduc_array = MSM_reduc_array * MSMDiagReduc_value;
    
    % Create IDU reduction arrary
    IDU_reduc_array = diagnosed_vector * Params.popIndicator(:,Params.pop_IDU)';
    IDU_reduc_array = IDU_reduc_array * IDUDiagReduc_value;
    
    % Combine
    ALL_reduc_array = HET_reduc_array + MSM_reduc_array + IDU_reduc_array;
    ALL_reduc_array = repmat(ALL_reduc_array,[1,1,Params.numStrats]);
    ALL_reduc_array = ALL_reduc_array + repm_undiag_vector; %add back in 1s from undiagnosed comparts
      
    % Male-female infected partners
        
        % Number of partnerships
            % Apply overall partners
                NumSexPartnersByCohort_MaleFemale = NumSexPartnersByCohort;
            % Remove MSM-MSM partnerships
                NumSexPartnersByCohort_MaleFemale(Params.idx_MSMandMSM) = 0; 
            
        % Number of infected partnerships- Removed sections like this
        % 4/3/2020 to improve runtime- Bates
            %repm_numSexPartnerMaleFemale = repmat(NumSexPartnersByCohort_MaleFemale,[1,1,Params.numComparts]);     
            %repm_numSexPartnerMaleFemale = permute(repm_numSexPartnerMaleFemale,[3,1,2]);     
            %numInfectedSexPartersMaleFemale = repm_numSexPartnerMaleFemale .* repm_PctInfPerStage .* ALL_reduc_array;
            
            
        % Number of infected partnerships
             numInfectedSexPartersMaleFemale = zeros(Params.numComparts, Params.numStrats, Params.numStrats);
            for i = 1: Params.numComparts
                numInfectedSexPartersMaleFemale(i,:,:) = PctOfInfectedPopulationInStage(i,:) .* NumSexPartnersByCohort_MaleFemale;
            end 
             numInfectedSexPartersMaleFemale = numInfectedSexPartersMaleFemale .* ALL_reduc_array;
            
     % Male-male infected partners
        
        % Number of partnerships
            % Apply overall partners
                NumSexPartnersByCohort_MaleMale = NumSexPartnersByCohort;
            % Remove non-MSM-MSM partnerships
                NumSexPartnersByCohort_MaleMale(Params.idx_allMaleAndFemaleMixing) = 0;

        % Number of infected partnerships
        numInfectedSexPartersMaleMale = zeros(Params.numComparts, Params.numStrats, Params.numStrats);
            for i = 1: Params.numComparts
                numInfectedSexPartersMaleMale(i,:,:) = PctOfInfectedPopulationInStage(i,:) .* NumSexPartnersByCohort_MaleMale;
            end
        numInfectedSexPartersMaleMale = numInfectedSexPartersMaleMale .* ALL_reduc_array;   

    end

    % 3.ii. Needle partners
    for sectionNeedlePartner = 1:1
        
    % read in reduction in needle sharing and subtract from one to get reduction percentage    
    NeedleDiagReduc_value = 1 - Params.pctReductDiagNeedlePartners; 
    
    % Create vector for Needle partner reduction due to diagnosis   
    NeedleDiagReduc_vector = (diagnosed_vector * NeedleDiagReduc_value) + undiag_vector;

    % Expand vector to 30x273x273 so it can be multiplied against the
    repm_NeedleDiagReduc_vector = repmat(NeedleDiagReduc_vector,[1,Params.numStrats,Params.numStrats]);       
        
    % Number of partnerships
            % [PartnersSex] * [Sexual Mixing]
            % [273x273]
        NumNeedlePartnersByCohort = Params.behav_partnersNeedle .* Params.behav_NeedleMixing;

    % Number of infected partners
            numInfectedNeedlePartners = zeros(Params.numComparts, Params.numStrats, Params.numStrats);
            for i = 1: Params.numComparts
                numInfectedNeedlePartners(i,:,:) = PctOfInfectedPopulationInStage(i,:) .* NumNeedlePartnersByCohort;
            end
            numInfectedNeedlePartners = numInfectedNeedlePartners .* repm_NeedleDiagReduc_vector;
            
    end
      
%% 4. Calculate the Lambdas (annual probability of an individual being infected)

    % 4.i. LambdaAonly_MM: Male-male partnership, AI transmission risk
    for sectionLambdaA = 1:1

    % Formula:
        % [lambdaAonly] = 1 - (1-betaAonly) ^ [# male-male infected partners) 
        
        % Calculate [30x273x273] lambda
        lambdaAonly_MM_store = oneMinusBetaAonly .^ numInfectedSexPartersMaleMale;
        InfRate.lambdaAonly_store = lambdaAonly_MM_store;
        
        % Multiply over all compartments [273x273]
        lambdaAonly_MM_store2 = prod(lambdaAonly_MM_store,1);
        InfRate.lambdaAonly_store2 = lambdaAonly_MM_store2;
        
        % Multiply over all partners [273x1]
            % Probability of not getting infected by all partners
        lambdaAonly_MM_store3 = permute(lambdaAonly_MM_store2,[2 3 1]);
        lambdaAonly_MM_store4 = prod(lambdaAonly_MM_store3,2);
        
        % Overall probability of getting infected 
            % [273x1]
        InfProb_lambdaAonly_MM = 1 - lambdaAonly_MM_store4;
       
    end
    
    % 4.ii. LambdaVonly_MF: Male-female partnership, VI transmission risk
    for sectionlambdaV = 1:1
    
    % Formula: 
        % [lambdaVonly] = [Percent people only V in M-F partnerships] * (1 - (1-betaVonly) ^ [# male-female infected partners) ) 
        
        %  Calculate [30x273x273] lambda
        lambdaVonly_MF_store1 = oneMinusBetaVonly .^ numInfectedSexPartersMaleFemale;         
            InfRate.lambdaVonly_store = lambdaVonly_MF_store1;
        
        % Multiply over all compartments [273x273]
        lambdaVonly_MF_store2 = prod(lambdaVonly_MF_store1, 1);
            InfRate.lambdaAonly_store2 = lambdaVonly_MF_store2;
        
        % Multiply over all partners [273x1]
        lambdaVonly_MF_store3 = permute(lambdaVonly_MF_store2,[2 3 1]);
        lambdaVonly_MF_store4 = prod(lambdaVonly_MF_store3,2);
        
        % Overall prob of getting infected [273x1]
        InfProb_lambdaVonly_MF = pctPeopleVonly_InMF .* (1 - lambdaVonly_MF_store4);
    
    end
        
    % 4.iii. LambdaVsomeA_MF and LambdaAsomeV_MF: Male-female Partnerships:
    % VI and AI transmission risks
    for sectionLamdaAsomeVI = 1:1
        

    % Infection from A: LambdaAsomeVI
    
        % Formula:
            % [% of people with A in M-F partnerships] * (1 - (1- [betaAsomeV] ^ ([# infected male-female partners] *[If person has A in M-F, % of partners who do A]) 

            %numinf(Params.numStrats,Params.numStrats)=0;
            %numinf(:,:) = numInfectedSexPartersMaleFemale(6,:,:);
            %numinfadj = numinf .* pctMFPartnersA_IfSomeAInMF';
            
            % Calculate [30x273x273] lambda
            lambdaAsomeVI_MF_store1 = oneMinusBetaAsomeVI .^ (numInfectedSexPartersMaleFemale .* pctMFPartnersA_IfSomeAInMF);

            % Multiply over all compartments [273x273]
            lambdaAsomeVI_MF_store2 = prod(lambdaAsomeVI_MF_store1,1);

            % Multiply over all partners [273x1]
                % Probability of not getting infected by all partners
            lambdaAsomeVI_MF_store3 = permute(lambdaAsomeVI_MF_store2,[2 3 1]);
            lambdaAsomeVI_MF_store4 = prod(lambdaAsomeVI_MF_store3,2);
            
            % Overall probability of getting infected 
                % [273x1]
            InfProb_lambdaAsomeVI_MF = pctPeopleWithA_InMF .* (1 - lambdaAsomeVI_MF_store4);


    % Infection from V
        % LambdaVsomeAI
  
        % Formula
        % [% people who have A in m-f relationships] * [ 
        %   (1 - (1-betaVonly) ^ ([# infected M-F partners] * [% partnerships that don't include A, if person has A in M-F partnerships]) )
        % + (1 - (1-betaVsomeA) ^ ([# infected M-F partners] * [% partnerships that include A, if person has A in M-F partnerships]) ) ]
    

        % Calculate [30x273x273] lambda
            % Male-female partnerships that don't include A, given person has some A in male-female partnerships
            lambdaVPersonBoth_MF_PartnerV__store1 = ...
                oneMinusBetaVonly .^ (numInfectedSexPartersMaleFemale .* (1 - pctMFPartnersA_IfSomeAInMF));
            % Male-female partnerships that include A
            lamVPersonBoth_MF_PartnerBoth_store1 = ...
                oneMinusBetaVsomeAI .^ (numInfectedSexPartersMaleFemale .* pctMFPartnersA_IfSomeAInMF);
                
        %  Multiply across compartments and cohorts until [273x1]
        	% Partnership only V
            lambdaVPersonBoth_MF_PartnerV_store2 = prod(lambdaVPersonBoth_MF_PartnerV__store1,1);
            lambdaVPersonBoth_MF_PartnerV_store3 = permute(lambdaVPersonBoth_MF_PartnerV_store2,[2 3 1]);
            lambdaVPersonBoth_MF_PartnerV_store4 = prod(lambdaVPersonBoth_MF_PartnerV_store3,2);       
           % Partnership includes A
            lamVPersonBoth_MF_PartnerBoth_store2 = prod(lamVPersonBoth_MF_PartnerBoth_store1,1);
            lamVPersonBoth_MF_PartnerBoth_store3 = permute(lamVPersonBoth_MF_PartnerBoth_store2,[2 3 1]);
            lamVPersonBoth_MF_PartnerBoth_store4 = prod(lamVPersonBoth_MF_PartnerBoth_store3,2);            
                
                
        % Calculate actual lambda
            InfProb_lambdaVsomeAI_MF = pctPeopleWithA_InMF .* ...
                 ( (1 - lambdaVPersonBoth_MF_PartnerV_store4 ) ...
                +  (1 - lamVPersonBoth_MF_PartnerBoth_store4));
                
        
                
    end

    % 4.iv. LambdaN: Needle partnerships, N transmission risk
    for sectionLambdaN = 1:1

    % Infection rate from needle transmission
    
        % Calculate [30x273x273] lambda
            lambdaN_store = oneMinusBetaN .^ numInfectedNeedlePartners;

        %  Multiply across compartments and cohorts until [273x1]
            lambdaN_store2 = prod(lambdaN_store,1);
            lambdaN_store3 = permute(lambdaN_store2,[2 3 1]);
            lambdaN_store4 = prod(lambdaN_store3,2);
            
        % Calculate actual lambda
            InfProb_lambdaN = 1 - lambdaN_store4;

    end

%% 5. Calculate the total infection rate and percentage of infections from each other subpop

    % Convert from probabilities to rates
    
    InfRate.lambdaVonly = -log(1-InfProb_lambdaVonly_MF);
    InfRate.lambdaAonly = -log(1-InfProb_lambdaAonly_MM);
    InfRate.lambdaVsomeAI = -log(1-InfProb_lambdaVsomeAI_MF);
    InfRate.lambdaAsomeVI = -log(1-InfProb_lambdaAsomeVI_MF);
    InfRate.lambdaN = -log(1-InfProb_lambdaN);

    % Reduction due to SSP
    InfRate.lambdaN = InfRate.lambdaN .* (1 - (TTProg.PctActivePWIDServedbySEP * Params.SEPTransmissionReduct));
    
    % Reduction due to PrEP % JCPrEPUpdate: modified to calculate reduced lambdas for oral and injectable PrEP
    if lambdaType == Params.ind_PrEP
        
        if PrEPType == Params.ind_OralPrEP_High
            
            InfRate.lambdaVonly = InfRate.lambdaVonly .* (1 - Params.inf_prepIncidenceReduct_Oral_HighAdherence);
            InfRate.lambdaAonly = InfRate.lambdaAonly .* (1 - Params.inf_prepIncidenceReduct_Oral_HighAdherence);
            InfRate.lambdaVsomeAI = InfRate.lambdaVsomeAI .* (1 - Params.inf_prepIncidenceReduct_Oral_HighAdherence);
            InfRate.lambdaAsomeVI = InfRate.lambdaAsomeVI .* (1 - Params.inf_prepIncidenceReduct_Oral_HighAdherence);
            InfRate.lambdaN = InfRate.lambdaN .* (1 - Params.inf_prepIncidenceReduct_Oral_HighAdherence);
            
        elseif PrEPType == Params.ind_OralPrEP_Low
            InfRate.lambdaVonly = InfRate.lambdaVonly .* (1 - Params.inf_prepIncidenceReduct_Oral_LowAdherence);
            InfRate.lambdaAonly = InfRate.lambdaAonly .* (1 - Params.inf_prepIncidenceReduct_Oral_LowAdherence);
            InfRate.lambdaVsomeAI = InfRate.lambdaVsomeAI .* (1 - Params.inf_prepIncidenceReduct_Oral_LowAdherence);
            InfRate.lambdaAsomeVI = InfRate.lambdaAsomeVI .* (1 - Params.inf_prepIncidenceReduct_Oral_LowAdherence);
            InfRate.lambdaN = InfRate.lambdaN .* (1 - Params.inf_prepIncidenceReduct_Oral_LowAdherence);
            
        elseif PrEPType == Params.ind_InjectPrEP_High    
            
            InfRate.lambdaVonly_Inject = InfRate.lambdaVonly .* (1 - Params.inf_prepIncidenceReduct_Inject_HighAdherence);
            InfRate.lambdaAonly_Inject = InfRate.lambdaAonly .* (1 - Params.inf_prepIncidenceReduct_Inject_HighAdherence);
            InfRate.lambdaVsomeAI_Inject = InfRate.lambdaVsomeAI .* (1 - Params.inf_prepIncidenceReduct_Inject_HighAdherence);
            InfRate.lambdaAsomeVI_Inject = InfRate.lambdaAsomeVI .* (1 - Params.inf_prepIncidenceReduct_Inject_HighAdherence);
            InfRate.lambdaN_Inject = InfRate.lambdaN .* (1 - Params.inf_prepIncidenceReduct_Inject_HighAdherence);
            
        else   
            
            InfRate.lambdaVonly_Inject = InfRate.lambdaVonly .* (1 - Params.inf_prepIncidenceReduct_Inject_LowAdherence);
            InfRate.lambdaAonly_Inject = InfRate.lambdaAonly .* (1 - Params.inf_prepIncidenceReduct_Inject_LowAdherence);
            InfRate.lambdaVsomeAI_Inject = InfRate.lambdaVsomeAI .* (1 - Params.inf_prepIncidenceReduct_Inject_LowAdherence);
            InfRate.lambdaAsomeVI_Inject = InfRate.lambdaAsomeVI .* (1 - Params.inf_prepIncidenceReduct_Inject_LowAdherence);
            InfRate.lambdaN_Inject = InfRate.lambdaN .* (1 - Params.inf_prepIncidenceReduct_Inject_LowAdherence);
            
        end
    end
    
    
    % Total infection rate 
        % sum across 5 forces of infection
    InfRate.totalInfRate = max(InfRate.lambdaVsomeAI + InfRate.lambdaVonly + ...
        InfRate.lambdaAonly + InfRate.lambdaAsomeVI + InfRate.lambdaN,0.0000000000000000000001);
    
    % [Proportion of total infection rate that is due to force of infection]
    InfRate.RelativeInfRates(:,Params.force_Vonly) = InfRate.lambdaVonly./InfRate.totalInfRate;
    InfRate.RelativeInfRates(:,Params.force_Aonly) = InfRate.lambdaAonly./InfRate.totalInfRate;
    InfRate.RelativeInfRates(:,Params.force_VsomeA) = InfRate.lambdaVsomeAI./InfRate.totalInfRate;
    InfRate.RelativeInfRates(:,Params.force_AsomeV) = InfRate.lambdaAsomeVI./InfRate.totalInfRate;
    InfRate.RelativeInfRates(:,Params.force_N) = InfRate.lambdaN./InfRate.totalInfRate;
    
    % Calculate percentage of infections from each subpopulation (273 x
    % 273) if selected to collect infections by source (on 'Model Settings'
    % sheet)
    if Params.CollectInfbySource == 1
        pctInfSourceBySubpop_Vonly = CalcPctInfSourceBySubpop(lambdaVonly_MF_store3);
        pctInfSourceBySubpop_Aonly = CalcPctInfSourceBySubpop(lambdaAonly_MM_store3);
        pctInfSourceBySubpop_AsomeV = CalcPctInfSourceBySubpop(lambdaAsomeVI_MF_store3);
        pctInfSourceBySubpop_VsomeA = CalcPctInfSourceBySubpop(1-((1 - lambdaVPersonBoth_MF_PartnerV_store3) + (1 - lamVPersonBoth_MF_PartnerBoth_store3)));
        pctInfSourceBySubpop_N = CalcPctInfSourceBySubpop(lambdaN_store3);
        for subpopk = 1:Params.numStrats
            pctInfSourceBySubpop(subpopk,:) = ...
                InfRate.RelativeInfRates(subpopk,Params.force_Vonly) .* ...
                    pctInfSourceBySubpop_Vonly(subpopk,:) + ...
                InfRate.RelativeInfRates(subpopk,Params.force_Aonly) .* ...
                    pctInfSourceBySubpop_Aonly(subpopk,:) + ...
                InfRate.RelativeInfRates(subpopk,Params.force_AsomeV) .* ...
                    pctInfSourceBySubpop_AsomeV(subpopk,:) + ...
                InfRate.RelativeInfRates(subpopk,Params.force_VsomeA) .* ...
                    pctInfSourceBySubpop_VsomeA(subpopk,:) + ...
                InfRate.RelativeInfRates(subpopk,Params.force_N) .* ...
                    pctInfSourceBySubpop_N(subpopk,:);
        end
        InfRate.pctInfSourceBySubpop = pctInfSourceBySubpop;
    end
    
    function PctInfSourceBySubpop = CalcPctInfSourceBySubpop(store3variable)
        
        % Percentage of infections from each subpopulation (273 x 273)
        PctInfSourceBySubpop(Params.numStrats, Params.numStrats)=0;
        probInfBySubpop = 1 - store3variable;
        for subpopi = 1:Params.numStrats
            for subpopj = 1:Params.numStrats
                if sum(probInfBySubpop(subpopi,:),2) > 0
                    PctInfSourceBySubpop(subpopi, subpopj) = ...
                        probInfBySubpop(subpopi, subpopj) ./ sum(probInfBySubpop(subpopi,:),2);
                end
            end    
        end

    end
    
end