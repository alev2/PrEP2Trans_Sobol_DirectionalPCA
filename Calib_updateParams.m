function [Params, CalibParams] = Calib_updateParams(Params, CalibParams, i) 
%% Purpose: A function which updates the model's parameters to reflect the
% sampled values

% Called from: HIVEpiModel - when running LHS calibration
% Called from: ______ - when running calibration by optimization
    % Note: This m file will only be called when a calibration is run.

%% 1. Set 2nd Time Period Values to 1st Time Period Values if Not Varying in Calibration

    if Params.Calib2ndPeriod == 0
    
        CalibParams.tt_testRefCase_2to4.paramValue(i) = CalibParams.tt_testRefCase_1.paramValue(i);
        CalibParams.tt_relRiskPop_2to4_MSM.paramValue(i) = CalibParams.tt_relRiskPop_1_MSM.paramValue(i);
        CalibParams.tt_relRiskPop_2to4_IDU.paramValue(i) = CalibParams.tt_relRiskPop_1_IDU.paramValue(i);
        CalibParams.tt_relRiskRace_2to4_H.paramValue(i) = CalibParams.tt_relRiskRace_1_H.paramValue(i);
        CalibParams.tt_relRiskRace_2to4_O.paramValue(i) = CalibParams.tt_relRiskRace_1_O.paramValue(i);
        CalibParams.tt_relRiskHIVstage_2to4_B.paramValue(i) = CalibParams.tt_relRiskHIVstage_1_B.paramValue(i);
        CalibParams.tt_relRiskHIVstage_2to4_D.paramValue(i) = CalibParams.tt_relRiskHIVstage_1_D.paramValue(i);
        CalibParams.tt_relRiskHIVstage_2to4_E.paramValue(i) = CalibParams.tt_relRiskHIVstage_1_E.paramValue(i);
        CalibParams.tt_relRiskHIVstage_2to4_F.paramValue(i) = CalibParams.tt_relRiskHIVstage_1_F.paramValue(i);
        CalibParams.tt_relRiskAge_2to4_1824.paramValue(i) = CalibParams.tt_relRiskAge_1_1824.paramValue(i);
        CalibParams.tt_relRiskAge_2to4_2534.paramValue(i) = CalibParams.tt_relRiskAge_1_2534.paramValue(i);
        CalibParams.tt_relRiskAge_2to4_3544.paramValue(i) = CalibParams.tt_relRiskAge_1_3544.paramValue(i);
        CalibParams.tt_relRiskAge_2to4_4554.paramValue(i) = CalibParams.tt_relRiskAge_1_4554.paramValue(i);
        CalibParams.tt_relRiskAge_2to4_5564.paramValue(i) = CalibParams.tt_relRiskAge_1_5564.paramValue(i);
        CalibParams.tt_relRiskAge_2to4_65.paramValue(i) = CalibParams.tt_relRiskAge_1_65.paramValue(i);
        CalibParams.tt_linkage_r_2to4_B.paramValue(i) = CalibParams.tt_linkage_r_1_B.paramValue(i);
        CalibParams.tt_linkage_r_2to4_H.paramValue(i) = CalibParams.tt_linkage_r_1_H.paramValue(i);
        CalibParams.tt_linkage_r_2to4_O.paramValue(i) = CalibParams.tt_linkage_r_1_O.paramValue(i);
        CalibParams.tt_RelRiskLTC_2to4_E.paramValue(i) = CalibParams.tt_RelRiskLTC_1_E.paramValue(i);
        CalibParams.tt_RelRiskLTC_2to4_F.paramValue(i) = CalibParams.tt_RelRiskLTC_1_F.paramValue(i);
        CalibParams.tt_dropOutProb_CareToAware_2to4_B.paramValue(i) = CalibParams.tt_dropOutProb_CareToAware_1_B.paramValue(i);
        CalibParams.tt_dropOutProb_CareToAware_2to4_H.paramValue(i) = CalibParams.tt_dropOutProb_CareToAware_1_H.paramValue(i);
        CalibParams.tt_dropOutProb_CareToAware_2to4_O.paramValue(i) = CalibParams.tt_dropOutProb_CareToAware_1_O.paramValue(i);
        CalibParams.tt_dropOutProb_ANVToCare_2to4_B.paramValue(i) = CalibParams.tt_dropOutProb_ANVToCare_1_B.paramValue(i);
        CalibParams.tt_dropOutProb_ANVToCare_2to4_H.paramValue(i) = CalibParams.tt_dropOutProb_ANVToCare_1_H.paramValue(i);
        CalibParams.tt_dropOutProb_ANVToCare_2to4_O.paramValue(i) = CalibParams.tt_dropOutProb_ANVToCare_1_O.paramValue(i);
        CalibParams.tt_dropOutProb_VLSToANV_2to4_B.paramValue(i) = CalibParams.tt_dropOutProb_VLSToANV_1_B.paramValue(i);
        CalibParams.tt_dropOutProb_VLSToANV_2to4_H.paramValue(i) = CalibParams.tt_dropOutProb_VLSToANV_1_H.paramValue(i);
        CalibParams.tt_dropOutProb_VLSToANV_2to4_O.paramValue(i) = CalibParams.tt_dropOutProb_VLSToANV_1_O.paramValue(i);
        CalibParams.tt_relRiskPop_VLSToANV_2to4_MSM.paramValue(i) = CalibParams.tt_relRiskPop_VLSToANV_1_MSM.paramValue(i);
        CalibParams.tt_relRiskPop_VLSToANV_2to4_PWID.paramValue(i) = CalibParams.tt_relRiskPop_VLSToANV_1_PWID.paramValue(i);
        CalibParams.tt_relRiskAge_VLSToANV_2to4_1824.paramValue(i) = CalibParams.tt_relRiskAge_VLSToANV_1_1824.paramValue(i);
        CalibParams.tt_relRiskAge_VLSToANV_2to4_2534.paramValue(i) = CalibParams.tt_relRiskAge_VLSToANV_1_2534.paramValue(i);
        CalibParams.tt_relRiskAge_VLSToANV_2to4_3544.paramValue(i) = CalibParams.tt_relRiskAge_VLSToANV_1_3544.paramValue(i);
        CalibParams.tt_relRiskAge_VLSToANV_2to4_4554.paramValue(i) = CalibParams.tt_relRiskAge_VLSToANV_1_4554.paramValue(i);
        CalibParams.tt_relRiskAge_VLSToANV_2to4_5564.paramValue(i) = CalibParams.tt_relRiskAge_VLSToANV_1_5564.paramValue(i);
        CalibParams.tt_relRiskAge_VLSToANV_2to4_65.paramValue(i) = CalibParams.tt_relRiskAge_VLSToANV_1_65.paramValue(i);
        CalibParams.tt_RelRiskARTInit_r_2to4_B.paramValue(i) = CalibParams.tt_RelRiskARTInit_r_1_B.paramValue(i);
        CalibParams.tt_RelRiskARTInit_r_2to4_H.paramValue(i) = CalibParams.tt_RelRiskARTInit_r_1_H.paramValue(i);
        CalibParams.tt_RelRiskARTInit_r_2to4_O.paramValue(i) = CalibParams.tt_RelRiskARTInit_r_1_O.paramValue(i);
        CalibParams.tt_BecomeVLSfromANV_r_2to4_B.paramValue(i) = CalibParams.tt_BecomeVLSfromANV_r_1_B.paramValue(i);
        CalibParams.tt_BecomeVLSfromANV_r_2to4_H.paramValue(i) = CalibParams.tt_BecomeVLSfromANV_r_1_H.paramValue(i);
        CalibParams.tt_BecomeVLSfromANV_r_2to4_O.paramValue(i) = CalibParams.tt_BecomeVLSfromANV_r_1_O.paramValue(i);
        CalibParams.tt_relRiskPop_ANVToVLS_2to4_MSM.paramValue(i) = CalibParams.tt_relRiskPop_ANVToVLS_1_MSM.paramValue(i);
        CalibParams.tt_relRiskPop_ANVToVLS_2to4_PWID.paramValue(i) = CalibParams.tt_relRiskPop_ANVToVLS_1_PWID.paramValue(i);
        CalibParams.tt_relRiskAge_ANVToVLS_2to4_1824.paramValue(i) = CalibParams.tt_relRiskAge_ANVToVLS_1_1824.paramValue(i);
        CalibParams.tt_relRiskAge_ANVToVLS_2to4_2534.paramValue(i) = CalibParams.tt_relRiskAge_ANVToVLS_1_2534.paramValue(i);
        CalibParams.tt_relRiskAge_ANVToVLS_2to4_3544.paramValue(i) = CalibParams.tt_relRiskAge_ANVToVLS_1_3544.paramValue(i);
        CalibParams.tt_relRiskAge_ANVToVLS_2to4_4554.paramValue(i) = CalibParams.tt_relRiskAge_ANVToVLS_1_4554.paramValue(i);
        CalibParams.tt_relRiskAge_ANVToVLS_2to4_5564.paramValue(i) = CalibParams.tt_relRiskAge_ANVToVLS_1_5564.paramValue(i);
        CalibParams.tt_relRiskAge_ANVToVLS_2to4_65.paramValue(i) = CalibParams.tt_relRiskAge_ANVToVLS_1_65.paramValue(i);
        CalibParams.pctPartnershipsVandA_2to4_13.paramValue(i) = CalibParams.pctPartnershipsVandA_1_13.paramValue(i);
        CalibParams.pctPartnershipsVandA_2to4_18.paramValue(i) = CalibParams.pctPartnershipsVandA_1_18.paramValue(i);
        CalibParams.pctPartnershipsVandA_2to4_25.paramValue(i) = CalibParams.pctPartnershipsVandA_1_25.paramValue(i);
        CalibParams.pctPartnershipsVandA_2to4_35.paramValue(i) = CalibParams.pctPartnershipsVandA_1_35.paramValue(i);
        CalibParams.pctPartnershipsVandA_2to4_45.paramValue(i) = CalibParams.pctPartnershipsVandA_1_45.paramValue(i);
        CalibParams.pctPartnershipsVandA_2to4_55.paramValue(i) = CalibParams.pctPartnershipsVandA_1_55.paramValue(i);
        CalibParams.pctPartnershipsVandA_2to4_65.paramValue(i) = CalibParams.pctPartnershipsVandA_1_65.paramValue(i);
                            
    end    
    
%% 2. Update Test and Treat Continuum Parameters
 
    % 2a. Testing
    for testSection = 1:1   
            
        store_tt_testRefCase_1 = CalibParams.tt_testRefCase_1.paramValue(i);
        store_tt_testRefCase_2to4 = CalibParams.tt_testRefCase_2to4.paramValue(i);
        
%Update the relative risks of being tested
            
     % Period 1 
        
        % Transmission group
            
            % HET (reference case) - not calibrated
            store_tt_relRiskPop_1(Params.pop_HET,1) = Params.tt_relRiskPop_1(Params.pop_HET);
            
            % MSM
            store_tt_relRiskPop_1(Params.pop_MSM,1) = CalibParams.tt_relRiskPop_1_MSM.paramValue(i);
            
            % IDU
            store_tt_relRiskPop_1(Params.pop_IDU,1) = CalibParams.tt_relRiskPop_1_IDU.paramValue(i);
            
       % Race/ethnicity
      
            % Black (reference case) - not calibrated
            store_tt_relRiskRace_1(Params.race_B,1) = Params.tt_relRiskRace_1(Params.race_B);
            
            % Hispanic
            store_tt_relRiskRace_1(Params.race_H,1) = CalibParams.tt_relRiskRace_1_H.paramValue(i);
            
            % Other
            store_tt_relRiskRace_1(Params.race_O,1) = CalibParams.tt_relRiskRace_1_O.paramValue(i);
            
        % Disease stage
            
            % Acute
            store_tt_relRiskHIVstage_1(Params.stage_Acute,1) = CalibParams.tt_relRiskHIVstage_1_B.paramValue(i);
            
            % CD4>500 (reference case) - not calibrated
            store_tt_relRiskHIVstage_1(Params.stage_LatentA,1) = Params.tt_relRiskHIVstage_1(Params.stage_LatentA);
            
            % CD4 350-500
            store_tt_relRiskHIVstage_1(Params.stage_LatentB,1) = CalibParams.tt_relRiskHIVstage_1_D.paramValue(i);
            
            % CD4 200-350
            store_tt_relRiskHIVstage_1(Params.stage_Late,1) = CalibParams.tt_relRiskHIVstage_1_E.paramValue(i);
            
            % CD4 < 200
            store_tt_relRiskHIVstage_1(Params.stage_AIDS,1) = CalibParams.tt_relRiskHIVstage_1_F.paramValue(i);
            
            % Unaware (added 2018/04/16 by JC)
            store_tt_relRiskHIVstage_1(Params.numHIVstages+1,1) = Params.tt_relRiskHIVstage_1(Params.numHIVstages+1);
            
        % Age group
            
            % 13-17
            store_tt_relRiskAge_1(Params.age_1,1) = Params.tt_relRiskAge_1(Params.age_1);
            
            % 18-24
            store_tt_relRiskAge_1(Params.age_2,1) = CalibParams.tt_relRiskAge_1_1824.paramValue(i);
            
            % 25-34
            store_tt_relRiskAge_1(Params.age_3,1) = CalibParams.tt_relRiskAge_1_2534.paramValue(i);
            
            % 35-44
            store_tt_relRiskAge_1(Params.age_4,1) = CalibParams.tt_relRiskAge_1_3544.paramValue(i);
            
            % 45-54
            store_tt_relRiskAge_1(Params.age_5,1) = CalibParams.tt_relRiskAge_1_4554.paramValue(i);
            
            % 55-64
            store_tt_relRiskAge_1(Params.age_6,1) = CalibParams.tt_relRiskAge_1_5564.paramValue(i);
            
            % 65+
            store_tt_relRiskAge_1(Params.age_7,1) = CalibParams.tt_relRiskAge_1_65.paramValue(i);
            
    %Period 2
        
        % Transmission group
        
            % HET (reference case) - not calibrated
            store_tt_relRiskPop_2to4(Params.pop_HET,1) = Params.tt_relRiskPop_2to4(Params.pop_HET);
            
            % MSM
            store_tt_relRiskPop_2to4(Params.pop_MSM,1) = CalibParams.tt_relRiskPop_2to4_MSM.paramValue(i);
            
            % IDU
            store_tt_relRiskPop_2to4(Params.pop_IDU,1) = CalibParams.tt_relRiskPop_2to4_IDU.paramValue(i);
            
        % Race/ethnicity
      
            % Black (reference case) - not calibrated
            store_tt_relRiskRace_2to4(Params.race_B,1) = Params.tt_relRiskRace_2to4(Params.race_B);
            
            % Hispanic
            store_tt_relRiskRace_2to4(Params.race_H,1) = CalibParams.tt_relRiskRace_2to4_H.paramValue(i);
            
            % Other
            store_tt_relRiskRace_2to4(Params.race_O,1) = CalibParams.tt_relRiskRace_2to4_O.paramValue(i);
            
        % Disease stage
            
            % Acute
            store_tt_relRiskHIVstage_2to4(Params.stage_Acute,1) = CalibParams.tt_relRiskHIVstage_2to4_B.paramValue(i);

            % CD4>500 (reference case) - not calibrated
            store_tt_relRiskHIVstage_2to4(Params.stage_LatentA,1) = Params.tt_relRiskHIVstage_2to4(Params.stage_LatentA);
            
            % CD4 350-500
            store_tt_relRiskHIVstage_2to4(Params.stage_LatentB,1) = CalibParams.tt_relRiskHIVstage_2to4_D.paramValue(i);
            
            % CD4 200-350
            store_tt_relRiskHIVstage_2to4(Params.stage_Late,1) = CalibParams.tt_relRiskHIVstage_2to4_E.paramValue(i);
            
            % CD4 < 200
            store_tt_relRiskHIVstage_2to4(Params.stage_AIDS,1) = CalibParams.tt_relRiskHIVstage_2to4_F.paramValue(i);
            
            % Unaware (added 2018/04/16 by JC)
            store_tt_relRiskHIVstage_2to4(Params.numHIVstages+1,1) = Params.tt_relRiskHIVstage_2to4(Params.numHIVstages+1);
            
        % Age group
            
            % 13-17
            store_tt_relRiskAge_2to4(Params.age_1,1) = Params.tt_relRiskAge_2to4(Params.age_1);
            
            % 18-24
            store_tt_relRiskAge_2to4(Params.age_2,1) = CalibParams.tt_relRiskAge_2to4_1824.paramValue(i);
            
            % 25-34
            store_tt_relRiskAge_2to4(Params.age_3,1) = CalibParams.tt_relRiskAge_2to4_2534.paramValue(i);
            
            % 35-44
            store_tt_relRiskAge_2to4(Params.age_4,1) = CalibParams.tt_relRiskAge_2to4_3544.paramValue(i);
            
            % 45-54
            store_tt_relRiskAge_2to4(Params.age_5,1) = CalibParams.tt_relRiskAge_2to4_4554.paramValue(i);
            
            % 55-64
            store_tt_relRiskAge_2to4(Params.age_6,1) = CalibParams.tt_relRiskAge_2to4_5564.paramValue(i);
            
            % 65+
            store_tt_relRiskAge_2to4(Params.age_7,1) = CalibParams.tt_relRiskAge_2to4_65.paramValue(i);

%Calculate the new annual testing rates using the
    %calibrated relative risks
             
    % Period 1
             
        testingStore=zeros(Params.numStrats,1);
             
        for s = 1:Params.numHIVstages + 1
            for r = 1:Params.numRace
              for p = 1:Params.numPop
                for a = 1:Params.numAge

                        pTest = Params.raceIndicator(:,r).* ...
                            Params.popIndicator(:,p) .* Params.ageIndicator(:,a) * (store_tt_relRiskRace_1(r)...
                            * store_tt_relRiskPop_1(p) * store_tt_relRiskAge_1(a) * store_tt_relRiskHIVstage_1(s) ...
                            * store_tt_testRefCase_1);
                        
                        testingStore = pTest+testingStore;
                end
              end
            end
               switch s
                   case 1
                        Params.tt_testRateAcute_rp_1 = testingStore;
                        clear testingStore
                        testingStore=zeros(Params.numStrats,1);
                   case 2
                        Params.tt_testRateLatentA_rp_1 = testingStore;
                        clear testingStore 
                        testingStore=zeros(Params.numStrats,1);
                   case 3
                        Params.tt_testRateLatentB_rp_1 = testingStore;
                        clear testingStore
                        testingStore=zeros(Params.numStrats,1);
                   case 4
                        Params.tt_testRateLate_rp_1 = testingStore;
                        clear testingStore 
                        testingStore=zeros(Params.numStrats,1);
                   case 5
                        Params.tt_testRateAIDS_rp_1 = testingStore;
                        clear testingStore 
                        testingStore=zeros(Params.numStrats,1);
                   case 6
                        Params.tt_testRateUninfected_rp_1 = testingStore;
                        clear testingStore 
                        testingStore = zeros(Params.numStrats,1);     
               end

        end
        
      % Period 2
        for s = 1:Params.numHIVstages + 1
            for r = 1:Params.numRace
              for p = 1:Params.numPop
                for a = 1:Params.numAge

                        pTest = Params.raceIndicator(:,r).* ...
                            Params.popIndicator(:,p) .* Params.ageIndicator(:,a) * (store_tt_relRiskRace_2to4(r)...
                            * store_tt_relRiskPop_2to4(p) * store_tt_relRiskAge_2to4(a) * store_tt_relRiskHIVstage_2to4(s) ...
                            * store_tt_testRefCase_2to4);

                        testingStore = pTest+testingStore;
                end
              end 
            end

               switch s
                   case 1
                        Params.tt_testRateAcute_rp_2to4 = testingStore;
                        clear testingStore
                        testingStore=zeros(Params.numStrats,1);
                   case 2
                        Params.tt_testRateLatentA_rp_2to4 = testingStore;
                        clear testingStore 
                        testingStore=zeros(Params.numStrats,1);
                   case 3
                        Params.tt_testRateLatentB_rp_2to4 = testingStore;
                        clear testingStore 
                        testingStore=zeros(Params.numStrats,1);
                   case 4
                        Params.tt_testRateLate_rp_2to4 = testingStore;
                        clear testingStore
                        testingStore=zeros(Params.numStrats,1);
                   case 5
                        Params.tt_testRateAIDS_rp_2to4 = testingStore;
                        clear testingStore
                        testingStore=zeros(Params.numStrats,1);
                   case 6
                        Params.tt_testRateUninfected_rp_2to4 = testingStore;
                        clear testingStore
                        testingStore = zeros(Params.numStrats,1);
               end
        end             
    end

    % 2b. Linked to Care First (removed on 2018/01/25 by JC)
    for linkageFirstSection = 1:1 
    
        % Calculate the probability of getting linked to care immediately
    
        %1st period
        %tt_linkageFirst_r_1(Params.race_B,1)=CalibParams.tt_linkageFirst_r_1_B.paramValue(i);
        %tt_linkageFirst_r_1(Params.race_H,1)=CalibParams.tt_linkageFirst_r_1_H.paramValue(i);
        %tt_linkageFirst_r_1(Params.race_O,1)=CalibParams.tt_linkageFirst_r_1_O.paramValue(i);
        
        %2nd period
        %tt_linkageFirst_r_2to4(Params.race_B,1)=CalibParams.tt_linkageFirst_r_2to4_B.paramValue(i);
        %tt_linkageFirst_r_2to4(Params.race_H,1)=CalibParams.tt_linkageFirst_r_2to4_H.paramValue(i);
        %tt_linkageFirst_r_2to4(Params.race_O,1)=CalibParams.tt_linkageFirst_r_2to4_O.paramValue(i);
    
    % Calculate new linkageFirst parameter
        %Params.tt_linkageFirst_r_1 = Params.raceIndicator*tt_linkageFirst_r_1;
        %Params.tt_linkageFirst_r_2to4 = Params.raceIndicator*tt_linkageFirst_r_2to4;
    
    end

    % 2c. Linked to Care
    for linkageSection = 1:1
        
% Calculate the probability of getting linked to care after 1 year

    % Apply calibrated values
    
        %1st period
        prob_baselinkage_r_1(Params.race_B,1)=CalibParams.tt_linkage_r_1_B.paramValue(i);
        prob_baselinkage_r_1(Params.race_H,1)=CalibParams.tt_linkage_r_1_H.paramValue(i);
        prob_baselinkage_r_1(Params.race_O,1)=CalibParams.tt_linkage_r_1_O.paramValue(i);
        
        %2nd period
        prob_baselinkage_r_2to4(Params.race_B,1)=CalibParams.tt_linkage_r_2to4_B.paramValue(i);
        prob_baselinkage_r_2to4(Params.race_H,1)=CalibParams.tt_linkage_r_2to4_H.paramValue(i);
        prob_baselinkage_r_2to4(Params.race_O,1)=CalibParams.tt_linkage_r_2to4_O.paramValue(i);
    
        
   % Convert to rates
        store_rate_baseLinkageAfter_1 = -log(1-prob_baselinkage_r_1);
        store_rate_baseLinkageAfter_2to4 = -log(1-prob_baselinkage_r_2to4);
        
   % Apply to full subpopulation matrix [273x1]
        rate_baseLinkageAfter_1 = Params.raceIndicator * store_rate_baseLinkageAfter_1;
        rate_baseLinkageAfter_2to4 = Params.raceIndicator * store_rate_baseLinkageAfter_2to4;   
        
       
   % Read in new relative risks
            % Note: only the Late and AIDS stages are calibrated
        Params.store_relRiskLinkageAfter_1(Params.stage_Late) = CalibParams.tt_RelRiskLTC_1_E.paramValue(i);
        Params.store_relRiskLinkageAfter_1(Params.stage_AIDS) = CalibParams.tt_RelRiskLTC_1_F.paramValue(i);
        
        Params.store_relRiskLinkageAfter_2to4(Params.stage_Late) = CalibParams.tt_RelRiskLTC_2to4_E.paramValue(i);
        Params.store_relRiskLinkageAfter_2to4(Params.stage_AIDS) = CalibParams.tt_RelRiskLTC_2to4_F.paramValue(i);
        
   % Apply relative risks to calculate new parameters
   
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
          
        
    % Period 2

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

           
   
    
    end
    
    % 2d. Drop out
    for dropOutSection = 1:1
        
    % [Care] -> [Aware]

        %1st period
        tt_dropOutProbCareToAware_r_1(Params.race_B,1) = CalibParams.tt_dropOutProb_CareToAware_1_B.paramValue(i);
        tt_dropOutProbCareToAware_r_1(Params.race_H,1) = CalibParams.tt_dropOutProb_CareToAware_1_H.paramValue(i);
        tt_dropOutProbCareToAware_r_1(Params.race_O,1) = CalibParams.tt_dropOutProb_CareToAware_1_O.paramValue(i);
        
        %2nd period
        tt_dropOutProbCareToAware_r_2to4(Params.race_B,1) = CalibParams.tt_dropOutProb_CareToAware_2to4_B.paramValue(i);
        tt_dropOutProbCareToAware_r_2to4(Params.race_H,1) = CalibParams.tt_dropOutProb_CareToAware_2to4_H.paramValue(i);
        tt_dropOutProbCareToAware_r_2to4(Params.race_O,1) = CalibParams.tt_dropOutProb_CareToAware_2to4_O.paramValue(i);

        
        % Convert to Rates
        dropOutRateCareToAware_1 = -log(1-tt_dropOutProbCareToAware_r_1);
        dropOutRateCareToAware_2to4 = -log(1-tt_dropOutProbCareToAware_r_2to4);
        
        % Calculate new parameter
        Params.tt_dropOutRate_CareToAware_1 = Params.raceIndicator * dropOutRateCareToAware_1;
        Params.tt_dropOutRate_CareToAware_2to4 = Params.raceIndicator * dropOutRateCareToAware_2to4;    
        
        
   % [ANV] -> [Aware] (removed on 2018/01/25 by JC)
        
        %1st period
        %tt_dropOutProbANVToAware_r_1(Params.race_B,1) = CalibParams.tt_dropOutProb_ANVToAware_1_B.paramValue(i);
        %tt_dropOutProbANVToAware_r_1(Params.race_H,1) = CalibParams.tt_dropOutProb_ANVToAware_1_H.paramValue(i);
        %tt_dropOutProbANVToAware_r_1(Params.race_O,1) = CalibParams.tt_dropOutProb_ANVToAware_1_O.paramValue(i);
        
        %2nd period
        %tt_dropOutProbANVToAware_r_2to4(Params.race_B,1) = CalibParams.tt_dropOutProb_ANVToAware_2to4_B.paramValue(i);
        %tt_dropOutProbANVToAware_r_2to4(Params.race_H,1) = CalibParams.tt_dropOutProb_ANVToAware_2to4_H.paramValue(i);
        %tt_dropOutProbANVToAware_r_2to4(Params.race_O,1) = CalibParams.tt_dropOutProb_ANVToAware_2to4_O.paramValue(i);

        
        % Convert to Rates
        %dropOutRateANVToAware_1 = -log(1-tt_dropOutProbANVToAware_r_1);
        %dropOutRateANVToAware_2to4 = -log(1-tt_dropOutProbANVToAware_r_2to4);
        
        % Calculate new parameter
        %Params.tt_dropOutRate_ANVToAware_1 = Params.raceIndicator * dropOutRateANVToAware_1;
        %Params.tt_dropOutRate_ANVToAware_2to4 = Params.raceIndicator * dropOutRateANVToAware_2to4;    
   
   
   % [ANV] -> [Care]
   
        %1st period
        tt_dropOutProbANVToCare_r_1(Params.race_B,1) = CalibParams.tt_dropOutProb_ANVToCare_1_B.paramValue(i);
        tt_dropOutProbANVToCare_r_1(Params.race_H,1) = CalibParams.tt_dropOutProb_ANVToCare_1_H.paramValue(i);
        tt_dropOutProbANVToCare_r_1(Params.race_O,1) = CalibParams.tt_dropOutProb_ANVToCare_1_O.paramValue(i);
        
        %2nd period
        tt_dropOutProbANVToCare_r_2to4(Params.race_B,1) = CalibParams.tt_dropOutProb_ANVToCare_2to4_B.paramValue(i);
        tt_dropOutProbANVToCare_r_2to4(Params.race_H,1) = CalibParams.tt_dropOutProb_ANVToCare_2to4_H.paramValue(i);
        tt_dropOutProbANVToCare_r_2to4(Params.race_O,1) = CalibParams.tt_dropOutProb_ANVToCare_2to4_O.paramValue(i);

        
        %Calculate relative risk of dropping off ANV to care by time period
        %1st period
        store_ttprobdropOutANVToCare_Overall_1 = CalibParams.tt_dropOutProb_ANVToCare_1_B.paramValue(i) * (Params.InitHIVPrevalence_Race(1,1)/sum(Params.InitHIVPrevalence_Race)) ...
            + CalibParams.tt_dropOutProb_ANVToCare_1_H.paramValue(i) * (Params.InitHIVPrevalence_Race(2,1)/sum(Params.InitHIVPrevalence_Race)) ...
            + CalibParams.tt_dropOutProb_ANVToCare_1_O.paramValue(i) * (Params.InitHIVPrevalence_Race(3,1)/sum(Params.InitHIVPrevalence_Race));
        
        store_relRiskDropOutANVAge_1(1,1) = Params.oddsRatio_DropOutANVAge(1,1) / (1 - store_ttprobdropOutANVToCare_Overall_1 ...
            + (store_ttprobdropOutANVToCare_Overall_1 * Params.oddsRatio_DropOutANVAge(1,1)));
        store_relRiskDropOutANVAge_1(2,1) = Params.oddsRatio_DropOutANVAge(2,1) / (1 - store_ttprobdropOutANVToCare_Overall_1 ...
            + (store_ttprobdropOutANVToCare_Overall_1 * Params.oddsRatio_DropOutANVAge(2,1)));
        store_relRiskDropOutANVAge_1(3,1) = Params.oddsRatio_DropOutANVAge(3,1) / (1 - store_ttprobdropOutANVToCare_Overall_1 ...
            + (store_ttprobdropOutANVToCare_Overall_1 * Params.oddsRatio_DropOutANVAge(3,1)));
        store_relRiskDropOutANVAge_1(4,1) = Params.oddsRatio_DropOutANVAge(4,1) / (1 - store_ttprobdropOutANVToCare_Overall_1 ...
            + (store_ttprobdropOutANVToCare_Overall_1 * Params.oddsRatio_DropOutANVAge(4,1)));
        store_relRiskDropOutANVAge_1(5,1) = Params.oddsRatio_DropOutANVAge(5,1) / (1 - store_ttprobdropOutANVToCare_Overall_1 ...
            + (store_ttprobdropOutANVToCare_Overall_1 * Params.oddsRatio_DropOutANVAge(5,1)));
        store_relRiskDropOutANVAge_1(6,1) = Params.oddsRatio_DropOutANVAge(6,1) / (1 - store_ttprobdropOutANVToCare_Overall_1 ...
            + (store_ttprobdropOutANVToCare_Overall_1 * Params.oddsRatio_DropOutANVAge(6,1)));
        store_relRiskDropOutANVAge_1(7,1) = Params.oddsRatio_DropOutANVAge(7,1) / (1 - store_ttprobdropOutANVToCare_Overall_1 ...
            + (store_ttprobdropOutANVToCare_Overall_1 * Params.oddsRatio_DropOutANVAge(7,1)));
        
        %2nd period
        store_ttprobdropOutANVToCare_Overall_2to4 = CalibParams.tt_dropOutProb_ANVToCare_2to4_B.paramValue(i) * (Params.InitHIVPrevalence_Race(1,1)/sum(Params.InitHIVPrevalence_Race)) ...
            + CalibParams.tt_dropOutProb_ANVToCare_2to4_H.paramValue(i) * (Params.InitHIVPrevalence_Race(2,1)/sum(Params.InitHIVPrevalence_Race)) ...
            + CalibParams.tt_dropOutProb_ANVToCare_2to4_O.paramValue(i) * (Params.InitHIVPrevalence_Race(3,1)/sum(Params.InitHIVPrevalence_Race));
        
        store_relRiskDropOutANVAge_2to4(1,1) = Params.oddsRatio_DropOutANVAge(1,1) / (1 - store_ttprobdropOutANVToCare_Overall_2to4 ...
            + (store_ttprobdropOutANVToCare_Overall_2to4 * Params.oddsRatio_DropOutANVAge(1,1)));
        store_relRiskDropOutANVAge_2to4(2,1) = Params.oddsRatio_DropOutANVAge(2,1) / (1 - store_ttprobdropOutANVToCare_Overall_2to4 ...
            + (store_ttprobdropOutANVToCare_Overall_2to4 * Params.oddsRatio_DropOutANVAge(2,1)));
        store_relRiskDropOutANVAge_2to4(3,1) = Params.oddsRatio_DropOutANVAge(3,1) / (1 - store_ttprobdropOutANVToCare_Overall_2to4 ...
            + (store_ttprobdropOutANVToCare_Overall_2to4 * Params.oddsRatio_DropOutANVAge(3,1)));
        store_relRiskDropOutANVAge_2to4(4,1) = Params.oddsRatio_DropOutANVAge(4,1) / (1 - store_ttprobdropOutANVToCare_Overall_2to4 ...
            + (store_ttprobdropOutANVToCare_Overall_2to4 * Params.oddsRatio_DropOutANVAge(4,1)));
        store_relRiskDropOutANVAge_2to4(5,1) = Params.oddsRatio_DropOutANVAge(5,1) / (1 - store_ttprobdropOutANVToCare_Overall_2to4 ...
            + (store_ttprobdropOutANVToCare_Overall_2to4 * Params.oddsRatio_DropOutANVAge(5,1)));
        store_relRiskDropOutANVAge_2to4(6,1) = Params.oddsRatio_DropOutANVAge(6,1) / (1 - store_ttprobdropOutANVToCare_Overall_2to4 ...
            + (store_ttprobdropOutANVToCare_Overall_2to4 * Params.oddsRatio_DropOutANVAge(6,1)));
        store_relRiskDropOutANVAge_2to4(7,1) = Params.oddsRatio_DropOutANVAge(7,1) / (1 - store_ttprobdropOutANVToCare_Overall_2to4 ...
            + (store_ttprobdropOutANVToCare_Overall_2to4 * Params.oddsRatio_DropOutANVAge(7,1)));
        
        % 
        tt_relRiskDropOutANVAge_1 = Params.ageIndicator * store_relRiskDropOutANVAge_1;
        tt_relRiskDropOutANVAge_2to4 = Params.ageIndicator * store_relRiskDropOutANVAge_2to4;
               
        % Convert to Rates
        dropOutRateANVToCare_1 = -log(1-tt_dropOutProbANVToCare_r_1);
        dropOutRateANVToCare_2to4 = -log(1-tt_dropOutProbANVToCare_r_2to4);
        
        % Calculate new parameter
        tt_dropOutRate_ANVToCare_1 = Params.raceIndicator * dropOutRateANVToCare_1;
        tt_dropOutRate_ANVToCare_2to4 = Params.raceIndicator * dropOutRateANVToCare_2to4;    
   
        Params.tt_dropOutRate_ANVToCare_1 = tt_dropOutRate_ANVToCare_1 .* tt_relRiskDropOutANVAge_1;
        Params.tt_dropOutRate_ANVToCare_2to4 = tt_dropOutRate_ANVToCare_2to4 .* tt_relRiskDropOutANVAge_2to4;
    
   % [VLS] -> [ANV] (edited by JC on 11/13/2017)
   
        %1st period
        tt_dropOutProbVLSToANV_r_1(Params.race_B,1) = CalibParams.tt_dropOutProb_VLSToANV_1_B.paramValue(i);
        tt_dropOutProbVLSToANV_r_1(Params.race_H,1) = CalibParams.tt_dropOutProb_VLSToANV_1_H.paramValue(i);
        tt_dropOutProbVLSToANV_r_1(Params.race_O,1) = CalibParams.tt_dropOutProb_VLSToANV_1_O.paramValue(i);
        
        %2nd period
        tt_dropOutProbVLSToANV_r_2to4(Params.race_B,1) = CalibParams.tt_dropOutProb_VLSToANV_2to4_B.paramValue(i);
        tt_dropOutProbVLSToANV_r_2to4(Params.race_H,1) = CalibParams.tt_dropOutProb_VLSToANV_2to4_H.paramValue(i);
        tt_dropOutProbVLSToANV_r_2to4(Params.race_O,1) = CalibParams.tt_dropOutProb_VLSToANV_2to4_O.paramValue(i);

        
        % Convert to Rates
        dropOutRateVLSToANV_1 = -log(1-tt_dropOutProbVLSToANV_r_1);
        dropOutRateVLSToANV_2to4 = -log(1-tt_dropOutProbVLSToANV_r_2to4);
        
        % Calculate new parameter
        tt_dropOutRate_VLSToANV_1 = Params.raceIndicator * dropOutRateVLSToANV_1;
        tt_dropOutRate_VLSToANV_2to4 = Params.raceIndicator * dropOutRateVLSToANV_2to4;  
        
        % Update RR by pop
        store_relRiskLoseVLSPop_1(Params.pop_HET,1) = Params.store_relRiskLoseVLSPop_1(Params.pop_HET,1);
        store_relRiskLoseVLSPop_1(Params.pop_MSM,1) = CalibParams.tt_relRiskPop_VLSToANV_1_MSM.paramValue(i);
        store_relRiskLoseVLSPop_1(Params.pop_IDU,1) = CalibParams.tt_relRiskPop_VLSToANV_1_PWID.paramValue(i);
        
        store_relRiskLoseVLSPop_2to4(Params.pop_HET,1) = Params.store_relRiskLoseVLSPop_2to4(Params.pop_HET,1);
        store_relRiskLoseVLSPop_2to4(Params.pop_MSM,1) = CalibParams.tt_relRiskPop_VLSToANV_2to4_MSM.paramValue(i);
        store_relRiskLoseVLSPop_2to4(Params.pop_IDU,1) = CalibParams.tt_relRiskPop_VLSToANV_2to4_PWID.paramValue(i);
        
        tt_relRiskLoseVLSPop_1 = Params.popIndicator * store_relRiskLoseVLSPop_1;
        tt_relRiskLoseVLSPop_2to4 = Params.popIndicator * store_relRiskLoseVLSPop_2to4;
        
        % Update RR by age
        store_relRiskLoseVLSAge_1(Params.age_1,1) = Params.store_relRiskLoseVLSAge_1(Params.age_1,1);
        store_relRiskLoseVLSAge_1(Params.age_2,1) = CalibParams.tt_relRiskAge_VLSToANV_1_1824.paramValue(i);
        store_relRiskLoseVLSAge_1(Params.age_3,1) = CalibParams.tt_relRiskAge_VLSToANV_1_2534.paramValue(i);
        store_relRiskLoseVLSAge_1(Params.age_4,1) = CalibParams.tt_relRiskAge_VLSToANV_1_3544.paramValue(i);
        store_relRiskLoseVLSAge_1(Params.age_5,1) = CalibParams.tt_relRiskAge_VLSToANV_1_4554.paramValue(i);
        store_relRiskLoseVLSAge_1(Params.age_6,1) = CalibParams.tt_relRiskAge_VLSToANV_1_5564.paramValue(i);
        store_relRiskLoseVLSAge_1(Params.age_7,1) = CalibParams.tt_relRiskAge_VLSToANV_1_65.paramValue(i);
        
        store_relRiskLoseVLSAge_2to4(Params.age_1,1) = Params.store_relRiskLoseVLSAge_2to4(Params.age_1,1);
        store_relRiskLoseVLSAge_2to4(Params.age_2,1) = CalibParams.tt_relRiskAge_VLSToANV_2to4_1824.paramValue(i);
        store_relRiskLoseVLSAge_2to4(Params.age_3,1) = CalibParams.tt_relRiskAge_VLSToANV_2to4_2534.paramValue(i);
        store_relRiskLoseVLSAge_2to4(Params.age_4,1) = CalibParams.tt_relRiskAge_VLSToANV_2to4_3544.paramValue(i);
        store_relRiskLoseVLSAge_2to4(Params.age_5,1) = CalibParams.tt_relRiskAge_VLSToANV_2to4_4554.paramValue(i);
        store_relRiskLoseVLSAge_2to4(Params.age_6,1) = CalibParams.tt_relRiskAge_VLSToANV_2to4_5564.paramValue(i);
        store_relRiskLoseVLSAge_2to4(Params.age_7,1) = CalibParams.tt_relRiskAge_VLSToANV_2to4_65.paramValue(i);
        
        tt_relRiskLoseVLSAge_1 = Params.ageIndicator * store_relRiskLoseVLSAge_1;
        tt_relRiskLoseVLSAge_2to4 = Params.ageIndicator * store_relRiskLoseVLSAge_2to4;
        
        % Update VLS-to-ANV drop-out rates
        Params.tt_dropOutRate_VLSToANV_1 = tt_dropOutRate_VLSToANV_1 .* tt_relRiskLoseVLSPop_1 .* tt_relRiskLoseVLSAge_1;
        Params.tt_dropOutRate_VLSToANV_2to4 = tt_dropOutRate_VLSToANV_2to4 .* tt_relRiskLoseVLSPop_2to4 .* tt_relRiskLoseVLSAge_2to4;
        
    end

    % 2e. ART Initiation
    for artInitSection = 1:1
        
    % Update the ART initiation (removed on 2018/01/25 by JC)
    
        % Convert to rates and remove ineligible stages
        %rate_ARTInitRate_B_1 = -log(1-CalibParams.tt_ARTInitiation_d_1_B.paramValue(i)) * Params.tt_ARTElig_1(Params.stage_Acute);
        %rate_ARTInitRate_C_1 = -log(1-CalibParams.tt_ARTInitiation_d_1_C.paramValue(i)) * Params.tt_ARTElig_1(Params.stage_LatentA);
        %rate_ARTInitRate_D_1 = -log(1-CalibParams.tt_ARTInitiation_d_1_D.paramValue(i)) * Params.tt_ARTElig_1(Params.stage_LatentB);
        %rate_ARTInitRate_E_1 = -log(1-CalibParams.tt_ARTInitiation_d_1_E.paramValue(i)) * Params.tt_ARTElig_1(Params.stage_Late);
        %rate_ARTInitRate_F_1 = -log(1-CalibParams.tt_ARTInitiation_d_1_F.paramValue(i)) * Params.tt_ARTElig_1(Params.stage_AIDS);
    
        %rate_ARTInitRate_B_2to4 = -log(1-CalibParams.tt_ARTInitiation_d_2to4_B.paramValue(i)) * Params.tt_ARTElig_2to4(Params.stage_Acute);
        %rate_ARTInitRate_C_2to4 = -log(1-CalibParams.tt_ARTInitiation_d_2to4_C.paramValue(i)) * Params.tt_ARTElig_2to4(Params.stage_LatentA);
        %rate_ARTInitRate_D_2to4 = -log(1-CalibParams.tt_ARTInitiation_d_2to4_D.paramValue(i)) * Params.tt_ARTElig_2to4(Params.stage_LatentB);
        %rate_ARTInitRate_E_2to4 = -log(1-CalibParams.tt_ARTInitiation_d_2to4_E.paramValue(i)) * Params.tt_ARTElig_2to4(Params.stage_Late);
        %rate_ARTInitRate_F_2to4 = -log(1-CalibParams.tt_ARTInitiation_d_2to4_F.paramValue(i)) * Params.tt_ARTElig_2to4(Params.stage_AIDS);
        
        
   % Update relative risks of ART initiation by race/ethnicity
    
        % Apply calibration value
        store_relRiskARTInitRace_1(Params.race_B,1) = CalibParams.tt_RelRiskARTInit_r_1_B.paramValue(i);
        store_relRiskARTInitRace_1(Params.race_H,1) = CalibParams.tt_RelRiskARTInit_r_1_H.paramValue(i);
        store_relRiskARTInitRace_1(Params.race_O,1) = CalibParams.tt_RelRiskARTInit_r_1_O.paramValue(i);
   
        store_relRiskARTInitRace_2to4(Params.race_B,1) = CalibParams.tt_RelRiskARTInit_r_2to4_B.paramValue(i);
        store_relRiskARTInitRace_2to4(Params.race_H,1) = CalibParams.tt_RelRiskARTInit_r_2to4_H.paramValue(i);
        store_relRiskARTInitRace_2to4(Params.race_O,1) = CalibParams.tt_RelRiskARTInit_r_2to4_O.paramValue(i);
        
        % Expand to full subpopulation matrix [273x1]
        relRiskARTInit_race_1 = Params.raceIndicator * store_relRiskARTInitRace_1;
        relRiskARTInit_race_2to4 = Params.raceIndicator * store_relRiskARTInitRace_2to4;
        
   % Apply to parameters and apply relative risk of initiation by race
   
       %1st period
       Params.tt_ARTInitRateAcute_1 =  Params.rate_ARTInitRate_B_1 * relRiskARTInit_race_1;
       Params.tt_ARTInitRateLatentA_1 = Params.rate_ARTInitRate_C_1 * relRiskARTInit_race_1;
       Params.tt_ARTInitRateLatentB_1 = Params.rate_ARTInitRate_D_1 * relRiskARTInit_race_1;
       Params.tt_ARTInitRateLate_1 = Params.rate_ARTInitRate_E_1 * relRiskARTInit_race_1;
       Params.tt_ARTInitRateAIDS_1 = Params.rate_ARTInitRate_F_1 * relRiskARTInit_race_1;

       %2nd period
       Params.tt_ARTInitRateAcute_2to4 =  Params.rate_ARTInitRate_B_2to4 * relRiskARTInit_race_2to4;
       Params.tt_ARTInitRateLatentA_2to4 = Params.rate_ARTInitRate_C_2to4 * relRiskARTInit_race_2to4;
       Params.tt_ARTInitRateLatentB_2to4 = Params.rate_ARTInitRate_D_2to4 * relRiskARTInit_race_2to4;
       Params.tt_ARTInitRateLate_2to4 = Params.rate_ARTInitRate_E_2to4 * relRiskARTInit_race_2to4;
       Params.tt_ARTInitRateAIDS_2to4 = Params.rate_ARTInitRate_F_2to4 * relRiskARTInit_race_2to4;
   
    end

    % 2f. Become VLS from ANV (edited by JC on 11/13/2017)
    
        %1st period
        tt_BecomeVLSfromANV_r_1(Params.race_B,1) = CalibParams.tt_BecomeVLSfromANV_r_1_B.paramValue(i);
        tt_BecomeVLSfromANV_r_1(Params.race_H,1) = CalibParams.tt_BecomeVLSfromANV_r_1_H.paramValue(i);
        tt_BecomeVLSfromANV_r_1(Params.race_O,1) = CalibParams.tt_BecomeVLSfromANV_r_1_O.paramValue(i);
        
        %2nd period
        tt_BecomeVLSfromANV_r_2to4(Params.race_B,1) = CalibParams.tt_BecomeVLSfromANV_r_2to4_B.paramValue(i);
        tt_BecomeVLSfromANV_r_2to4(Params.race_H,1) = CalibParams.tt_BecomeVLSfromANV_r_2to4_H.paramValue(i);
        tt_BecomeVLSfromANV_r_2to4(Params.race_O,1) = CalibParams.tt_BecomeVLSfromANV_r_2to4_O.paramValue(i);

        
        % Convert to Rates
        BecomeVLSfromANV_1 = -log(1-tt_BecomeVLSfromANV_r_1);
        BecomeVLSfromANV_2to4 = -log(1-tt_BecomeVLSfromANV_r_2to4);
        
        % Calculate new parameter
        tt_ARTInitRateFromANV_1 = Params.raceIndicator * BecomeVLSfromANV_1;
        tt_ARTInitRateFromANV_2to4 = Params.raceIndicator * BecomeVLSfromANV_2to4;
        
        % Update RR by pop
        store_relRiskBecomeVLSPop_1(Params.pop_HET,1) = Params.store_relRiskBecomeVLSPop_1(Params.pop_HET,1);
        store_relRiskBecomeVLSPop_1(Params.pop_MSM,1) = CalibParams.tt_relRiskPop_ANVToVLS_1_MSM.paramValue(i);
        store_relRiskBecomeVLSPop_1(Params.pop_IDU,1) = CalibParams.tt_relRiskPop_ANVToVLS_1_PWID.paramValue(i);
        
        store_relRiskBecomeVLSPop_2to4(Params.pop_HET,1) = Params.store_relRiskBecomeVLSPop_2to4(Params.pop_HET,1);
        store_relRiskBecomeVLSPop_2to4(Params.pop_MSM,1) = CalibParams.tt_relRiskPop_ANVToVLS_2to4_MSM.paramValue(i);
        store_relRiskBecomeVLSPop_2to4(Params.pop_IDU,1) = CalibParams.tt_relRiskPop_ANVToVLS_2to4_PWID.paramValue(i);
        
        tt_relRiskBecomeVLSPop_1 = Params.popIndicator * store_relRiskBecomeVLSPop_1;
        tt_relRiskBecomeVLSPop_2to4 = Params.popIndicator * store_relRiskBecomeVLSPop_2to4;
        
        % Update RR by age
        store_relRiskBecomeVLSAge_1(Params.age_1,1) = Params.store_relRiskBecomeVLSAge_1(Params.age_1,1);
        store_relRiskBecomeVLSAge_1(Params.age_2,1) = CalibParams.tt_relRiskAge_ANVToVLS_1_1824.paramValue(i);
        store_relRiskBecomeVLSAge_1(Params.age_3,1) = CalibParams.tt_relRiskAge_ANVToVLS_1_2534.paramValue(i);
        store_relRiskBecomeVLSAge_1(Params.age_4,1) = CalibParams.tt_relRiskAge_ANVToVLS_1_3544.paramValue(i);
        store_relRiskBecomeVLSAge_1(Params.age_5,1) = CalibParams.tt_relRiskAge_ANVToVLS_1_4554.paramValue(i);
        store_relRiskBecomeVLSAge_1(Params.age_6,1) = CalibParams.tt_relRiskAge_ANVToVLS_1_5564.paramValue(i);
        store_relRiskBecomeVLSAge_1(Params.age_7,1) = CalibParams.tt_relRiskAge_ANVToVLS_1_65.paramValue(i);
        
        store_relRiskBecomeVLSAge_2to4(Params.age_1,1) = Params.store_relRiskBecomeVLSAge_2to4(Params.age_1,1);
        store_relRiskBecomeVLSAge_2to4(Params.age_2,1) = CalibParams.tt_relRiskAge_ANVToVLS_2to4_1824.paramValue(i);
        store_relRiskBecomeVLSAge_2to4(Params.age_3,1) = CalibParams.tt_relRiskAge_ANVToVLS_2to4_2534.paramValue(i);
        store_relRiskBecomeVLSAge_2to4(Params.age_4,1) = CalibParams.tt_relRiskAge_ANVToVLS_2to4_3544.paramValue(i);
        store_relRiskBecomeVLSAge_2to4(Params.age_5,1) = CalibParams.tt_relRiskAge_ANVToVLS_2to4_4554.paramValue(i);
        store_relRiskBecomeVLSAge_2to4(Params.age_6,1) = CalibParams.tt_relRiskAge_ANVToVLS_2to4_5564.paramValue(i);
        store_relRiskBecomeVLSAge_2to4(Params.age_7,1) = CalibParams.tt_relRiskAge_ANVToVLS_2to4_65.paramValue(i);
        
        tt_relRiskBecomeVLSAge_1 = Params.ageIndicator * store_relRiskBecomeVLSAge_1;
        tt_relRiskBecomeVLSAge_2to4 = Params.ageIndicator * store_relRiskBecomeVLSAge_2to4;
        
        % Update VLS-to-ANV drop-out rates
        Params.tt_ARTInitRateFromANV_1 = tt_ARTInitRateFromANV_1 .* tt_relRiskBecomeVLSPop_1 .* tt_relRiskBecomeVLSAge_1;
        Params.tt_ARTInitRateFromANV_2to4 = tt_ARTInitRateFromANV_2to4 .* tt_relRiskBecomeVLSPop_2to4 .* tt_relRiskBecomeVLSAge_2to4;
                    
    
    % 2g. Distribution of ANV (removed on 2018/01/25 by JC)
    for anvDistSection = 1:1
          
        % Percent of ANV in care (ART + not on ART)
       %Params.tt_PctANVWhoAreInCareOrART_1 = CalibParams.tt_PctANVWhoAreInCareOrART_1.paramValue(i);
       %Params.tt_PctANVWhoAreInCareOrART_2to4 = CalibParams.tt_PctANVWhoAreInCareOrART_2to4.paramValue(i);

       % Percent of in-care ANV who are on ART
       %Params.tt_PctInCareANVWhoAreOnART_1 = CalibParams.tt_PctInCareANVWhoAreOnART_1.paramValue(i);
       %Params.tt_PctInCareANVWhoAreOnART_2to4 = CalibParams.tt_PctInCareANVWhoAreOnART_2to4.paramValue(i);
        
    end

%% 3. Update Disease Progression Parameters
   
    % Length of time in ART-not-VLS (modified 08/03/2018)
        % Applied to PLWH in [ART-not-VLS] compartments
        
        % Note: Acute individuals not eligible for ART
        % Note: AIDS/ANV receive ART mortality (from NA-ACCORD)
        
        hiv_durHIVStage_ANV_LatentA(1,1:4) = CalibParams.hiv_durStage_ANV_LatentA_1344.paramValue(i);
        hiv_durHIVStage_ANV_LatentA(1,5:Params.numAge-1) = CalibParams.hiv_durStage_ANV_LatentA_4564.paramValue(i);
        hiv_durHIVStage_ANV_LatentA(1,Params.numAge) = CalibParams.hiv_durStage_ANV_LatentA_65.paramValue(i);
        
        hiv_durHIVStage_ANV_LatentB(1,1:4) = CalibParams.hiv_durStage_ANV_LatentB_1344.paramValue(i);
        hiv_durHIVStage_ANV_LatentB(1,5:Params.numAge-1) = CalibParams.hiv_durStage_ANV_LatentB_4564.paramValue(i);
        hiv_durHIVStage_ANV_LatentB(1,Params.numAge) = CalibParams.hiv_durStage_ANV_LatentB_65.paramValue(i);
        
        hiv_durHIVStage_ANV_Late(1,1:4) = CalibParams.hiv_durStage_ANV_Late_1344.paramValue(i);
        hiv_durHIVStage_ANV_Late(1,5:Params.numAge-1) = CalibParams.hiv_durStage_ANV_Late_4564.paramValue(i);
        hiv_durHIVStage_ANV_Late(1,Params.numAge) = CalibParams.hiv_durStage_ANV_Late_65.paramValue(i);
        
        Params.hiv_durHIVStage_ANV_LatentA = hiv_durHIVStage_ANV_LatentA * Params.ageIndicator';
        Params.hiv_durHIVStage_ANV_LatentB = hiv_durHIVStage_ANV_LatentB * Params.ageIndicator';
        Params.hiv_durHIVStage_ANV_Late = hiv_durHIVStage_ANV_Late * Params.ageIndicator';        
       
        
% DISEASE STAGE TRANSITION RATES WHILE VLS

    % Decrease in CD4 count 
        % Rate of transitioning to a lower HIV stage while VLS

        % Apply to parameters
        Params.hiv_rateVLSCD4decr_LatentA = CalibParams.hiv_rateVLSCD4decr_LatentA.paramValue(i);
        Params.hiv_rateVLSCD4decr_LatentB = CalibParams.hiv_rateVLSCD4decr_LatentB.paramValue(i);
        Params.hiv_rateVLSCD4decr_Late = CalibParams.hiv_rateVLSCD4decr_Late.paramValue(i);
       
        
    % Increase in CD4 count (modified 08/03/2018)
        % Rate of transitioning to a higher HIV stage while VLS
        hiv_rateVLSCD4incr_LatentB(1,1:4) = CalibParams.hiv_rateVLSCD4incr_LatentB_1344.paramValue(i);
        hiv_rateVLSCD4incr_LatentB(1,5:Params.numAge-1) = CalibParams.hiv_rateVLSCD4incr_LatentB_4564.paramValue(i);
        hiv_rateVLSCD4incr_LatentB(1,Params.numAge) = CalibParams.hiv_rateVLSCD4incr_LatentB_65.paramValue(i);
        
        hiv_rateVLSCD4incr_Late(1,1:4) = CalibParams.hiv_rateVLSCD4incr_Late_1344.paramValue(i);
        hiv_rateVLSCD4incr_Late(1,5:Params.numAge-1) = CalibParams.hiv_rateVLSCD4incr_Late_4564.paramValue(i);
        hiv_rateVLSCD4incr_Late(1,Params.numAge) = CalibParams.hiv_rateVLSCD4incr_Late_65.paramValue(i);
        
        hiv_rateVLSCD4incr_AIDS(1,1:4) = CalibParams.hiv_rateVLSCD4incr_AIDS_1344.paramValue(i);
        hiv_rateVLSCD4incr_AIDS(1,5:Params.numAge-1) = CalibParams.hiv_rateVLSCD4incr_AIDS_4564.paramValue(i);
        hiv_rateVLSCD4incr_AIDS(1,Params.numAge) = CalibParams.hiv_rateVLSCD4incr_AIDS_65.paramValue(i);
        
        % Apply to parameters
        Params.hiv_rateVLSCD4incr_LatentB = hiv_rateVLSCD4incr_LatentB * Params.ageIndicator';
        Params.hiv_rateVLSCD4incr_Late = hiv_rateVLSCD4incr_Late * Params.ageIndicator';
        Params.hiv_rateVLSCD4incr_AIDS = hiv_rateVLSCD4incr_AIDS * Params.ageIndicator';        

%% 4. Update the Sexual Mixing Matrix  
    
    % TRANSMISSION GROUP / SEX MIXING
    
        % Calculate partnerships that are dependent on calibrated values
         % KH updated on 20Sep23 to make HET mixing for HRHETs and LRHETs

            behav_SexualMixing_LRHETM_IDUF = 1 - ...
                CalibParams.behav_SexualMixing_LRHETM_HETF.paramValue(i);

            behav_SexualMixing_LRHETF_MSM = 1 - ...
                 CalibParams.behav_SexualMixing_LRHETF_HETM.paramValue(i) - ...
                 CalibParams.behav_SexualMixing_LRHETF_IDUM.paramValue(i);
            
            behav_SexualMixing_HRHETM_IDUF = 1 - ...
                CalibParams.behav_SexualMixing_HRHETM_HETF.paramValue(i);

            behav_SexualMixing_HRHETF_MSM = 1 - ...
                 CalibParams.behav_SexualMixing_HRHETF_HETM.paramValue(i) - ...
                 CalibParams.behav_SexualMixing_HRHETF_IDUM.paramValue(i);

            behav_SexualMixing_MSM_MSM = 1 - ...
                CalibParams.behav_SexualMixing_MSM_HETF.paramValue(i) - ...
                CalibParams.behav_SexualMixing_MSM_IDUF.paramValue(i);

            behav_SexualMixing_IDUM_HETF = 1 - ...
                CalibParams.behav_SexualMixing_IDUM_IDUF.paramValue(i);

            behav_SexualMixing_IDUF_HETM = 1 -...
                CalibParams.behav_SexualMixing_IDUF_IDUM.paramValue(i) - ...
                Params.mix_IDUF_MSM_value;


        % Calculate the 273x273 matrix for each combination of trans group/sex
        % for sexual mixing
         % KH updated on 20Sep23 to make HET mixing for HRHETs and LRHETs

            % Applying individual factors to a 273x273 matrix
            mix_LRHETM_HETF = Params.idx_LRHETM_HETF * CalibParams.behav_SexualMixing_LRHETM_HETF.paramValue(i);
            mix_LRHETM_IDUF = Params.idx_LRHETM_IDUF * behav_SexualMixing_LRHETM_IDUF;
            mix_LRHETF_HETM = Params.idx_LRHETF_HETM * CalibParams.behav_SexualMixing_LRHETF_HETM.paramValue(i) ;
            mix_LRHETF_MSM = Params.idx_LRHETF_MSM * behav_SexualMixing_LRHETF_MSM; 
            mix_LRHETF_IDUM = Params.idx_LRHETF_IDUM * CalibParams.behav_SexualMixing_LRHETF_IDUM.paramValue(i); 
            mix_HRHETM_HETF = Params.idx_HRHETM_HETF * CalibParams.behav_SexualMixing_HRHETM_HETF.paramValue(i);
            mix_HRHETM_IDUF = Params.idx_HRHETM_IDUF * behav_SexualMixing_HRHETM_IDUF;
            mix_HRHETF_HETM = Params.idx_HRHETF_HETM * CalibParams.behav_SexualMixing_HRHETF_HETM.paramValue(i) ;
            mix_HRHETF_MSM = Params.idx_HRHETF_MSM * behav_SexualMixing_HRHETF_MSM; 
            mix_HRHETF_IDUM = Params.idx_HRHETF_IDUM * CalibParams.behav_SexualMixing_HRHETF_IDUM.paramValue(i); 
            mix_MSM_HETF = Params.idx_MSM_HETF * CalibParams.behav_SexualMixing_MSM_HETF.paramValue(i);
            mix_MSM_MSM = Params.idx_MSM_MSM * behav_SexualMixing_MSM_MSM;
            mix_MSM_IDUF = Params.idx_MSM_IDUF * CalibParams.behav_SexualMixing_MSM_IDUF.paramValue(i);
            mix_IDUM_HETF = Params.idx_IDUM_HETF * behav_SexualMixing_IDUM_HETF;
            mix_IDUM_IDUF = Params.idx_IDUM_IDUF * CalibParams.behav_SexualMixing_IDUM_IDUF.paramValue(i);
            mix_IDUF_HETM = Params.idx_IDUF_HETM * behav_SexualMixing_IDUF_HETM ;
            % IDUF_MSM isn't calibrated and is initialized in InitializeParameters.m
            mix_IDUF_IDUM = Params.idx_IDUF_IDUM * CalibParams.behav_SexualMixing_IDUF_IDUM.paramValue(i);
        
            % Combine each partnership into an overall matrix
                % Note: this is equivalent to the "sexualMix_SexRiskType"
                % matrix in InitializeParameters
            mixMatrix_popSex = ...
                  mix_LRHETM_HETF + mix_LRHETM_IDUF ...
                + mix_LRHETF_HETM + mix_LRHETF_MSM + mix_LRHETF_IDUM ...
                + mix_HRHETM_HETF + mix_HRHETM_IDUF ...
                + mix_HRHETF_HETM + mix_HRHETF_MSM + mix_HRHETF_IDUM ...
                + mix_MSM_HETF + mix_MSM_MSM + mix_MSM_IDUF ...
                + mix_IDUM_HETF + mix_IDUM_IDUF ...
                + mix_IDUF_HETM + Params.mix_IDUF_MSM_matrix + mix_IDUF_IDUM;

        
    % RISK LEVEL MIXING
    
        % Read in values and format to parallel code imported from intializeParameters  
        sexMix_riskLevel_HET(1,1) = CalibParams.behav_SexualMixing_HET_LOW_LOW.paramValue(i);
        sexMix_riskLevel_HET(1,2) = (1-CalibParams.behav_SexualMixing_HET_LOW_LOW.paramValue(i));
        
        sexMix_riskLevel_HET(2,2) = CalibParams.behav_SexualMixing_HET_HIGH_HIGH.paramValue(i);
        sexMix_riskLevel_HET(2,1) = (1-CalibParams.behav_SexualMixing_HET_HIGH_HIGH.paramValue(i));
        
        sexMix_riskLevel_MSM(1,1) = CalibParams.behav_SexualMixing_MSM_LOW_LOW.paramValue(i);
        sexMix_riskLevel_MSM(1,2) = (1-CalibParams.behav_SexualMixing_MSM_LOW_LOW.paramValue(i));
        
        sexMix_riskLevel_MSM(2,2) = CalibParams.behav_SexualMixing_MSM_HIGH_HIGH.paramValue(i);
        sexMix_riskLevel_MSM(2,1) = (1-CalibParams.behav_SexualMixing_MSM_HIGH_HIGH.paramValue(i));
        
        sexMix_riskLevel_IDU = sexMix_riskLevel_HET;
                       
        % Define indicators
            % Combine risk level and transmission group
            % [273x2]: dim 1 is low risk, dim 2 is high risk
            % If the subpopulation isn't specified, both dim are 0
            
            numRiskLevel = 2;
            
            for nRL = 1:numRiskLevel
               riskLevelIndicator_HET(:,nRL) = Params.riskLevelIndicator(:,nRL) .* Params.popIndicator(:,Params.pop_HET);
            end

            for nRL = 1:numRiskLevel
               riskLevelIndicator_MSM(:,nRL) = Params.riskLevelIndicator(:,nRL) .* Params.popIndicator(:,Params.pop_MSM);
            end

            for nRL = 1:numRiskLevel
               riskLevelIndicator_IDU(:,nRL) = Params.riskLevelIndicator(:,nRL) .* Params.popIndicator(:,Params.pop_IDU);
            end

        
        % Combining the mixing factors to get the riskLevel Mixing Matrix:
            % Note: this is equivalent to the "behav_riskLevel_S" matrix in
            % InitializeParameters.m
        mixMatrix_riskLevel = ...
                riskLevelIndicator_HET * sexMix_riskLevel_HET * Params.riskLevelIndicator' ...
              + riskLevelIndicator_MSM * sexMix_riskLevel_MSM * Params.riskLevelIndicator' ...
              + riskLevelIndicator_IDU * sexMix_riskLevel_IDU * Params.riskLevelIndicator';
        
        % Adjust when IDU is partner #2
            % Whenever a population has IDU as a partner, the risk level
            % mix is 1.
            
            % Record index when IDU is partner 2
            idx_IDUisPartner2 = Params.popIndicator(:,Params.pop_IDU) == 1;
        
            % Adjust in full matrix
            mixMatrix_riskLevel(:,idx_IDUisPartner2) = 1;
        
            
    % RACE MIXING
    
        
        
        % Calculated factors for race mixing (scalars) - HET
            behav_SexualMixingHET_black_black = 1 - CalibParams.behav_SexualMixingHET_BLK_HISP.paramValue(i) - CalibParams.behav_SexualMixingHET_BLK_OTH.paramValue(i);
            behav_SexualMixingHET_hispanic_hispanic = 1 - CalibParams.behav_SexualMixingHET_HISP_BLK.paramValue(i) - CalibParams.behav_SexualMixingHET_HISP_OTH.paramValue(i);
            behav_SexualMixingHET_other_other =  1 - CalibParams.behav_SexualMixingHET_OTH_BLK.paramValue(i) - CalibParams.behav_SexualMixingHET_OTH_HISP.paramValue(i);
            
        % Apply scalars to a full 273x273 matrices - HET    
        mixHET_black_black = Params.idx_black_black .* Params.HET_IDU_indicator * behav_SexualMixingHET_black_black;
        mixHET_black_hispanic = Params.idx_black_hispanic .* Params.HET_IDU_indicator * CalibParams.behav_SexualMixingHET_BLK_HISP.paramValue(i);
        mixHET_black_other = Params.idx_black_other .* Params.HET_IDU_indicator * CalibParams.behav_SexualMixingHET_BLK_OTH.paramValue(i);
        mixHET_hispanic_black = Params.idx_hispanic_black .* Params.HET_IDU_indicator * CalibParams.behav_SexualMixingHET_HISP_BLK.paramValue(i);
        mixHET_hispanic_hispanic = Params.idx_hispanic_hispanic .* Params.HET_IDU_indicator * behav_SexualMixingHET_hispanic_hispanic;
        mixHET_hispanic_other = Params.idx_hispanic_other .* Params.HET_IDU_indicator * CalibParams.behav_SexualMixingHET_HISP_OTH.paramValue(i);
        mixHET_other_black = Params.idx_other_black .* Params.HET_IDU_indicator * CalibParams.behav_SexualMixingHET_OTH_BLK.paramValue(i);
        mixHET_other_hispanic = Params.idx_other_hispanic .* Params.HET_IDU_indicator * CalibParams.behav_SexualMixingHET_OTH_HISP.paramValue(i);
        mixHET_other_other = Params.idx_other_other .* Params.HET_IDU_indicator * behav_SexualMixingHET_other_other;
        
        mixMatrixHET_race = ...
              mixHET_black_black + mixHET_black_hispanic + mixHET_black_other ...
            + mixHET_hispanic_black + mixHET_hispanic_hispanic + mixHET_hispanic_other ...
            + mixHET_other_black + mixHET_other_hispanic + mixHET_other_other;
        
        % Calculated factors for race mixing (scalars) - MSM
            behav_SexualMixingMSM_black_black = 1 - CalibParams.behav_SexualMixingMSM_BLK_HISP.paramValue(i) - CalibParams.behav_SexualMixingMSM_BLK_OTH.paramValue(i);
            behav_SexualMixingMSM_hispanic_hispanic = 1 - CalibParams.behav_SexualMixingMSM_HISP_BLK.paramValue(i) - CalibParams.behav_SexualMixingMSM_HISP_OTH.paramValue(i);
            behav_SexualMixingMSM_other_other =  1 - CalibParams.behav_SexualMixingMSM_OTH_BLK.paramValue(i) - CalibParams.behav_SexualMixingMSM_OTH_HISP.paramValue(i);
            
        % Apply scalars to a full 273x273 matrices - MSM    
        mixMSM_black_black = Params.idx_black_black .* Params.MSM_indicator * behav_SexualMixingMSM_black_black;
        mixMSM_black_hispanic = Params.idx_black_hispanic .* Params.MSM_indicator * CalibParams.behav_SexualMixingMSM_BLK_HISP.paramValue(i);
        mixMSM_black_other = Params.idx_black_other .* Params.MSM_indicator * CalibParams.behav_SexualMixingMSM_BLK_OTH.paramValue(i);
        mixMSM_hispanic_black = Params.idx_hispanic_black .* Params.MSM_indicator * CalibParams.behav_SexualMixingMSM_HISP_BLK.paramValue(i);
        mixMSM_hispanic_hispanic = Params.idx_hispanic_hispanic .* Params.MSM_indicator * behav_SexualMixingMSM_hispanic_hispanic;
        mixMSM_hispanic_other = Params.idx_hispanic_other .* Params.MSM_indicator * CalibParams.behav_SexualMixingMSM_HISP_OTH.paramValue(i);
        mixMSM_other_black = Params.idx_other_black .* Params.MSM_indicator * CalibParams.behav_SexualMixingMSM_OTH_BLK.paramValue(i);
        mixMSM_other_hispanic = Params.idx_other_hispanic .* Params.MSM_indicator * CalibParams.behav_SexualMixingMSM_OTH_HISP.paramValue(i);
        mixMSM_other_other = Params.idx_other_other .* Params.MSM_indicator * behav_SexualMixingMSM_other_other;
        
        mixMatrixMSM_race = ...
              mixMSM_black_black + mixMSM_black_hispanic + mixMSM_black_other ...
            + mixMSM_hispanic_black + mixMSM_hispanic_hispanic + mixMSM_hispanic_other ...
            + mixMSM_other_black + mixMSM_other_hispanic + mixMSM_other_other;
        
        % Combine race-mixing matrices
        mixMatrix_race = mixMatrixHET_race + mixMatrixMSM_race;
        
    % AGE MIXING (edited by JC on 08/03/2018)
    
        % Pre-allocate
            store_mix_age_S_MSM(Params.numAge,Params.numAge)=0;
                
                store_mix_age_S_MSM(1,1) = 1 - CalibParams.behav_SexualMixingMSM_1317_1824.paramValue(i) - CalibParams.behav_SexualMixingMSM_1317_2534.paramValue(i) ...
                    - CalibParams.behav_SexualMixingMSM_1317_3544.paramValue(i) - CalibParams.behav_SexualMixingMSM_1317_4554.paramValue(i) ...
                    - CalibParams.behav_SexualMixingMSM_1317_5564.paramValue(i) - CalibParams.behav_SexualMixingMSM_1317_65.paramValue(i);
                store_mix_age_S_MSM(1,2) = CalibParams.behav_SexualMixingMSM_1317_1824.paramValue(i);
                store_mix_age_S_MSM(1,3) = CalibParams.behav_SexualMixingMSM_1317_2534.paramValue(i);
                store_mix_age_S_MSM(1,4) = CalibParams.behav_SexualMixingMSM_1317_3544.paramValue(i);
                store_mix_age_S_MSM(1,5) = CalibParams.behav_SexualMixingMSM_1317_4554.paramValue(i);
                store_mix_age_S_MSM(1,6) = CalibParams.behav_SexualMixingMSM_1317_5564.paramValue(i);
                store_mix_age_S_MSM(1,7) = CalibParams.behav_SexualMixingMSM_1317_65.paramValue(i);
                
                store_mix_age_S_MSM(2,1) = CalibParams.behav_SexualMixingMSM_1824_1317.paramValue(i);
                store_mix_age_S_MSM(2,2) = 1 - CalibParams.behav_SexualMixingMSM_1824_1317.paramValue(i) - CalibParams.behav_SexualMixingMSM_1824_2534.paramValue(i) ...
                    - CalibParams.behav_SexualMixingMSM_1824_3544.paramValue(i) - CalibParams.behav_SexualMixingMSM_1824_4554.paramValue(i) ...
                    - CalibParams.behav_SexualMixingMSM_1824_5564.paramValue(i) - CalibParams.behav_SexualMixingMSM_1824_65.paramValue(i);
                store_mix_age_S_MSM(2,3) = CalibParams.behav_SexualMixingMSM_1824_2534.paramValue(i);
                store_mix_age_S_MSM(2,4) = CalibParams.behav_SexualMixingMSM_1824_3544.paramValue(i);
                store_mix_age_S_MSM(2,5) = CalibParams.behav_SexualMixingMSM_1824_4554.paramValue(i);
                store_mix_age_S_MSM(2,6) = CalibParams.behav_SexualMixingMSM_1824_5564.paramValue(i);
                store_mix_age_S_MSM(2,7) = CalibParams.behav_SexualMixingMSM_1824_65.paramValue(i);
                
                store_mix_age_S_MSM(3,1) = CalibParams.behav_SexualMixingMSM_2534_1317.paramValue(i);
                store_mix_age_S_MSM(3,2) = CalibParams.behav_SexualMixingMSM_2534_1824.paramValue(i);
                store_mix_age_S_MSM(3,3) = 1 - CalibParams.behav_SexualMixingMSM_2534_1317.paramValue(i) - CalibParams.behav_SexualMixingMSM_2534_1824.paramValue(i) ...
                    - CalibParams.behav_SexualMixingMSM_2534_3544.paramValue(i) - CalibParams.behav_SexualMixingMSM_2534_4554.paramValue(i) ...
                    - CalibParams.behav_SexualMixingMSM_2534_5564.paramValue(i) - CalibParams.behav_SexualMixingMSM_2534_65.paramValue(i);
                store_mix_age_S_MSM(3,4) = CalibParams.behav_SexualMixingMSM_2534_3544.paramValue(i);
                store_mix_age_S_MSM(3,5) = CalibParams.behav_SexualMixingMSM_2534_4554.paramValue(i);
                store_mix_age_S_MSM(3,6) = CalibParams.behav_SexualMixingMSM_2534_5564.paramValue(i);
                store_mix_age_S_MSM(3,7) = CalibParams.behav_SexualMixingMSM_2534_65.paramValue(i);
                
                store_mix_age_S_MSM(4,1) = CalibParams.behav_SexualMixingMSM_3544_1317.paramValue(i);
                store_mix_age_S_MSM(4,2) = CalibParams.behav_SexualMixingMSM_3544_1824.paramValue(i);
                store_mix_age_S_MSM(4,3) = CalibParams.behav_SexualMixingMSM_3544_2534.paramValue(i);
                store_mix_age_S_MSM(4,4) = 1 - CalibParams.behav_SexualMixingMSM_3544_1317.paramValue(i) - CalibParams.behav_SexualMixingMSM_3544_1824.paramValue(i) ...
                    - CalibParams.behav_SexualMixingMSM_3544_2534.paramValue(i) - CalibParams.behav_SexualMixingMSM_3544_4554.paramValue(i) ...
                    - CalibParams.behav_SexualMixingMSM_3544_5564.paramValue(i) - CalibParams.behav_SexualMixingMSM_3544_65.paramValue(i);
                store_mix_age_S_MSM(4,5) = CalibParams.behav_SexualMixingMSM_3544_4554.paramValue(i);
                store_mix_age_S_MSM(4,6) = CalibParams.behav_SexualMixingMSM_3544_5564.paramValue(i);
                store_mix_age_S_MSM(4,7) = CalibParams.behav_SexualMixingMSM_3544_65.paramValue(i);
                
                store_mix_age_S_MSM(5,1) = CalibParams.behav_SexualMixingMSM_4554_1317.paramValue(i);
                store_mix_age_S_MSM(5,2) = CalibParams.behav_SexualMixingMSM_4554_1824.paramValue(i);
                store_mix_age_S_MSM(5,3) = CalibParams.behav_SexualMixingMSM_4554_2534.paramValue(i);
                store_mix_age_S_MSM(5,4) = CalibParams.behav_SexualMixingMSM_4554_3544.paramValue(i);
                store_mix_age_S_MSM(5,5) = 1 - CalibParams.behav_SexualMixingMSM_4554_1317.paramValue(i) - CalibParams.behav_SexualMixingMSM_4554_1824.paramValue(i) ...
                    - CalibParams.behav_SexualMixingMSM_4554_2534.paramValue(i) - CalibParams.behav_SexualMixingMSM_4554_3544.paramValue(i) ...
                    - CalibParams.behav_SexualMixingMSM_4554_5564.paramValue(i) - CalibParams.behav_SexualMixingMSM_4554_65.paramValue(i);
                store_mix_age_S_MSM(5,6) = CalibParams.behav_SexualMixingMSM_4554_5564.paramValue(i);
                store_mix_age_S_MSM(5,7) = CalibParams.behav_SexualMixingMSM_4554_65.paramValue(i);
                
                store_mix_age_S_MSM(6,1) = CalibParams.behav_SexualMixingMSM_5564_1317.paramValue(i);
                store_mix_age_S_MSM(6,2) = CalibParams.behav_SexualMixingMSM_5564_1824.paramValue(i);
                store_mix_age_S_MSM(6,3) = CalibParams.behav_SexualMixingMSM_5564_2534.paramValue(i);
                store_mix_age_S_MSM(6,4) = CalibParams.behav_SexualMixingMSM_5564_3544.paramValue(i);
                store_mix_age_S_MSM(6,5) = CalibParams.behav_SexualMixingMSM_5564_4554.paramValue(i);
                store_mix_age_S_MSM(6,6) = 1 - CalibParams.behav_SexualMixingMSM_5564_1317.paramValue(i) - CalibParams.behav_SexualMixingMSM_5564_1824.paramValue(i) ...
                    - CalibParams.behav_SexualMixingMSM_5564_2534.paramValue(i) - CalibParams.behav_SexualMixingMSM_5564_3544.paramValue(i) ...
                    - CalibParams.behav_SexualMixingMSM_5564_4554.paramValue(i) - CalibParams.behav_SexualMixingMSM_5564_65.paramValue(i);
                store_mix_age_S_MSM(6,7) = CalibParams.behav_SexualMixingMSM_5564_65.paramValue(i);
                
                store_mix_age_S_MSM(7,1) = CalibParams.behav_SexualMixingMSM_65_1317.paramValue(i);
                store_mix_age_S_MSM(7,2) = CalibParams.behav_SexualMixingMSM_65_1824.paramValue(i);
                store_mix_age_S_MSM(7,3) = CalibParams.behav_SexualMixingMSM_65_2534.paramValue(i);
                store_mix_age_S_MSM(7,4) = CalibParams.behav_SexualMixingMSM_65_3544.paramValue(i);
                store_mix_age_S_MSM(7,5) = CalibParams.behav_SexualMixingMSM_65_4554.paramValue(i);
                store_mix_age_S_MSM(7,6) = CalibParams.behav_SexualMixingMSM_65_5564.paramValue(i);
                store_mix_age_S_MSM(7,7) = 1 - CalibParams.behav_SexualMixingMSM_65_1317.paramValue(i) - CalibParams.behav_SexualMixingMSM_65_1824.paramValue(i) ...
                    - CalibParams.behav_SexualMixingMSM_65_2534.paramValue(i) - CalibParams.behav_SexualMixingMSM_65_3544.paramValue(i) ...
                    - CalibParams.behav_SexualMixingMSM_65_4554.paramValue(i) - CalibParams.behav_SexualMixingMSM_65_5564.paramValue(i);
                
                Params.behav_mix_age_S_MSM = ...
                      Params.ageIndicator ...
                    * store_mix_age_S_MSM ...
                    * Params.ageIndicator';  

                % Zero out HET / IDU Strat
                Params.behav_mix_age_S_MSM = Params.behav_mix_age_S_MSM .* Params.MSM_indicator;

                % Combined age mixing
                Params.behav_mix_age_S = Params.behav_mix_age_S_HET + Params.behav_mix_age_S_MSM;
       

    % Combine entire Sexual Mixing Matrix 
        % (without the race-specific circumcision factor applied)
        %And without MSM-specific matrix
        
        behav_mixing_store_S = ...
               mixMatrix_riskLevel .* mixMatrix_popSex ...
            .* mixMatrix_race .* Params.behav_mix_age_S;

    % Apply the circumcision factor
        % Loop through all the rows of the mixing matrix to apply
        % race-specific circ
        
        % Note: circ/uncirc multipliers are intialized in
        % InitializeParameters
        
        for j = 1:Params.numStrats
            Params.behav_SexualMixing(j,:) = ...
                Params.maleUncircMultiplier' .* Params.maleCircMultiplier' .*behav_mixing_store_S(j,:);
        end

        
    % Check: sum(Params.behav_SexualMixing,2) = 1 for all rows

%% 5. Update the Needle Mixing Matrix

    % Calculate values dependent on calibrated values
        % Scalars
        
        behav_NeedleMixing_IDUM_IDUM = 1 - ...
                CalibParams.behav_NeedleMixing_IDUM_IDUF.paramValue(i);

        behav_NeedleMixing_IDUF_IDUF = 1 - ...
                CalibParams.behav_NeedleMixing_IDUF_IDUM.paramValue(i);
    
    % TRANSMISSION GROUP / SEX MIXING        

        % Applying individual factors to a 273x273 matrix
        mixNeedle_IDUM_IDUM = Params.idx_IDUM_IDUM * behav_NeedleMixing_IDUM_IDUM;
        mixNeedle_IDUM_IDUF = Params.idx_IDUM_IDUF * CalibParams.behav_NeedleMixing_IDUM_IDUF.paramValue(i);
        mixNeedle_IDUF_IDUM = Params.idx_IDUF_IDUM * CalibParams.behav_NeedleMixing_IDUF_IDUM.paramValue(i);
        mixNeedle_IDUF_IDUF = Params.idx_IDUF_IDUF * behav_NeedleMixing_IDUF_IDUF;
    
        % Combine to one single 273x273 matrix
        needleMix_SexRiskType = ...
              mixNeedle_IDUM_IDUM + mixNeedle_IDUM_IDUF ...
            + mixNeedle_IDUF_IDUM + mixNeedle_IDUF_IDUF;
    
    % Calculate full mixing matrix 
        % (without circumcision factor applied)
    behav_mixing_store_N = needleMix_SexRiskType .* ...
        Params.behav_mix_race_N .* Params.behav_mix_age_N;
    
    % Apply circumcision facto
    for stratcount = 1:Params.numStrats
        Params.behav_NeedleMixing(stratcount,:)= ...
            Params.maleUncircMultiplier' .* Params.maleCircMultiplier' .* behav_mixing_store_N(stratcount,:);
    end

    % Check: sum(Params.behav_NeedleMixing,2) = 1 for all rows with IDUs
    % (non-IDU rows should be 0)
    
%% 6. Update the Annual Number of Sexual Contacts (added by JC on 12/05/2017)

    % 6.i. Calculate the annual number of sexual contacts for LR MSM
    for sectionLRMSMContacts = 1:1
        
        store_LRMSMContacts_Black_13 = (Params.store_behav_MSMOverallNumContacts(1) - (Params.MSM_PercentHR * ((Params.MSMInitPopSize_Race(1)/sum(Params.MSMInitPopSize_Race(:)))*CalibParams.behav_HRMSMNumContacts_Black_13.paramValue(i) ...
            + (Params.MSMInitPopSize_Race(2)/sum(Params.MSMInitPopSize_Race(:))) * CalibParams.behav_HRMSMNumContacts_Hispanic_13.paramValue(i) ...
            + (Params.MSMInitPopSize_Race(3)/sum(Params.MSMInitPopSize_Race(:))) * CalibParams.behav_HRMSMNumContacts_Other_13.paramValue(i)))) ...
            / ((1 - Params.MSM_PercentHR)*((Params.MSMInitPopSize_Race(1)/sum(Params.MSMInitPopSize_Race(:))) ...
            + (Params.MSMInitPopSize_Race(2)/sum(Params.MSMInitPopSize_Race(:))) * (CalibParams.behav_HRMSMNumContacts_Hispanic_13.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_13.paramValue(i)) ...
            + (Params.MSMInitPopSize_Race(3)/sum(Params.MSMInitPopSize_Race(:))) * (CalibParams.behav_HRMSMNumContacts_Other_13.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_13.paramValue(i))));
        store_LRMSMContacts_Black_18 = (Params.store_behav_MSMOverallNumContacts(2) - (Params.MSM_PercentHR * ((Params.MSMInitPopSize_Race(1)/sum(Params.MSMInitPopSize_Race(:)))*CalibParams.behav_HRMSMNumContacts_Black_18.paramValue(i) ...
            + (Params.MSMInitPopSize_Race(2)/sum(Params.MSMInitPopSize_Race(:))) * CalibParams.behav_HRMSMNumContacts_Hispanic_18.paramValue(i) ...
            + (Params.MSMInitPopSize_Race(3)/sum(Params.MSMInitPopSize_Race(:))) * CalibParams.behav_HRMSMNumContacts_Other_18.paramValue(i)))) ...
            / ((1 - Params.MSM_PercentHR)*((Params.MSMInitPopSize_Race(1)/sum(Params.MSMInitPopSize_Race(:))) ...
            + (Params.MSMInitPopSize_Race(2)/sum(Params.MSMInitPopSize_Race(:))) * (CalibParams.behav_HRMSMNumContacts_Hispanic_18.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_18.paramValue(i)) ...
            + (Params.MSMInitPopSize_Race(3)/sum(Params.MSMInitPopSize_Race(:))) * (CalibParams.behav_HRMSMNumContacts_Other_18.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_18.paramValue(i))));
        store_LRMSMContacts_Black_25 = (Params.store_behav_MSMOverallNumContacts(3) - (Params.MSM_PercentHR * ((Params.MSMInitPopSize_Race(1)/sum(Params.MSMInitPopSize_Race(:)))*CalibParams.behav_HRMSMNumContacts_Black_25.paramValue(i) ...
            + (Params.MSMInitPopSize_Race(2)/sum(Params.MSMInitPopSize_Race(:))) * CalibParams.behav_HRMSMNumContacts_Hispanic_25.paramValue(i) ...
            + (Params.MSMInitPopSize_Race(3)/sum(Params.MSMInitPopSize_Race(:))) * CalibParams.behav_HRMSMNumContacts_Other_25.paramValue(i)))) ...
            / ((1 - Params.MSM_PercentHR)*((Params.MSMInitPopSize_Race(1)/sum(Params.MSMInitPopSize_Race(:))) ...
            + (Params.MSMInitPopSize_Race(2)/sum(Params.MSMInitPopSize_Race(:))) * (CalibParams.behav_HRMSMNumContacts_Hispanic_25.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_25.paramValue(i)) ...
            + (Params.MSMInitPopSize_Race(3)/sum(Params.MSMInitPopSize_Race(:))) * (CalibParams.behav_HRMSMNumContacts_Other_25.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_25.paramValue(i))));
        store_LRMSMContacts_Black_35 = (Params.store_behav_MSMOverallNumContacts(4) - (Params.MSM_PercentHR * ((Params.MSMInitPopSize_Race(1)/sum(Params.MSMInitPopSize_Race(:)))*CalibParams.behav_HRMSMNumContacts_Black_35.paramValue(i) ...
            + (Params.MSMInitPopSize_Race(2)/sum(Params.MSMInitPopSize_Race(:))) * CalibParams.behav_HRMSMNumContacts_Hispanic_35.paramValue(i) ...
            + (Params.MSMInitPopSize_Race(3)/sum(Params.MSMInitPopSize_Race(:))) * CalibParams.behav_HRMSMNumContacts_Other_35.paramValue(i)))) ...
            / ((1 - Params.MSM_PercentHR)*((Params.MSMInitPopSize_Race(1)/sum(Params.MSMInitPopSize_Race(:))) ...
            + (Params.MSMInitPopSize_Race(2)/sum(Params.MSMInitPopSize_Race(:))) * (CalibParams.behav_HRMSMNumContacts_Hispanic_35.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_35.paramValue(i)) ...
            + (Params.MSMInitPopSize_Race(3)/sum(Params.MSMInitPopSize_Race(:))) * (CalibParams.behav_HRMSMNumContacts_Other_35.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_35.paramValue(i))));
        store_LRMSMContacts_Black_45 = (Params.store_behav_MSMOverallNumContacts(5) - (Params.MSM_PercentHR * ((Params.MSMInitPopSize_Race(1)/sum(Params.MSMInitPopSize_Race(:)))*CalibParams.behav_HRMSMNumContacts_Black_45.paramValue(i) ...
            + (Params.MSMInitPopSize_Race(2)/sum(Params.MSMInitPopSize_Race(:))) * CalibParams.behav_HRMSMNumContacts_Hispanic_45.paramValue(i) ...
            + (Params.MSMInitPopSize_Race(3)/sum(Params.MSMInitPopSize_Race(:))) * CalibParams.behav_HRMSMNumContacts_Other_45.paramValue(i)))) ...
            / ((1 - Params.MSM_PercentHR)*((Params.MSMInitPopSize_Race(1)/sum(Params.MSMInitPopSize_Race(:))) ...
            + (Params.MSMInitPopSize_Race(2)/sum(Params.MSMInitPopSize_Race(:))) * (CalibParams.behav_HRMSMNumContacts_Hispanic_45.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_45.paramValue(i)) ...
            + (Params.MSMInitPopSize_Race(3)/sum(Params.MSMInitPopSize_Race(:))) * (CalibParams.behav_HRMSMNumContacts_Other_45.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_45.paramValue(i))));
        store_LRMSMContacts_Black_55 = (Params.store_behav_MSMOverallNumContacts(6) - (Params.MSM_PercentHR * ((Params.MSMInitPopSize_Race(1)/sum(Params.MSMInitPopSize_Race(:)))*CalibParams.behav_HRMSMNumContacts_Black_55.paramValue(i) ...
            + (Params.MSMInitPopSize_Race(2)/sum(Params.MSMInitPopSize_Race(:))) * CalibParams.behav_HRMSMNumContacts_Hispanic_55.paramValue(i) ...
            + (Params.MSMInitPopSize_Race(3)/sum(Params.MSMInitPopSize_Race(:))) * CalibParams.behav_HRMSMNumContacts_Other_55.paramValue(i)))) ...
            / ((1 - Params.MSM_PercentHR)*((Params.MSMInitPopSize_Race(1)/sum(Params.MSMInitPopSize_Race(:))) ...
            + (Params.MSMInitPopSize_Race(2)/sum(Params.MSMInitPopSize_Race(:))) * (CalibParams.behav_HRMSMNumContacts_Hispanic_55.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_55.paramValue(i)) ...
            + (Params.MSMInitPopSize_Race(3)/sum(Params.MSMInitPopSize_Race(:))) * (CalibParams.behav_HRMSMNumContacts_Other_55.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_55.paramValue(i))));
        store_LRMSMContacts_Black_65 = (Params.store_behav_MSMOverallNumContacts(7) - (Params.MSM_PercentHR * ((Params.MSMInitPopSize_Race(1)/sum(Params.MSMInitPopSize_Race(:)))*CalibParams.behav_HRMSMNumContacts_Black_65.paramValue(i) ...
            + (Params.MSMInitPopSize_Race(2)/sum(Params.MSMInitPopSize_Race(:))) * CalibParams.behav_HRMSMNumContacts_Hispanic_65.paramValue(i) ...
            + (Params.MSMInitPopSize_Race(3)/sum(Params.MSMInitPopSize_Race(:))) * CalibParams.behav_HRMSMNumContacts_Other_65.paramValue(i)))) ...
            / ((1 - Params.MSM_PercentHR)*((Params.MSMInitPopSize_Race(1)/sum(Params.MSMInitPopSize_Race(:))) ...
            + (Params.MSMInitPopSize_Race(2)/sum(Params.MSMInitPopSize_Race(:))) * (CalibParams.behav_HRMSMNumContacts_Hispanic_65.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_65.paramValue(i)) ...
            + (Params.MSMInitPopSize_Race(3)/sum(Params.MSMInitPopSize_Race(:))) * (CalibParams.behav_HRMSMNumContacts_Other_65.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_65.paramValue(i))));
        
        store_LRMSMContacts_Hispanic_13 = store_LRMSMContacts_Black_13 * (CalibParams.behav_HRMSMNumContacts_Hispanic_13.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_13.paramValue(i));
        store_LRMSMContacts_Hispanic_18 = store_LRMSMContacts_Black_18 * (CalibParams.behav_HRMSMNumContacts_Hispanic_18.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_18.paramValue(i));
        store_LRMSMContacts_Hispanic_25 = store_LRMSMContacts_Black_25 * (CalibParams.behav_HRMSMNumContacts_Hispanic_25.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_25.paramValue(i));
        store_LRMSMContacts_Hispanic_35 = store_LRMSMContacts_Black_35 * (CalibParams.behav_HRMSMNumContacts_Hispanic_35.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_35.paramValue(i));
        store_LRMSMContacts_Hispanic_45 = store_LRMSMContacts_Black_45 * (CalibParams.behav_HRMSMNumContacts_Hispanic_45.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_45.paramValue(i));
        store_LRMSMContacts_Hispanic_55 = store_LRMSMContacts_Black_55 * (CalibParams.behav_HRMSMNumContacts_Hispanic_55.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_55.paramValue(i));
        store_LRMSMContacts_Hispanic_65 = store_LRMSMContacts_Black_65 * (CalibParams.behav_HRMSMNumContacts_Hispanic_65.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_65.paramValue(i));
        
        store_LRMSMContacts_Other_13 = store_LRMSMContacts_Black_13 * (CalibParams.behav_HRMSMNumContacts_Other_13.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_13.paramValue(i));
        store_LRMSMContacts_Other_18 = store_LRMSMContacts_Black_18 * (CalibParams.behav_HRMSMNumContacts_Other_18.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_18.paramValue(i));
        store_LRMSMContacts_Other_25 = store_LRMSMContacts_Black_25 * (CalibParams.behav_HRMSMNumContacts_Other_25.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_25.paramValue(i));
        store_LRMSMContacts_Other_35 = store_LRMSMContacts_Black_35 * (CalibParams.behav_HRMSMNumContacts_Other_35.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_35.paramValue(i));
        store_LRMSMContacts_Other_45 = store_LRMSMContacts_Black_45 * (CalibParams.behav_HRMSMNumContacts_Other_45.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_45.paramValue(i));
        store_LRMSMContacts_Other_55 = store_LRMSMContacts_Black_55 * (CalibParams.behav_HRMSMNumContacts_Other_55.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_55.paramValue(i));
        store_LRMSMContacts_Other_65 = store_LRMSMContacts_Black_65 * (CalibParams.behav_HRMSMNumContacts_Other_65.paramValue(i)/CalibParams.behav_HRMSMNumContacts_Black_65.paramValue(i));
        
    end
    
    % 6.ii. Calculate the annual # of contacts per partner
    for sectionMSMContactsPerPartner = 1:1
        store_MSMLow_Black(1,1) = store_LRMSMContacts_Black_13 / Params.behav_partnersSex_MSM(1,1);
        store_MSMLow_Black(2,1) = store_LRMSMContacts_Black_18 / Params.behav_partnersSex_MSM(1,1);
        store_MSMLow_Black(3,1) = store_LRMSMContacts_Black_25 / Params.behav_partnersSex_MSM(1,1);
        store_MSMLow_Black(4,1) = store_LRMSMContacts_Black_35 / Params.behav_partnersSex_MSM(1,1);
        store_MSMLow_Black(5,1) = store_LRMSMContacts_Black_45 / Params.behav_partnersSex_MSM(1,1);
        store_MSMLow_Black(6,1) = store_LRMSMContacts_Black_55 / Params.behav_partnersSex_MSM(1,1);
        store_MSMLow_Black(7,1) = store_LRMSMContacts_Black_65 / Params.behav_partnersSex_MSM(1,1);
        
        store_MSMHigh_Black(1,1) = CalibParams.behav_HRMSMNumContacts_Black_13.paramValue(i) / Params.behav_partnersSex_MSM(1,2);
        store_MSMHigh_Black(2,1) = CalibParams.behav_HRMSMNumContacts_Black_18.paramValue(i) / Params.behav_partnersSex_MSM(1,2);
        store_MSMHigh_Black(3,1) = CalibParams.behav_HRMSMNumContacts_Black_25.paramValue(i) / Params.behav_partnersSex_MSM(1,2);
        store_MSMHigh_Black(4,1) = CalibParams.behav_HRMSMNumContacts_Black_35.paramValue(i) / Params.behav_partnersSex_MSM(1,2);
        store_MSMHigh_Black(5,1) = CalibParams.behav_HRMSMNumContacts_Black_45.paramValue(i) / Params.behav_partnersSex_MSM(1,2);
        store_MSMHigh_Black(6,1) = CalibParams.behav_HRMSMNumContacts_Black_55.paramValue(i) / Params.behav_partnersSex_MSM(1,2);
        store_MSMHigh_Black(7,1) = CalibParams.behav_HRMSMNumContacts_Black_65.paramValue(i) / Params.behav_partnersSex_MSM(1,2);
        
        store_MSMLow_Hisp(1,1) = store_LRMSMContacts_Hispanic_13 / Params.behav_partnersSex_MSM(2,1);
        store_MSMLow_Hisp(2,1) = store_LRMSMContacts_Hispanic_18 / Params.behav_partnersSex_MSM(2,1);
        store_MSMLow_Hisp(3,1) = store_LRMSMContacts_Hispanic_25 / Params.behav_partnersSex_MSM(2,1);
        store_MSMLow_Hisp(4,1) = store_LRMSMContacts_Hispanic_35 / Params.behav_partnersSex_MSM(2,1);
        store_MSMLow_Hisp(5,1) = store_LRMSMContacts_Hispanic_45 / Params.behav_partnersSex_MSM(2,1);
        store_MSMLow_Hisp(6,1) = store_LRMSMContacts_Hispanic_55 / Params.behav_partnersSex_MSM(2,1);
        store_MSMLow_Hisp(7,1) = store_LRMSMContacts_Hispanic_65 / Params.behav_partnersSex_MSM(2,1);
        
        store_MSMHigh_Hisp(1,1) = CalibParams.behav_HRMSMNumContacts_Hispanic_13.paramValue(i) / Params.behav_partnersSex_MSM(2,2);
        store_MSMHigh_Hisp(2,1) = CalibParams.behav_HRMSMNumContacts_Hispanic_18.paramValue(i) / Params.behav_partnersSex_MSM(2,2);
        store_MSMHigh_Hisp(3,1) = CalibParams.behav_HRMSMNumContacts_Hispanic_25.paramValue(i) / Params.behav_partnersSex_MSM(2,2);
        store_MSMHigh_Hisp(4,1) = CalibParams.behav_HRMSMNumContacts_Hispanic_35.paramValue(i) / Params.behav_partnersSex_MSM(2,2);
        store_MSMHigh_Hisp(5,1) = CalibParams.behav_HRMSMNumContacts_Hispanic_45.paramValue(i) / Params.behav_partnersSex_MSM(2,2);
        store_MSMHigh_Hisp(6,1) = CalibParams.behav_HRMSMNumContacts_Hispanic_55.paramValue(i) / Params.behav_partnersSex_MSM(2,2);
        store_MSMHigh_Hisp(7,1) = CalibParams.behav_HRMSMNumContacts_Hispanic_65.paramValue(i) / Params.behav_partnersSex_MSM(2,2);
        
        store_MSMLow_Other(1,1) = store_LRMSMContacts_Other_13 / Params.behav_partnersSex_MSM(3,1);
        store_MSMLow_Other(2,1) = store_LRMSMContacts_Other_18 / Params.behav_partnersSex_MSM(3,1);
        store_MSMLow_Other(3,1) = store_LRMSMContacts_Other_25 / Params.behav_partnersSex_MSM(3,1);
        store_MSMLow_Other(4,1) = store_LRMSMContacts_Other_35 / Params.behav_partnersSex_MSM(3,1);
        store_MSMLow_Other(5,1) = store_LRMSMContacts_Other_45 / Params.behav_partnersSex_MSM(3,1);
        store_MSMLow_Other(6,1) = store_LRMSMContacts_Other_55 / Params.behav_partnersSex_MSM(3,1);
        store_MSMLow_Other(7,1) = store_LRMSMContacts_Other_65 / Params.behav_partnersSex_MSM(3,1);
        
        store_MSMHigh_Other(1,1) = CalibParams.behav_HRMSMNumContacts_Other_13.paramValue(i) / Params.behav_partnersSex_MSM(3,2);
        store_MSMHigh_Other(2,1) = CalibParams.behav_HRMSMNumContacts_Other_18.paramValue(i) / Params.behav_partnersSex_MSM(3,2);
        store_MSMHigh_Other(3,1) = CalibParams.behav_HRMSMNumContacts_Other_25.paramValue(i) / Params.behav_partnersSex_MSM(3,2);
        store_MSMHigh_Other(4,1) = CalibParams.behav_HRMSMNumContacts_Other_35.paramValue(i) / Params.behav_partnersSex_MSM(3,2);
        store_MSMHigh_Other(5,1) = CalibParams.behav_HRMSMNumContacts_Other_45.paramValue(i) / Params.behav_partnersSex_MSM(3,2);
        store_MSMHigh_Other(6,1) = CalibParams.behav_HRMSMNumContacts_Other_55.paramValue(i) / Params.behav_partnersSex_MSM(3,2);
        store_MSMHigh_Other(7,1) = CalibParams.behav_HRMSMNumContacts_Other_65.paramValue(i) / Params.behav_partnersSex_MSM(3,2);
        
    end
    
    % 6.iii. Calculate the annual # of contacts per partner
    for sectionMSMContactsPerPartner = 1:1
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
                factor_adjNumSexContact_NoPrEP(:,:,Params.AwareComparts) = ...
                    factor_adjNumSexContact_NoPrEP(:,:,Params.AwareComparts) * (1 - Params.hiv_reductSexActsDueToAwareness);

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
        
    end
    
%% 7. Update Duration in Age Group Params

    % durAge_adjusted takes multiplier from init params or opti algorithm
    % and multiplies by the base duration in age group. Then sets Params
    % variable as it is in `InitializeParameters'
    durAge_adjusted = zeros(7,1);
    durAge_adjusted(1) = Params.durAge_normal(1) * CalibParams.durInAgeMultiplier_1317.paramValue(i);
    durAge_adjusted(2) = Params.durAge_normal(2) * CalibParams.durInAgeMultiplier_1824.paramValue(i);
    durAge_adjusted(3) = Params.durAge_normal(3) * CalibParams.durInAgeMultiplier_2534.paramValue(i);
    durAge_adjusted(4) = Params.durAge_normal(4) * CalibParams.durInAgeMultiplier_3544.paramValue(i);
    durAge_adjusted(5) = Params.durAge_normal(5) * CalibParams.durInAgeMultiplier_4554.paramValue(i);
    durAge_adjusted(6) = Params.durAge_normal(6) * CalibParams.durInAgeMultiplier_5564.paramValue(i);
    durAge_adjusted(7) = Params.durAge_normal(7);

    Params.hiv_durAge_a = Params.ageIndicator * durAge_adjusted;
    
%% 8. Update the AI Parameters

    % 8.i. Update percent partnerships
    for sectionPctPartnerships = 1:1
       
    % Update pct partnerships V and A (modified 08/03/2017)                 
        
        store_pctPartnershipsAIfSomeA_1(1,1) = CalibParams.pctPartnershipsVandA_1_13.paramValue(i);
        store_pctPartnershipsAIfSomeA_1(2,1) = CalibParams.pctPartnershipsVandA_1_18.paramValue(i);
        store_pctPartnershipsAIfSomeA_1(3,1) = CalibParams.pctPartnershipsVandA_1_25.paramValue(i);
        store_pctPartnershipsAIfSomeA_1(4,1) = CalibParams.pctPartnershipsVandA_1_35.paramValue(i);
        store_pctPartnershipsAIfSomeA_1(5,1) = CalibParams.pctPartnershipsVandA_1_45.paramValue(i);
        store_pctPartnershipsAIfSomeA_1(6,1) = CalibParams.pctPartnershipsVandA_1_55.paramValue(i);
        store_pctPartnershipsAIfSomeA_1(7,1) = CalibParams.pctPartnershipsVandA_1_65.paramValue(i);
        
        store_pctPartnershipsAIfSomeA_2to4(1,1) = CalibParams.pctPartnershipsVandA_2to4_13.paramValue(i);
        store_pctPartnershipsAIfSomeA_2to4(2,1) = CalibParams.pctPartnershipsVandA_2to4_18.paramValue(i);
        store_pctPartnershipsAIfSomeA_2to4(3,1) = CalibParams.pctPartnershipsVandA_2to4_25.paramValue(i);
        store_pctPartnershipsAIfSomeA_2to4(4,1) = CalibParams.pctPartnershipsVandA_2to4_35.paramValue(i);
        store_pctPartnershipsAIfSomeA_2to4(5,1) = CalibParams.pctPartnershipsVandA_2to4_45.paramValue(i);
        store_pctPartnershipsAIfSomeA_2to4(6,1) = CalibParams.pctPartnershipsVandA_2to4_55.paramValue(i);
        store_pctPartnershipsAIfSomeA_2to4(7,1) = CalibParams.pctPartnershipsVandA_2to4_65.paramValue(i);
        
        % Apply to 273x1
        pctMFPartnersA_IfSomeAInMF_1 = Params.ageIndicator * store_pctPartnershipsAIfSomeA_1;
        pctMFPartnersA_IfSomeAInMF_2to4 = Params.ageIndicator * store_pctPartnershipsAIfSomeA_2to4;
        
        pctMFPartnersA_IfSomeAInMF_1_2D(Params.numStrats,Params.numStrats) = 0;
        pctMFPartnersA_IfSomeAInMF_2to4_2D(Params.numStrats,Params.numStrats) = 0;
        % Figure out where in the mixing matrix people don't mix with each other
        Params.idx_doesnotmix = Params.behav_SexualMixing == 0;
        % Initially - everyone can mix together
        Params.allMixCombos = ones(Params.numStrats, Params.numStrats);
        % Remove all the people who don't mix with each other
        Params.allMixCombos(Params.idx_doesnotmix) = 0;
        allNonMSM_MSM_MixCombos = Params.allMixCombos;
        allNonMSM_MSM_MixCombos(Params.idx_MSMandMSM)=0;
        Params.idx_allMaleAndFemaleMixing = (allNonMSM_MSM_MixCombos == 1);
        
        % Loop over partner cohorts
        for s = 1:Params.numStrats     
           pctMFPartnersA_IfSomeAInMF_1_2D(:,s) = allNonMSM_MSM_MixCombos(:,s) .* ...
               pctMFPartnersA_IfSomeAInMF_1;
           pctMFPartnersA_IfSomeAInMF_2to4_2D(:,s) = allNonMSM_MSM_MixCombos(:,s) .* ...
               pctMFPartnersA_IfSomeAInMF_2to4;                  
        end
        
        % Ceate 3D matrix
        pctMFPartnersA_IfSomeAInMF_1_3D = repmat(pctMFPartnersA_IfSomeAInMF_1_2D,[1,1,Params.numComparts]);
        pctMFPartnersA_IfSomeAInMF_2to4_3D = repmat(pctMFPartnersA_IfSomeAInMF_2to4_2D,[1,1,Params.numComparts]);
        
        % Re-arrange matrix for Params struct to be used in CalcInfectionRates (30x273x273)
        Params.pctMFPartnersA_IfSomeAInMF_1 = permute(pctMFPartnersA_IfSomeAInMF_1_3D,[3,1,2]);
        Params.pctMFPartnersA_IfSomeAInMF_2to4 = permute(pctMFPartnersA_IfSomeAInMF_2to4_3D,[3,1,2]);
        
        

    end
    
    % 8.ii. Update percent of contacts that are VI (this needs to be
    % recalculated since AllMixCombos variable is updated - changed on July
    % 3, 2019)
    for sectionPctV = 1:1
      
        % New way to update pctVaginal (updated 22May2014 based on new
        % stratifications)
       
            % Percent of contacts that are anal
                % Given individual is in a partnership with both V and A
            %store_AgeMale(1:Params.numAge,1)= [...
            %    CalibParams.factor_pctContactsA_VandA_M_13.paramValue(i), ...
            %    CalibParams.factor_pctContactsA_VandA_M_18.paramValue(i), ...
            %    CalibParams.factor_pctContactsA_VandA_M_25.paramValue(i), ...
            %    CalibParams.factor_pctContactsA_VandA_M_35.paramValue(i), ...
            %    CalibParams.factor_pctContactsA_VandA_M_45.paramValue(i)];
            
            %store_AgeFemale(1:Params.numAge,1) = [...
            %    CalibParams.factor_pctContactsA_VandA_F_13.paramValue(i), ...
            %    CalibParams.factor_pctContactsA_VandA_F_18.paramValue(i), ...
            %    CalibParams.factor_pctContactsA_VandA_F_25.paramValue(i), ...
            %    CalibParams.factor_pctContactsA_VandA_F_35.paramValue(i), ...
            %    CalibParams.factor_pctContactsA_VandA_F_45.paramValue(i)];
            
            % Populate a [273x1] vector of pct contacts anal
            %    pctContactsAnal = ...
            %           Params.ageIndicator * store_AgeMale .* Params.sexIndicator(:,Params.sex_Male) ...
            %        +  Params.ageIndicator * store_AgeFemale .* Params.sexIndicator(:,Params.sex_Female);
            
            % Calculate percent of contacts that are vaginal
            %    pctContactsVaginal = 1 - pctContactsAnal;
            
        % Loop through each period
        for RatePeriod = 1:(Params.numRateInputPeriod - 1)
            % Apply the pctContactsVaginal and PctContactsAnal from InitParams            
            switch RatePeriod
                case 1
                    pctContactsAnal = Params.pctContactsAnal_1;
                    pctContactsVaginal = Params.pctContactsVaginal_1;
                case 2
                    pctContactsAnal = Params.pctContactsAnal_2to4;
                    pctContactsVaginal = Params.pctContactsVaginal_2to4;                             
            end
            
            
            %Preallocate
            pctVaginal_2D(Params.numStrats,Params.numStrats) = 0;
            pctAnal_2D = zeros(Params.numStrats,Params.numStrats);
                
                % Loop over all possible partner cohorts
                for s = 1:Params.numStrats 
                    pctVaginal_2D(:,s) = allNonMSM_MSM_MixCombos(:,s) .* ...
                        pctContactsVaginal;
                    pctAnal_2D(:,s) = allNonMSM_MSM_MixCombos(:,s) .* ...
                        pctContactsAnal;
                end

                % Expand to a 3d matrix so that it can be applied to the
                % betas (apply the 2D mixing matrix to every Compartment)
                
                pctVaginal_3D = repmat(pctVaginal_2D,[1,1,Params.numComparts]);
                pctAnal_3D = repmat(pctAnal_2D,[1,1,Params.numComparts]);
                switch RatePeriod
                    case 1
                        BetaFactor.pctContactsVaginal_InMFWithA_1 = pctVaginal_3D;
                        BetaFactor.pctContactsAnal_InMFWithA_1 = pctAnal_3D;
                    case 2      
                        BetaFactor.pctContactsVaginal_InMFWithA_2to4 = pctVaginal_3D;
                        BetaFactor.pctContactsAnal_InMFWithA_2to4 = pctAnal_3D;                    
                end
            
            % Read in beta factors established in InitializeParams (they
            % are not saved when that m file is run)
            %BetaFactor.pctContactsVaginal_InMFWithA_1 = Params.pctContactsVaginal_InMFWithA_1;
            %BetaFactor.pctContactsAnal_InMFWithA_1 = Params.pctContactsAnal_InMFWithA_1;
            %BetaFactor.pctContactsVaginal_InMFWithA_2 = Params.pctContactsVaginal_InMFWithA_2;
            %BetaFactor.pctContactsAnal_InMFWithA_2 = Params.pctContactsAnal_InMFWithA_2;
        end
    end
    
    % 8.iii. Update percent of MSM contacts that are insertive
    for sectionPctInsertive = 1:1
    
        store_behav_pctMSMRec_Black(1,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Black_13.paramValue(i);
        store_behav_pctMSMRec_Black(2,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Black_18.paramValue(i);
        store_behav_pctMSMRec_Black(3,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Black_25.paramValue(i);
        store_behav_pctMSMRec_Black(4,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Black_35.paramValue(i);
        store_behav_pctMSMRec_Black(5,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Black_45.paramValue(i);
        store_behav_pctMSMRec_Black(6,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Black_55.paramValue(i);
        store_behav_pctMSMRec_Black(7,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Black_65.paramValue(i);
        
        store_behav_pctMSMRec_Hisp(1,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Hispanic_13.paramValue(i);
        store_behav_pctMSMRec_Hisp(2,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Hispanic_18.paramValue(i);
        store_behav_pctMSMRec_Hisp(3,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Hispanic_25.paramValue(i);
        store_behav_pctMSMRec_Hisp(4,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Hispanic_35.paramValue(i);
        store_behav_pctMSMRec_Hisp(5,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Hispanic_45.paramValue(i);
        store_behav_pctMSMRec_Hisp(6,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Hispanic_55.paramValue(i);
        store_behav_pctMSMRec_Hisp(7,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Hispanic_65.paramValue(i);
        
        store_behav_pctMSMRec_Other(1,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Other_13.paramValue(i);
        store_behav_pctMSMRec_Other(2,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Other_18.paramValue(i);
        store_behav_pctMSMRec_Other(3,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Other_25.paramValue(i);
        store_behav_pctMSMRec_Other(4,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Other_35.paramValue(i);
        store_behav_pctMSMRec_Other(5,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Other_45.paramValue(i);
        store_behav_pctMSMRec_Other(6,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Other_55.paramValue(i);
        store_behav_pctMSMRec_Other(7,1) = 1 - CalibParams.behav_MSMpctContacts_Insertive_Other_65.paramValue(i);
        
        
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
    
    end
    % 8.iv. Update percent of MSM contacts that are insertive
    for sectionMSMPctCondomUse = 1:1
        
        store_basePctCondomA_MSM_Black(1,1) = CalibParams.behav_MSMpctCondomUse_Black_13.paramValue(i);
        store_basePctCondomA_MSM_Black(2,1) = CalibParams.behav_MSMpctCondomUse_Black_18.paramValue(i);
        store_basePctCondomA_MSM_Black(3,1) = CalibParams.behav_MSMpctCondomUse_Black_25.paramValue(i);
        store_basePctCondomA_MSM_Black(4,1) = CalibParams.behav_MSMpctCondomUse_Black_35.paramValue(i);
        store_basePctCondomA_MSM_Black(5,1) = CalibParams.behav_MSMpctCondomUse_Black_45.paramValue(i);
        store_basePctCondomA_MSM_Black(6,1) = CalibParams.behav_MSMpctCondomUse_Black_55.paramValue(i);
        store_basePctCondomA_MSM_Black(7,1) = CalibParams.behav_MSMpctCondomUse_Black_65.paramValue(i);
        
        store_basePctCondomA_MSM_Hisp(1,1) = CalibParams.behav_MSMpctCondomUse_Hispanic_13.paramValue(i);
        store_basePctCondomA_MSM_Hisp(2,1) = CalibParams.behav_MSMpctCondomUse_Hispanic_18.paramValue(i);
        store_basePctCondomA_MSM_Hisp(3,1) = CalibParams.behav_MSMpctCondomUse_Hispanic_25.paramValue(i);
        store_basePctCondomA_MSM_Hisp(4,1) = CalibParams.behav_MSMpctCondomUse_Hispanic_35.paramValue(i);
        store_basePctCondomA_MSM_Hisp(5,1) = CalibParams.behav_MSMpctCondomUse_Hispanic_45.paramValue(i);
        store_basePctCondomA_MSM_Hisp(6,1) = CalibParams.behav_MSMpctCondomUse_Hispanic_55.paramValue(i);
        store_basePctCondomA_MSM_Hisp(7,1) = CalibParams.behav_MSMpctCondomUse_Hispanic_65.paramValue(i);
        
        store_basePctCondomA_MSM_Other(1,1) = CalibParams.behav_MSMpctCondomUse_Other_13.paramValue(i);
        store_basePctCondomA_MSM_Other(2,1) = CalibParams.behav_MSMpctCondomUse_Other_18.paramValue(i);
        store_basePctCondomA_MSM_Other(3,1) = CalibParams.behav_MSMpctCondomUse_Other_25.paramValue(i);
        store_basePctCondomA_MSM_Other(4,1) = CalibParams.behav_MSMpctCondomUse_Other_35.paramValue(i);
        store_basePctCondomA_MSM_Other(5,1) = CalibParams.behav_MSMpctCondomUse_Other_45.paramValue(i);
        store_basePctCondomA_MSM_Other(6,1) = CalibParams.behav_MSMpctCondomUse_Other_55.paramValue(i);
        store_basePctCondomA_MSM_Other(7,1) = CalibParams.behav_MSMpctCondomUse_Other_65.paramValue(i);
        
        % Apply M-F and F-M to [273x273]
            % Male-Female and Female-Male AI condom use
            % Assumes every contact is between 2 people of opposite sexes
        
        storeHET_basePctCondomA = Params.storeHET_basePctCondomA;
        storeMSM_basePctCondomA = ((Params.ageIndicator * store_basePctCondomA_MSM_Black ...
            .* Params.raceIndicator(:,Params.race_B) ...
            + Params.ageIndicator * store_basePctCondomA_MSM_Hisp ...
            .* Params.raceIndicator(:,Params.race_H) ...
            + Params.ageIndicator * store_basePctCondomA_MSM_Other ...
            .* Params.raceIndicator(:,Params.race_O))...
            .* Params.popIndicator(:, Params.pop_MSM))...
            * transpose(Params.popIndicator(:, Params.pop_MSM));
        
        
        % Apply value for MSM-MSM partnerships [273x273]
        store3_basePctCondomA = storeHET_basePctCondomA + storeMSM_basePctCondomA;

        % Expand to [273x273x30]
        behav_basePctCondomA = repmat(store3_basePctCondomA,[1,1,Params.numComparts]);
        
        
        % PERCENT OF CONTACTS PROTECTED WITH A CONDOM (ADJUSTED BASED ON AWARENESS AND PrEP)
        
        % AI Contacts
        
            % set aware compartments variable
            inf_awareComparts_273x273x30=zeros(Params.numStrats,Params.numStrats,Params.numComparts);
            inf_awareComparts_273x273x30(:,:,Params.AwareComparts)=1;
        
            % Apply increase due to awareness
            store_factor_adjPctCondomUseA_NoPrEP_1to4 = behav_basePctCondomA .* ...
                (1 + inf_awareComparts_273x273x30 * Params.behav_incrCondomUseDueToDiag_1to4);
            store_factor_adjPctCondomUseA_NoPrEP_5 = behav_basePctCondomA .* ...
                (1 + inf_awareComparts_273x273x30 * Params.behav_incrCondomUseDueToDiag_5);

            % Apply decrease due to PrEP
            store_factor_adjPctCondomUseA_PrEP_1to4 = ...
                store_factor_adjPctCondomUseA_NoPrEP_1to4 * (1 - Params.behav_decrCondomUseDueToPrEP);
            store_factor_adjPctCondomUseA_PrEP_5 = ...
                store_factor_adjPctCondomUseA_NoPrEP_5 * (1 - Params.behav_decrCondomUseDueToPrEP);
        
            % Limit maximum pct condom usage to 1 (these are used in
            % CalcBetas later)
            Params.factor_adjPctCondomUseA_NoPrEP_1to4 = min(store_factor_adjPctCondomUseA_NoPrEP_1to4, 1);
            Params.factor_adjPctCondomUseA_PrEP_1to4 = min(store_factor_adjPctCondomUseA_PrEP_1to4, 1);
            Params.factor_adjPctCondomUseA_NoPrEP_5 = min(store_factor_adjPctCondomUseA_NoPrEP_5, 1);
            Params.factor_adjPctCondomUseA_PrEP_5 = min(store_factor_adjPctCondomUseA_PrEP_5, 1);
        
        
        
    end

%% 9. Update the Factors That Go Into Betas
   
    % Vaginal Insertive
        % Pre-allocate
        inf_basePerActProb_VaginalIns = zeros(Params.numStrats, Params.numStrats, Params.numComparts);

        inf_basePerActProb_VaginalIns(:,:,Params.AcuteComparts)=CalibParams.inf_vaginalInsRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_Acute);
        inf_basePerActProb_VaginalIns(:,:,Params.LatentAComparts)=CalibParams.inf_vaginalInsRisk.paramValue(i)  * Params.store_relRiskSexTransmission(Params.stage_LatentA);
        inf_basePerActProb_VaginalIns(:,:,Params.LatentBComparts)=CalibParams.inf_vaginalInsRisk.paramValue(i)  * Params.store_relRiskSexTransmission(Params.stage_LatentB);
        inf_basePerActProb_VaginalIns(:,:,Params.LateComparts)=CalibParams.inf_vaginalInsRisk.paramValue(i)  * Params.store_relRiskSexTransmission(Params.stage_Late);
        inf_basePerActProb_VaginalIns(:,:,Params.AIDSComparts)=CalibParams.inf_vaginalInsRisk.paramValue(i)  * Params.store_relRiskSexTransmission(Params.stage_AIDS);

        
    % Anal Insertive
        % Pre-allocate
        inf_basePerActProb_AnalIns = zeros(Params.numStrats, Params.numStrats, Params.numComparts);


        inf_basePerActProb_AnalIns(:,:,Params.AcuteComparts)=CalibParams.inf_analInsRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_Acute);
        inf_basePerActProb_AnalIns(:,:,Params.LatentAComparts)=CalibParams.inf_analInsRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_LatentA);
        inf_basePerActProb_AnalIns(:,:,Params.LatentBComparts)=CalibParams.inf_analInsRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_LatentB);
        inf_basePerActProb_AnalIns(:,:,Params.LateComparts)=CalibParams.inf_analInsRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_Late);
        inf_basePerActProb_AnalIns(:,:,Params.AIDSComparts)=CalibParams.inf_analInsRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_AIDS);

        
    % Vaginal Receptive
        % Pre-allocate
        inf_basePerActProb_VaginalRec = zeros(Params.numStrats, Params.numStrats, Params.numComparts);


        inf_basePerActProb_VaginalRec(:,:,Params.AcuteComparts)=CalibParams.inf_vaginalRecRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_Acute);
        inf_basePerActProb_VaginalRec(:,:,Params.LatentAComparts)=CalibParams.inf_vaginalRecRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_LatentA);
        inf_basePerActProb_VaginalRec(:,:,Params.LatentBComparts)=CalibParams.inf_vaginalRecRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_LatentB);
        inf_basePerActProb_VaginalRec(:,:,Params.LateComparts)=CalibParams.inf_vaginalRecRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_Late);
        inf_basePerActProb_VaginalRec(:,:,Params.AIDSComparts)=CalibParams.inf_vaginalRecRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_AIDS);

        
    % Anal Receptive
        % Pre-allocate
        inf_basePerActProb_AnalRec = zeros(Params.numStrats, Params.numStrats, Params.numComparts);

        inf_basePerActProb_AnalRec(:,:,Params.AcuteComparts)=CalibParams.inf_analRecRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_Acute);
        inf_basePerActProb_AnalRec(:,:,Params.LatentAComparts)=CalibParams.inf_analRecRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_LatentA);
        inf_basePerActProb_AnalRec(:,:,Params.LatentBComparts)=CalibParams.inf_analRecRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_LatentB);
        inf_basePerActProb_AnalRec(:,:,Params.LateComparts)=CalibParams.inf_analRecRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_Late);
        inf_basePerActProb_AnalRec(:,:,Params.AIDSComparts)=CalibParams.inf_analRecRisk.paramValue(i) * Params.store_relRiskSexTransmission(Params.stage_AIDS);


    % Needle Risk
        % Pre-allocate
        inf_basePerNeedleProb = zeros(Params.numComparts,Params.numStrats,Params.numStrats);
        
        inf_basePerNeedleProb(Params.HIVComparts,:,:) = ...
            Params.idx_3d_IDU_IDU(Params.HIVComparts,:,:) * CalibParams.inf_sharedNeedleRisk.paramValue(i);
        

% Update BetaFactors (read into CalcBetas.m)

    % Vaginal Insertive Risk
    
        % Reductions due to
            % Circumcision
            % VLS
        
        % Apply circumcision reduction
        BetaFactor.adjPerActVaginalInsProb = ...
            inf_basePerActProb_VaginalIns .* ...
            (1-Params.inf_circReductVI_MF);
        
        % Apply VLS reduction
            for c = Params.VLSComparts
                BetaFactor.adjPerActVaginalInsProb(:,:,c)= ...
                    BetaFactor.adjPerActVaginalInsProb(:,:,c) .*...
                    (1 - Params.inf_vlsReductS_273x273);
            end
            
    
    % Vaginal Receptive Risk
        % Adjusted per-unprotected-receptive-act sexual risk 

        % Reduction due to 
            % VLS
   
        % Preallocate: apply base risk to full 273x273x31 matrix
        BetaFactor.adjPerActVaginalRecProb = inf_basePerActProb_VaginalRec;
        
        % Apply VLS reduction
        for c = Params.VLSComparts
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
        for c = Params.VLSComparts
           BetaFactor.adjPerActAnalInsProb(:,:,c) = ...     
                BetaFactor.adjPerActAnalInsProb(:,:,c) .* ...
                (1 - Params.inf_vlsReductS_273x273);
        end
        
        
        
   % Anal Receptive Risk
    
        % Reduction due to 
            % VLS
            
        % Preallocate: apply base risk to full 273x273x31 matrix
        BetaFactor.adjPerActAnalRecProb = inf_basePerActProb_AnalRec;
        
         % Apply VLS Reduction
        for c = Params.VLSComparts
            BetaFactor.adjPerActAnalRecProb(:,:,c) = ...
                inf_basePerActProb_AnalRec(:,:,c) .* ...
                (1 - Params.inf_vlsReductS_273x273);
        end
    
%% 10. Recalculate Betas

    % Note: The zero refers to not being called from InitParams. This is
    % included because the period 3 betas do not have to be updated when
    % a calibration is run.
    
    % Sexual betas
    Params = CalcBetas(Params,BetaFactor,0); 

    % Needle Betas
    for betaNSection = 1:1
    

% CALCULATED: NUMBER OF NEEDLES SHARED PER PARTNER
    % Will be adjusted later adjusted based on awareness
    % Calculated as: 
        % [Number injections per year] * [Pct inj shared] / [Number of partners]
    % Limited to 

    if Params.behav_defaultpartnersNeedle > 0

        behav_needlesPerPartner = ...
               Params.behav_numInjectionsPerYear ...
            .* CalibParams.inf_pctInjectionsShared.paramValue(i) ...
            / Params.behav_defaultpartnersNeedle;

    else
        Params.behav_defaultpartnersNeedle = 0;
    end
    
% NUMBER OF NEEDLES SHARED PER PARTNER (ADJUSTED)
    % Adjusted based on awareness

   %Pre-allocate
    inf_adjNeedlesShared_NoPrEP_1to4 = zeros(Params.numComparts,Params.numStrats,Params.numStrats);
    inf_adjNeedlesShared_NoPrEP_5 = zeros(Params.numComparts,Params.numStrats,Params.numStrats);
    
    % Calculate number of needles shared (unaware/aware)
        behav_needlesSharedUnaware = behav_needlesPerPartner;
        behav_needlesSharedAware_1to4 = behav_needlesPerPartner * (1-Params.behav_needleReductDueToAware_1to4);
        behav_needlesSharedAware_5 = behav_needlesPerPartner * (1-Params.behav_needleReductDueToAware_5);
    
    % Apply to parameter
            % Note: It's ok that they're all filled because the betas will remove the
            % compartments and people who won't be infected via needle.
        inf_adjNeedlesShared_NoPrEP_1to4(:,:,:) = behav_needlesSharedUnaware;
        inf_adjNeedlesShared_NoPrEP_5(:,:,:) = behav_needlesSharedUnaware;
            %Overwrite for the people who are aware of the their HIV+ status
        inf_adjNeedlesShared_NoPrEP_1to4(Params.AwareComparts,:,:) = behav_needlesSharedAware_1to4;
        inf_adjNeedlesShared_NoPrEP_5(Params.AwareComparts,:,:) = behav_needlesSharedAware_5;
        
        % Adjust for PrEP
        inf_adjNeedlesShared_PrEP_1to4 = ...
            inf_adjNeedlesShared_NoPrEP_1to4 .* (1 + Params.behav_incrNeedleShareDueToPrEP);
        inf_adjNeedlesShared_PrEP_5 = ...
            inf_adjNeedlesShared_NoPrEP_5 .* (1 + Params.behav_incrNeedleShareDueToPrEP);
    
        
% CALCULATE: BETA N

    % Pre-allocate adjusted needle risk (set equal to base risk)
    inf_adjPerNeedleRisk = inf_basePerNeedleProb;

    % Reduce transmission risk for PLWH who are VLS
    inf_adjPerNeedleRisk(Params.VLSComparts,:,:) = ...
        inf_basePerNeedleProb(Params.VLSComparts,:,:) * (1 - CalibParams.inf_vlsNeedleReduct.paramValue(i));
    
            % No PrEP
        inf_betaN_NoPrEP_1to4 = 1 - (1 - inf_adjPerNeedleRisk) .^ inf_adjNeedlesShared_NoPrEP_1to4;
        Params.oneMinusBetaN_NoPrEP_1to4 = 1 - inf_betaN_NoPrEP_1to4;
        
        inf_betaN_NoPrEP_5 = 1 - (1 - inf_adjPerNeedleRisk) .^ inf_adjNeedlesShared_NoPrEP_5;
        Params.oneMinusBetaN_NoPrEP_5 = 1 - inf_betaN_NoPrEP_5;

        % With PrEP
            % Only calculate when PrEP is run
        if  Params.CalcPrEPSpecificBetas == Params.ind_PrEP            
            inf_betaN_PrEP_1to4 = 1 - (1 - inf_adjPerNeedleRisk) .^ inf_adjNeedlesShared_PrEP_1to4;
            Params.oneMinusBetaN_PrEP_1to4 = 1 - inf_betaN_PrEP_1to4;
            inf_betaN_PrEP_5 = 1 - (1 - inf_adjPerNeedleRisk) .^ inf_adjNeedlesShared_PrEP_5;
            Params.oneMinusBetaN_PrEP_5 = 1 - inf_betaN_PrEP_5;
        end

        
    end

end