function Params = CalcBetas(Params,BetaFactor,calledFromInitParams)
%% Purpose: Calculate the per partnership risk of infection

%% 1. Calculate the betas

    % Initialize counter (increments each time a beta is calculated)
    BetaCounter = 1;

    % Initialize force of infection variables
        % 1 is [V only], 2 is [A only], 3 is [V some A], 4 is [A some V]
    Vonly = Params.force_Vonly;
    Aonly = Params.force_Aonly;
    VsomeA = Params.force_VsomeA;
    AsomeV = Params.force_AsomeV;

    % PrEP indicators

        % Only calculate betas with PrEP if PrEP behavior changes applied
        if Params.CalcPrEPSpecificBetas == Params.ind_PrEP
            betasWithPrEP = Params.ind_PrEP;
        else
            betasWithPrEP = Params.ind_NoPrEP; % This will set the loop to 1:1 and it won't calculate the PrEP betas
        end
    
% Loop through the non-PrEP/PrEP betas    
for decisionPrEP = Params.ind_NoPrEP:betasWithPrEP % Either 1:1 or 1:2
    
% Loop through forces of infection
for forceOfInf = Vonly : AsomeV % 1:4
    
    % Assign parameter values as appropriate - JC to move numSexContacts
    % within period loop
    
        % VI betas
        if (forceOfInf == Vonly || forceOfInf == VsomeA)

            CondomEfficacy = Params.behav_condomEfficacy';

            TransProbRec = BetaFactor.adjPerActVaginalRecProb;
            TransProbIns = BetaFactor.adjPerActVaginalInsProb; 

            pctCondomUse_NoPrEP_1to4 = Params.factor_adjPctCondomUseV_NoPrEP_1to4;
            pctCondomUse_NoPrEP_5 = Params.factor_adjPctCondomUseV_NoPrEP_5;
            pctCondomUse_PrEP_1to4 = Params.factor_adjPctCondomUseV_PrEP_1to4;
            pctCondomUse_PrEP_5 = Params.factor_adjPctCondomUseV_PrEP_5;
            
        % AI betas
        else

            CondomEfficacy = Params.behav_condomEfficacy';

            TransProbRec = BetaFactor.adjPerActAnalRecProb;
            TransProbIns = BetaFactor.adjPerActAnalInsProb;           

            pctCondomUse_NoPrEP_1to4 = Params.factor_adjPctCondomUseA_NoPrEP_1to4;
            pctCondomUse_NoPrEP_5 = Params.factor_adjPctCondomUseA_NoPrEP_5;
            pctCondomUse_PrEP_1to4 = Params.factor_adjPctCondomUseA_PrEP_1to4;
            pctCondomUse_PrEP_5 = Params.factor_adjPctCondomUseA_PrEP_5;

        end                            
                                          
    % Calculate the betas
        % Loops through each model time period
       for Period = 1:Params.numPeriod 
          
          % Pull num sex contacts and pct condom used for each time period 
          switch Period 
               case 1
                   % No PrEP
                   if decisionPrEP == Params.ind_NoPrEP
                       numSexContacts = Params.factor_adjNumSexContact_NoPrEP_1and5;
                       pctCondomUse = pctCondomUse_NoPrEP_1to4; 
                   else % PrEP
                       numSexContacts = Params.factor_adjNumSexContact_PrEP_1and5;
                       pctCondomUse = pctCondomUse_PrEP_1to4; 
                   end                                                  
               case 2
                   % No PrEP
                   if decisionPrEP == Params.ind_NoPrEP
                       numSexContacts = Params.factor_adjNumSexContact_NoPrEP_2;
                       pctCondomUse = pctCondomUse_NoPrEP_1to4; 
                   else % PrEP
                       numSexContacts = Params.factor_adjNumSexContact_PrEP_2;
                       pctCondomUse = pctCondomUse_PrEP_1to4; 
                   end
               case 3
                   % No PrEP
                   if decisionPrEP == Params.ind_NoPrEP
                       numSexContacts = Params.factor_adjNumSexContact_NoPrEP_3;
                       pctCondomUse = pctCondomUse_NoPrEP_1to4; 
                   else % PrEP
                       numSexContacts = Params.factor_adjNumSexContact_PrEP_3;
                       pctCondomUse = pctCondomUse_PrEP_1to4; 
                   end
               case 4
                   % No PrEP
                   if decisionPrEP == Params.ind_NoPrEP
                       numSexContacts = Params.factor_adjNumSexContact_NoPrEP_4;
                       pctCondomUse = pctCondomUse_NoPrEP_1to4; 
                   else % PrEP
                       numSexContacts = Params.factor_adjNumSexContact_PrEP_4;
                       pctCondomUse = pctCondomUse_PrEP_1to4; 
                   end
               case 5
                   % No PrEP
                   if decisionPrEP == Params.ind_NoPrEP
                       numSexContacts = Params.factor_adjNumSexContact_NoPrEP_1and5;
                       pctCondomUse = pctCondomUse_NoPrEP_5; 
                   else % PrEP
                       numSexContacts = Params.factor_adjNumSexContact_PrEP_1and5;
                       pctCondomUse = pctCondomUse_PrEP_5; 
                   end
          end           
           
           % Calculate the exponents 
           % Note: don't include pctContactsThisType 
           ExponentRecUnpro = numSexContacts .* (1 - pctCondomUse) .* ...
               Params.factor_proportionRec;          
           ExponentRecPro = numSexContacts .* pctCondomUse .* ...
               Params.factor_proportionRec;         
           ExponentInsUnpro = numSexContacts .* (1 - pctCondomUse) .* ...
               (1-Params.factor_proportionRec);
           ExponentInsPro = numSexContacts .* pctCondomUse .* ...
               (1-Params.factor_proportionRec);          
          
           % This if statement prevents the code from calculating the betas
           % for periods 3, 4, 5 if CalcBetas wasn't called from InitParams.
           % The calibration code only needs to update the betas for
           % periods 1 and 2 (updated by MClinkscales on 5/17/2022 for T4
           % and T5)
           % MC had set second condition in If statement to be "Period <3" but 
           % KH updated on 12/27/22 to also update betas for periods 3 and
           % 4 since the values for periods 2-4 are set by the second set of calibrated
           % input values (versus the first set that is only for time period 1). 
           % Condition is now "Period <5." 
           if calledFromInitParams == 1 || Period <5 % Not called from within initParams or not period 5
               
               % Determine what percent of contacts are this type
               switch forceOfInf
                   case Vonly
                       PercentContactsThisType = 1;
                   case Aonly
                       PercentContactsThisType = 1;
                   case VsomeA
                       switch Period
                           case 1
                                PercentContactsThisType = BetaFactor.pctContactsVaginal_InMFWithA_1;
                           case {2,3,4}
                                PercentContactsThisType = BetaFactor.pctContactsVaginal_InMFWithA_2to4;
                           case 5
                                PercentContactsThisType = BetaFactor.pctContactsVaginal_InMFWithA_5;                           
                       end
                   case AsomeV
                        switch Period
                           case 1
                                PercentContactsThisType = BetaFactor.pctContactsAnal_InMFWithA_1;
                           case {2,3,4}
                                PercentContactsThisType = BetaFactor.pctContactsAnal_InMFWithA_2to4;
                           case 5
                                PercentContactsThisType = BetaFactor.pctContactsAnal_InMFWithA_5;                           
                        end
               end

               % Calculate the beta [30x273x273]
               store = (1 - ...
                    (1- TransProbRec) .^ (ExponentRecUnpro .* PercentContactsThisType) ...
                    .* (1 - TransProbRec .* (1- repmat(CondomEfficacy,1,Params.numStrats,Params.numComparts))) .^ (ExponentRecPro .* PercentContactsThisType) ...
                    .* (1 - TransProbIns) .^ (ExponentInsUnpro .* PercentContactsThisType)...
                    .* (1 - TransProbIns .* (1- repmat(CondomEfficacy,1,Params.numStrats,Params.numComparts))) .^ (ExponentInsPro .* PercentContactsThisType));

               % Put the calculated beta in the proper form
                store_2 = permute(store,[3 1 2]);
                store_3 = min(1,store_2);
                store_oneMinusBeta = 1 - store_3;   
           
           else
                store_oneMinusBeta = 0;
           end
           
            % Apply the calculated beta to the proper parameter
            switch BetaCounter % JC change; Updated by MClinkscales on 5/17/2022
                case 1
                    Params.oneMinusBetaVonly_NoPrEP_period1 = store_oneMinusBeta;
                case 2
                    Params.oneMinusBetaVonly_NoPrEP_period2 = store_oneMinusBeta;
                case 3
                    Params.oneMinusBetaVonly_NoPrEP_period3 = store_oneMinusBeta;
                case 4
                    Params.oneMinusBetaVonly_NoPrEP_period4 = store_oneMinusBeta;
                case 5
                    Params.oneMinusBetaVonly_NoPrEP_period5 = store_oneMinusBeta;
                case 6
                    Params.oneMinusBetaAonly_NoPrEP_period1 = store_oneMinusBeta; 
                case 7
                    Params.oneMinusBetaAonly_NoPrEP_period2 = store_oneMinusBeta;
                case 8
                    Params.oneMinusBetaAonly_NoPrEP_period3 = store_oneMinusBeta;    
                case 9
                    Params.oneMinusBetaAonly_NoPrEP_period4 = store_oneMinusBeta;    
                case 10
                    Params.oneMinusBetaAonly_NoPrEP_period5 = store_oneMinusBeta;  
                case 11
                    Params.oneMinusBetaVsomeAI_NoPrEP_period1 = store_oneMinusBeta; 
                case 12
                    Params.oneMinusBetaVsomeAI_NoPrEP_period2 = store_oneMinusBeta; 
                case 13
                    Params.oneMinusBetaVsomeAI_NoPrEP_period3 = store_oneMinusBeta;
                case 14
                    Params.oneMinusBetaVsomeAI_NoPrEP_period4 = store_oneMinusBeta;
                case 15
                    Params.oneMinusBetaVsomeAI_NoPrEP_period5 = store_oneMinusBeta;
                case 16
                    Params.oneMinusBetaAsomeVI_NoPrEP_period1 = store_oneMinusBeta;              
                case 17
                    Params.oneMinusBetaAsomeVI_NoPrEP_period2 = store_oneMinusBeta;
                case 18
                    Params.oneMinusBetaAsomeVI_NoPrEP_period3 = store_oneMinusBeta;
                case 19
                    Params.oneMinusBetaAsomeVI_NoPrEP_period4 = store_oneMinusBeta;
                case 20
                    Params.oneMinusBetaAsomeVI_NoPrEP_period5 = store_oneMinusBeta;
                case 21
                    Params.oneMinusBetaVonly_PrEP_period1 = store_oneMinusBeta; % cases 21-40 only pulled if PrEP is applied and behavior changes related to PrEP are used
                case 22
                    Params.oneMinusBetaVonly_PrEP_period2 = store_oneMinusBeta;
                case 23
                    Params.oneMinusBetaVonly_PrEP_period3 = store_oneMinusBeta;
                case 24
                    Params.oneMinusBetaVonly_PrEP_period4 = store_oneMinusBeta;
                case 25
                    Params.oneMinusBetaVonly_PrEP_period5 = store_oneMinusBeta;
                case 26
                    Params.oneMinusBetaAonly_PrEP_period1 = store_oneMinusBeta;
                case 27
                    Params.oneMinusBetaAonly_PrEP_period2 = store_oneMinusBeta;
                case 28
                    Params.oneMinusBetaAonly_PrEP_period3 = store_oneMinusBeta;
                case 29
                    Params.oneMinusBetaAonly_PrEP_period4 = store_oneMinusBeta;
                case 30
                    Params.oneMinusBetaAonly_PrEP_period5 = store_oneMinusBeta;
                case 31
                    Params.oneMinusBetaVsomeAI_PrEP_period1 = store_oneMinusBeta;
                case 32
                    Params.oneMinusBetaVsomeAI_PrEP_period2 = store_oneMinusBeta;
                case 33
                    Params.oneMinusBetaVsomeAI_PrEP_period3 = store_oneMinusBeta;
                case 34
                    Params.oneMinusBetaVsomeAI_PrEP_period4 = store_oneMinusBeta;
                case 35
                    Params.oneMinusBetaVsomeAI_PrEP_period5 = store_oneMinusBeta;
                case 36
                    Params.oneMinusBetaAsomeVI_PrEP_period1 = store_oneMinusBeta;
                case 37
                    Params.oneMinusBetaAsomeVI_PrEP_period2 = store_oneMinusBeta;
                case 38
                    Params.oneMinusBetaAsomeVI_PrEP_period3 = store_oneMinusBeta;
                case 39
                    Params.oneMinusBetaAsomeVI_PrEP_period4 = store_oneMinusBeta;
                case 40
                    Params.oneMinusBetaAsomeVI_PrEP_period5 = store_oneMinusBeta;
            end
                
            % Update counter
            BetaCounter = BetaCounter + 1;
       end
end
end

%% 3Oct2022 KH commented out because it was wrong (missing two PrEP calcs and didn't actually make code much more efficient
%% 2. Apply period 5 betas (if appropriate) % Changes will need to be reflected in calibration code for changes 
% 
%     % If this was called from within InitParams apply period 5 betas.
%     % Updated by Mclinkscales on 5/16/2022
%     if calledFromInitParams == 1
%         Params.oneMinusBetaVonly_NoPrEP_period5 = oneMinusBetaVonly_NoPrEP_period5;
%         Params.oneMinusBetaAonly_NoPrEP_period5 = oneMinusBetaAonly_NoPrEP_period5;
%         Params.oneMinusBetaVsomeAI_NoPrEP_period5 = oneMinusBetaVsomeAI_NoPrEP_period5;
%         Params.oneMinusBetaAsomeVI_NoPrEP_period5 = oneMinusBetaAsomeVI_NoPrEP_period5;
%         
%         if Params.CalcPrEPSpecificBetas == Params.ind_PrEP
%             Params.oneMinusBetaVsomeAI_PrEP_period5 = oneMinusBetaVsomeAI_PrEP_period5;
%             Params.oneMinusBetaAsomeVI_PrEP_period5 = oneMinusBetaAsomeVI_PrEP_period5;
%         end
%     end

end