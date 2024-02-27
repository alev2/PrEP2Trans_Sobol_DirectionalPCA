function delta_pop =  RatesODE(Y, Params, TransRates, TTProg, nFevalsPerTStep, Year)
%% Purpose: A function that calculates the number of people flowing 
% between compartments.
% Called from: StepStateCont

%% 1. Initialize Compartments and Global Variables

    % Apply ODE-compatible compartment sizes [5460x1] to the 
    % TransRates-compatible "Compartment" variable [30x273]
    Compartments = reshape(Y', Params.numComparts,Params.numStrats);

% Declare global variables
    % These are included to pull the ODE-calculated rates out of the
    % function "RatesODE" - used to calculate the "Transitions" variable.
    
    % You can't pull them directly out of ODE45 - the only outcomes are Y
    % and T (i.e., compartments and microstep)- but to calculate our outcomes
    % we need these rates too.
    
    % Note: global variables are marked with a g_
    global g_TransRatesAnnual;
    global g_RelativeInfRates;
    global g_lambdaVonly;
    global g_lambdaAonly;
    global g_lambdaVsomeAI;
    global g_lambdaAsomeVI;
    global g_lambdaN;
    global g_odeStepCounter;
    
    global g_RelativeInfRates_PrEP_Oral_High;
    global g_lambdaVonly_PrEP_Oral_High;
    global g_lambdaAonly_PrEP_Oral_High;
    global g_lambdaVsomeAI_PrEP_Oral_High;
    global g_lambdaAsomeVI_PrEP_Oral_High;
    global g_lambdaN_PrEP_Oral_High;
    
    global g_RelativeInfRates_PrEP_Oral_Low;
    global g_lambdaVonly_PrEP_Oral_Low;
    global g_lambdaAonly_PrEP_Oral_Low;
    global g_lambdaVsomeAI_PrEP_Oral_Low;
    global g_lambdaAsomeVI_PrEP_Oral_Low;
    global g_lambdaN_PrEP_Oral_Low;
    
    global g_RelativeInfRates_PrEP_Inject_High;
    global g_lambdaVonly_PrEP_Inject_High;
    global g_lambdaAonly_PrEP_Inject_High;
    global g_lambdaVsomeAI_PrEP_Inject_High;
    global g_lambdaAsomeVI_PrEP_Inject_High;
    global g_lambdaN_PrEP_Inject_High;
    
    global g_RelativeInfRates_PrEP_Inject_Low;
    global g_lambdaVonly_PrEP_Inject_Low;
    global g_lambdaAonly_PrEP_Inject_Low;
    global g_lambdaVsomeAI_PrEP_Inject_Low;
    global g_lambdaAsomeVI_PrEP_Inject_Low;
    global g_lambdaN_PrEP_Inject_Low;
    
%% 2. Calculate infection rates
  
    % Calculate infection rates
        % Note: this code is duplicated in CalcTransRates

        % Indicators    
        lambdasNoPrEP = Params.ind_NoPrEP; 
        lambdasPrEP = Params.ind_PrEP;
        
        % Normal infection rate
        InfRate = CalcInfectionRates(Params, Compartments, TTProg, lambdasNoPrEP, Params.ind_OralPrEP_High); 

        % PrEP infection rate
        
            % No PrEP behavior differences or direct entry and no PrEP
            % initiation
            if Params.CalcPrEPSpecificBetas == Params.ind_NoPrEP
                
                % Reduce normal lambdas by PrEP reduction % JCPrEPUpdate: modified to calculated reduced lambdas
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

            else
                % Calculate PrEP-specific lambdas
                    % (Incorporates behavior change) % JCPrEPUpdate: modified to call CalcInfectionRates separately for oral and injectable PrEP
                    InfRate_PrEP_Oral_High = CalcInfectionRates(Params, Compartments, TTProg, lambdasPrEP, Params.ind_OralPrEP_High);
                    InfRate_PrEP_Oral_Low = CalcInfectionRates(Params, Compartments, TTProg, lambdasPrEP, Params.ind_OralPrEP_Low);
                    InfRate_PrEP_Inject_High = CalcInfectionRates(Params, Compartments, TTProg, lambdasPrEP, Params.ind_InjectPrEP_High);
                    InfRate_PrEP_Inject_Low = CalcInfectionRates(Params, Compartments, TTProg, lambdasPrEP, Params.ind_InjectPrEP_Low);
            end
            

    % Apply total infection rate to TransRates struct % JCPrEPUpdate: modified to include
        % transitions from all PrEP states
     TransRates.GetInfected(Params.A1,Params.B1,:) = InfRate.totalInfRate;
     TransRates.GetInfected(Params.A6,Params.B3,:) = InfRate_PrEP_Oral_High.totalInfRate;
     TransRates.GetInfected(Params.A7,Params.B3,:) = InfRate_PrEP_Oral_Low.totalInfRate;
     TransRates.GetInfected(Params.A8,Params.B3,:) = InfRate_PrEP_Inject_High.totalInfRate;
     TransRates.GetInfected(Params.A9,Params.B3,:) = InfRate_PrEP_Inject_Low.totalInfRate;
     
     
%% 3. Calculate Full TransRates matrix

    % Error checking
        chknonzero = (TransRates.BecomeAware>0) + (TransRates.GetInfected>0) ...
               + (TransRates.ProgressHIV>0) + (TransRates.LinkToCare>0) ...
               + (TransRates.DropOut>0) + (TransRates.StartART>0) ...
               + (TransRates.DieNormal>0) ...
               + (TransRates.DieAIDS>0) + (TransRates.StayDeadAged>0);            
           
        if max(max(max(chknonzero)))>1
            error('ERROR! Elements of TransRates matrix conflict! Nonzero elems in more than 1 Trans')
        end
           clear chknonzero

    % Calculate full TransRates matrix
    TransRates.AnnualRates = ...
            TransRates.BecomeAware + TransRates.GetInfected ...
           + TransRates.ProgressHIV + TransRates.LinkToCare ...
           + TransRates.DropOut + TransRates.StartART ...
           + TransRates.DieNormal ...
           + TransRates.DieAIDS;           
       
     % Make Diagonal equal to the negative sum of the rows
     for i = 1:Params.numStrats;
         TransRates.AnnualRates(:,:,i) = TransRates.AnnualRates(:,:,i)+...
                diag(-sum(TransRates.AnnualRates(:,:,i),2));
     end

%% 4. Assign Global Variables

% Assign global variables their values
    if isempty(g_lambdaVonly) == 1 % First run 
        
        if mod(g_odeStepCounter,nFevalsPerTStep) == 2 % starts at 2nd call
            
            g_TransRatesAnnual = TransRates.AnnualRates;    
            g_RelativeInfRates = InfRate.RelativeInfRates;
            g_lambdaVonly = InfRate.lambdaVonly;
            g_lambdaAonly = InfRate.lambdaAonly;
            g_lambdaVsomeAI = InfRate.lambdaVsomeAI;
            g_lambdaAsomeVI = InfRate.lambdaAsomeVI;
            g_lambdaN = InfRate.lambdaN;
                        
            g_RelativeInfRates_PrEP_Oral_High = InfRate_PrEP_Oral_High.RelativeInfRates;
            g_lambdaVonly_PrEP_Oral_High = InfRate_PrEP_Oral_High.lambdaVonly;
            g_lambdaAonly_PrEP_Oral_High = InfRate_PrEP_Oral_High.lambdaAonly;
            g_lambdaVsomeAI_PrEP_Oral_High = InfRate_PrEP_Oral_High.lambdaVsomeAI;
            g_lambdaAsomeVI_PrEP_Oral_High = InfRate_PrEP_Oral_High.lambdaAsomeVI;
            g_lambdaN_PrEP_Oral_High = InfRate_PrEP_Oral_High.lambdaN;
            
            g_RelativeInfRates_PrEP_Oral_Low = InfRate_PrEP_Oral_Low.RelativeInfRates;
            g_lambdaVonly_PrEP_Oral_Low = InfRate_PrEP_Oral_Low.lambdaVonly;
            g_lambdaAonly_PrEP_Oral_Low = InfRate_PrEP_Oral_Low.lambdaAonly;
            g_lambdaVsomeAI_PrEP_Oral_Low = InfRate_PrEP_Oral_Low.lambdaVsomeAI;
            g_lambdaAsomeVI_PrEP_Oral_Low = InfRate_PrEP_Oral_Low.lambdaAsomeVI;
            g_lambdaN_PrEP_Oral_Low = InfRate_PrEP_Oral_Low.lambdaN;
            
            g_RelativeInfRates_PrEP_Inject_High = InfRate_PrEP_Inject_High.RelativeInfRates;
            g_lambdaVonly_PrEP_Inject_High = InfRate_PrEP_Inject_High.lambdaVonly;
            g_lambdaAonly_PrEP_Inject_High = InfRate_PrEP_Inject_High.lambdaAonly;
            g_lambdaVsomeAI_PrEP_Inject_High = InfRate_PrEP_Inject_High.lambdaVsomeAI;
            g_lambdaAsomeVI_PrEP_Inject_High = InfRate_PrEP_Inject_High.lambdaAsomeVI;
            g_lambdaN_PrEP_Inject_High = InfRate_PrEP_Inject_High.lambdaN;
            
            g_RelativeInfRates_PrEP_Inject_Low = InfRate_PrEP_Inject_Low.RelativeInfRates;
            g_lambdaVonly_PrEP_Inject_Low = InfRate_PrEP_Inject_Low.lambdaVonly;
            g_lambdaAonly_PrEP_Inject_Low = InfRate_PrEP_Inject_Low.lambdaAonly;
            g_lambdaVsomeAI_PrEP_Inject_Low = InfRate_PrEP_Inject_Low.lambdaVsomeAI;
            g_lambdaAsomeVI_PrEP_Inject_Low = InfRate_PrEP_Inject_Low.lambdaAsomeVI;
            g_lambdaN_PrEP_Inject_Low = InfRate_PrEP_Inject_Low.lambdaN;
        end
        
    else % All other runs, concatenate the current values to the existing variable
         if mod(g_odeStepCounter,nFevalsPerTStep) == 2 % grabs every 6th one starting at number 2
             
            g_TransRatesAnnual = cat(4,g_TransRatesAnnual,TransRates.AnnualRates);    
            g_RelativeInfRates = cat(3,g_RelativeInfRates, InfRate.RelativeInfRates);
            g_lambdaVonly = cat(2,g_lambdaVonly, InfRate.lambdaVonly);
            g_lambdaAonly = cat(2,g_lambdaAonly, InfRate.lambdaAonly);
            g_lambdaVsomeAI = cat(2,g_lambdaVsomeAI, InfRate.lambdaVsomeAI);
            g_lambdaAsomeVI = cat(2,g_lambdaAsomeVI, InfRate.lambdaAsomeVI);
            g_lambdaN = cat(2,g_lambdaN,InfRate.lambdaN);
            
            g_RelativeInfRates_PrEP_Oral_High = cat(3,g_RelativeInfRates_PrEP_Oral_High, InfRate_PrEP_Oral_High.RelativeInfRates);
            g_lambdaVonly_PrEP_Oral_High = cat(2,g_lambdaVonly_PrEP_Oral_High, InfRate_PrEP_Oral_High.lambdaVonly);
            g_lambdaAonly_PrEP_Oral_High = cat(2,g_lambdaAonly_PrEP_Oral_High, InfRate_PrEP_Oral_High.lambdaAonly);
            g_lambdaVsomeAI_PrEP_Oral_High = cat(2,g_lambdaVsomeAI_PrEP_Oral_High, InfRate_PrEP_Oral_High.lambdaVsomeAI);
            g_lambdaAsomeVI_PrEP_Oral_High = cat(2,g_lambdaAsomeVI_PrEP_Oral_High, InfRate_PrEP_Oral_High.lambdaAsomeVI);
            g_lambdaN_PrEP_Oral_High = cat(2,g_lambdaN_PrEP_Oral_High,InfRate_PrEP_Oral_High.lambdaN);
            
            g_RelativeInfRates_PrEP_Oral_Low = cat(3,g_RelativeInfRates_PrEP_Oral_Low, InfRate_PrEP_Oral_Low.RelativeInfRates);
            g_lambdaVonly_PrEP_Oral_Low = cat(2,g_lambdaVonly_PrEP_Oral_Low, InfRate_PrEP_Oral_Low.lambdaVonly);
            g_lambdaAonly_PrEP_Oral_Low = cat(2,g_lambdaAonly_PrEP_Oral_Low, InfRate_PrEP_Oral_Low.lambdaAonly);
            g_lambdaVsomeAI_PrEP_Oral_Low = cat(2,g_lambdaVsomeAI_PrEP_Oral_Low, InfRate_PrEP_Oral_Low.lambdaVsomeAI);
            g_lambdaAsomeVI_PrEP_Oral_Low = cat(2,g_lambdaAsomeVI_PrEP_Oral_Low, InfRate_PrEP_Oral_Low.lambdaAsomeVI);
            g_lambdaN_PrEP_Oral_Low = cat(2,g_lambdaN_PrEP_Oral_Low,InfRate_PrEP_Oral_Low.lambdaN);
            
            g_RelativeInfRates_PrEP_Inject_High = cat(3,g_RelativeInfRates_PrEP_Inject_High, InfRate_PrEP_Inject_High.RelativeInfRates);
            g_lambdaVonly_PrEP_Inject_High = cat(2,g_lambdaVonly_PrEP_Inject_High, InfRate_PrEP_Inject_High.lambdaVonly);
            g_lambdaAonly_PrEP_Inject_High = cat(2,g_lambdaAonly_PrEP_Inject_High, InfRate_PrEP_Inject_High.lambdaAonly);
            g_lambdaVsomeAI_PrEP_Inject_High = cat(2,g_lambdaVsomeAI_PrEP_Inject_High, InfRate_PrEP_Inject_High.lambdaVsomeAI);
            g_lambdaAsomeVI_PrEP_Inject_High = cat(2,g_lambdaAsomeVI_PrEP_Inject_High, InfRate_PrEP_Inject_High.lambdaAsomeVI);
            g_lambdaN_PrEP_Inject_High = cat(2,g_lambdaN_PrEP_Inject_High,InfRate_PrEP_Inject_High.lambdaN);
            
            g_RelativeInfRates_PrEP_Inject_Low = cat(3,g_RelativeInfRates_PrEP_Inject_Low, InfRate_PrEP_Inject_Low.RelativeInfRates);
            g_lambdaVonly_PrEP_Inject_Low = cat(2,g_lambdaVonly_PrEP_Inject_Low, InfRate_PrEP_Inject_Low.lambdaVonly);
            g_lambdaAonly_PrEP_Inject_Low = cat(2,g_lambdaAonly_PrEP_Inject_Low, InfRate_PrEP_Inject_Low.lambdaAonly);
            g_lambdaVsomeAI_PrEP_Inject_Low = cat(2,g_lambdaVsomeAI_PrEP_Inject_Low, InfRate_PrEP_Inject_Low.lambdaVsomeAI);
            g_lambdaAsomeVI_PrEP_Inject_Low = cat(2,g_lambdaAsomeVI_PrEP_Inject_Low, InfRate_PrEP_Inject_Low.lambdaAsomeVI);
            g_lambdaN_PrEP_Inject_Low = cat(2,g_lambdaN_PrEP_Inject_Low,InfRate_PrEP_Inject_Low.lambdaN);
         end
    end

clear diagsum diagones chkdiag Transitions

%% 5. Calculate number of people flowing between compartments
    
    % This code replicates the compartments vector to make it the same size
    % as the transrates matrix (B)
    CompartRep = permute(repmat(Compartments,[1, 1, Params.numComparts]), [1 3 2]);

    % CompartRep transpose (Bp)
    CompartRepT = permute(CompartRep, [2 1 3]);

    % Big matrix for TransRates
    Rateslong = TransRates.AnnualRates;

    
    % TransRatesLong transpose 
    RateslongT = permute(TransRates.AnnualRates, [2 1 3]);

    Rates_30x273 = squeeze(sum((CompartRep.*Rateslong - CompartRepT.*RateslongT),1));
   
    % Calculate new population (based on rates)
    Compartments_AfterRates = Compartments + Rates_30x273;
    
    % Age the new population
    [Compartments_AfterRatesAndAging, ~] = Aging(Params,Compartments_AfterRates, TTProg);
    
    % Net number of people who move
    delta_pop_30x273 = Compartments_AfterRatesAndAging - Compartments;
    
    % Convert to ODE-compatible variable size [5460x1]
    delta_pop = reshape(delta_pop_30x273,1,Params.numComparts*Params.numStrats)';

    % Increment the counter
        % It will add 1 everytime ratesODE is called.
        % If RefineVal is 1, only the results from every 6 function call
        % will be used.
    g_odeStepCounter = g_odeStepCounter + 1;

end