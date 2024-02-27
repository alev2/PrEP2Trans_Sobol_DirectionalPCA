function [Compartments, Results] = StepStateDisc(Params, TransRates, Compartments, InfRate, InfRate_PrEP_Oral_High, InfRate_PrEP_Oral_Low, InfRate_PrEP_Inject_High, InfRate_PrEP_Inject_Low, Year, k, Results, TTProg)
%% Purpose: For the discrete version of the model: a function which 
% advances the state of the system.
% Called from: HIVEpiModel

%% 1. Calculate the rates of change for all health states

        % Initialize variables
            delta_pop = zeros(Params.numComparts, Params.numStrats);
            onemat = ones(Params.numComparts, 1);
            Transitions(Params.numComparts,Params.numComparts,Params.numStrats) = 0;

        % Calculate Rates
            for i = 1:Params.numStrats;

                % This is a 31 x 31 matrix of the population
                C = repmat(Compartments(:,i),[1,Params.numComparts]);

                % Peel off strat
                Tr = TransRates.TSRates(:,:,i);

                % Flows in from each compartment
                Transitions(:,:,i) = C .* Tr;    
                delta_pop(:,i) = (onemat'*(Transitions(:,:,i)))' - (Transitions(:,:,i))*onemat;
            end

%% 2. Update Compartments

        % Update Compartments
        Compartments = max(0,Compartments + delta_pop);

        % Age population
        [Compartments, TTProg] = Aging(Params, Compartments, TTProg); 

%% 3. Collect Results

    % Call Collect Results
    Results = CollectResults(Params, Compartments, Year, InfRate, InfRate_PrEP_Oral_High, InfRate_PrEP_Oral_Low, InfRate_PrEP_Inject_High, InfRate_PrEP_Inject_Low, k, Transitions, Results, TTProg);

end