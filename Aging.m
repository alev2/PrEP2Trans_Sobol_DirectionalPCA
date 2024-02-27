function [Compartments, TTProg] = Aging(Params,Compartments,TTProg)
%% Purpose: A function which ages the population each time step

%% 1. Calculate number of people aging

    % Percent aging
    pctCohortAgingAnnual = 1./Params.hiv_durAge_a;
    
    % Zero out aging for 65+ age group (added 08/02/2018)
    indicator_65plus_zeros = (Params.ageIndicator(:,1) + Params.ageIndicator(:,2) ...
            + Params.ageIndicator(:,3) + Params.ageIndicator(:,4) ...
            + Params.ageIndicator(:,5) + Params.ageIndicator(:,6));
    pctCohortAgingAnnual = pctCohortAgingAnnual .* indicator_65plus_zeros;
    
    pctCohortAgingEachStep = pctCohortAgingAnnual *  Params.tt_tstep;

    % Calculating the number of people aging

    % numPeopleAging [30x273]
    numPeopleAging = zeros(Params.numComparts,Params.numStrats);
 
%   Original aging code    
%     for c = Params.NonAbsorbingComparts
%     % for c = 1:Params.numComparts - Params.numAbsorbingStates % c = 1:24
%        numPeopleAging(c,:) = Compartments(c,:) .* pctCohortAgingEachStep';
%     end
    
%     % Revised aging code by CG
%     for s = 1:Params.numStrats
%         if mod(s,Params.numAge) == 0 % 7th group is 65+ and doesn't age
%             continue
%         end
%         factor = sum(Compartments(:,s)) ./ sum(Compartments(:,s+1)); % compensates for strats of varying size
%         for c = Params.NonAbsorbingComparts
%             numPeopleAging(c,s) = (Compartments(c,s) + (Compartments(c,s+1) .* factor))/2  .* pctCohortAgingEachStep(s);
%         end
%     end
    
    % Revised aging code by KH so that distribution across the care 
    % continuum of aging group is more like the distribution of the next
    % age group
    pctSubpopPLWHInC = zeros(Params.numComparts,Params.numStrats);
    NumPLWHInStrat = zeros(Params.numStrats,1);
    %Can reprogram with vectors / matrix mult to be more efficient
    for s = 1:Params.numStrats
        for c = Params.HIVComparts
            NumPLWHInStrat(s) = NumPLWHInStrat(s) + sum(Compartments(c,s),1);
        end
    end    

    for s = 1:Params.numStrats
        for c = Params.HIVComparts
            pctSubpopPLWHInC(c,s) = Compartments(c,s) ./ NumPLWHInStrat(s);
        end
    end
        
    for s = 1:Params.numStrats
        if mod(s,Params.numAge) == 0 % Last group doesn't age
            continue
        end
        for c = Params.HIVComparts
            numPeopleAging(c,s) = Params.tt_tstep * ...
                (NumPLWHInStrat(s)./Params.hiv_durAge_a(s)   .*   ...
                    (pctSubpopPLWHInC(c,s) + (pctSubpopPLWHInC(c,s+1) - pctSubpopPLWHInC(c,s)) .* ...
                    Params.hiv_durAge_a(s) ./ (Params.hiv_durAge_a(s+1) + Params.hiv_durAge_a(s))) ...
                );
        end
    end
    
    for c = Params.UninfectedComparts
       numPeopleAging(c,:) = Compartments(c,:) .* pctCohortAgingEachStep';
    end
    
    numPeopleAging = min(numPeopleAging , Compartments);

%% 2. Initialize Indices

    % Gives the main 24 model states
    NonAbsorbingStates = Params.A1:Params.F5;

    % Index which includes that cohorts who are age 1
        store_age1 = Params.ageIndicator(:,1) == 1;
        idx_age1 = store_age1';
    
    % Index which includes the cohorts who are ages 2, 3, 4, 5, 6 and 7
        store234567 = (Params.ageIndicator(:,2) + Params.ageIndicator(:,3) ...
            + Params.ageIndicator(:,4) + Params.ageIndicator(:,5) ...
            + Params.ageIndicator(:,6) + Params.ageIndicator(:,7))';
        idx_age234567 = store234567 == 1;

    % Index which includes the cohorts who are ages 1, 2, 3, 4, 5 and 6
        store123456 = (Params.ageIndicator(:,1)+ Params.ageIndicator(:,2)...
            + Params.ageIndicator(:,3) + Params.ageIndicator(:,4) ...
            + Params.ageIndicator(:,5) + Params.ageIndicator(:,6))';
        idx_age123456 = store123456 == 1;
        
        
%% 3. Age people between cohorts (age stratifications)

        %This ensures the indices are correct when the people are added to
        %Compartments
        PeopleGoingInto234567 = zeros(Params.numComparts,Params.numStrats);
        PeopleGoingInto234567(NonAbsorbingStates, idx_age234567) = ...
            numPeopleAging(NonAbsorbingStates, idx_age123456);

    % Aging - add people from states (1,2,3,4,5,6) to (2,3,4,5,6,7), subtract from
    % (2,3,4,5,6,7)
        Compartments(NonAbsorbingStates,idx_age234567) = ...
            Compartments(NonAbsorbingStates,idx_age234567) ...
            + PeopleGoingInto234567(NonAbsorbingStates,idx_age234567) ...
            - numPeopleAging(NonAbsorbingStates,idx_age234567);

        %Subtract from state 1:
        Compartments(NonAbsorbingStates, idx_age1)= ...
            Compartments(NonAbsorbingStates, idx_age1) - ...
            numPeopleAging(NonAbsorbingStates,idx_age1);

%% 4. Age into the population

    % "Age in" = add new pops
           % not using indices allows the model to be flexible on where the
           % new pops flow into (may start at 13 or at 18)
                % Excel auto populates the correct age (either 1 or 2) for
                % new pops to start
   %% Example indicator code DELETE
%         Results.ann_DropOff_PrEP_HighRiskIDUs(yearCount,:) = ...
%                            Results.ann_DropOff_PrEP(yearCount,:) ...
%                         .* Params.popIndicator(:,Params.pop_IDU)' ...
%                         .* Params.riskLevelIndicator(:,Params.risk_Casual)';
%repmat(Params.popIndicator(:,Params.pop_IDU)',[30,1])
    
    % create 30x273 indicators for each trans group / sex / race combo
    counter = 0;
    for i = 1:Params.numPop
        for j = 1:2
            for k = 1:Params.numRace
                counter = counter + 1;
                indicator_popSexRace(:,:,counter) = repmat( ...
                    Params.popIndicator(:,i)'.* ...
                    Params.sexIndicator(:,j)' .* ...
                    Params.raceIndicator(:,k)' ...
                    ,[Params.numComparts,1]);
            end
        end
    end
    % check that all indicators added together makes full matrix of 1s
    % CONFIRMED, seems to be correct
    
    % this loop multiplies the compartments by the pop / sex / race
    % indicator and then the pop / sex / race growth rate
    counter = 0;
    for i = 1:Params.numPop
        for j = 1:2
            for k = 1:Params.numRace
                counter = counter + 1;
                newPopStrat(:,:,counter) = Compartments .* ...
                    indicator_popSexRace(:,:,counter) * ...
                    Params.growthRate(counter);
            end
        end
    end
    
    % sums all new pops together (each of the 15 combos), then puts all new people into A1, then
    % expands back to 30x273 (where only A1 should have non-zeros)
    newPopStrat = sum(sum(newPopStrat,3),1); %first sum dim 3 is 15 combos to 1, then dim 1 is 30 comparts to 1
    newPop = zeros(Params.numComparts,Params.numStrats);
    newPop(Params.A1,:) = newPopStrat;
    newPop(Params.AbsorbingComparts,:) = 0;
    
    % Applies entry rate by age group to new pop (normally 100% for 13-17
    % or 18-24)
    for i = 1:Params.numAge
        newPopAge(:,:,i) = repmat(Params.ageIndicator(:,i)',[Params.numComparts,1]) .* newPop * Params.entryRateApplied(i);
    end
    
    % sums dim 3 (5 age groups) to get back to 30x273 from 30x273x5
    newPop = sum(newPopAge,3);
    
    Compartments = Compartments + newPop * Params.tt_tstep;

end