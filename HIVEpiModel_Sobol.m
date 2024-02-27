function [Outcome,Params] = HIVEpiModel_Sobol(ExcelValues_AllParameters,ExcelValues_Populations)
%% Purpose: Main file for conducting analyses using the HIV epidemic model

%% 1. Import inputs and set up key runtime variables   

    % 1.i. Reset seed
     mystream = RandStream('mt19937ar','Seed',sum(100*clock));
     RandStream.setGlobalStream(mystream);

    % 1.ii. Populate parameter set.
    [Params] = InitializeParameters_Modified(ExcelValues_AllParameters,ExcelValues_Populations);

    % 1.iii. Set up key runtime variables
    Epi1CalLHS2CalOpt3RAOpt4SA5EE6UA7 = Params.ModelRunType; 
    runtime = Params.tt_timeHorizon; % Length of run time in years
    tstep =  Params.tt_tstep; % Time step in years

%% 2. Code specific to model run type
   
    % 2.i. Single Epi Model
    % 2.ii. Calibration using LHS
    % 2.iii. Calibration using optimization
    % 2.iv. Resource allocation optimization
    % 2.v. One-way SA (single epi model)    
    % 2.vi. Morris SA
    % 2.vii. Uncertainty
    
    switch Epi1CalLHS2CalOpt3RAOpt4SA5EE6UA7
        
        % 2.i. Single model run
        case 1 
            
            % Run Epi Model
            InitOutcome = CalcEpiOutcomes();
            
            % If collecting calibration outcomes
            if Params.OutcomesToCollect == 3
                
                % Initialize calib targets
                [CalibTargets] = Calib_initializeTargets(ExcelFileName);
                
                CalibOutputStruct = [];
                CalibOutputStruct = Calib_collectResults(InitOutcome,CalibOutputStruct, 1);
                
                [~, OOBPenaltyObjValue_Calib, ~] = ...
                    Calib_calcObjValue(CalibOutputStruct, CalibTargets, 1);
                
                InitOutcome.OOBobjvalue = OOBPenaltyObjValue_Calib;                
                
                %Calculate number of targets out of bounds
                outcomefield = fieldnames(CalibOutputStruct);
                TargetField = fieldnames(CalibTargets);
                nOutcomes = length(outcomefield);              
                nTargets = length(TargetField);
                
                % Pre-allocate
                LB = 1;
                UB = 2;
                targets_LB(nTargets,1)=0;
                targets_UB(nTargets,1)=0;
                modelOutcomes(nOutcomes,1)=0;
                
                
                for nTarget = 1:nTargets
                    targets_LB(nTarget,1)=CalibTargets.(TargetField{nTarget}).range(LB);
                    targets_UB(nTarget,1)=CalibTargets.(TargetField{nTarget}).range(UB);
                end
                
                for nOutcome=1:nOutcomes
                    modelOutcomes(nOutcome,1)=CalibOutputStruct.(outcomefield{nOutcome}).paramValue;
                end
                
                nOutsideRange = 0;
                for nTarget = 1:nTargets
                    if and(CalibTargets.(TargetField{nTarget}).weight > 0, ...
                        or(modelOutcomes(nTarget,1)<targets_LB(nTarget,1), ...
                        modelOutcomes(nTarget,1)>targets_UB(nTarget,1)))
                            nOutsideRange = nOutsideRange + 1;
                    end
                end
                
                InitOutcome.numTargetsOOB = nOutsideRange;
                
            end
            
            % Pull outcome list based on user selection
            Outcome = PickOutcomeList(InitOutcome, Params.OutcomesToCollect);

        %-----------------------------------------------------------------
        % 2.ii. Calibration using LHS
        case 2
            global LHSrun_counter calibLHS_results
            LHSrun_counter = 0;
            calibLHS_results = [];
            % Preallocate
            CalibOutputStruct = [];
                
            % 2.ii.1 Initialize Calibration Parameters
            [CalibParams,nSamples, ~,calibObjValueType] = Calib_initializeParams(ExcelFileName, Params);
            
            % 2.ii.2 Generate the hypercube sets
            [CalibParams] = Calib_genHypercubeSets(nSamples, CalibParams);
            
            % 2.ii.3 Initialize Targets
            [CalibTargets] = Calib_initializeTargets(ExcelFileName);
            
            % 2.ii.4 Loop: run the model using each generated parameter set                       
            
            for i = 1:nSamples
                
                % 2.ii.3.a Update parameters
                [Params, CalibParams] = Calib_updateParams(Params, CalibParams, i);  
                
                % 2.ii.3.b Run Epi model
                  
                % Run the model
                Outcome = CalcEpiOutcomes();                

                % 2.ii.3.c Store key outputs for calibration    
                CalibOutputStruct = Calib_collectResults(Outcome,CalibOutputStruct, 1);
                
                % 2.ii.3.d Collect Objective Value
                [ObjValue_Calib, OOBPenaltyObjValue_Calib, TargetErrorObjValue_Calib] = ...
                    Calib_calcObjValue(CalibOutputStruct, CalibTargets, calibObjValueType);
                
                % 2.ii.3.e Reset the system and back up results
                
                field = fieldnames(CalibParams);
                outcomefield = fieldnames(CalibOutputStruct);
                TargetField = fieldnames(CalibTargets);
                nCalibParams = length(field);
                nOutcomes = length(outcomefield);              
                nTargets = length(TargetField);
                
                % Pre-allocate
                LB = 1;
                UB = 2;
                targets_LB(nTargets,1)=0;
                targets_UB(nTargets,1)=0;
                modelParams(nCalibParams,1)=0;
                modelOutcomes(nOutcomes,1)=0;

                
                for nParam=1:nCalibParams
                    modelParams(nParam,1)=CalibParams.(field{nParam}).paramValue(i);
                end
                
                %Calculate number of targets out of bounds
                
                for nTarget = 1:nTargets
                    targets_LB(nTarget,1)=CalibTargets.(TargetField{nTarget}).range(LB);
                    targets_UB(nTarget,1)=CalibTargets.(TargetField{nTarget}).range(UB);
                end
                
                for nOutcome=1:nOutcomes
                    modelOutcomes(nOutcome,1)=CalibOutputStruct.(outcomefield{nOutcome}).paramValue;
                end
                
                nOutsideRange = 0;
                for nTarget = 1:nTargets
                    if and(CalibTargets.(TargetField{nTarget}).weight > 0, ...
                        or(modelOutcomes(nTarget,1)<targets_LB(nTarget,1), ...
                        modelOutcomes(nTarget,1)>targets_UB(nTarget,1)))
                            nOutsideRange = nOutsideRange + 1;
                    end
                end
                
                
                LHSrun_counter = LHSrun_counter + 1;
                temp1 = vertcat(ObjValue_Calib,OOBPenaltyObjValue_Calib,TargetErrorObjValue_Calib,nOutsideRange,modelParams,modelOutcomes);
                calibLHS_results = horzcat(calibLHS_results, temp1);
                save('calibLHS_results','calibLHS_results');
                

                fprintf('\nRun: %.0d', LHSrun_counter)
                if ObjValue_Calib == OOBPenaltyObjValue_Calib
                    fprintf('\nOOB obj value: %.10f', OOBPenaltyObjValue_Calib)
                    fprintf('\nTE obj value: %.10f', TargetErrorObjValue_Calib)
                else
                    fprintf('\nTE obj value: %.10f', TargetErrorObjValue_Calib)
                    fprintf('\nOOB obj value: %.10f', OOBPenaltyObjValue_Calib)
                end
                fprintf('\nTargets OOB: %.0d\n\n', nOutsideRange)
               
                %Reset
                clear Outcome
                
                
            end   

            
            Outcome.inputs = CalibParams;

            
        %-----------------------------------------------------------------
        % 2.iii. Calibration using Optimization
        case 3 
            
%             %temp to delete CSV file where results are stored
%             tic
%             if exist('calib_results.dat') >= 1
%                     delete('calib_results.dat');
%             end

            % alternate save method and global variable intilization
            % tic is for time. run counter and calib_results are for saving
            % outcomes.
            tic
            global run_counter calib_results calibOVStop runs_below_threshold calibNumBelowStop
            run_counter = 0;
            calib_results = [];
            calibOVStop = 0;
            runs_below_threshold = 0;
            calibNumBelowStop = 0;
            
            % 2.iii.1 Initialize Calibration Parameters
            [CalibParams, ~, calibMethod, calibObjValueType,initValueSource, calibStopCondition, calibTimeToRun, nRandomStartRuns, calibOVStop, calibNumBelowStop, CalibParams_reduced] = Calib_initializeParams(ExcelFileName, Params); 
            [CalibTargets] = Calib_initializeTargets(ExcelFileName);
            ParamFields = fieldnames(CalibParams);
            nCalibratedParameters = length(ParamFields);            
                        
            % 2.iii.2 Generate inputs for optimization model

                % fun: function being minimized by fmincon. fun is a function that accepts a vector x and returns a scalar f, the objective function evaluated at x. 
                %        fun can be specified as a function handle for a file: x = fmincon(@myfun,x0,A,b), where myfun is a MATLAB® function such as 
                %    function ObjValue = CalcObjValue(DVAllocation)

                % x0: initial solution; can be a scalar, vector, or matrix.
                % A / b: defines linear inequalities A*x <= b. If no inequalities exist, set A = [] and b = []. If A large w/ few nonzero entries, make sparse.
                % Aeq / beq: defines linear equalities Aeq*x = beq. If no equalities exist, set Aeq = [] and beq = []
                % lb / ub: lower and upper bounds on the design variables in x, so that the solution is always in the range lb <= x <= ub. If no bounds exist, set lb = [] and/or ub = [].If x(i) is unbounded below, set lb(i) = -Inf, and if x(i) is unbounded above, set ub(i) = Inf.
                % nonlcon: defines nonlinear inequalities c(x) or equalities ceq(x). fmincon optimizes such that c(x) <= 0 and ceq(x) = 0. accepts a vector x and returns the two vectors c and ceq. c is a vector that contains the nonlinear inequalities evaluated at x, and ceq is a vector that contains the nonlinear equalities evaluated at x. nonlcon should be specified as a function handle to a file or to an anonymous function, such as mycon: x = fmincon(@myfun,x0,A,b,Aeq,beq,lb,ub,@mycon), where mycon is a MATLAB function such as 
                %       function [c,ceq] = mycon(x)
                %       c = ...     % Compute nonlinear inequalities at x.
                %       ceq = ...   % Compute nonlinear equalities at x.

            % 2.iii.3 Run optimization
            %[OptAllocation,fval,exitflag,output] = fmincon(@CalcObjValue,x0,A,b,Aeq,beq,lb,ub,@nonlcon,options); 

            if initValueSource == 1
                nRandomStartRuns = 1;
            end
            for n=1:nRandomStartRuns    
                 %Fmincon
                [x0,A,b,Aeq,beq,lb,ub,options] = Calib_GenOptComponents(CalibParams_reduced, initValueSource, 3, calibStopCondition, calibTimeToRun);
                [OptParSet, fval, exitflag, output] = fmincon(@CalcObjValue_Calib,x0,A,b,Aeq,beq,lb,ub,@nonlcon,options);

                % matches OptParSet from a completed optimization to the full solution vector
                % that we use inside 'calib_results'. Then append to end of 'calib_results'
                % JC updated on 06/08/2020 so that
                % temp.allCalibInputVectors is calculated based on
                % total array of calib params (versus referring to
                % length of OptParSet, which may be less than the total
                % array if only calibrating values for 1 time period).
                % KH commented out code below b/c it was causing errors
%                 load('calib_results');
%                 temp.allCalibInputVectors = calib_results(Params.calib_num_outputVectorErrorValues:(Params.calib_num_outputVectorErrorValues + nCalibratedParameters - 1),:);
%                 [temp.NA, temp.index] = ismember(OptParSet',temp.allCalibInputVectors', 'rows'); % have to transpose to use 'rows' option
%                 calib_results =  horzcat(calib_results,calib_results(:,temp.index));
%                 save('calib_results', 'calib_results');
%                 fprintf('added best run to last column in calib_results')
%                 clearvars temp
            end

            % 2.iii.3 Assign key outcomes

                % OptParSet: value of decision variable
                % fval: value of the objective function fun at the solution x.
                % exitflag: Integer identifying the reason the algorithm terminated.
                % output: Structure containing information about the optimization

            Outcome.OptParSet = OptParSet;
            Outcome.fval = fval;
            Outcome.exitflag = exitflag;
            Outcome.output = output;
            
        %-----------------------------------------------------------------
        % 2.iv. Resource allocation optimization run
        case 4 
            
            % tic is for time. run counter and calib_results are for saving
            % outcomes.
            tic
            global run_counter RA_results
            run_counter = 0;
            RA_results = [];
            
            if or(Params.SolveCont~=0,Params.tt_progressionSource~=3)
                Outcome.output = ...
                    ['Optimization runs require allocation-based progression '...
                    'and discrete difference equations. Please reset in the '...
                    'Settings sheet and try again.'];
            else
    
                % 2.iv.1 Generate inputs for optimization model

                    % fun: function being minimized by fmincon. fun is a function that accepts a vector x and returns a scalar f, the objective function evaluated at x. 
                    %        fun can be specified as a function handle for a file: x = fmincon(@myfun,x0,A,b), where myfun is a MATLAB® function such as 
                    %    function ObjValue_RA = CalcObjValue_RA(DVAllocation)

                    % x0: initial solution; can be a scalar, vector, or matrix.
                    % A / b: defines linear inequalities A*x <= b. If no inequalities exist, set A = [] and b = []. If A large w/ few nonzero entries, make sparse.
                    % Aeq / beq: defines linear equalities Aeq*x = beq. If no equalities exist, set Aeq = [] and beq = []
                    % lb / ub: lower and upper bounds on the design variables in x, so that the solution is always in the range lb <= x <= ub. If no bounds exist, set lb = [] and/or ub = [].If x(i) is unbounded below, set lb(i) = -Inf, and if x(i) is unbounded above, set ub(i) = Inf.
                    % nonlcon: defines nonlinear inequalities c(x) or equalities ceq(x). fmincon optimizes such that c(x) <= 0 and ceq(x) = 0. accepts a vector x and returns the two vectors c and ceq. c is a vector that contains the nonlinear inequalities evaluated at x, and ceq is a vector that contains the nonlinear equalities evaluated at x. nonlcon should be specified as a function handle to a file or to an anonymous function, such as mycon: x = fmincon(@myfun,x0,A,b,Aeq,beq,lb,ub,@mycon), where mycon is a MATLAB function such as 
                    %       function [c,ceq] = mycon(x)
                    %       c = ...     % Compute nonlinear inequalities at x.
                    %       ceq = ...   % Compute nonlinear equalities at x.

                [x0,A,b,Aeq,beq,lb,ub,options] = GenOptComponents(Params); 
                format long;

                % 2.iv.2 Run optimization
                if Params.optObjNum == 4 %Minimize intervention cost to hit incidence targets

                    FirstOutcomeYr = max(Params.tt_modelStartYear,Params.outcomeCollectionStartYr);
                    LastOutcomeYr = min(Params.tt_modelStartYear+Params.tt_timeHorizon-1,Params.outcomeCollectionEndYr);
                    TargetBaseYr = Params.tt_periodFiveStartYear - ...
                            Params.OutcomesIndex_AllocationStartYr(1)+1;
                    SecondTargetYr = Params.tt_periodFiveStartYear + ...
                            Params.minCost_Tgt2NumYears - 1;
  
%                     KH added this if statement on 30Sept2021. It needs debugging before implementing.
%                     if (FirstOutcomeYr <= TargetBaseYr) && (LastOutcomeYr >= SecondTargetYr)
% 
%                         h = msgbox('You only collected outcomes from ' && ...
%                             num2str(FirstOutcomeYr) && ...
%                             ' to ' && num2str(LastOutcomeYr) && ...
%                             '. Please revise to include all years that are used to calculate the incidence targets (' &&...
%                             num2str(TargetBaseYr) && ' to ' &&  ...
%                             num2str(SecondTargetYr) && ').');
%                         
%                     else
                        [OptAllocation, fval, exitflag, output] = ...
                            fmincon(@CalcObjValue_RA, x0, A, b, Aeq, beq, lb, ub, @nonlcon_hitIncidenceTargets,options); 
%                     end                    
                    
                elseif Params.alloc_ConsiderARTAndCareCosts == 1
                    [OptAllocation, fval, exitflag, output] = ...
                        fmincon(@CalcObjValue_RA, x0, A, b, Aeq, beq, lb, ub, @nonlcon_withARTAndCareCosts,options); 
                    
                    if exist('DVAllocation','var')>0
                        OptAllocation(Params.numIntns+Params.alloc_ConsiderARTAndCareCosts)= ...
                            DVAllocation(Params.numIntns+Params.alloc_ConsiderARTAndCareCosts);
                    else
                        OptAllocation(Params.numIntns+Params.alloc_ConsiderARTAndCareCosts)= ...
                            Params.alloc_TotalIntnBudget - ...
                            sum(OptAllocation(1:Params.numIntns));
                    end
                    
                else
                     [OptAllocation, fval, exitflag, output] = ...
                        fmincon(@CalcObjValue_RA, x0, A, b, Aeq, beq, lb, ub, @nonlcon, options); 
                end


                % 2.iii.3 Assign key outcomes

                    % OptAllocation: value of decision variable
                    % fval: value of the objective function fun at the solution x.
                    % exitflag: Integer identifying the reason the algorithm terminated.
                    % output: Structure containing information about the optimization

                Outcome.OptAllocation = OptAllocation;
                Outcome.fval = fval;
                Outcome.exitflag = exitflag;
                Outcome.output = output;
            end
        %-----------------------------------------------------------------
        % 2.v.  One-way sensitivity analysis run
        case 5  
            
            % 2.v.1 Run Epi Model
            

                % Run model
                AllOutcomes = CalcEpiOutcomes();

                
            % 2.v.2 Set key outcome
                % Outcome is a the single number to compare via tornado
                % diagram
               
                    % Note the order of the outcome numbers should be
                    % consistent with the order of the drop-down list on
                    % the One-WAY SA sheet in Excel.
                
                % QALYs
                if Params.outcomeNumber == 1 % QALYS
                    Outcome = AllOutcomes.total_QALYs;
                    
                % Life-years
                elseif Params.outcomeNumber == 2
                    Outcome = AllOutcomes.total_LifeYears;
                
                % All Costs
                elseif Params.outcomeNumber == 3 
                    Outcome = AllOutcomes.total_ARTCarePrEPTransitionAndSEPCost_Disc;
                
                % QALYs and Costs
                    % Note: this is a special case where the QALYs are
                    % concatenated with the total costs. If you make
                    % changes to this, please also update the VBA code in
                    % the corresponding Excel file.
                elseif Params.outcomeNumber == 4 
                    Outcome = strcat(num2str(AllOutcomes.total_QALYs),'-', num2str(AllOutcomes.total_ARTCarePrEPTransitionAndSEPCost_Disc));
                
                % New infections in 2020
                elseif Params.outcomeNumber == 5
                    Outcome = sum(AllOutcomes.ann_TotalNewInfections(Params.year_2020,:),2);
                    
                % Black new infections in 2020
                elseif Params.outcomeNumber == 6
                    Outcome = sum(AllOutcomes.ann_NewInfections_Blk(Params.year_2020,:),2);
                    
                % Hispanic new infections in 2020
                elseif Params.outcomeNumber == 7
                    Outcome = sum(AllOutcomes.ann_NewInfections_Hisp(Params.year_2020,:),2);
                
                % Other new infections in 2020
                elseif Params.outcomeNumber == 8
                    Outcome = sum(AllOutcomes.ann_NewInfections_Oth(Params.year_2020,:),2);
                    
                % 2015 new infections among Black MSM 18-24
                elseif Params.outcomeNumber == 9
                    Outcome = AllOutcomes.calib_NewInfections_2015_MSM_18_24_B;
                    
                % 2015 new infections among Hispanic MSM 18-24
                elseif Params.outcomeNumber == 10
                    Outcome = AllOutcomes.calib_NewInfections_2015_MSM_18_24_H;
                    
                % 2015 new infections among Other MSM 18-24
                elseif Params.outcomeNumber == 11
                    Outcome = AllOutcomes.calib_NewInfections_2015_MSM_18_24_O;
                    
                % 2015 new infections among Black MSM 25-34
                elseif Params.outcomeNumber == 12
                    Outcome = AllOutcomes.calib_NewInfections_2015_MSM_25_34_B;
                    
                % 2015 new infections among Hispanic MSM 25-34
                elseif Params.outcomeNumber == 13
                    Outcome = AllOutcomes.calib_NewInfections_2015_MSM_25_34_H;
                    
                % 2015 new infections among Other MSM 25-34
                elseif Params.outcomeNumber == 14
                    Outcome = AllOutcomes.calib_NewInfections_2015_MSM_25_34_O;
                    
                % MSM 18-34 new infections 2030
                elseif Params.outcomeNumber == 15
                    Outcome = sum(AllOutcomes.ann_TotalNewInfections(Params.year_2015+15, [16 17 23 24 51 52 58 59 107 108 114 115 142 143 149 150 198 199 205 206 233 234 240 241]));
                    
                % MSM 18-34 cumulative new infections 2018-2030
                elseif Params.outcomeNumber == 16
                    Outcome = sum(sum(AllOutcomes.ann_TotalNewInfections(Params.year_2015+3:Params.year_2015+15, [16 17 23 24 51 52 58 59 107 108 114 115 142 143 149 150 198 199 205 206 233 234 240 241])'));                   
                    
                % QALYs and Costs (for ICER)
                % Note: this is a special case where the QALYs are
                    % concatenated with the total costs. If you make
                    % changes to this, please also update the VBA code in
                    % the corresponding Excel file.
               
                % Cumulative HIV incidence - Clinkscales 7/27/2021 (Edited
                % by Laurel)
                elseif Params.outcomeNumber == 17
                    Outcome = AllOutcomes.TotalNewInfections_Disc;
                else
                    Outcome = AllOutcomes;
                
                end            
            
        %-----------------------------------------------------------------
        % 2.vi. Morris/Elementary Effects SA run
        case 6
        
            % 2.vi.1 Run Epi Model
            
                % Run model
                AllOutcomes = CalcEpiOutcomes();

                
            % 2.vi.2 Set key outcome
               
                    % Note the order of the outcome numbers should be
                    % consistent with the order of the drop-down list on
                    % the Morris SA sheet in Excel.
                
                % QALYs
                if Params.outcomeNumber == 1 % QALYS
                    Outcome = AllOutcomes.total_QALYs;
                    
                % Life-years
                elseif Params.outcomeNumber == 2
                    Outcome = AllOutcomes.total_LifeYears;
                
                % All Costs
                elseif Params.outcomeNumber == 3 
                    Outcome = AllOutcomes.total_ARTCarePrEPTransitionAndSEPCost_Disc;
                
                % QALYs and Costs
                    % Note: this is a special case where the QALYs are
                    % concatenated with the total costs. If you make
                    % changes to this, please also update the VBA code in
                    % the corresponding Excel file.
                elseif Params.outcomeNumber == 4 
                    Outcome = strcat(num2str(AllOutcomes.total_QALYs),'-', num2str(AllOutcomes.total_ARTCarePrEPTransitionAndSEPCost_Disc));
                
                % New infections in 2020
                elseif Params.outcomeNumber == 5
                    Outcome = sum(AllOutcomes.ann_TotalNewInfections(Params.year_2020,:),2);
                    
                % Black new infections in 2020
                elseif Params.outcomeNumber == 6
                    Outcome = sum(AllOutcomes.ann_NewInfections_Blk(Params.year_2020,:),2);
                
                % Hispanic new infections in 2020
                elseif Params.outcomeNumber == 7
                    Outcome = sum(AllOutcomes.ann_NewInfections_Hisp(Params.year_2020,:),2);
                
                % Other new infections in 2020
                elseif Params.outcomeNumber == 8
                    Outcome = sum(AllOutcomes.ann_NewInfections_Oth(Params.year_2020,:),2);
                
                % 2015 new infections among Black MSM 18-24
                elseif Params.outcomeNumber == 9
                    Outcome = AllOutcomes.calib_NewInfections_2015_MSM_18_24_B;
                    
                % 2015 new infections among Hispanic MSM 18-24
                elseif Params.outcomeNumber == 10
                    Outcome = AllOutcomes.calib_NewInfections_2015_MSM_18_24_H;
                    
                % 2015 new infections among Other MSM 18-24
                elseif Params.outcomeNumber == 11
                    Outcome = AllOutcomes.calib_NewInfections_2015_MSM_18_24_O;
                    
                % 2015 new infections among Black MSM 25-34
                elseif Params.outcomeNumber == 12
                    Outcome = AllOutcomes.calib_NewInfections_2015_MSM_25_34_B;
                    
                % 2015 new infections among Hispanic MSM 25-34
                elseif Params.outcomeNumber == 13
                    Outcome = AllOutcomes.calib_NewInfections_2015_MSM_25_34_H;
                    
                % 2015 new infections among Other MSM 25-34
                elseif Params.outcomeNumber == 14
                    Outcome = AllOutcomes.calib_NewInfections_2015_MSM_25_34_O;
                    
                 % MSM 18-34 new infections 2030
                elseif Params.outcomeNumber == 15
                    Outcome = sum(AllOutcomes.ann_TotalNewInfections(Params.year_2015+15, [16 17 23 24 51 52 58 59 107 108 114 115 142 143 149 150 198 199 205 206 233 234 240 241]));
                    
                % MSM 18-34 cumulative new infections 2018-2030
                elseif Params.outcomeNumber == 16
                    Outcome = sum(sum(AllOutcomes.ann_TotalNewInfections(Params.year_2015+3:Params.year_2015+15, [16 17 23 24 51 52 58 59 107 108 114 115 142 143 149 150 198 199 205 206 233 234 240 241])'));                   
                                       
                % QALYs and Costs (for ICER)
                % Note: this is a special case where the QALYs are
                    % concatenated with the total costs. If you make
                    % changes to this, please also update the VBA code in
                    % the corresponding Excel file.
                    
                % Cumulative HIV incidence - Clinkscales 7/27/2021. Edited
                % by Laurel
                elseif Params.outcomeNumber == 17
                    Outcome = AllOutcomes.TotalNewInfections_Disc;
                    
                else
                    Outcome = AllOutcomes;
                
                end
        
        
    end

%% 3. CalcEpiOutcomes function: Run Epi model 
    
    % Note: this is a nested function within the main HIVEpiModel.m
    % function
      
    function EpiOutcomes = CalcEpiOutcomes(Opt_DV)

        % If this isn't an optimization, no decision variable passed in
         if ~exist('Opt_DV','var')
              Opt_DV = [];
         end

        % 3.i. Initialize key parameters

        % Clear out any existing global variables besides params
            if exist('Results','var') 
                clear Results;
            end
            if exist('Compartments','var') 
                clear Compartments;
            end  
            if exist('TransRates','var') 
                clear TransRates;
            end    
            if exist('Transitions','var') 
                clear Transitions;
            end    
            if exist('TTProg','var') 
                clear TTProg;
            end     
            if exist('InfRate','var') 
                clear InfRate;
            end                

            Results_p = intializeResults(Params);

        % Pre-allocate

            TransRates.AnnualRates(Params.numComparts,Params.numComparts,Params.numStrats) = 0.000;
            Transitions(Params.numComparts,Params.numComparts,Params.numStrats)=0;
            TTProg.All = 0;

        % Initialize variables at the beginning of the model
            Compartments = Params.initPop;

            Year = Params.tt_modelStartYear;


        % Assign optimization DVs to correct variables in Params
        switch Epi1CalLHS2CalOpt3RAOpt4SA5EE6UA7

            case 3 % Calib optimization
                % when I (CG) added code to only run non-zero range calib
                % params, during an optimization calibration, 'Opt_DV'
                % is actually the reduced version and needs to be rebuilt
                % here            
                fNames_reduced = fieldnames(CalibParams_reduced);                          

                for j = 1:length(fNames_reduced)
                   update_params.(fNames_reduced{j}) = Opt_DV(j); 
                end
                for j = 1:length(fNames_reduced)
                    CalibParams.(fNames_reduced{j}).paramValue = update_params.(fNames_reduced{j});
                end
                [Params, CalibParams] = Calib_updateParams(Params, CalibParams, 1);

            case 4 % RA optimization
                % Establish number of DVs per allocation period
                if Params.TargetIntnsToYMSM == 1
                    
                    nNumDVsPerAlloc = Params.numIntns_YMSM+Params.alloc_ConsiderARTAndCareCosts;
                
                    
                    Params.intn_Testing_Investment_MSM_Low(1) = Opt_DV(1);
                    Params.intn_Testing_Investment_MSM_High(1) = Opt_DV(2);                    
                    Params.intn_LTCatDiag_Investment(1) = Opt_DV(3);
                    Params.intn_LTCafterDiag_Investment(1) = Opt_DV(4);
                    Params.intn_ARTInitiation_Investment(1) = Opt_DV(5);
                    Params.intn_ARTAdher5to4_Investment(1) = Opt_DV(6);
                    Params.intn_ARTAdher4to5_Investment(1) = Opt_DV(7);                   
                    Params.intn_PrEP_Oral_Investment_MSM_B(1) = Opt_DV(8);
                    Params.intn_PrEP_Oral_Investment_MSM_H(1) = Opt_DV(9);
                    Params.intn_PrEP_Oral_Investment_MSM_O(1) = Opt_DV(10);
                    Params.intn_PrEP_Inject_Investment_MSM_B(1) = Opt_DV(11);
                    Params.intn_PrEP_Inject_Investment_MSM_H(1) = Opt_DV(12);
                    Params.intn_PrEP_Inject_Investment_MSM_O(1) = Opt_DV(13);
                    % Don't need to assign Opt_DV(15) since it's the
                    % treatment and care costs output
                    if Params.numAllocationPeriods >=2                        
                        Params.intn_Testing_Investment_MSM_Low(2) = Opt_DV(nNumDVsPerAlloc+1);
                        Params.intn_Testing_Investment_MSM_High(2) = Opt_DV(nNumDVsPerAlloc+2);                        
                        Params.intn_LTCatDiag_Investment(2) = Opt_DV(nNumDVsPerAlloc+3);
                        Params.intn_LTCafterDiag_Investment(2) = Opt_DV(nNumDVsPerAlloc+4);
                        Params.intn_ARTInitiation_Investment(2) = Opt_DV(nNumDVsPerAlloc+5);
                        Params.intn_ARTAdher5to4_Investment(2) = Opt_DV(nNumDVsPerAlloc+6);
                        Params.intn_ARTAdher4to5_Investment(2) = Opt_DV(nNumDVsPerAlloc+7);                        
                        Params.intn_PrEP_Oral_Investment_MSM_B(2) = Opt_DV(nNumDVsPerAlloc+8);
                        Params.intn_PrEP_Oral_Investment_MSM_H(2) = Opt_DV(nNumDVsPerAlloc+9);
                        Params.intn_PrEP_Oral_Investment_MSM_O(2) = Opt_DV(nNumDVsPerAlloc+10);
                        Params.intn_PrEP_Inject_Investment_MSM_B(2) = Opt_DV(nNumDVsPerAlloc+11);
                        Params.intn_PrEP_Inject_Investment_MSM_H(2) = Opt_DV(nNumDVsPerAlloc+12);
                        Params.intn_PrEP_Inject_Investment_MSM_O(2) = Opt_DV(nNumDVsPerAlloc+13);
                    end
                    if Params.numAllocationPeriods ==3                        
                        Params.intn_Testing_Investment_MSM_Low(3) = Opt_DV(2*(nNumDVsPerAlloc)+1);
                        Params.intn_Testing_Investment_MSM_High(3) = Opt_DV(2*(nNumDVsPerAlloc)+2);                        
                        Params.intn_LTCatDiag_Investment(3) = Opt_DV(2*(nNumDVsPerAlloc)+3);
                        Params.intn_LTCafterDiag_Investment(3) = Opt_DV(2*(nNumDVsPerAlloc)+4);
                        Params.intn_ARTInitiation_Investment(3) = Opt_DV(2*(nNumDVsPerAlloc)+5);
                        Params.intn_ARTAdher5to4_Investment(3) = Opt_DV(2*(nNumDVsPerAlloc)+6);
                        Params.intn_ARTAdher4to5_Investment(3) = Opt_DV(2*(nNumDVsPerAlloc)+7);                       
                        Params.intn_PrEP_Oral_Investment_MSM_B(3) = Opt_DV(2*(nNumDVsPerAlloc)+8);
                        Params.intn_PrEP_Oral_Investment_MSM_H(3) = Opt_DV(2*(nNumDVsPerAlloc)+9);
                        Params.intn_PrEP_Oral_Investment_MSM_O(3) = Opt_DV(2*(nNumDVsPerAlloc)+10);
                        Params.intn_PrEP_Inject_Investment_MSM_B(3) = Opt_DV(2*(nNumDVsPerAlloc)+11);
                        Params.intn_PrEP_Inject_Investment_MSM_H(3) = Opt_DV(2*(nNumDVsPerAlloc)+12);
                        Params.intn_PrEP_Inject_Investment_MSM_O(3) = Opt_DV(2*(nNumDVsPerAlloc)+13);
                    end
                    
                else
                    
                    nNumDVsPerAlloc = Params.numIntns+Params.alloc_ConsiderARTAndCareCosts;

                    Params.intn_Testing_Investment_HET_Low(1) = Opt_DV(1);
                    Params.intn_Testing_Investment_HET_High(1) = Opt_DV(2);
                    Params.intn_Testing_Investment_MSM_Low(1) = Opt_DV(3);
                    Params.intn_Testing_Investment_MSM_High(1) = Opt_DV(4);
                    Params.intn_Testing_Investment_IDU(1) = Opt_DV(5);
                    Params.intn_LTCatDiag_Investment(1) = Opt_DV(6);
                    Params.intn_LTCafterDiag_Investment(1) = Opt_DV(7);
                    Params.intn_ARTInitiation_Investment(1) = Opt_DV(8);
                    Params.intn_ARTAdher5to4_Investment(1) = Opt_DV(9);
                    Params.intn_ARTAdher4to5_Investment(1) = Opt_DV(10);
                    Params.intn_SEP_Investment_B(1) = Opt_DV(11);
                    Params.intn_SEP_Investment_H(1) = Opt_DV(12);
                    Params.intn_SEP_Investment_O(1) = Opt_DV(13);
                    Params.intn_PrEP_Oral_Investment_HETM_B(1) = Opt_DV(14);
                    Params.intn_PrEP_Oral_Investment_HETM_H(1) = Opt_DV(15);
                    Params.intn_PrEP_Oral_Investment_HETM_O(1) = Opt_DV(16);
                    Params.intn_PrEP_Oral_Investment_HETF_B(1) = Opt_DV(17);
                    Params.intn_PrEP_Oral_Investment_HETF_H(1) = Opt_DV(18);
                    Params.intn_PrEP_Oral_Investment_HETF_O(1) = Opt_DV(19);
                    Params.intn_PrEP_Oral_Investment_MSM_B(1) = Opt_DV(20);
                    Params.intn_PrEP_Oral_Investment_MSM_H(1) = Opt_DV(21);
                    Params.intn_PrEP_Oral_Investment_MSM_O(1) = Opt_DV(22);
                    Params.intn_PrEP_Oral_Investment_IDU_B(1) = Opt_DV(23);
                    Params.intn_PrEP_Oral_Investment_IDU_H(1) = Opt_DV(24);
                    Params.intn_PrEP_Oral_Investment_IDU_O(1) = Opt_DV(25);
                    Params.intn_PrEP_Inject_Investment_HETM_B(1) = Opt_DV(26);
                    Params.intn_PrEP_Inject_Investment_HETM_H(1) = Opt_DV(27);
                    Params.intn_PrEP_Inject_Investment_HETM_O(1) = Opt_DV(28);
                    Params.intn_PrEP_Inject_Investment_HETF_B(1) = Opt_DV(29);
                    Params.intn_PrEP_Inject_Investment_HETF_H(1) = Opt_DV(30);
                    Params.intn_PrEP_Inject_Investment_HETF_O(1) = Opt_DV(31);
                    Params.intn_PrEP_Inject_Investment_MSM_B(1) = Opt_DV(32);
                    Params.intn_PrEP_Inject_Investment_MSM_H(1) = Opt_DV(33);
                    Params.intn_PrEP_Inject_Investment_MSM_O(1) = Opt_DV(34);
                    Params.intn_PrEP_Inject_Investment_IDU_B(1) = Opt_DV(35);
                    Params.intn_PrEP_Inject_Investment_IDU_H(1) = Opt_DV(36);
                    Params.intn_PrEP_Inject_Investment_IDU_O(1) = Opt_DV(37);
                    % Don't need to assign Opt_DV(15) since it's the
                    % treatment and care costs output
                    if Params.numAllocationPeriods >=2
                        Params.intn_Testing_Investment_HET_Low(2) = Opt_DV(nNumDVsPerAlloc+1);
                        Params.intn_Testing_Investment_HET_High(2) = Opt_DV(nNumDVsPerAlloc+2);
                        Params.intn_Testing_Investment_MSM_Low(2) = Opt_DV(nNumDVsPerAlloc+3);
                        Params.intn_Testing_Investment_MSM_High(2) = Opt_DV(nNumDVsPerAlloc+4);
                        Params.intn_Testing_Investment_IDU(2) = Opt_DV(nNumDVsPerAlloc+5);
                        Params.intn_LTCatDiag_Investment(2) = Opt_DV(nNumDVsPerAlloc+6);
                        Params.intn_LTCafterDiag_Investment(2) = Opt_DV(nNumDVsPerAlloc+7);
                        Params.intn_ARTInitiation_Investment(2) = Opt_DV(nNumDVsPerAlloc+8);
                        Params.intn_ARTAdher5to4_Investment(2) = Opt_DV(nNumDVsPerAlloc+9);
                        Params.intn_ARTAdher4to5_Investment(2) = Opt_DV(nNumDVsPerAlloc+10);
                        Params.intn_SEP_Investment_B(2) = Opt_DV(nNumDVsPerAlloc+11);
                        Params.intn_SEP_Investment_H(2) = Opt_DV(nNumDVsPerAlloc+12);
                        Params.intn_SEP_Investment_O(2) = Opt_DV(nNumDVsPerAlloc+13);
                        Params.intn_PrEP_Oral_Investment_HETM_B(2) = Opt_DV(nNumDVsPerAlloc+14);
                        Params.intn_PrEP_Oral_Investment_HETM_H(2) = Opt_DV(nNumDVsPerAlloc+15);
                        Params.intn_PrEP_Oral_Investment_HETM_O(2) = Opt_DV(nNumDVsPerAlloc+16);
                        Params.intn_PrEP_Oral_Investment_HETF_B(2) = Opt_DV(nNumDVsPerAlloc+17);
                        Params.intn_PrEP_Oral_Investment_HETF_H(2) = Opt_DV(nNumDVsPerAlloc+18);
                        Params.intn_PrEP_Oral_Investment_HETF_O(2) = Opt_DV(nNumDVsPerAlloc+19);
                        Params.intn_PrEP_Oral_Investment_MSM_B(2) = Opt_DV(nNumDVsPerAlloc+20);
                        Params.intn_PrEP_Oral_Investment_MSM_H(2) = Opt_DV(nNumDVsPerAlloc+21);
                        Params.intn_PrEP_Oral_Investment_MSM_O(2) = Opt_DV(nNumDVsPerAlloc+22);
                        Params.intn_PrEP_Oral_Investment_IDU_B(2) = Opt_DV(nNumDVsPerAlloc+23);
                        Params.intn_PrEP_Oral_Investment_IDU_H(2) = Opt_DV(nNumDVsPerAlloc+24);
                        Params.intn_PrEP_Oral_Investment_IDU_O(2) = Opt_DV(nNumDVsPerAlloc+25);
                        Params.intn_PrEP_Inject_Investment_HETM_B(2) = Opt_DV(nNumDVsPerAlloc+26);
                        Params.intn_PrEP_Inject_Investment_HETM_H(2) = Opt_DV(nNumDVsPerAlloc+27);
                        Params.intn_PrEP_Inject_Investment_HETM_O(2) = Opt_DV(nNumDVsPerAlloc+28);
                        Params.intn_PrEP_Inject_Investment_HETF_B(2) = Opt_DV(nNumDVsPerAlloc+29);
                        Params.intn_PrEP_Inject_Investment_HETF_H(2) = Opt_DV(nNumDVsPerAlloc+30);
                        Params.intn_PrEP_Inject_Investment_HETF_O(2) = Opt_DV(nNumDVsPerAlloc+31);
                        Params.intn_PrEP_Inject_Investment_MSM_B(2) = Opt_DV(nNumDVsPerAlloc+32);
                        Params.intn_PrEP_Inject_Investment_MSM_H(2) = Opt_DV(nNumDVsPerAlloc+33);
                        Params.intn_PrEP_Inject_Investment_MSM_O(2) = Opt_DV(nNumDVsPerAlloc+34);
                        Params.intn_PrEP_Inject_Investment_IDU_B(2) = Opt_DV(nNumDVsPerAlloc+35);
                        Params.intn_PrEP_Inject_Investment_IDU_H(2) = Opt_DV(nNumDVsPerAlloc+36);
                        Params.intn_PrEP_Inject_Investment_IDU_O(2) = Opt_DV(nNumDVsPerAlloc+37);
                    end
                    if Params.numAllocationPeriods ==3
                        Params.intn_Testing_Investment_HET_Low(3) = Opt_DV(2*(nNumDVsPerAlloc)+1);
                        Params.intn_Testing_Investment_HET_High(3) = Opt_DV(2*(nNumDVsPerAlloc)+2);
                        Params.intn_Testing_Investment_MSM_Low(3) = Opt_DV(2*(nNumDVsPerAlloc)+3);
                        Params.intn_Testing_Investment_MSM_High(3) = Opt_DV(2*(nNumDVsPerAlloc)+4);
                        Params.intn_Testing_Investment_IDU(3) = Opt_DV(2*(nNumDVsPerAlloc)+5);
                        Params.intn_LTCatDiag_Investment(3) = Opt_DV(2*(nNumDVsPerAlloc)+6);
                        Params.intn_LTCafterDiag_Investment(3) = Opt_DV(2*(nNumDVsPerAlloc)+7);
                        Params.intn_ARTInitiation_Investment(3) = Opt_DV(2*(nNumDVsPerAlloc)+8);
                        Params.intn_ARTAdher5to4_Investment(3) = Opt_DV(2*(nNumDVsPerAlloc)+9);
                        Params.intn_ARTAdher4to5_Investment(3) = Opt_DV(2*(nNumDVsPerAlloc)+10);
                        Params.intn_SEP_Investment_B(3) = Opt_DV(2*(nNumDVsPerAlloc)+11);
                        Params.intn_SEP_Investment_H(3) = Opt_DV(2*(nNumDVsPerAlloc)+12);
                        Params.intn_SEP_Investment_O(3) = Opt_DV(2*(nNumDVsPerAlloc)+13);
                        Params.intn_PrEP_Oral_Investment_HETM_B(3) = Opt_DV(2*(nNumDVsPerAlloc)+14);
                        Params.intn_PrEP_Oral_Investment_HETM_H(3) = Opt_DV(2*(nNumDVsPerAlloc)+15);
                        Params.intn_PrEP_Oral_Investment_HETM_O(3) = Opt_DV(2*(nNumDVsPerAlloc)+16);
                        Params.intn_PrEP_Oral_Investment_HETF_B(3) = Opt_DV(2*(nNumDVsPerAlloc)+17);
                        Params.intn_PrEP_Oral_Investment_HETF_H(3) = Opt_DV(2*(nNumDVsPerAlloc)+18);
                        Params.intn_PrEP_Oral_Investment_HETF_O(3) = Opt_DV(2*(nNumDVsPerAlloc)+19);
                        Params.intn_PrEP_Oral_Investment_MSM_B(3) = Opt_DV(2*(nNumDVsPerAlloc)+20);
                        Params.intn_PrEP_Oral_Investment_MSM_H(3) = Opt_DV(2*(nNumDVsPerAlloc)+21);
                        Params.intn_PrEP_Oral_Investment_MSM_O(3) = Opt_DV(2*(nNumDVsPerAlloc)+22);
                        Params.intn_PrEP_Oral_Investment_IDU_B(3) = Opt_DV(2*(nNumDVsPerAlloc)+23);
                        Params.intn_PrEP_Oral_Investment_IDU_H(3) = Opt_DV(2*(nNumDVsPerAlloc)+24);
                        Params.intn_PrEP_Oral_Investment_IDU_O(3) = Opt_DV(2*(nNumDVsPerAlloc)+25);
                        Params.intn_PrEP_Inject_Investment_HETM_B(3) = Opt_DV(2*(nNumDVsPerAlloc)+26);
                        Params.intn_PrEP_Inject_Investment_HETM_H(3) = Opt_DV(2*(nNumDVsPerAlloc)+27);
                        Params.intn_PrEP_Inject_Investment_HETM_O(3) = Opt_DV(2*(nNumDVsPerAlloc)+28);
                        Params.intn_PrEP_Inject_Investment_HETF_B(3) = Opt_DV(2*(nNumDVsPerAlloc)+29);
                        Params.intn_PrEP_Inject_Investment_HETF_H(3) = Opt_DV(2*(nNumDVsPerAlloc)+30);
                        Params.intn_PrEP_Inject_Investment_HETF_O(3) = Opt_DV(2*(nNumDVsPerAlloc)+31);
                        Params.intn_PrEP_Inject_Investment_MSM_B(3) = Opt_DV(2*(nNumDVsPerAlloc)+32);
                        Params.intn_PrEP_Inject_Investment_MSM_H(3) = Opt_DV(2*(nNumDVsPerAlloc)+33);
                        Params.intn_PrEP_Inject_Investment_MSM_O(3) = Opt_DV(2*(nNumDVsPerAlloc)+34);
                        Params.intn_PrEP_Inject_Investment_IDU_B(3) = Opt_DV(2*(nNumDVsPerAlloc)+35);
                        Params.intn_PrEP_Inject_Investment_IDU_H(3) = Opt_DV(2*(nNumDVsPerAlloc)+36);
                        Params.intn_PrEP_Inject_Investment_IDU_O(3) = Opt_DV(2*(nNumDVsPerAlloc)+37);
                    end
                    
                end


        end

        % 3.ii. Progress the state of the system 1 model time step
                % For discrete 1 model time step is < 1 year
                % For continuous 1 model time step is 1 year
            % This will loop over each time step in the time horizon,
            % progressing the system each loop.                    
            
        for HIVEpiStep = 1:(runtime/tstep)

            % 3.ii.1 Calculate transition rates
            [TransRates, TTProg, Params, InfRate, InfRate_PrEP_Oral_High, InfRate_PrEP_Oral_Low, InfRate_PrEP_Inject_High, InfRate_PrEP_Inject_Low] = ...
                    CalcTransRates(Params, TransRates, Year,  Compartments, TTProg, HIVEpiStep);
            %InfRate_all_ts{HIVEpiStep} = InfRate;
            %InfRate_PREP_all_ts{HIVEpiStep} = InfRate_PrEP;
            %TTProg_all_ts{HIVEpiStep} = TTProg;

            % 3.ii.2 Calculate the progression of the population

                % 3.ii.2.a Solve based on continuous differential equations
                if Params.SolveCont == 1
                    % Note: TTProg needs to be an argument out because
                    % it is used to record the number of ODE
                    % steps per year (which is used in CollectResults)

                    [Compartments, Results_p, TTProg]  = StepStateCont(Params, TransRates, Compartments, Year, HIVEpiStep, Results_p,TTProg );

                else  

                % 3.ii.2.b Solve based on discretized difference equations

                    [Compartments, Results_p] = StepStateDisc(Params, TransRates, Compartments, InfRate, InfRate_PrEP_Oral_High, InfRate_PrEP_Oral_Low, InfRate_PrEP_Inject_High, InfRate_PrEP_Inject_Low, Year, HIVEpiStep, Results_p,TTProg);
                    %Compartments_all_ts(:,:,HIVEpiStep) = Compartments;
                    %Transitions_all_ts(:,:,:,HIVEpiStep) = Transitions;

                end

            % 3.ii.3 Reset system at end of time step

                % Update year
                Year = Params.tt_modelStartYear + HIVEpiStep*tstep;

                % Show countdown on screen if it's a regular model run
                    % Note: if this is shown during the one-way SA run,
                    % the VBA code will throw an error
                if Epi1CalLHS2CalOpt3RAOpt4SA5EE6UA7 == 1          
                    disp(runtime/tstep - HIVEpiStep)
                end

               % Assign outcomes 
               EpiOutcomes=Results_p;
        end
    end
        
%% 4. CalcObjValue_RA function: Calculate value of optimization's objective function for resource allocation
    
    % Note: this is a nested function within HIVEpiModel.m
    
   function ObjValue_RA = CalcObjValue_RA(DV_Allocation)

       Million = 1000000;

        % 4.i. Run Epi Model
           EpiOutcomesForRA = CalcEpiOutcomes(DV_Allocation);
           
           if Params.TargetIntnsToYMSM == 1
               numIntns = Params.numIntns_YMSM;
           else
               numIntns = Params.numIntns;
           end
           
           if and(Params.optObjNum ~= 4, Params.alloc_ConsiderARTAndCareCosts == 1)
                DV_Allocation(numIntns+1)=...
                    mean(sum(EpiOutcomesForRA.ann_TotalARTAndCareCost_Undisc...
                    (Params.OutcomesIndex_AllocationStartYr(1):Params.OutcomesIndex_AllocationEndYr(1)),2))/1000000*...
                    Params.alloc_PctARTAndCareCoveredByBudget;
                if Params.numAllocationPeriods >= 2
                    DV_Allocation(2*(numIntns+1))=...
                        mean(sum(EpiOutcomesForRA.ann_TotalARTAndCareCost_Undisc...
                        (Params.OutcomesIndex_AllocationStartYr(2):Params.OutcomesIndex_AllocationEndYr(2)),2))/1000000*...
                        Params.alloc_PctARTAndCareCoveredByBudget;
                end                
                if Params.numAllocationPeriods == 3
                    DV_Allocation(3*(numIntns+1))=...
                        mean(sum(EpiOutcomesForRA.ann_TotalARTAndCareCost_Undisc...
                        (Params.OutcomesIndex_AllocationStartYr(3):Params.OutcomesIndex_AllocationEndYr(3)),2))/1000000*...
                        Params.alloc_PctARTAndCareCoveredByBudget;
                end                
           end           
           
           run_counter = run_counter + 1;
           time = toc;
           time = time / 60;
           fprintf('\nRun: %.0d', run_counter)
           fprintf('\nTime elapsed (minutes): %.2f', time)

           DV_Allocation
%            %Compare allocation to constraints
%             %Budget constraint
%             AllocVsBudget = Params.alloc_TotalIntnBudget-sum(DV_Allocation);
%             if Params.alloc_ConsiderARTAndCareCosts == 1
%                 AllocVsMinByType = ...
%                     Params.ind_IntnsByIntnType*DV_Allocation - b(1:Params.numIntnTypes,1); 
%                 AllocVsMaxByType = b(1+Params.numIntnTypes:2*Params.numIntnTypes,1) - ...
%                     Params.ind_IntnsByIntnType*DV_Allocation;
%             else
%                 AllocVsMinByType = ...
%                     Params.ind_IntnsByIntnType*DV_Allocation - b(2:1+Params.numIntnTypes,1); 
%                 AllocVsMaxByType = b(2+Params.numIntnTypes:1+2*Params.numIntnTypes,1) - ...
%                     Params.ind_IntnsByIntnType*DV_Allocation;
%             end

        % 4.2 Pull objective function value from epi model outcomes
            if Params.optObjNum == 1 %Minimize discounted new infections
                ObjValue_RA = EpiOutcomesForRA.TotalNewInfections_Disc
            elseif Params.optObjNum == 2 % Minimize discounted new infections among YMSM
                ObjValue_RA = EpiOutcomesForRA.TotalNewInfections_YMSM_Disc    
            elseif Params.optObjNum == 3 % Maximize QALYs (= minimize negative QALYs)
                ObjValue_RA = -sum(EpiOutcomesForRA.total_QALYs,2) / Million
            else %Params.optObjNum == 4 %Minimize intervention spending to hit incidence targets
                
                NumYrsInTotalAllocPeriod = max(Params.outcomeCollectionEndYr - Params.tt_periodFiveStartYear + 1,0);
                
                NumYrsInAllocPeriod1 = Params.OutcomesIndex_AllocationEndYr(1) - Params.OutcomesIndex_AllocationStartYr(1) + 1;
                TotalSpending = sum(DV_Allocation(1:numIntns,1)) * NumYrsInAllocPeriod1 * Million;
               
                if Params.numAllocationPeriods >= 2
                   
                    NumYrsInAllocPeriod2 = Params.OutcomesIndex_AllocationEndYr(2) - Params.OutcomesIndex_AllocationStartYr(2) + 1;
                    TotalSpending = TotalSpending + (sum(DV_Allocation(numIntns + 1:(2 * numIntns),1)) * NumYrsInAllocPeriod2 * Million);
                   
                end                
                if Params.numAllocationPeriods == 3
                    NumYrsInAllocPeriod3 = Params.OutcomesIndex_AllocationEndYr(3) - Params.OutcomesIndex_AllocationStartYr(3) + 1;
                    TotalSpending = TotalSpending + (sum(DV_Allocation((2 * numIntns + 1):(3 * numIntns),1)) * NumYrsInAllocPeriod3 * Million);
                end
                
                ObjValue_RA = TotalSpending / NumYrsInTotalAllocPeriod
                               
            end

        % 4.3 Save outcomes

        %Verticaly combines
        %outcomes of interest. Then combines newest run with all
        %older runs and saves. If running cost minimization, also combines
        %includes incidence reduction targets
        if Params.optObjNum == 4
            
            IncidenceReductionTarget1 = EpiOutcomesForRA.IncidenceReductionTarget1
            IncidenceReductionTarget2 = EpiOutcomesForRA.IncidenceReductionTarget2
            AnnualUndiscInfections = sum(EpiOutcomesForRA.ann_TotalNewInfections,2);
            AnnualUndiscTxandCareCost_inM = sum(EpiOutcomesForRA.ann_TotalARTAndCareCost_Undisc,2) / 1000000;
            
            temp1 = vertcat(time,ObjValue_RA, IncidenceReductionTarget1, ...
                IncidenceReductionTarget2, DV_Allocation, AnnualUndiscInfections, AnnualUndiscTxandCareCost_inM);
            
        else
            
            temp1 = vertcat(time,ObjValue_RA, ...
                DV_Allocation);
            
        end
        
        RA_results = horzcat(RA_results, temp1);
        save('RA_results','RA_results');




   end

%% 5. CalcObjValue_Calib function: Calculate value of optimization's objective function for calibration
    
    % Note: this is a nested function within HIVEpiModel.m
    
   function ObjValue_Calib = CalcObjValue_Calib(DV_InputSet)

       % 5.i. Run Epi Model
       EpiOutcomesForCalib = CalcEpiOutcomes(DV_InputSet);
        DV_InputSet;

        % 5.ii. Pull key outputs for calibration    
        CalibOutputStruct = Calib_collectResults(EpiOutcomesForCalib,[], 1);

        % 5.iii Calculate ObjValue_Calib = aggregated error measure

        [ObjValue_Calib, OOBPenaltyObjValue_Calib, TargetErrorObjValue_Calib] = ...
            Calib_calcObjValue(CalibOutputStruct, CalibTargets, calibObjValueType);

        % 5.iv  Calculate number of outcomes outside of bounds

            % Collect names and count number of fields in CalibTargets
            TargetFields = fieldnames(CalibTargets);
            nTargetOutcomes = length(TargetFields);
            ResultsFields = fieldnames(CalibOutputStruct);
            nResults = length(ResultsFields);
            assert(nTargetOutcomes == nResults, ...
                'CalibOutputStruct and CalibTargets are not of equal size');

            % Pre-allocate / initiate variables
            LB = 1;
            UB = 2;
            target_LB(nTargetOutcomes,1)=0;
            target_UB(nTargetOutcomes,1)=0;
            modelOutcomeValues(nResults,1)=0;

            % Calculate number of targets (with non-zero weights) that
            %   fall outside their corresponding range
            for nTargetOutcome = 1:nTargetOutcomes
                target_LB(nTargetOutcome,1)=CalibTargets.(TargetFields{nTargetOutcome}).range(LB);
                target_UB(nTargetOutcome,1)=CalibTargets.(TargetFields{nTargetOutcome}).range(UB);
            end

            for nResult = 1:nResults
                modelOutcomeValues(nResult,1)=CalibOutputStruct.(TargetFields{nResult}).paramValue;
            end   

            % If all non-zero-weighted targets fall within target ranges,
            % save current input set and results
            nOutsideRange = 0;
            for nTargetOutcome = 1:nTargetOutcomes
                if and(CalibTargets.(TargetFields{nTargetOutcome}).weight > 0, ...
                    or(modelOutcomeValues(nTargetOutcome,1)<target_LB(nTargetOutcome,1), ...
                    modelOutcomeValues(nTargetOutcome,1)>target_UB(nTargetOutcome,1)))
                        nOutsideRange = nOutsideRange + 1;
                end
            end
            nOutsideRange;
            
        % 5.v Set variable for full set of calibrated input values
        CalibFields = fieldnames(CalibParams);
        nCalibratedParams = length(CalibFields);
        CalibParamValues(nCalibratedParams,1) = 0;
        for j = 1:nCalibratedParams
            CalibParamValues(j,1) = CalibParams.(CalibFields{j}).paramValue(1);
        end
        
        % 5.vi Save outcomes

        %increments run counter for saving then verticaly combines
        %outcomes of interest. Then combines newest run with all
        %older runs and saves.
        run_counter = run_counter + 1;
        temp1 = vertcat(ObjValue_Calib, ...
            OOBPenaltyObjValue_Calib, TargetErrorObjValue_Calib,  ...
            nOutsideRange, CalibParamValues, modelOutcomeValues);
        calib_results = horzcat(calib_results, temp1);
        save('calib_results','calib_results');

        time = toc;
        time = time / 60;
        fprintf('\nRun: %.0d', run_counter)
        fprintf('\nTime elapsed (minutes): %.2f', time)
        if ObjValue_Calib == OOBPenaltyObjValue_Calib
            fprintf('\nOOB obj value: %.10f', OOBPenaltyObjValue_Calib)
            fprintf('\nTE obj value: %.10f', TargetErrorObjValue_Calib)
        else
            fprintf('\nTE obj value: %.10f', TargetErrorObjValue_Calib)
            fprintf('\nOOB obj value: %.10f', OOBPenaltyObjValue_Calib)
        end
        fprintf('\nTargets OOB: %.0d\n\n', nOutsideRange)

        %stop condition based on number of runs below objective value
        %threshold.  Increments a counter for runs that have occured
        %with OV below threshold, then stops if that number is equal
        %the user input for runs to stop at

        if ObjValue_Calib < calibOVStop
            runs_below_threshold = runs_below_threshold + 1;
        end

        if runs_below_threshold == calibNumBelowStop
            error('This is not an error. Runs below objective value threshold have been reached. Results can now be viewed.')
        end

   end
 
%% 6. CalcObjValue_Cont function: Calculate value of optimization's objective function for reaching continuum targets
    
    % Note: this is a nested function within HIVEpiModel.m
    
    function ObjValue_Cont = CalcObjValue_Cont(DV_InputValue)

        % 4.i. Run Epi Model
           EpiOutcomesForNHAS = CalcEpiOutcomes(DV_InputValue);
            DV_InputValue

        %JC, set objective function to be the absolute difference
        %between the value of the selected target outcome in the model
        %and the target value for that outcome. I would guess that
        %you'd just have a different if statement for each of the
        %different target outcomes so that it pulls the right one.

        % 4.ii. Pull selected outcome
        TargetOutcome = Params.nhas_TargetOutcome;

        switch TargetOutcome

            case 1 % Percentage of PLWH who know their HIV status in 2020: Black
                OutcomeVal = EpiOutcomesForNHAS.nhas_pctOfPLWH_Aware_2020_Blk;
            case 2 % Percentage of PLWH who know their HIV status in 2020: Hispanic/Latino
                OutcomeVal = EpiOutcomesForNHAS.nhas_pctOfPLWH_Aware_2020_Hisp;
            case 3 % Percentage of PLWH who know their HIV status in 2020: Other
                OutcomeVal = EpiOutcomesForNHAS.nhas_pctOfPLWH_Aware_2020_Oth;
            case 4 % Percentage of diagnosed PLWH who are VLS in 2020: Black
                OutcomeVal = EpiOutcomesForNHAS.nhas_pctOfAware_VLS_2020_Blk;
            case 5 % Percentage of diagnosed PLWH who are VLS in 2020: Hispanic/Latino
                OutcomeVal = EpiOutcomesForNHAS.nhas_pctOfAware_VLS_2020_Hisp;
            case 6 % Percentage of diagnosed PLWH who are VLS in 2020: Other
                OutcomeVal = EpiOutcomesForNHAS.nhas_pctOfAware_VLS_2020_Oth;
        end    

        % 4.iii. Calculate objective value (absolute difference of outcome and
        % target outcome value)
        ObjValue_Cont = abs(OutcomeVal - Params.nhas_TargetOutcomeVal);


   end
   
%% 7. Non-linear constraints for optimization algorithms
   
    function [c,ceq] = nonlcon_withARTAndCareCosts(DVAllocation)
    % PER KH in June 2019: THIS CODE NEEDS REVISING. FIRST, IT IS CURRENTLY RUNNING THE MODEL ONCE
    % FOR EACH ALLOCATION PERIOD. IT INSTEAD NEEDS TO RUN IT ONCE AND PULL THE OUTCOMES
    % CORRESPONDING TO EACH PERIOD. ALSO, FMINCON IS UNHAPPY WITH C BEING A
    % VECTOR. FOR NOW, I AM JUST FORCING THE USER NOT TO ALLOW THE
    % CONSIDERATION OF TREATMENT AND CARE COSTS SO THAT THIS CODE DOESN'T GET
    % CALLED.
        
        c = sum(DVAllocation(1:Params.numIntns))+ ...
            CalcTxAndCareCosts_RA(DVAllocation, ...
            Params.OutcomesIndex_AllocationStartYr(1), Params.OutcomesIndex_AllocationEndYr(1)) ...
            -Params.alloc_TotalIntnBudget(1);     % Compute nonlinear inequalities at x.
        if Params.numAllocationPeriods >= 2
            c = [c;
                sum(DVAllocation((Params.numIntns+1+1):2*(Params.numIntns+1)-1))+ ...
                CalcTxAndCareCosts_RA(DVAllocation, ...
                Params.OutcomesIndex_AllocationStartYr(2), Params.OutcomesIndex_AllocationEndYr(2)) ...
                -Params.alloc_TotalIntnBudget(2)];     % Compute nonlinear inequalities at x.
        end
        if Params.numAllocationPeriods == 3
            c = [c;
                sum(DVAllocation((2*(Params.numIntns+1)+1):3*(Params.numIntns+1)-1))+ ...
                CalcTxAndCareCosts_RA(DVAllocation, ...
                Params.OutcomesIndex_AllocationStartYr(3), Params.OutcomesIndex_AllocationEndYr(3)) ...
                -Params.alloc_TotalIntnBudget(3)];     % Compute nonlinear inequalities at x.
        end
        ceq = [];   % Compute nonlinear equalities at x.
   end

   function [c,ceq] = nonlcon_hitIncidenceTargets(DV_Allocation)
                       
        %Calculate c(x) so that constraint is c(x)<=0
        c = CalcRelativeIncidence_RA(DV_Allocation);
        ceq = [];   % Compute nonlinear equalities at x.
        
   end

   function [c,ceq] = nonlcon(DVAllocation)
        c = [];     % Compute nonlinear inequalities at x.
        ceq = [];   % Compute nonlinear equalities at x.
   end

%% 8. Functions for enabling non-linear constraints in resource allocation optimization
% Note: these are nested functions within HIVEpiModel.m

% 8a. CalcTxAndCareCosts_RA function: Calculate treatment and care costs resulting from specific allocation of funding
    
   function MeanTxAndCareCosts_RA_InMil = CalcTxAndCareCosts_RA(DV_Allocation, IndexStart, IndexEnd)

        % 8a.i. Run Epi Model
       AllEpiOutcomesForTxAndCareCostCalc = CalcEpiOutcomes(DV_Allocation);

        % 8a.ii Pull treatment and care costs from epi model outcomes (in
        % millions)
       MeanTxAndCareCosts_RA_InMil = ...
           mean(sum(AllEpiOutcomesForTxAndCareCostCalc.ann_TotalARTAndCareCost_Undisc(IndexStart:IndexEnd,:),2))/1000000*...
           Params.alloc_PctARTAndCareCoveredByBudget;

   end

% 8b. CalcRelativeIncidence_RA function: Calculate difference between 
% target incidence values and HOPE reduction in incidence from year before 
% time period 3 starts to each of 2target years (must be <=0) 

   function RelativeIncidence_RA = CalcRelativeIncidence_RA(DV_Allocation)
        %Method: 
        % -(Incidence reduction outcome) + (Target % reduction)
        
        % 8.i. Run Epi Model
        AllEpiOutcomesForRelIncidCalc = CalcEpiOutcomes(DV_Allocation);
               
        
        % 8.ii Calc outcomes
        RelativeIncidence_RA(1,1) = -AllEpiOutcomesForRelIncidCalc.IncidenceReductionTarget1 + Params.minCost_Tgt1PctReduce;
        RelativeIncidence_RA(2,1) = -AllEpiOutcomesForRelIncidCalc.IncidenceReductionTarget2 + Params.minCost_Tgt2PctReduce;
 
%         if nargout > 1
%              RelativeIncidence_RA_2 = abs(dif);
%         end


   end
   

%% 9. PickOutcomeList function: Pull outcomes based on the outcome list selected by the user in the Model Settings sheet
    
    function [SelectedOutcomes] = PickOutcomeList(AllOutcomes, ListChoice)
        
        FirstOutcomeYr = max(Params.tt_modelStartYear,Params.outcomeCollectionStartYr);
        LastOutcomeYr = min(Params.tt_modelStartYear+Params.tt_timeHorizon-1,Params.outcomeCollectionEndYr);
        LastOutcomeYrIndex = max(floor(LastOutcomeYr - FirstOutcomeYr + 1),0);
        OutcomesIndexTarget0=max(Params.OutcomesIndex_AllocationStartYr(1)-1,1);
        OutcomesIndexTarget1=Params.OutcomesIndex_AllocationStartYr(1)+Params.minCost_Tgt1NumYears-1;
        OutcomesIndexTarget2=Params.OutcomesIndex_AllocationStartYr(1)+Params.minCost_Tgt2NumYears-1;
        rateYr = Params.diagRateYr;
        if FirstOutcomeYr <= rateYr && (LastOutcomeYr >= rateYr)
            rateYrIndex = max(floor(rateYr - FirstOutcomeYr + 1),0);
        else
            rateYrIndex = -1;
        end
              
               
        switch ListChoice
                
            case 1 %All outcomes

                SelectedOutcomes = AllOutcomes;

            case 2 %Key outcomes only

                SelectedOutcomes.ann_PopulationSize = sum(AllOutcomes.ann_PopulationSize,2);
                SelectedOutcomes.ann_HIVPrevalence_HET = sum(AllOutcomes.ann_HIVPrevalence_HET,2);
                SelectedOutcomes.ann_HIVPrevalence_MSM = sum(AllOutcomes.ann_HIVPrevalence_MSM,2);
                SelectedOutcomes.ann_HIVPrevalence_IDU = sum(AllOutcomes.ann_HIVPrevalence_IDU,2);
                SelectedOutcomes.ann_TotalNewInfections = sum(AllOutcomes.ann_TotalNewInfections,2);
                SelectedOutcomes.ann_LifeYears = sum(AllOutcomes.ann_LifeYears,2);
                SelectedOutcomes.ann_QALYs = sum(AllOutcomes.ann_QALYs,2);
                SelectedOutcomes.ann_PctAware = AllOutcomes.ann_PctAware;
                SelectedOutcomes.ann_PctVLSamongdiag = AllOutcomes.ann_PctVLSamongdiag;
                SelectedOutcomes.total_ARTCarePrEPTransitionAndSEPCost_Disc = sum(AllOutcomes.total_ARTCarePrEPTransitionAndSEPCost_Disc,2);

            case 3 %Calibration outcomes            
            
                nCalibYrsIncluded = 0;
                %Removed 2010 targets. Bates 6/11/19. Updated 4/21/20 for
                %PrEP analysis.
                % 2019 targets added. Clinkscales 6/29/2021
                % 2019 NewDiag targets added (lines 1275-76). Clinkscales 04/21/2022
                % InCare targets added. Bates 11/4/2022
                % Additional targets updated to 2019. Bates 11/22/22
                % Updated # on PrEP targets. Bates 02/10/2023
                if AllOutcomes.set_FirstOutcomeYr <= 2010 &&(AllOutcomes.set_LastOutcomeYr >= 2018) == 1
                    
                    SelectedOutcomes.OOBobjvalue = AllOutcomes.OOBobjvalue;
                    SelectedOutcomes.numTargetsOOB = AllOutcomes.numTargetsOOB;
                 
                    SelectedOutcomes.calib_HIVPrevalence2019 = AllOutcomes.calib_HIVPrevalence2019;
                    SelectedOutcomes.calib_HIVPrevalence2019_B = AllOutcomes.calib_HIVPrevalence2019_B;
                    SelectedOutcomes.calib_HIVPrevalence2019_H = AllOutcomes.calib_HIVPrevalence2019_H;
                    SelectedOutcomes.calib_HIVPrevalence2019_O = AllOutcomes.calib_HIVPrevalence2019_O;
                    SelectedOutcomes.calib_HIVPrevalence2019_HETM = AllOutcomes.calib_HIVPrevalence2019_HETM;
                    SelectedOutcomes.calib_HIVPrevalence2019_HETF = AllOutcomes.calib_HIVPrevalence2019_HETF;
                    SelectedOutcomes.calib_HIVPrevalence2019_MSM = AllOutcomes.calib_HIVPrevalence2019_MSM;
                    SelectedOutcomes.calib_HIVPrevalence2019_IDUM = AllOutcomes.calib_HIVPrevalence2019_IDUM;
                    SelectedOutcomes.calib_HIVPrevalence2019_IDUF = AllOutcomes.calib_HIVPrevalence2019_IDUF;
                    SelectedOutcomes.calib_HIVPrevalence2016_13_24 = AllOutcomes.calib_HIVPrevalence2016_13_24;
                    SelectedOutcomes.calib_HIVPrevalence2016_25_34 = AllOutcomes.calib_HIVPrevalence2016_25_34;
                    SelectedOutcomes.calib_HIVPrevalence2016_35_44 = AllOutcomes.calib_HIVPrevalence2016_35_44;
                    SelectedOutcomes.calib_HIVPrevalence2016_45_54 = AllOutcomes.calib_HIVPrevalence2016_45_54;
                    SelectedOutcomes.calib_HIVPrevalence2016_55_64 = AllOutcomes.calib_HIVPrevalence2016_55_64;
                    SelectedOutcomes.calib_HIVPrevalence2016_65 = AllOutcomes.calib_HIVPrevalence2016_65;
                    SelectedOutcomes.calib_HIVPrevalence2016_B_13_24 = AllOutcomes.calib_HIVPrevalence2016_B_13_24;
                    SelectedOutcomes.calib_HIVPrevalence2016_H_13_24 = AllOutcomes.calib_HIVPrevalence2016_H_13_24;
                    SelectedOutcomes.calib_HIVPrevalence2016_O_13_24 = AllOutcomes.calib_HIVPrevalence2016_O_13_24;                    
                    SelectedOutcomes.calib_HIVPrevalence2016_B_25_34 = AllOutcomes.calib_HIVPrevalence2016_B_25_34;
                    SelectedOutcomes.calib_HIVPrevalence2016_H_25_34 = AllOutcomes.calib_HIVPrevalence2016_H_25_34;
                    SelectedOutcomes.calib_HIVPrevalence2016_O_25_34 = AllOutcomes.calib_HIVPrevalence2016_O_25_34;                    
                    SelectedOutcomes.calib_TTdist2019_Diagnosed_B = AllOutcomes.calib_TTdist2019_Diagnosed_B;
                    SelectedOutcomes.calib_TTdist2019_Diagnosed_H = AllOutcomes.calib_TTdist2019_Diagnosed_H;
                    SelectedOutcomes.calib_TTdist2019_Diagnosed_O = AllOutcomes.calib_TTdist2019_Diagnosed_O;
                    SelectedOutcomes.calib_TTdist2016_Diagnosed_13_24 = AllOutcomes.calib_TTdist2016_Diagnosed_13_24;
                    SelectedOutcomes.calib_TTdist2016_Diagnosed_25_34 = AllOutcomes.calib_TTdist2016_Diagnosed_25_34;
                    SelectedOutcomes.calib_TTdist2016_Diagnosed_35_44 = AllOutcomes.calib_TTdist2016_Diagnosed_35_44;
                    SelectedOutcomes.calib_TTdist2016_Diagnosed_45_54 = AllOutcomes.calib_TTdist2016_Diagnosed_45_54;
                    SelectedOutcomes.calib_TTdist2016_Diagnosed_55_64 = AllOutcomes.calib_TTdist2016_Diagnosed_55_64;
                    SelectedOutcomes.calib_TTdist2016_Diagnosed_65 = AllOutcomes.calib_TTdist2016_Diagnosed_65;
                    SelectedOutcomes.calib_TTdist2016_Diagnosed_MSM_B_13_24 = AllOutcomes.calib_TTdist2016_Diagnosed_MSM_B_13_24;
                    SelectedOutcomes.calib_TTdist2016_Diagnosed_MSM_H_13_24 = AllOutcomes.calib_TTdist2016_Diagnosed_MSM_H_13_24;
                    SelectedOutcomes.calib_TTdist2016_Diagnosed_MSM_O_13_24 = AllOutcomes.calib_TTdist2016_Diagnosed_MSM_O_13_24;                    
                    SelectedOutcomes.calib_TTdist2016_Diagnosed_MSM_B_25_34 = AllOutcomes.calib_TTdist2016_Diagnosed_MSM_B_25_34;
                    SelectedOutcomes.calib_TTdist2016_Diagnosed_MSM_H_25_34 = AllOutcomes.calib_TTdist2016_Diagnosed_MSM_H_25_34;
                    SelectedOutcomes.calib_TTdist2016_Diagnosed_MSM_O_25_34 = AllOutcomes.calib_TTdist2016_Diagnosed_MSM_O_25_34; 
                    SelectedOutcomes.calib_TTdist2019_Diagnosed_HET = AllOutcomes.calib_TTdist2019_Diagnosed_HET;
                    SelectedOutcomes.calib_TTdist2019_Diagnosed_MSM = AllOutcomes.calib_TTdist2019_Diagnosed_MSM;
                    SelectedOutcomes.calib_TTdist2019_Diagnosed_PWID = AllOutcomes.calib_TTdist2019_Diagnosed_PWID;
                    SelectedOutcomes.calib_TTdist2019_Diagnosed_Total = AllOutcomes.calib_TTdist2019_Diagnosed_Total;
                    SelectedOutcomes.calib_TTdist2019_InCareAmongDiag_B = AllOutcomes.calib_TTdist2019_InCareAmongDiag_B;
                    SelectedOutcomes.calib_TTdist2019_InCareAmongDiag_H = AllOutcomes.calib_TTdist2019_InCareAmongDiag_H;
                    SelectedOutcomes.calib_TTdist2019_InCareAmongDiag_O = AllOutcomes.calib_TTdist2019_InCareAmongDiag_O;
                    SelectedOutcomes.calib_TTdist2019_InCareAmongDiag_Total = AllOutcomes.calib_TTdist2019_InCareAmongDiag_Total;                    
                    SelectedOutcomes.calib_TTdist2019_VLSamongdiag_B = AllOutcomes.calib_TTdist2019_VLSamongdiag_B;
                    SelectedOutcomes.calib_TTdist2019_VLSamongdiag_H = AllOutcomes.calib_TTdist2019_VLSamongdiag_H;
                    SelectedOutcomes.calib_TTdist2019_VLSamongdiag_O = AllOutcomes.calib_TTdist2019_VLSamongdiag_O;
                    SelectedOutcomes.calib_TTdist2018_VLSamongdiag_HET = AllOutcomes.calib_TTdist2018_VLSamongdiag_HET;
                    SelectedOutcomes.calib_TTdist2018_VLSamongdiag_MSM = AllOutcomes.calib_TTdist2018_VLSamongdiag_MSM;
                    SelectedOutcomes.calib_TTdist2018_VLSamongdiag_PWID = AllOutcomes.calib_TTdist2018_VLSamongdiag_PWID;
                    SelectedOutcomes.calib_TTdist2015_VLSamongdiag_13_24 = AllOutcomes.calib_TTdist2015_VLSamongdiag_13_24;
                    SelectedOutcomes.calib_TTdist2015_VLSamongdiag_25_34 = AllOutcomes.calib_TTdist2015_VLSamongdiag_25_34;
                    SelectedOutcomes.calib_TTdist2015_VLSamongdiag_35_44 = AllOutcomes.calib_TTdist2015_VLSamongdiag_35_44;
                    SelectedOutcomes.calib_TTdist2015_VLSamongdiag_45_54 = AllOutcomes.calib_TTdist2015_VLSamongdiag_45_54;
                    SelectedOutcomes.calib_TTdist2015_VLSamongdiag_55_64 = AllOutcomes.calib_TTdist2015_VLSamongdiag_55_64;
                    SelectedOutcomes.calib_TTdist2015_VLSamongdiag_65 = AllOutcomes.calib_TTdist2015_VLSamongdiag_65;
                    SelectedOutcomes.calib_TTdist2016_VLSamongdiag_13_24_MSM_B = AllOutcomes.calib_TTdist2016_VLSamongdiag_13_24_MSM_B;
                    SelectedOutcomes.calib_TTdist2016_VLSamongdiag_13_24_MSM_H = AllOutcomes.calib_TTdist2016_VLSamongdiag_13_24_MSM_H;
                    SelectedOutcomes.calib_TTdist2016_VLSamongdiag_13_24_MSM_O = AllOutcomes.calib_TTdist2016_VLSamongdiag_13_24_MSM_O;                  
                    SelectedOutcomes.calib_TTdist2016_VLSamongdiag_25_34_MSM_B = AllOutcomes.calib_TTdist2016_VLSamongdiag_25_34_MSM_B;
                    SelectedOutcomes.calib_TTdist2016_VLSamongdiag_25_34_MSM_H = AllOutcomes.calib_TTdist2016_VLSamongdiag_25_34_MSM_H;
                    SelectedOutcomes.calib_TTdist2016_VLSamongdiag_25_34_MSM_O = AllOutcomes.calib_TTdist2016_VLSamongdiag_25_34_MSM_O;                  
                    SelectedOutcomes.calib_TTdist2019_VLSamongdiag_Total = AllOutcomes.calib_TTdist2019_VLSamongdiag_Total;
                    SelectedOutcomes.calib_AwarePLWHDeaths2019 = AllOutcomes.calib_AwarePLWHDeaths2019;
                    SelectedOutcomes.calib_AwarePLWHDeaths2016_45_54 = AllOutcomes.calib_AwarePLWHDeaths2016_45_54;
                    SelectedOutcomes.calib_AwarePLWHDeaths2016_55_64 = AllOutcomes.calib_AwarePLWHDeaths2016_55_64;
                    SelectedOutcomes.calib_AwarePLWHDeaths2016_65 = AllOutcomes.calib_AwarePLWHDeaths2016_65;
                    SelectedOutcomes.calib_AwarePLWHAIDSDeaths2019 = AllOutcomes.calib_AwarePLWHAIDSDeaths2019;
                    SelectedOutcomes.calib_UnawarePLWHAIDSDeaths2019 = AllOutcomes.calib_UnawarePLWHAIDSDeaths2019;
                    SelectedOutcomes.calib_NewDiagnoses2016_HETM = AllOutcomes.calib_NewDiagnoses2016_HETM;
                    SelectedOutcomes.calib_NewDiagnoses2016_HETF = AllOutcomes.calib_NewDiagnoses2016_HETF;
                    SelectedOutcomes.calib_NewDiagnoses2016_MSM = AllOutcomes.calib_NewDiagnoses2016_MSM;
                    SelectedOutcomes.calib_NewDiagnoses2016_IDUM = AllOutcomes.calib_NewDiagnoses2016_IDUM;
                    SelectedOutcomes.calib_NewDiagnoses2016_IDUF = AllOutcomes.calib_NewDiagnoses2016_IDUF;
                    SelectedOutcomes.calib_NewDiagnoses2016_Total = AllOutcomes.calib_NewDiagnoses2016_Total;
                    SelectedOutcomes.calib_NewDiagnoses2019_PctB1C1 = AllOutcomes.calib_NewDiagnoses2019_PctB1C1;      
                    SelectedOutcomes.calib_NewDiagnoses2019_Total = AllOutcomes.calib_NewDiagnoses2019_Total;
                    SelectedOutcomes.calib_NewInfections2019_B = AllOutcomes.calib_NewInfections2019_B;
                    SelectedOutcomes.calib_NewInfections2019_H = AllOutcomes.calib_NewInfections2019_H;
                    SelectedOutcomes.calib_NewInfections2019_O = AllOutcomes.calib_NewInfections2019_O;
                    SelectedOutcomes.calib_NewInfections2019_HETM = AllOutcomes.calib_NewInfections2019_HETM;
                    SelectedOutcomes.calib_NewInfections2019_HETF = AllOutcomes.calib_NewInfections2019_HETF;
                    SelectedOutcomes.calib_NewInfections2019_MSM = AllOutcomes.calib_NewInfections2019_MSM;
                    SelectedOutcomes.calib_NewInfections2019_IDUM = AllOutcomes.calib_NewInfections2019_IDUM;
                    SelectedOutcomes.calib_NewInfections2019_IDUF = AllOutcomes.calib_NewInfections2019_IDUF;
                    SelectedOutcomes.calib_NewInfections_2016_MSM_13_24 = AllOutcomes.calib_NewInfections_2016_MSM_13_24;
                    SelectedOutcomes.calib_NewInfections_2016_MSM_13_24_B = AllOutcomes.calib_NewInfections_2016_MSM_13_24_B;
                    SelectedOutcomes.calib_NewInfections_2016_MSM_13_24_H = AllOutcomes.calib_NewInfections_2016_MSM_13_24_H;
                    SelectedOutcomes.calib_NewInfections_2016_MSM_13_24_O = AllOutcomes.calib_NewInfections_2016_MSM_13_24_O;
                    SelectedOutcomes.calib_NewInfections_2016_MSM_25_34 = AllOutcomes.calib_NewInfections_2016_MSM_25_34;
                    SelectedOutcomes.calib_NewInfections_2016_MSM_25_34_B = AllOutcomes.calib_NewInfections_2016_MSM_25_34_B;
                    SelectedOutcomes.calib_NewInfections_2016_MSM_25_34_H = AllOutcomes.calib_NewInfections_2016_MSM_25_34_H;
                    SelectedOutcomes.calib_NewInfections_2016_MSM_25_34_O = AllOutcomes.calib_NewInfections_2016_MSM_25_34_O;
                    SelectedOutcomes.calib_NewInfections2019_Total = AllOutcomes.calib_NewInfections2019_Total;
                    SelectedOutcomes.calib_NewInfections2016_45_54 = AllOutcomes.calib_NewInfections2016_45_54;
                    SelectedOutcomes.calib_NewInfections2016_55_64 = AllOutcomes.calib_NewInfections2016_55_64;
                    SelectedOutcomes.calib_NewInfections2016_65 = AllOutcomes.calib_NewInfections2016_65;
                    SelectedOutcomes.calib_NumOnPrEP2019_Male = AllOutcomes.calib_NumOnPrEP2019_Male;           
                    SelectedOutcomes.calib_NumOnPrEP2019_Female = AllOutcomes.calib_NumOnPrEP2019_Female;
                    SelectedOutcomes.calib_NumOnPrEP2019_Black = AllOutcomes.calib_NumOnPrEP2019_Black;
                    SelectedOutcomes.calib_NumOnPrEP2019_Hispanic = AllOutcomes.calib_NumOnPrEP2019_Hispanic;
                    SelectedOutcomes.calib_NumOnPrEP2019_Other = AllOutcomes.calib_NumOnPrEP2019_Other;
                    SelectedOutcomes.calib_NumOnPrEP2019_Total = AllOutcomes.calib_NumOnPrEP2019_Total;
                    SelectedOutcomes.calib_NumOnPrEP2021_Male = AllOutcomes.calib_NumOnPrEP2021_Male;
                    SelectedOutcomes.calib_NumOnPrEP2021_Female = AllOutcomes.calib_NumOnPrEP2021_Female;
                    SelectedOutcomes.calib_NumOnPrEP2021_Black = AllOutcomes.calib_NumOnPrEP2021_Black;
                    SelectedOutcomes.calib_NumOnPrEP2021_Hispanic = AllOutcomes.calib_NumOnPrEP2021_Hispanic;
                    SelectedOutcomes.calib_NumOnPrEP2021_Other = AllOutcomes.calib_NumOnPrEP2021_Other;
                    SelectedOutcomes.calib_NumOnPrEP2021_Total = AllOutcomes.calib_NumOnPrEP2021_Total;
                    SelectedOutcomes.calib_OverallPrev2019v2010 = AllOutcomes.calib_OverallPrev2019v2010;
                    SelectedOutcomes.calib_HETPrev2019v2010_LR = AllOutcomes.calib_HETPrev2019v2010_LR;
                    SelectedOutcomes.calib_HETPrev2019v2010_HR = AllOutcomes.calib_HETPrev2019v2010_HR;
                    SelectedOutcomes.calib_populationMaleBlk2019 = AllOutcomes.calib_populationMaleBlk2019;
                    SelectedOutcomes.calib_populationMaleHisp2019 = AllOutcomes.calib_populationMaleHisp2019;
                    SelectedOutcomes.calib_populationMaleOth2019 = AllOutcomes.calib_populationMaleOth2019;
                    SelectedOutcomes.calib_populationFemaleBlk2019 = AllOutcomes.calib_populationFemaleBlk2019;
                    SelectedOutcomes.calib_populationFemaleHisp2019 = AllOutcomes.calib_populationFemaleHisp2019;
                    SelectedOutcomes.calib_populationFemaleOth2019 = AllOutcomes.calib_populationFemaleOth2019;
                    SelectedOutcomes.calib_population_1334_2019 = AllOutcomes.calib_population_1334_2019;
                    SelectedOutcomes.calib_population_3564_2019 = AllOutcomes.calib_population_3564_2019;
                    SelectedOutcomes.calib_population_65_2019 = AllOutcomes.calib_population_65_2019;

                    OutcomeFields = fieldnames(SelectedOutcomes);

                    nOutcomes = length(OutcomeFields);

                    for OCount=1:nOutcomes
                        SelectedOutcomes.summary(OCount,1) = SelectedOutcomes.(OutcomeFields{OCount});
                    end
                    
                else
                    h = msgbox(['You only collected outcomes from ' num2str(AllOutcomes.set_FirstOutcomeYr ) ...
                        ' to ' num2str(AllOutcomes.set_LastOutcomeYr) ...
                        '. Please revise to include all years for which we have calibration outcomes.'] , 'Error');
                    
                    SelectedOutcomes = [];
                    
                end

            case 4 % Resource allocation outcomes  
                                               
                Million = 1000000;
                
                MaxReach.Testing = Params.intn_Testing_MaxReach;
                MaxReach.LTCatDiag = Params.intn_LTCatDiag_MaxReach;
                MaxReach.LTCPostDiag = Params.intn_LTCafterDiag_MaxReach;
                MaxReach.ARTPrescription = Params.intn_ARTInitiation_MaxReach;
                MaxReach.ARTAdherance = Params.intn_TxAdherence_MaxReach;
                MaxReach.SEP = Params.intn_SEP_MaxReach;
                MaxReach.PrEP_Oral_HET = Params.intn_PrEP_Oral_HET_MaxReach;
                MaxReach.PrEP_Oral_MSM = Params.intn_PrEP_Oral_MSM_MaxReach;
                MaxReach.PrEP_Oral_IDU = Params.intn_PrEP_Oral_IDU_MaxReach;
                MaxReach.PrEP_Inject_HET = Params.intn_PrEP_Inject_HET_MaxReach;
                MaxReach.PrEP_Inject_MSM = Params.intn_PrEP_Inject_MSM_MaxReach;
                MaxReach.PrEP_Inject_IDU = Params.intn_PrEP_Inject_IDU_MaxReach;
                
                Allocation.alloc_Testing_HET_Low = AllOutcomes.alloc_Testing_HET_Low;
                Allocation.alloc_Testing_HET_High = AllOutcomes.alloc_Testing_HET_High;
                Allocation.alloc_Testing_MSM_Low = AllOutcomes.alloc_Testing_MSM_Low;
                Allocation.alloc_Testing_MSM_High = AllOutcomes.alloc_Testing_MSM_High;
                Allocation.alloc_Testing_IDU = AllOutcomes.alloc_Testing_IDU;
                Allocation.alloc_LTCatDiag = AllOutcomes.alloc_LTCatDiag;
                Allocation.alloc_LTCafterDiag = AllOutcomes.alloc_LTCafterDiag;
                Allocation.alloc_ARTInitiation = AllOutcomes.alloc_ARTInitiation;
                Allocation.alloc_ARTAdher5to4 = AllOutcomes.alloc_ARTAdher5to4;
                Allocation.alloc_ARTAdher4to5 = AllOutcomes.alloc_ARTAdher4to5;              
                Allocation.alloc_SEP_B = AllOutcomes.alloc_SEP_B;
                Allocation.alloc_SEP_H = AllOutcomes.alloc_SEP_H;
                Allocation.alloc_SEP_O = AllOutcomes.alloc_SEP_O;
                Allocation.alloc_PrEP_Oral_HETM_B = AllOutcomes.alloc_PrEP_Oral_HETM_B;
                Allocation.alloc_PrEP_Oral_HETM_H = AllOutcomes.alloc_PrEP_Oral_HETM_H;
                Allocation.alloc_PrEP_Oral_HETM_O = AllOutcomes.alloc_PrEP_Oral_HETM_O;
                Allocation.alloc_PrEP_Oral_HETF_B = AllOutcomes.alloc_PrEP_Oral_HETF_B;
                Allocation.alloc_PrEP_Oral_HETF_H = AllOutcomes.alloc_PrEP_Oral_HETF_H;
                Allocation.alloc_PrEP_Oral_HETF_O = AllOutcomes.alloc_PrEP_Oral_HETF_O;
                Allocation.alloc_PrEP_Oral_MSM_B = AllOutcomes.alloc_PrEP_Oral_MSM_B;
                Allocation.alloc_PrEP_Oral_MSM_H = AllOutcomes.alloc_PrEP_Oral_MSM_H;
                Allocation.alloc_PrEP_Oral_MSM_O = AllOutcomes.alloc_PrEP_Oral_MSM_O;
                Allocation.alloc_PrEP_Oral_IDU_B = AllOutcomes.alloc_PrEP_Oral_IDU_B;
                Allocation.alloc_PrEP_Oral_IDU_H = AllOutcomes.alloc_PrEP_Oral_IDU_H;
                Allocation.alloc_PrEP_Oral_IDU_O = AllOutcomes.alloc_PrEP_Oral_IDU_O;
                Allocation.alloc_PrEP_Inject_HETM_B = AllOutcomes.alloc_PrEP_Inject_HETM_B;
                Allocation.alloc_PrEP_Inject_HETM_H = AllOutcomes.alloc_PrEP_Inject_HETM_H;
                Allocation.alloc_PrEP_Inject_HETM_O = AllOutcomes.alloc_PrEP_Inject_HETM_O;
                Allocation.alloc_PrEP_Inject_HETF_B = AllOutcomes.alloc_PrEP_Inject_HETF_B;
                Allocation.alloc_PrEP_Inject_HETF_H = AllOutcomes.alloc_PrEP_Inject_HETF_H;
                Allocation.alloc_PrEP_Inject_HETF_O = AllOutcomes.alloc_PrEP_Inject_HETF_O;
                Allocation.alloc_PrEP_Inject_MSM_B = AllOutcomes.alloc_PrEP_Inject_MSM_B;
                Allocation.alloc_PrEP_Inject_MSM_H = AllOutcomes.alloc_PrEP_Inject_MSM_H;
                Allocation.alloc_PrEP_Inject_MSM_O = AllOutcomes.alloc_PrEP_Inject_MSM_O;
                Allocation.alloc_PrEP_Inject_IDU_B = AllOutcomes.alloc_PrEP_Inject_IDU_B;
                Allocation.alloc_PrEP_Inject_IDU_H = AllOutcomes.alloc_PrEP_Inject_IDU_H;
                Allocation.alloc_PrEP_Inject_IDU_O = AllOutcomes.alloc_PrEP_Inject_IDU_O;                
                % Collect treatment and care costs in each allocation period
                sum_TotalARTAndCareCost_Undisc(Params.numAllocationPeriods)=0;
                Intn_Allocation_EndYr(Params.numAllocationPeriods)=0;
                for nAllocPeriod = 1:Params.numAllocationPeriods
                    if nAllocPeriod < Params.numAllocationPeriods
                        Intn_Allocation_EndYr(nAllocPeriod)=Params.Intn_Allocation_StartYr(nAllocPeriod+1);
                    else
                        Intn_Allocation_EndYr(nAllocPeriod)=min(Params.outcomeCollectionEndYr, ...
                            Params.tt_modelStartYear + Params.tt_timeHorizon - 1);
                    end
                end
                for nOutcomeYear = Params.outcomeCollectionStartYr:Params.outcomeCollectionEndYr
                    for nAllocPeriod = 1:Params.numAllocationPeriods
                        if and(nOutcomeYear>=Params.Intn_Allocation_StartYr(nAllocPeriod),...
                                nOutcomeYear<Intn_Allocation_EndYr(nAllocPeriod))
                            sum_TotalARTAndCareCost_Undisc(nAllocPeriod) = ...
                                sum_TotalARTAndCareCost_Undisc(nAllocPeriod) + ...
                                sum(AllOutcomes.ann_TotalARTAndCareCost_Undisc(nOutcomeYear-Params.outcomeCollectionStartYr+1,:),2);
                        end
                    end
                end  
                for nAllocPeriod = 1:Params.numAllocationPeriods
                    Allocation.mean_annARTAndCareCost_Undisc_PaidByThisBudget_inM(nAllocPeriod,1) = ...
                    	sum_TotalARTAndCareCost_Undisc(nAllocPeriod) / ...
                        (Intn_Allocation_EndYr(nAllocPeriod)-Params.Intn_Allocation_StartYr(nAllocPeriod)+1) ...
                        * Params.alloc_PctARTAndCareCoveredByBudget/Million;
                end

                AnnOutcomes.ann_undiscNewInfectionsHET = AllOutcomes.ann_undiscNewInfectionsHET;
                SumOutcomes.sum_undiscNewInfectionsHET = sum(AllOutcomes.ann_undiscNewInfectionsHET);
                AnnOutcomes.ann_undiscNewInfectionsMSM = AllOutcomes.ann_undiscNewInfectionsMSM;
                SumOutcomes.sum_undiscNewInfectionsMSM = sum(AllOutcomes.ann_undiscNewInfectionsMSM);
                AnnOutcomes.ann_undiscNewInfectionsIDU = AllOutcomes.ann_undiscNewInfectionsIDU;
                SumOutcomes.sum_undiscNewInfectionsIDU = sum(AllOutcomes.ann_undiscNewInfectionsIDU);
                AnnOutcomes.ann_undiscNewInfections = AllOutcomes.ann_undiscNewInfections;
                SumOutcomes.sum_undiscNewInfections = sum(AllOutcomes.ann_undiscNewInfections);
                SumOutcomes.sum_discNewInfections = AllOutcomes.TotalNewInfections_Disc;
                
                if FirstOutcomeYr <= Params.tt_periodFiveStartYear - OutcomesIndexTarget0 && (LastOutcomeYr >= Params.tt_periodFiveStartYear + OutcomesIndexTarget1 - OutcomesIndexTarget0 - 1)
                   
                    SumOutcomes.IncidenceReductionTarget1 = AllOutcomes.IncidenceReductionTarget1;
                    
                end
                
                if FirstOutcomeYr <= Params.tt_periodFiveStartYear - OutcomesIndexTarget0 && (LastOutcomeYr >= Params.tt_periodFiveStartYear + OutcomesIndexTarget2 - OutcomesIndexTarget0 - 1)
                    
                    SumOutcomes.IncidenceReductionTarget2 = AllOutcomes.IncidenceReductionTarget2;
                    
                end
                
                SumOutcomes.sum_discQALYsHET_inM = AllOutcomes.total_discQALYs_pop(Params.pop_HET) / Million;
                SumOutcomes.sum_discQALYsMSM_inM = AllOutcomes.total_discQALYs_pop(Params.pop_MSM) / Million;
                SumOutcomes.sum_discQALYsIDU_inM = AllOutcomes.total_discQALYs_pop(Params.pop_IDU) / Million;
                AnnOutcomes.ann_discQALYs_inM = AllOutcomes.ann_discQALYs / Million;
                SumOutcomes.sum_discQALYs_inM = sum(AllOutcomes.ann_discQALYs) / Million;
                AnnOutcomes.ann_TotalARTAndCareCost_inM = sum(AllOutcomes.ann_TotalARTAndCareCost_Disc,2)/Million;
                SumOutcomes.sum_TotalARTAndCareCost_inM = sum(sum(AllOutcomes.ann_TotalARTAndCareCost_Disc,2))/Million;
                SumOutcomes.sum_TotalARTAndCareCost_PaidByThisBudget_inM = ...
                    SumOutcomes.sum_TotalARTAndCareCost_inM*Params.alloc_PctARTAndCareCoveredByBudget;
                AnnOutcomes.ann_TotalPrEPCost_inclInAlloc_inM = sum(AllOutcomes.ann_HealthStateCost_PrEP_Disc,2)/Million;
                SumOutcomes.sum_TotalPrEPCost_inclInAlloc_inM = sum(sum(AllOutcomes.ann_HealthStateCost_PrEP_Disc,2))/Million;

                if FirstOutcomeYr <= rateYr && (LastOutcomeYr >= rateYr)
                   
                    SumOutcomes.tr_LRHAllEligTestProb = AllOutcomes.tr_LRHAllEligTestProb_AllYrs(rateYrIndex);
                    SumOutcomes.tr_HRHAllEligTestProb = AllOutcomes.tr_HRHAllEligTestProb_AllYrs(rateYrIndex);
                    SumOutcomes.tr_LRMSMAllEligTestProb = AllOutcomes.tr_LRMSMAllEligTestProb_AllYrs(rateYrIndex);
                    SumOutcomes.tr_HRMSMAllEligTestProb = AllOutcomes.tr_HRMSMAllEligTestProb_AllYrs(rateYrIndex);
                    SumOutcomes.tr_IDUAllEligTestProb = AllOutcomes.tr_IDUAllEligTestProb_AllYrs(rateYrIndex);
                    SumOutcomes.ann_PctLinkToCareFirst = AllOutcomes.ann_PctLinkToCareFirst(rateYrIndex);
                    SumOutcomes.ann_LinkageAfterProb = AllOutcomes.ann_LinkageAfterProb(rateYrIndex);
                    SumOutcomes.ann_ARTPrescrProb = AllOutcomes.ann_ARTPrescrProb(rateYrIndex);
                    SumOutcomes.ann_VLSToANVProb = AllOutcomes.ann_VLSToANVProb(rateYrIndex);
                    SumOutcomes.ann_ANVToVLSProb = AllOutcomes.ann_ANVToVLSProb(rateYrIndex);
                    SumOutcomes.ann_PctUninfPWIDServedbySEP = AllOutcomes.ann_PctUninfPWIDServedbySEP(rateYrIndex);
                    SumOutcomes.ann_PrEPInitProb_HighRiskHETs = AllOutcomes.ann_PrEPInitProb_HighRiskHETs(rateYrIndex);
                    SumOutcomes.ann_PrEPInitProb_HighRiskMSM = AllOutcomes.ann_PrEPInitProb_HighRiskMSM(rateYrIndex);
                    SumOutcomes.ann_PrEPInitProb_IDUs = AllOutcomes.ann_PrEPInitProb_IDUs(rateYrIndex);
                end
                
                SumOutcomes.ann_lastyrDiagnosed = AllOutcomes.ann_PctAware(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrDiagnosed_HRHET = AllOutcomes.ann_PctAware_HRHET(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrDiagnosed_OHET = AllOutcomes.ann_PctAware_OHET(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrDiagnosed_HRMSM = AllOutcomes.ann_PctAware_HRMSM(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrDiagnosed_OMSM = AllOutcomes.ann_PctAware_OMSM(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrDiagnosed_PWID = AllOutcomes.ann_PctAware_PWID(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrLinkedtoCare = AllOutcomes.ann_PctInCare(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrOnART = AllOutcomes.ann_PctOnART(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrVLS = AllOutcomes.ann_PctVLS(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrVLSamongdiag = AllOutcomes.ann_PctVLSamongdiag(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrPctUninfPWIDServedbySEP = AllOutcomes.ann_PctUninfPWIDServedbySEP(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrPctEligOnPrEP_HRHET = AllOutcomes.ann_PctEligOnPrEP_HRHET(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrPctEligOnPrEP_HRMSM = AllOutcomes.ann_PctEligOnPrEP_HRMSM(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrPctEligOnPrEP_PWID = AllOutcomes.ann_PctEligOnPrEP_PWID(LastOutcomeYrIndex);
                AnnOutcomes.ann_TotalARTAndCareCost_Undisc = sum(AllOutcomes.ann_TotalARTAndCareCost_Undisc,2);
                SumOutcomeFields = fieldnames(SumOutcomes);

                SelectedOutcomes = SumOutcomes;
                
                for nAllocPeriod = 1:Params.numAllocationPeriods
                    SelectedOutcomes.TotalIntnBudget_inM(nAllocPeriod,1) = Params.alloc_TotalIntnBudget(nAllocPeriod);
                    SelectedOutcomes.summary(nAllocPeriod,1) = Params.alloc_TotalIntnBudget(nAllocPeriod);
                end
                SummaryRows = Params.numAllocationPeriods;
                
                SelectedOutcomes.summary(SummaryRows+1,1) = ...
                    Params.alloc_PctARTAndCareCoveredByBudget;
                SummaryRows = SummaryRows + 1;
                
                MaxReachFields = fieldnames(MaxReach);
                nMaxReachFields = length(MaxReachFields);
                for nMRCount=1:nMaxReachFields
                    SelectedOutcomes.MaxReach(nMRCount,1) = MaxReach.(MaxReachFields{nMRCount});
                    SelectedOutcomes.summary(SummaryRows+nMRCount,1) = MaxReach.(MaxReachFields{nMRCount});
                end
                SummaryRows = SummaryRows + nMaxReachFields;
                
                SelectedOutcomes.periodStartYrs = Params.Intn_Allocation_StartYr(1:Params.numAllocationPeriods)';
                SelectedOutcomes.summary ...
                    ((SummaryRows+1):(SummaryRows+Params.numAllocationPeriods),1)= ...
                    SelectedOutcomes.periodStartYrs;
                SummaryRows = SummaryRows + Params.numAllocationPeriods;
                
                AllocFields = fieldnames(Allocation);
                nAllocFields = length(AllocFields);
                for nAFCount=1:nAllocFields
                    for nAllocPeriod = 1:Params.numAllocationPeriods
                        SelectedOutcomes.Allocation_inM(nAllocFields*(nAllocPeriod-1)+nAFCount,1) = Allocation.(AllocFields{nAFCount})(nAllocPeriod);
                        SelectedOutcomes.summary(SummaryRows+nAllocFields*(nAllocPeriod-1)+nAFCount,1) = Allocation.(AllocFields{nAFCount})(nAllocPeriod);
                    end
                end
                SummaryRows = SummaryRows + nAllocFields*Params.numAllocationPeriods;
                
                for SOFCount=1:length(SumOutcomeFields)
                    SelectedOutcomes.summary(SummaryRows + SOFCount,1) = ...
                        SumOutcomes.(SumOutcomeFields{SOFCount});
                end
                SummaryRows = SummaryRows + length(SumOutcomeFields);
                
                SelectedOutcomes.summary(SummaryRows+1:SummaryRows + length(AnnOutcomes.ann_TotalARTAndCareCost_Undisc)) = ...
                    AnnOutcomes.ann_TotalARTAndCareCost_Undisc/Million;
                
                SelectedOutcomes.AnnualOutcomes = AnnOutcomes;
                
            case 5 % NHAS analysis outcomes 
                
                SelectedOutcomes.TotalNewInfections = sum(AllOutcomes.ann_TotalNewInfections,2);
                SelectedOutcomes.NewInfections_Blk = sum(AllOutcomes.ann_NewInfections_Blk,2);
                SelectedOutcomes.NewInfections_Hisp = sum(AllOutcomes.ann_NewInfections_Hisp,2);
                SelectedOutcomes.NewInfections_Oth = sum(AllOutcomes.ann_NewInfections_Oth,2);
                SelectedOutcomes.HIVPrevalence = sum(AllOutcomes.ann_HIVPrevalence,2);
                SelectedOutcomes.HIVPrevalence_Blk = sum(AllOutcomes.ann_HIVPrevalence_Blk,2);
                SelectedOutcomes.HIVPrevalence_Hisp = sum(AllOutcomes.ann_HIVPrevalence_Hisp,2);
                SelectedOutcomes.HIVPrevalence_Oth = sum(AllOutcomes.ann_HIVPrevalence_Oth,2);
                SelectedOutcomes.pctOfPLWH_Undiagnosed_2020_Blk = AllOutcomes.nhas_pctOfPLWH_Undiagnosed_2020_Blk;
                SelectedOutcomes.pctOfPLWH_Aware_2020_Blk = AllOutcomes.nhas_pctOfPLWH_Aware_2020_Blk;
                SelectedOutcomes.pctOfPLWH_InCare_2020_Blk = AllOutcomes.nhas_pctOfPLWH_InCare_2020_Blk;
                SelectedOutcomes.pctOfPLWH_OnART_2020_Blk = AllOutcomes.nhas_pctOfPLWH_OnART_2020_Blk;
                SelectedOutcomes.pctOfPLWH_VLS_2020_Blk = AllOutcomes.nhas_pctOfPLWH_VLS_2020_Blk;
                SelectedOutcomes.pctOfPLWH_Undiagnosed_2020_Hisp = AllOutcomes.nhas_pctOfPLWH_Undiagnosed_2020_Hisp;
                SelectedOutcomes.pctOfPLWH_Aware_2020_Hisp = AllOutcomes.nhas_pctOfPLWH_Aware_2020_Hisp;
                SelectedOutcomes.pctOfPLWH_InCare_2020_Hisp = AllOutcomes.nhas_pctOfPLWH_InCare_2020_Hisp;
                SelectedOutcomes.pctOfPLWH_OnART_2020_Hisp = AllOutcomes.nhas_pctOfPLWH_OnART_2020_Hisp;
                SelectedOutcomes.pctOfPLWH_VLS_2020_Hisp = AllOutcomes.nhas_pctOfPLWH_VLS_2020_Hisp;
                SelectedOutcomes.pctOfPLWH_Undiagnosed_2020_Oth = AllOutcomes.nhas_pctOfPLWH_Undiagnosed_2020_Oth;
                SelectedOutcomes.pctOfPLWH_Aware_2020_Oth = AllOutcomes.nhas_pctOfPLWH_Aware_2020_Oth;
                SelectedOutcomes.pctOfPLWH_InCare_2020_Oth = AllOutcomes.nhas_pctOfPLWH_InCare_2020_Oth;
                SelectedOutcomes.pctOfPLWH_OnART_2020_Oth = AllOutcomes.nhas_pctOfPLWH_OnART_2020_Oth;
                SelectedOutcomes.pctOfPLWH_VLS_2020_Oth = AllOutcomes.nhas_pctOfPLWH_VLS_2020_Oth;
                SelectedOutcomes.nhas_pctOfPLWH_Undiagnosed_2020_Total = AllOutcomes.nhas_pctOfPLWH_Undiagnosed_2020_Total;
                SelectedOutcomes.pctOfPLWH_Aware_2020_Total = AllOutcomes.nhas_pctOfPLWH_Aware_2020_Total;
                SelectedOutcomes.pctOfPLWH_InCare_2020_Total = AllOutcomes.nhas_pctOfPLWH_InCare_2020_Total;
                SelectedOutcomes.pctOfPLWH_OnART_2020_Total = AllOutcomes.nhas_pctOfPLWH_OnART_2020_Total;
                SelectedOutcomes.pctOfPLWH_VLS_2020_Total = AllOutcomes.nhas_pctOfPLWH_VLS_2020_Total;
                SelectedOutcomes.pctOfAware_InCare_2020_Blk = AllOutcomes.nhas_pctOfAware_InCare_2020_Blk;
                SelectedOutcomes.pctOfAware_OnART_2020_Blk = AllOutcomes.nhas_pctOfAware_OnART_2020_Blk;
                SelectedOutcomes.pctOfAware_VLS_2020_Blk = AllOutcomes.nhas_pctOfAware_VLS_2020_Blk;
                SelectedOutcomes.pctOfAware_InCare_2020_Hisp = AllOutcomes.nhas_pctOfAware_InCare_2020_Hisp;
                SelectedOutcomes.pctOfAware_OnART_2020_Hisp = AllOutcomes.nhas_pctOfAware_OnART_2020_Hisp;
                SelectedOutcomes.pctOfAware_VLS_2020_Hisp = AllOutcomes.nhas_pctOfAware_VLS_2020_Hisp;
                SelectedOutcomes.pctOfAware_InCare_2020_Oth = AllOutcomes.nhas_pctOfAware_InCare_2020_Oth;
                SelectedOutcomes.pctOfAware_OnART_2020_Oth = AllOutcomes.nhas_pctOfAware_OnART_2020_Oth;
                SelectedOutcomes.pctOfAware_VLS_2020_Oth = AllOutcomes.nhas_pctOfAware_VLS_2020_Oth;
                SelectedOutcomes.pctOfAware_InCare_2020_Total = AllOutcomes.nhas_pctOfAware_InCare_2020_Total;
                SelectedOutcomes.pctOfAware_OnART_2020_Total = AllOutcomes.nhas_pctOfAware_OnART_2020_Total;
                SelectedOutcomes.pctOfAware_VLS_2020_Total = AllOutcomes.nhas_pctOfAware_VLS_2020_Total;
                SelectedOutcomes.pctOfAffByART_VLS_2020_Blk = AllOutcomes.nhas_pctOfAffByART_VLS_2020_Blk;
                SelectedOutcomes.pctOfAffByART_VLS_2020_Hisp = AllOutcomes.nhas_pctOfAffByART_VLS_2020_Hisp;
                SelectedOutcomes.pctOfAffByART_VLS_2020_Oth = AllOutcomes.nhas_pctOfAffByART_VLS_2020_Oth;
                SelectedOutcomes.pctOfAffByART_VLS_2020_Total = AllOutcomes.nhas_pctOfAffByART_VLS_2020_Total;
                SelectedOutcomes.NumberUndiagnosed_2020_Blk = AllOutcomes.nhas_NumberUndiagnosed_2020_Blk;
                SelectedOutcomes.NumberAware_2020_Blk = AllOutcomes.nhas_NumberAware_2020_Blk ;
                SelectedOutcomes.NumberInCare_2020_Blk = AllOutcomes.nhas_NumberInCare_2020_Blk;
                SelectedOutcomes.NumberOnART_2020_Blk = AllOutcomes.nhas_NumberOnART_2020_Blk;
                SelectedOutcomes.NumberVLS_2020_Blk = AllOutcomes.nhas_NumberVLS_2020_Blk;
                SelectedOutcomes.NumberUndiagnosed_2020_Hisp = AllOutcomes.nhas_NumberUndiagnosed_2020_Hisp;
                SelectedOutcomes.NumberAware_2020_Hisp = AllOutcomes.nhas_NumberAware_2020_Hisp ;
                SelectedOutcomes.NumberInCare_2020_Hisp = AllOutcomes.nhas_NumberInCare_2020_Hisp;
                SelectedOutcomes.NumberOnART_2020_Hisp = AllOutcomes.nhas_NumberOnART_2020_Hisp;
                SelectedOutcomes.NumberVLS_2020_Hisp = AllOutcomes.nhas_NumberVLS_2020_Hisp;
                SelectedOutcomes.NumberUndiagnosed_2020_Oth = AllOutcomes.nhas_NumberUndiagnosed_2020_Oth;
                SelectedOutcomes.NumberAware_2020_Oth = AllOutcomes.nhas_NumberAware_2020_Oth;
                SelectedOutcomes.NumberInCare_2020_Oth = AllOutcomes.nhas_NumberInCare_2020_Oth;
                SelectedOutcomes.NumberOnART_2020_Oth = AllOutcomes.nhas_NumberOnART_2020_Oth;
                SelectedOutcomes.NumberVLS_2020_Oth = AllOutcomes.nhas_NumberVLS_2020_Oth;
                SelectedOutcomes.NumberUndiagnosed_2020_Total = AllOutcomes.nhas_NumberUndiagnosed_2020_Total;
                SelectedOutcomes.NumberAware_2020_Total = AllOutcomes.nhas_NumberAware_2020_Total;
                SelectedOutcomes.NumberInCare_2020_Total = AllOutcomes.nhas_NumberInCare_2020_Total;
                SelectedOutcomes.NumberOnART_2020_Total = AllOutcomes.nhas_NumberOnART_2020_Total;
                SelectedOutcomes.NumberVLS_2020_Total = AllOutcomes.nhas_NumberVLS_2020_Total;
                SelectedOutcomes.set_FirstOutcomeYr = AllOutcomes.set_FirstOutcomeYr; 
                SelectedOutcomes.Key_PctDiagByRace =  ...
                    [AllOutcomes.nhas_pctOfPLWH_Aware_2020_Blk;...
                    AllOutcomes.nhas_pctOfPLWH_Aware_2020_Hisp;...
                    AllOutcomes.nhas_pctOfPLWH_Aware_2020_Oth];
                SelectedOutcomes.Key_PctVLSAmongDiagByRace =  ...
                    [AllOutcomes.nhas_pctOfAware_VLS_2020_Blk;...
                    AllOutcomes.nhas_pctOfAware_VLS_2020_Hisp;...
                    AllOutcomes.nhas_pctOfAware_VLS_2020_Oth];
                
            case 6 %Testing frequency outcomes
                
                SelectedOutcomes.ann_LifeYears_HET = AllOutcomes.ann_LifeYears_HET;
                SelectedOutcomes.ann_QALYs_HET = AllOutcomes.ann_QALYs_HET;
                SelectedOutcomes.ann_CostTestAndNotify_HET_Disc = AllOutcomes.ann_CostTestAndNotify_HET_Disc;
                SelectedOutcomes.ann_ARTCarePrEPCost_Disc_HET = AllOutcomes.ann_ARTCarePrEPCost_Disc_HET; 
                SelectedOutcomes.ann_ARTCarePrEPTransitionAndSEPCost_Disc = AllOutcomes.ann_ARTCarePrEPTransitionAndSEPCost_Disc; 
                SelectedOutcomes.total_ARTCarePrEPTransitionAndSEPCost_Disc = AllOutcomes.total_ARTCarePrEPTransitionAndSEPCost_Disc;
                SelectedOutcomes.ann_DistHRHByComparts_LastYr = AllOutcomes.ann_DistHRHByComparts_LastYr;
                
                if rateYrIndex > -1 
                    SelectedOutcomes.tr_WeightedAllEligTestRateByCohort = AllOutcomes.tr_WeightedAllEligTestRateByCohort;
                    SelectedOutcomes.tr_raceAllEligTestRate = AllOutcomes.tr_raceAllEligTestRate; 
                    SelectedOutcomes.tr_riskGroupAllEligTestRate = AllOutcomes.tr_riskGroupAllEligTestRate; 
                    SelectedOutcomes.tr_HRHAllEligTestRate = AllOutcomes.tr_HRHAllEligTestRate; 
                    SelectedOutcomes.tr_LRHAllEligTestRate = AllOutcomes.tr_LRHAllEligTestRate; 
                    SelectedOutcomes.tr_HRMSMAllEligTestRate = AllOutcomes.tr_HRMSMAllEligTestRate; 
                    SelectedOutcomes.tr_LRMSMAllEligTestRate = AllOutcomes.tr_LRMSMAllEligTestRate; 
                    SelectedOutcomes.tr_IDUAllEligTestRate = AllOutcomes.tr_IDUAllEligTestRate; 
                    SelectedOutcomes.tr_AllEligTestRateByRiskGpLevel = AllOutcomes.tr_AllEligTestRateByRiskGpLevel;
                    SelectedOutcomes.tr_WeightedHIVPosTestRateByCohort = AllOutcomes.tr_WeightedHIVPosTestRateByCohort; 
                    SelectedOutcomes.tr_raceHIVPosTestRate = AllOutcomes.tr_raceHIVPosTestRate; 
                    SelectedOutcomes.tr_riskGroupHIVPosTestRate = AllOutcomes.tr_riskGroupHIVPosTestRate; 
                    SelectedOutcomes.tr_HRHHIVPosTestRate = AllOutcomes.tr_HRHHIVPosTestRate; 
                    SelectedOutcomes.tr_LRHHIVPosTestRate = AllOutcomes.tr_LRHHIVPosTestRate; 
                    SelectedOutcomes.tr_uninfectedTestRate = AllOutcomes.tr_uninfectedTestRate; 
                    SelectedOutcomes.tr_diseaseStageTestRate = AllOutcomes.tr_diseaseStageTestRate; 
                    SelectedOutcomes.tr_raceDiagRate = AllOutcomes.tr_raceDiagRate; 
                    SelectedOutcomes.tr_riskGroupDiagRate = AllOutcomes.tr_riskGroupDiagRate; 
                    SelectedOutcomes.tr_diseaseStageDiagRate = AllOutcomes.tr_diseaseStageDiagRate; 
                    SelectedOutcomes.tr_LRHDiagRate = AllOutcomes.tr_LRHDiagRate; 
                    SelectedOutcomes.tr_HRHDiagRate = AllOutcomes.tr_HRHDiagRate; 
                end
                SelectedOutcomes.ann_CostTestAndNotify = AllOutcomes.ann_TransCost_CostTestAndNotify_Disc;
                SelectedOutcomes.ann_ARTCarePrEPCost_Disc = AllOutcomes.ann_ARTCarePrEPCost_Disc;
                SelectedOutcomes.ann_NewInfectionsHET = AllOutcomes.ann_NewInfectionsHET;
                SelectedOutcomes.ann_NewInfectionsHET_LR = sum(AllOutcomes.ann_NewInfections_LowRiskHETs,2);
                SelectedOutcomes.ann_NewInfectionsHET_HR = sum(AllOutcomes.ann_NewInfections_HighRiskHETs,2);
                SelectedOutcomes.ann_TotalNewInfections = AllOutcomes.ann_TotalNewInfections;
                SelectedOutcomes.TotalNewInfections = sum(AllOutcomes.ann_TotalNewInfections,2);
                SelectedOutcomes.ann_LifeYears = AllOutcomes.ann_LifeYears;
                SelectedOutcomes.ann_QALYs = AllOutcomes.ann_QALYs;
                SelectedOutcomes.ann_CostTestAndNotify = AllOutcomes.ann_TransCost_CostTestAndNotify_Disc;
                SelectedOutcomes.ann_ARTCarePrEPCost_Disc = AllOutcomes.ann_ARTCarePrEPCost_Disc;
                SelectedOutcomes.ann_TotalNewDiagnosesHET = AllOutcomes.ann_TotalNewDiagnosesHET;
                SelectedOutcomes.ann_TotalNewDiagnoses = AllOutcomes.ann_TotalNewDiagnoses;
                SelectedOutcomes.HIVPrevalence = sum(AllOutcomes.ann_HIVPrevalence,2);
                SelectedOutcomes.HIVPrevalenceHET_LR = sum(AllOutcomes.ann_HIVPrevalence_LowRiskHETs,2);
                SelectedOutcomes.HIVPrevalenceHET_HR = sum(AllOutcomes.ann_HIVPrevalence_HighRiskHETs,2);                
                
                if AllOutcomes.set_FirstOutcomeYr <=2006 &&(AllOutcomes.set_LastOutcomeYr >= 2015)
                    SelectedOutcomes.calib_LRHETPrevalence2015 = AllOutcomes.calib_LRHETPrevalence2015;
                    SelectedOutcomes.calib_LRHETPrevalence2006 = AllOutcomes.calib_LRHETPrevalence2006;
                    SelectedOutcomes.calib_HETPrev2015v2006_LR = AllOutcomes.calib_HETPrev2015v2006_LR;
                    SelectedOutcomes.calib_HRHETPrevalence2015 = AllOutcomes.calib_HRHETPrevalence2015;
                    SelectedOutcomes.calib_HRHETPrevalence2006 = AllOutcomes.calib_HRHETPrevalence2006;
                    SelectedOutcomes.calib_HETPrev2015v2006_HR = AllOutcomes.calib_HETPrev2015v2006_HR;
                end
                
                % Infections    
                SelectedOutcomes.TotalNewInfections_LowRiskHETs = sum(sum(AllOutcomes.ann_TotalNewInfections) .* Params.LowRiskHETIndicator');
                SelectedOutcomes.TotalNewInfections_HighRiskHETs = sum(sum(AllOutcomes.ann_TotalNewInfections) .* Params.HRHIndicator');

                % Diagnoses
                SelectedOutcomes.TotalNewDiagnoses_LowRiskHETs = sum(sum(AllOutcomes.ann_TotalNewDiagnoses) .* Params.LowRiskHETIndicator');
                SelectedOutcomes.TotalNewDiagnoses_HighRiskHETs = sum(sum(AllOutcomes.ann_TotalNewDiagnoses) .* Params.HRHIndicator');

                % Life years
                SelectedOutcomes.TotalLifeYears_LowRiskHETs = sum(sum(AllOutcomes.ann_LifeYears_HET) .* Params.LowRiskHETIndicator');
                SelectedOutcomes.TotalLifeYears_HighRiskHETs = sum(sum(AllOutcomes.ann_LifeYears_HET) .* Params.HRHIndicator');

                % QALYs
                SelectedOutcomes.TotalQALYs_LowRiskHETs = sum(sum(AllOutcomes.ann_QALYs_HET).* Params.LowRiskHETIndicator');
                SelectedOutcomes.TotalQALYs_HighRiskHETs = sum(sum(AllOutcomes.ann_QALYs_HET) .* Params.HRHIndicator');

                % Testing Costs
                SelectedOutcomes.TotalCostTestNotify_LowRiskHETs = sum(sum(AllOutcomes.ann_TransCost_CostTestAndNotify_Disc).* Params.LowRiskHETIndicator');
                SelectedOutcomes.TotalCostTestNotify_HighRiskHETs = sum(sum(AllOutcomes.ann_TransCost_CostTestAndNotify_Disc) .* Params.HRHIndicator');

                % Health state costs
                SelectedOutcomes.TotalCostHealthState_LowRiskHETs = sum(sum(AllOutcomes.ann_ARTCarePrEPCost_Disc) .* Params.LowRiskHETIndicator');
                SelectedOutcomes.TotalCostHealthState_HighRiskHETs = sum(sum(AllOutcomes.ann_ARTCarePrEPCost_Disc) .* Params.HRHIndicator');

                % Total costs
                SelectedOutcomes.ARTCarePrEPTransitionAndSEPCost_LowRiskHETs = sum(sum(AllOutcomes.ann_ARTCarePrEPTransitionAndSEPCost_Disc_noSEP) .* Params.LowRiskHETIndicator');
                SelectedOutcomes.ARTCarePrEPTransitionAndSEPCost_HighRiskHETs = sum(sum(AllOutcomes.ann_ARTCarePrEPTransitionAndSEPCost_Disc_noSEP) .* Params.HRHIndicator');

            case 7 %Young MSM analysis outcomes            
                 
                %Investment
                SelectedOutcomes.alloc_Testing_HET_Low = AllOutcomes.alloc_Testing_HET_Low;
                SelectedOutcomes.alloc_Testing_HET_High = AllOutcomes.alloc_Testing_HET_High;
                SelectedOutcomes.alloc_Testing_MSM_Low = AllOutcomes.alloc_Testing_MSM_Low;
                SelectedOutcomes.alloc_Testing_MSM_High = AllOutcomes.alloc_Testing_MSM_High;
                SelectedOutcomes.alloc_Testing_IDU = AllOutcomes.alloc_Testing_IDU;
                SelectedOutcomes.alloc_LTCatDiag = AllOutcomes.alloc_LTCatDiag;
                SelectedOutcomes.alloc_LTCafterDiag = AllOutcomes.alloc_LTCafterDiag;
                SelectedOutcomes.alloc_ARTInitiation = AllOutcomes.alloc_ARTInitiation;
                SelectedOutcomes.alloc_ARTAdher5to4 = AllOutcomes.alloc_ARTAdher5to4;
                SelectedOutcomes.alloc_ARTAdher4to5 = AllOutcomes.alloc_ARTAdher4to5;
                SelectedOutcomes.alloc_SEP_B = AllOutcomes.alloc_SEP_B;
                SelectedOutcomes.alloc_SEP_H = AllOutcomes.alloc_SEP_H;
                SelectedOutcomes.alloc_SEP_O = AllOutcomes.alloc_SEP_O;
                SelectedOutcomes.alloc_PrEP_Oral_HETM_B = AllOutcomes.alloc_PrEP_Oral_HETM_B;
                SelectedOutcomes.alloc_PrEP_Oral_HETM_H = AllOutcomes.alloc_PrEP_Oral_HETM_H;
                SelectedOutcomes.alloc_PrEP_Oral_HETM_O = AllOutcomes.alloc_PrEP_Oral_HETM_O;
                SelectedOutcomes.alloc_PrEP_Oral_HETF_B = AllOutcomes.alloc_PrEP_Oral_HETF_B;
                SelectedOutcomes.alloc_PrEP_Oral_HETF_H = AllOutcomes.alloc_PrEP_Oral_HETF_H;
                SelectedOutcomes.alloc_PrEP_Oral_HETF_O = AllOutcomes.alloc_PrEP_Oral_HETF_O;
                SelectedOutcomes.alloc_PrEP_Oral_MSM_B = AllOutcomes.alloc_PrEP_Oral_MSM_B;
                SelectedOutcomes.alloc_PrEP_Oral_MSM_H = AllOutcomes.alloc_PrEP_Oral_MSM_H;
                SelectedOutcomes.alloc_PrEP_Oral_MSM_O = AllOutcomes.alloc_PrEP_Oral_MSM_O;
                SelectedOutcomes.alloc_PrEP_Oral_IDU_B = AllOutcomes.alloc_PrEP_Oral_IDU_B;
                SelectedOutcomes.alloc_PrEP_Oral_IDU_H = AllOutcomes.alloc_PrEP_Oral_IDU_H;
                SelectedOutcomes.alloc_PrEP_Oral_IDU_O = AllOutcomes.alloc_PrEP_Oral_IDU_O;
                SelectedOutcomes.alloc_PrEP_Inject_HETM_B = AllOutcomes.alloc_PrEP_Inject_HETM_B;
                SelectedOutcomes.alloc_PrEP_Inject_HETM_H = AllOutcomes.alloc_PrEP_Inject_HETM_H;
                SelectedOutcomes.alloc_PrEP_Inject_HETM_O = AllOutcomes.alloc_PrEP_Inject_HETM_O;
                SelectedOutcomes.alloc_PrEP_Inject_HETF_B = AllOutcomes.alloc_PrEP_Inject_HETF_B;
                SelectedOutcomes.alloc_PrEP_Inject_HETF_H = AllOutcomes.alloc_PrEP_Inject_HETF_H;
                SelectedOutcomes.alloc_PrEP_Inject_HETF_O = AllOutcomes.alloc_PrEP_Inject_HETF_O;
                SelectedOutcomes.alloc_PrEP_Inject_MSM_B = AllOutcomes.alloc_PrEP_Inject_MSM_B;
                SelectedOutcomes.alloc_PrEP_Inject_MSM_H = AllOutcomes.alloc_PrEP_Inject_MSM_H;
                SelectedOutcomes.alloc_PrEP_Inject_MSM_O = AllOutcomes.alloc_PrEP_Inject_MSM_O;
                SelectedOutcomes.alloc_PrEP_Inject_IDU_B = AllOutcomes.alloc_PrEP_Inject_IDU_B;
                SelectedOutcomes.alloc_PrEP_Inject_IDU_H = AllOutcomes.alloc_PrEP_Inject_IDU_H;
                SelectedOutcomes.alloc_PrEP_Inject_IDU_O = AllOutcomes.alloc_PrEP_Inject_IDU_O;                
                
                %Young MSM (18-24 and 25-34) incidence by age group & race/eth
                SelectedOutcomes.YMSM_NewInfections_MSM_18_24_B = AllOutcomes.YMSM_NewInfections_MSM_18_24_B;
                SelectedOutcomes.YMSM_NewInfections_MSM_18_24_H = AllOutcomes.YMSM_NewInfections_MSM_18_24_H;
                SelectedOutcomes.YMSM_NewInfections_MSM_18_24_O = AllOutcomes.YMSM_NewInfections_MSM_18_24_O;
                SelectedOutcomes.YMSM_NewInfections_MSM_25_34_B = AllOutcomes.YMSM_NewInfections_MSM_25_34_B;
                SelectedOutcomes.YMSM_NewInfections_MSM_25_34_H = AllOutcomes.YMSM_NewInfections_MSM_25_34_H;
                SelectedOutcomes.YMSM_NewInfections_MSM_25_34_O = AllOutcomes.YMSM_NewInfections_MSM_25_34_O;                
                %Total incidence by risk gp and overall
                SelectedOutcomes.YMSM_NewInfectionsHET = sum(sum(AllOutcomes.ann_NewInfectionsHET,2));
                SelectedOutcomes.YMSM_NewInfectionsMSM = sum(sum(AllOutcomes.ann_NewInfectionsMSM,2));
                SelectedOutcomes.YMSM_NewInfectionsIDU = sum(sum(AllOutcomes.ann_NewInfectionsIDU,2)); 
                SelectedOutcomes.YMSM_NewInfectionsTotal = sum(sum(AllOutcomes.ann_TotalNewInfections,2));
                %Young MSM percent diagnosed by age group & race/eth (final year)
                SelectedOutcomes.YMSM_PctAwareLastYr_18_24_MSM_B = AllOutcomes.YMSM_PctAwareLastYr_18_24_MSM_B;
                SelectedOutcomes.YMSM_PctAwareLastYr_18_24_MSM_H = AllOutcomes.YMSM_PctAwareLastYr_18_24_MSM_H;
                SelectedOutcomes.YMSM_PctAwareLastYr_18_24_MSM_O = AllOutcomes.YMSM_PctAwareLastYr_18_24_MSM_O; 
                SelectedOutcomes.YMSM_PctAwareLastYr_25_34_MSM_B = AllOutcomes.YMSM_PctAwareLastYr_25_34_MSM_B;
                SelectedOutcomes.YMSM_PctAwareLastYr_25_34_MSM_H = AllOutcomes.YMSM_PctAwareLastYr_25_34_MSM_H;
                SelectedOutcomes.YMSM_PctAwareLastYr_25_34_MSM_O = AllOutcomes.YMSM_PctAwareLastYr_25_34_MSM_O;
                SelectedOutcomes.YMSM_PctAwareLastYr_HET = AllOutcomes.YMSM_PctAwareLastYr_HET;
                SelectedOutcomes.YMSM_PctAwareLastYr_MSM = AllOutcomes.YMSM_PctAwareLastYr_MSM;
                SelectedOutcomes.YMSM_PctAwareLastYr_IDU = AllOutcomes.YMSM_PctAwareLastYr_IDU;
                %Young MSM percent VLS among diagnosed by age group & race/eth (final year)
                SelectedOutcomes.YMSM_VLSAmongDiagLastYr_18_24_MSM_B = AllOutcomes.YMSM_VLSAmongDiagLastYr_18_24_MSM_B;
                SelectedOutcomes.YMSM_VLSAmongDiagLastYr_18_24_MSM_H = AllOutcomes.YMSM_VLSAmongDiagLastYr_18_24_MSM_H;
                SelectedOutcomes.YMSM_VLSAmongDiagLastYr_18_24_MSM_O = AllOutcomes.YMSM_VLSAmongDiagLastYr_18_24_MSM_O; 
                SelectedOutcomes.YMSM_VLSAmongDiagLastYr_25_34_MSM_B = AllOutcomes.YMSM_VLSAmongDiagLastYr_25_34_MSM_B;
                SelectedOutcomes.YMSM_VLSAmongDiagLastYr_25_34_MSM_H = AllOutcomes.YMSM_VLSAmongDiagLastYr_25_34_MSM_H;
                SelectedOutcomes.YMSM_VLSAmongDiagLastYr_25_34_MSM_O = AllOutcomes.YMSM_VLSAmongDiagLastYr_25_34_MSM_O;                  
                SelectedOutcomes.YMSM_VLSAmongDiagLastYr_HET = AllOutcomes.YMSM_VLSAmongDiagLastYr_HET;
                SelectedOutcomes.YMSM_VLSAmongDiagLastYr_MSM = AllOutcomes.YMSM_VLSAmongDiagLastYr_MSM;
                SelectedOutcomes.YMSM_VLSAmongDiagLastYr_IDU = AllOutcomes.YMSM_VLSAmongDiagLastYr_IDU;
                
                OutcomeFields = fieldnames(SelectedOutcomes);

                nOutcomes = length(OutcomeFields);

                for OCount=1:nOutcomes
                    SelectedOutcomes.summary(OCount,1) = SelectedOutcomes.(OutcomeFields{OCount});
                end                                
                
                
                %Outcomes for calculating cost per person effectively diagnosed
                SelectedOutcomes.ann_CostTestAndNotify_LRMSM = AllOutcomes.ann_CostTestAndNotify_LRMSM_Disc;
                SelectedOutcomes.ann_CostTestAndNotify_HRMSM = AllOutcomes.ann_CostTestAndNotify_HRMSM_Disc;
                SelectedOutcomes.ann_TotalNewDiagnoses_LRMSM = AllOutcomes.ann_TotalNewDiagnoses_LRMSM;
                SelectedOutcomes.ann_TotalNewDiagnoses_HRMSM = AllOutcomes.ann_TotalNewDiagnoses_HRMSM;
        
          case 8 % YMSM resource allocation outcomes  
                
                % Create indicator to apply costs only to interventions
                % targeted by allocations (YMSM only)               
                AllocPopsTargetedIndicator = Params.YoungMSMIndicator;                                       
                
                Million = 1000000;
                
                MaxReach.Testing = Params.intn_Testing_MaxReach;
                MaxReach.LTCatDiag = Params.intn_LTCatDiag_MaxReach;
                MaxReach.LTCPostDiag = Params.intn_LTCafterDiag_MaxReach;
                MaxReach.ARTPrescription = Params.intn_ARTInitiation_MaxReach;
                MaxReach.ARTAdherance = Params.intn_TxAdherence_MaxReach;
                MaxReach.PrEP_Oral_HET = Params.intn_PrEP_Oral_HET_MaxReach;
                MaxReach.PrEP_Oral_MSM = Params.intn_PrEP_Oral_MSM_MaxReach;
                MaxReach.PrEP_Oral_IDU = Params.intn_PrEP_Oral_IDU_MaxReach;
                MaxReach.PrEP_Inject_HET = Params.intn_PrEP_Inject_HET_MaxReach;
                MaxReach.PrEP_Inject_MSM = Params.intn_PrEP_Inject_MSM_MaxReach;
                MaxReach.PrEP_Inject_IDU = Params.intn_PrEP_Inject_IDU_MaxReach;
                
                Allocation.alloc_Testing_MSM_Low_YMSM = AllOutcomes.alloc_Testing_MSM_Low;
                Allocation.alloc_Testing_MSM_High_YMSM = AllOutcomes.alloc_Testing_MSM_High;
                Allocation.alloc_LTCatDiag_YMSM = AllOutcomes.alloc_LTCatDiag;
                Allocation.alloc_LTCafterDiag_YMSM = AllOutcomes.alloc_LTCafterDiag;
                Allocation.alloc_ARTInitiation_YMSM = AllOutcomes.alloc_ARTInitiation;
                Allocation.alloc_ARTAdher5to4_YMSM = AllOutcomes.alloc_ARTAdher5to4;
                Allocation.alloc_ARTAdher4to5_YMSM = AllOutcomes.alloc_ARTAdher4to5;

                Allocation.alloc_PrEP_Oral_MSM_B_YMSM = AllOutcomes.alloc_PrEP_Oral_MSM_B;
                Allocation.alloc_PrEP_Oral_MSM_H_YMSM = AllOutcomes.alloc_PrEP_Oral_MSM_H;
                Allocation.alloc_PrEP_Oral_MSM_O_YMSM = AllOutcomes.alloc_PrEP_Oral_MSM_O;
                Allocation.alloc_PrEP_Inject_MSM_B_YMSM = AllOutcomes.alloc_PrEP_Inject_MSM_B;
                Allocation.alloc_PrEP_Inject_MSM_H_YMSM = AllOutcomes.alloc_PrEP_Inject_MSM_H;
                Allocation.alloc_PrEP_Inject_MSM_O_YMSM = AllOutcomes.alloc_PrEP_Inject_MSM_O;
                
                AnnOutcomes.ann_undiscNewInfections_YMSM = AllOutcomes.ann_undiscNewInfectionsYMSM;
                SumOutcomes.sum_undiscNewInfections_YMSM = sum(AllOutcomes.ann_undiscNewInfectionsYMSM);                
                AnnOutcomes.ann_undiscNewInfections_Overall = AllOutcomes.ann_undiscNewInfections;
                SumOutcomes.sum_undiscNewInfections_Overall = sum(AllOutcomes.ann_undiscNewInfections);                              
                SumOutcomes.sum_discQALYsYMSM_inM = sum(AllOutcomes.ann_discQALYs_YMSM)  / Million;                
                AnnOutcomes.ann_discQALYs_YMSM_inM = AllOutcomes.ann_discQALYs_YMSM / Million;                
                AnnOutcomes.ann_TotalARTAndCareCost_YMSM_inM = sum(AllOutcomes.ann_TotalARTAndCareCost_Disc,2)/Million;
                SumOutcomes.sum_TotalARTAndCareCost_YMSM_inM = sum(sum(AllOutcomes.ann_TotalARTAndCareCost_Disc,2))/Million;                
                AnnOutcomes.ann_TotalPrEPCost_YMSM_inclInAlloc_inM = sum(AllOutcomes.ann_HealthStateCost_PrEP_Disc,2)/Million;
                SumOutcomes.sum_TotalPrEPCost_YMSM_inclInAlloc_inM = sum(sum(AllOutcomes.ann_HealthStateCost_PrEP_Disc,2))/Million;
                AnnOutcomes.ann_TotalInterventionCosts_Overall_inM = AllOutcomes.ann_TotalTransCost_Undisc/Million;
                SumOutcomes.sum_TotalInterventionCosts_Overall_inM = sum(AllOutcomes.ann_TotalTransCost_Undisc)/Million;

                if FirstOutcomeYr <= rateYr && (LastOutcomeYr >= rateYr)
                                       
                    SumOutcomes.tr_LRYMSMAllEligTestProb = AllOutcomes.tr_LRYMSMAllEligTestProb_AllYrs(rateYrIndex);
                    SumOutcomes.tr_HRYMSMAllEligTestProb = AllOutcomes.tr_HRYMSMAllEligTestProb_AllYrs(rateYrIndex);                    
                    SumOutcomes.ann_PctLinkToCareFirst_YMSM = AllOutcomes.ann_PctLinkToCareFirst(rateYrIndex);
                    SumOutcomes.ann_LinkageAfterProb_YMSM = AllOutcomes.ann_LinkageAfterProb(rateYrIndex);
                    SumOutcomes.ann_ARTPrescrProb_YMSM = AllOutcomes.ann_ARTPrescrProb(rateYrIndex);
                    SumOutcomes.ann_VLSToANVProb_YMSM = AllOutcomes.ann_VLSToANVProb(rateYrIndex);
                    SumOutcomes.ann_ANVToVLSProb_YMSM = AllOutcomes.ann_ANVToVLSProb(rateYrIndex);                                     
                    SumOutcomes.ann_PrEPInitProb_HighRiskYMSM = AllOutcomes.ann_PrEPInitProb_HighRiskMSM(rateYrIndex);                   
                end
                
                SumOutcomes.ann_lastyrDiagnosed_YMSM = AllOutcomes.ann_PctAware(LastOutcomeYrIndex);                
                SumOutcomes.ann_lastyrDiagnosed_HRYMSM = AllOutcomes.ann_PctAware_HRMSM(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrDiagnosed_OYMSM = AllOutcomes.ann_PctAware_OMSM(LastOutcomeYrIndex);                
                SumOutcomes.ann_lastyrLinkedtoCare_YMSM = AllOutcomes.ann_PctInCare(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrOnART_YMSM = AllOutcomes.ann_PctOnART(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrVLS_YMSM = AllOutcomes.ann_PctVLS(LastOutcomeYrIndex);
                SumOutcomes.ann_lastyrVLSamongdiag_YMSM = AllOutcomes.ann_PctVLSamongdiag(LastOutcomeYrIndex);                
                SumOutcomes.ann_lastyrPctEligOnPrEP_HRYMSM = AllOutcomes.ann_PctEligOnPrEP_HRMSM(LastOutcomeYrIndex);               
                AnnOutcomes.ann_TotalARTAndCareCost_YMSM_Undisc = sum(AllOutcomes.ann_TotalARTAndCareCost_Undisc,2);
                SumOutcomeFields = fieldnames(SumOutcomes);

                SelectedOutcomes = SumOutcomes;
                
                for nAllocPeriod = 1:Params.numAllocationPeriods
                    SelectedOutcomes.TotalIntnBudget_inM(nAllocPeriod,1) = Params.alloc_TotalIntnBudget(nAllocPeriod);
                    SelectedOutcomes.summary(nAllocPeriod,1) = Params.alloc_TotalIntnBudget(nAllocPeriod);
                end
                SummaryRows = Params.numAllocationPeriods;
                
                SelectedOutcomes.summary(SummaryRows+1,1) = ...
                    Params.alloc_PctARTAndCareCoveredByBudget;
                SummaryRows = SummaryRows + 1;
                
                MaxReachFields = fieldnames(MaxReach);
                nMaxReachFields = length(MaxReachFields);
                for nMRCount=1:nMaxReachFields
                    SelectedOutcomes.MaxReach(nMRCount,1) = MaxReach.(MaxReachFields{nMRCount});
                    SelectedOutcomes.summary(SummaryRows+nMRCount,1) = MaxReach.(MaxReachFields{nMRCount});
                end
                SummaryRows = SummaryRows + nMaxReachFields;
                
                SelectedOutcomes.periodStartYrs = Params.Intn_Allocation_StartYr(1:Params.numAllocationPeriods)';
                SelectedOutcomes.summary ...
                    ((SummaryRows+1):(SummaryRows+Params.numAllocationPeriods),1)= ...
                    SelectedOutcomes.periodStartYrs;
                SummaryRows = SummaryRows + Params.numAllocationPeriods;
                
                AllocFields = fieldnames(Allocation);
                nAllocFields = length(AllocFields);
                for nAFCount=1:nAllocFields
                    for nAllocPeriod = 1:Params.numAllocationPeriods
                        SelectedOutcomes.Allocation_inM(nAllocFields*(nAllocPeriod-1)+nAFCount,1) = Allocation.(AllocFields{nAFCount})(nAllocPeriod);
                        SelectedOutcomes.summary(SummaryRows+nAllocFields*(nAllocPeriod-1)+nAFCount,1) = Allocation.(AllocFields{nAFCount})(nAllocPeriod);
                    end
                end
                SummaryRows = SummaryRows + nAllocFields*Params.numAllocationPeriods;
                
                for SOFCount=1:length(SumOutcomeFields)
                    SelectedOutcomes.summary(SummaryRows + SOFCount,1) = ...
                        SumOutcomes.(SumOutcomeFields{SOFCount});
                end
                SummaryRows = SummaryRows + length(SumOutcomeFields);
                
                SelectedOutcomes.summary(SummaryRows+1:SummaryRows + length(AnnOutcomes.ann_TotalARTAndCareCost_YMSM_inM)) = ...
                    AnnOutcomes.ann_TotalARTAndCareCost_YMSM_inM;
                
                SelectedOutcomes.AnnualOutcomes = AnnOutcomes;        

            case 9 %COVID analysis outcomes- Added 1/10/2023 Bates
                                      
                SelectedOutcomes.ann_PctAware = sum(AllOutcomes.ann_PctAware,2);
                SelectedOutcomes.ann_PctVLSamongdiag = sum(AllOutcomes.ann_PctVLSamongdiag,2);
                SelectedOutcomes.ann_NumberOnART = sum(AllOutcomes.ann_NumberOnART,2);
                SelectedOutcomes.ann_NumberVLS = sum(AllOutcomes.ann_NumberVLS,2);
                SelectedOutcomes.ann_NumberOnPrEP = sum(AllOutcomes.ann_NumberOnPrEP,2);
                SelectedOutcomes.ann_numTotalTests = sum(AllOutcomes.ann_numTotalTests,2);
                SelectedOutcomes.tr_WeightedAllEligTestRate_AllYrs = AllOutcomes.tr_WeightedAllEligTestRate_AllYrs; 
                SelectedOutcomes.tr_DiagRateAll_AllYrs = AllOutcomes.tr_DiagRateAll_AllYrs; 
                SelectedOutcomes.ann_TotalNewDiagnoses = sum(AllOutcomes.ann_TotalNewDiagnoses,2);                
                SelectedOutcomes.ann_TotalNewInfections = sum(AllOutcomes.ann_TotalNewInfections,2);

                                                                                           
                OutcomeFields = fieldnames(SelectedOutcomes);

                nOutcomes = length(OutcomeFields);
                                

            case 10 %Injectable PrEP impact - Added 12/1/2021 Clinkscales
                
                Million = 1000000;
                                                
                SelectedOutcomes.ann_NumberOnPrEP_Oral = sum(AllOutcomes.ann_NumberOnPrEP_Oral,2);
                SelectedOutcomes.ann_NumberOnPrEP_Inject = sum(AllOutcomes.ann_NumberOnPrEP_Inject,2);
                SelectedOutcomes.ann_NumberOnPrEP = sum(AllOutcomes.ann_NumberOnPrEP,2);
                SelectedOutcomes.ann_NumberOnPrEP_HighRiskHETs_M = sum(AllOutcomes.ann_NumberOnPrEP_HighRiskHETs_M,2);
                SelectedOutcomes.ann_NumberOnPrEP_HighRiskHETs_F = sum(AllOutcomes.ann_NumberOnPrEP_HighRiskHETs_F,2);
                SelectedOutcomes.ann_NumberOnPrEP_HighRiskMSM = sum(AllOutcomes.ann_NumberOnPrEP_HighRiskMSM,2);
                SelectedOutcomes.ann_NumberOnPrEP_HighRiskIDUs = sum(AllOutcomes.ann_NumberOnPrEP_HighRiskIDUs,2);  
                SelectedOutcomes.ann_PersonYears_onPrEP_Oral = sum(AllOutcomes.ann_PersonYears_onPrEP_Oral,2);
                SelectedOutcomes.ann_PersonYears_onPrEP_Inject = sum(AllOutcomes.ann_PersonYears_onPrEP_Inject,2);
                SelectedOutcomes.ann_PersonYears_onPrEP = sum(AllOutcomes.ann_PersonYears_onPrEP,2);
                SelectedOutcomes.ann_PersonYears_onPrEP_HETM = sum(AllOutcomes.ann_PersonYears_onPrEP_HETM,2);
                SelectedOutcomes.ann_PersonYears_onPrEP_HETF = sum(AllOutcomes.ann_PersonYears_onPrEP_HETF,2);
                SelectedOutcomes.ann_PersonYears_onPrEP_MSM = sum(AllOutcomes.ann_PersonYears_onPrEP_MSM,2);
                SelectedOutcomes.ann_PersonYears_onPrEP_IDU = sum(AllOutcomes.ann_PersonYears_onPrEP_IDU,2);
                SelectedOutcomes.ann_StartPrEP_Oral = sum(AllOutcomes.ann_StartPrEP_Oral,2);
                SelectedOutcomes.ann_StartPrEP_Inject = sum(AllOutcomes.ann_StartPrEP_Inject,2);
                SelectedOutcomes.ann_StartPrEP = sum(AllOutcomes.ann_StartPrEP,2);                
                SelectedOutcomes.ann_DropOff_PrEP_Oral = sum(AllOutcomes.ann_DropOff_PrEP_Oral,2);
                SelectedOutcomes.ann_DropOff_PrEP_Inject = sum(AllOutcomes.ann_DropOff_PrEP_Inject,2);
                SelectedOutcomes.ann_DropOff_PrEP = sum(AllOutcomes.ann_DropOff_PrEP,2);
                SelectedOutcomes.ann_NewInfectionsHETM = sum(AllOutcomes.ann_NewInfectionsHETM,2);
                SelectedOutcomes.ann_NewInfectionsHETF = sum(AllOutcomes.ann_NewInfectionsHETF,2);
                SelectedOutcomes.ann_NewInfectionsIDUM = sum(AllOutcomes.ann_NewInfectionsIDUM,2);
                SelectedOutcomes.ann_NewInfectionsIDUF = sum(AllOutcomes.ann_NewInfectionsIDUF,2);
                SelectedOutcomes.ann_NewInfectionsMSM = sum(AllOutcomes.ann_NewInfectionsMSM,2);
                SelectedOutcomes.ann_TotalNewInfections = sum(AllOutcomes.ann_TotalNewInfections,2);
                SelectedOutcomes.ann_HealthStateCost_PrEP_Oral_inM = sum(AllOutcomes.ann_healthStateCost_PrEP_Oral_Undisc,2) / Million; %discounted version: ann_HealthStateCost_PrEP_Oral_Disc
                SelectedOutcomes.ann_HealthStateCost_PrEP_Inject_inM = sum(AllOutcomes.ann_healthStateCost_PrEP_Inject_Undisc,2) / Million; %discounted version: ann_HealthStateCost_PrEP_Inject_Disc
                SelectedOutcomes.ann_HealthStateCost_PrEP_inM = sum(AllOutcomes.ann_healthStateCost_PrEP_Undisc,2) / Million; %discounted version: ann_HealthStateCost_PrEP
                SelectedOutcomes.ann_TotalARTAndCareCost_inM = sum(AllOutcomes.ann_TotalARTAndCareCost_Undisc,2) / Million; %discounted version: ann_TotalARTAndCareCost
                
                SelectedOutcomes.ann_PctEligOnPrEP_HRHET = AllOutcomes.ann_PctEligOnPrEP_HRHET;
                SelectedOutcomes.ann_PctEligOnPrEP_HRMSM = AllOutcomes.ann_PctEligOnPrEP_HRMSM;
                SelectedOutcomes.ann_PctEligOnPrEP_PWID = AllOutcomes.ann_PctEligOnPrEP_PWID;                           
                                                                                           
                OutcomeFields = fieldnames(SelectedOutcomes);

                nOutcomes = length(OutcomeFields);

                rowcount = 1;
                
                for OCount=1:nOutcomes
                    if contains(OutcomeFields(OCount),{'PctEligOnPrEP', 'NumberOnPrEP'}) == 0
                        SelectedOutcomes.CumulativeSummary(rowcount,1) = sum(SelectedOutcomes.(OutcomeFields{OCount}));
                        
                        rowcount = rowcount + 1;
                    end
                end                                
                
                rowcount = 1;
                
                for OCount=1:nOutcomes
                    if contains(OutcomeFields(OCount),'PersonYears_onPrEP') == 0
                        SelectedOutcomes.AnnualSummary(length(SelectedOutcomes.ann_NumberOnPrEP_Oral) *(rowcount - 1) + 1: length(SelectedOutcomes.ann_NumberOnPrEP_Oral) * rowcount,1) = SelectedOutcomes.(OutcomeFields{OCount});
                        
                        rowcount = rowcount + 1;
                    end
                end
             
                
%                 for PCount=1:nOutcomes
%                     for i=1:10
%                         SelectedOutcomes.Summary(OCount,i) = SelectedOutcomes.(OutcomeFields{OCount,i});
%                     end
%                 end


            case 11 %Disparities analysis outcomes
                SelectedOutcomes.ann_PctAware_Overall = sum(AllOutcomes.ann_PctAware,2);
                SelectedOutcomes.ann_PctAware_Blk = sum(AllOutcomes.ann_PctAware_Blk,2);
                SelectedOutcomes.ann_PctAware_Hisp = sum(AllOutcomes.ann_PctAware_Hisp,2);
                SelectedOutcomes.ann_PctAware_Oth = sum(AllOutcomes.ann_PctAware_Oth,2);
                SelectedOutcomes.ann_PctInCareAmongDiag_Overall = sum(AllOutcomes.ann_NumberInCare,2) ./ sum(AllOutcomes.ann_NumberAware,2);
                SelectedOutcomes.ann_PctInCareAmongDiag_Blk = sum(AllOutcomes.ann_NumberInCare_Blk,2) ./ sum(AllOutcomes.ann_NumberAware_Blk,2);
                SelectedOutcomes.ann_PctInCareAmongDiag_Hisp = sum(AllOutcomes.ann_NumberInCare_Hisp,2) ./ sum(AllOutcomes.ann_NumberAware_Hisp,2);
                SelectedOutcomes.ann_PctInCareAmongDiag_Oth = sum(AllOutcomes.ann_NumberInCare_Oth,2) ./ sum(AllOutcomes.ann_NumberAware_Oth,2);
                SelectedOutcomes.ann_PctVLSamongdiag_Overall = AllOutcomes.ann_PctVLSamongdiag;
                SelectedOutcomes.ann_PctVLSamongdiag_Blk = sum(AllOutcomes.ann_PctVLSamongdiag_Blk,2);
                SelectedOutcomes.ann_PctVLSamongdiag_Hisp = sum(AllOutcomes.ann_PctVLSamongdiag_Hisp,2);
                SelectedOutcomes.ann_PctVLSamongdiag_Oth = sum(AllOutcomes.ann_PctVLSamongdiag_Oth,2);
                SelectedOutcomes.ann_PctUninfPWIDServedbySEP_Overall = sum(AllOutcomes.ann_PctUninfPWIDServedbySEP,2);
                SelectedOutcomes.ann_PctUninfPWIDServedbySEP_Blk = sum(AllOutcomes.ann_PctUninfPWIDServedbySEP_Blk,2);
                SelectedOutcomes.ann_PctUninfPWIDServedbySEP_Hisp = sum(AllOutcomes.ann_PctUninfPWIDServedbySEP_Hisp,2);
                SelectedOutcomes.ann_PctUninfPWIDServedbySEP_Oth = sum(AllOutcomes.ann_PctUninfPWIDServedbySEP_Oth,2);
                SelectedOutcomes.ann_PctSusOnPrEP_OverallHR = sum(AllOutcomes.ann_NumberOnPrEP,2) ./ sum(AllOutcomes.ann_NumberEligForPrEPNotOnPrEP,2);
                SelectedOutcomes.ann_PctSusOnPrEP_HRBlk = sum(AllOutcomes.ann_PctSusOnPrEP_HighRiskBlk,2);
                SelectedOutcomes.ann_PctSusOnPrEP_HRHisp = sum(AllOutcomes.ann_PctSusOnPrEP_HighRiskHisp,2);
                SelectedOutcomes.ann_PctSusOnPrEP_HROth = sum(AllOutcomes.ann_PctSusOnPrEP_HighRiskOth,2);
                SelectedOutcomes.ann_TotalNewInfections_Overall = sum(AllOutcomes.ann_TotalNewInfections,2);
                SelectedOutcomes.ann_NewInfections_Blk = sum(AllOutcomes.ann_NewInfections_Blk,2);
                SelectedOutcomes.ann_NewInfections_Hisp = sum(AllOutcomes.ann_NewInfections_Hisp,2);
                SelectedOutcomes.ann_NewInfections_Oth = sum(AllOutcomes.ann_NewInfections_Oth,2);
                SelectedOutcomes.ann_TotalNewInfections_PerPerson_Overall = sum(AllOutcomes.ann_TotalNewInfections,2) ./ sum(AllOutcomes.ann_PopulationSize,2);
                SelectedOutcomes.ann_TotalNewInfections_PerPerson_Blk = sum(AllOutcomes.ann_NewInfections_Blk,2) ./ sum(AllOutcomes.ann_popSize_Blk,2);
                SelectedOutcomes.ann_TotalNewInfections_PerPerson_Hisp = sum(AllOutcomes.ann_NewInfections_Hisp,2) ./ sum(AllOutcomes.ann_popSize_Hisp,2);
                SelectedOutcomes.ann_TotalNewInfections_PerPerson_Oth = sum(AllOutcomes.ann_NewInfections_Oth,2) ./ sum(AllOutcomes.ann_popSize_Oth,2);
                SelectedOutcomes.ann_TotalNewInfections_PerPLWH_Overall = sum(AllOutcomes.ann_TotalNewInfections,2) ./ sum(AllOutcomes.ann_HIVPrevalence,2);
                SelectedOutcomes.ann_TotalNewInfections_PerPLWH_Blk = sum(AllOutcomes.ann_NewInfections_Blk,2) ./ sum(AllOutcomes.ann_HIVPrevalence_Blk,2);
                SelectedOutcomes.ann_TotalNewInfections_PerPLWH_Hisp = sum(AllOutcomes.ann_NewInfections_Hisp,2) ./ sum(AllOutcomes.ann_HIVPrevalence_Hisp,2);
                SelectedOutcomes.ann_TotalNewInfections_PerPLWH_Oth = sum(AllOutcomes.ann_NewInfections_Oth,2) ./ sum(AllOutcomes.ann_HIVPrevalence_Oth,2);
                SelectedOutcomes.ann_ARTCarePrEPTransitionAndSEPCost_Disc_Overall = AllOutcomes.ann_ARTCarePrEPTransitionAndSEPCost_Disc;
                SelectedOutcomes.ann_ARTCarePrEPTransitionAndSEPCost_Disc_Blk = AllOutcomes.ann_ARTCarePrEPTransitionAndSEPCost_Disc_Blk;
                SelectedOutcomes.ann_ARTCarePrEPTransitionAndSEPCost_Disc_Hisp = AllOutcomes.ann_ARTCarePrEPTransitionAndSEPCost_Disc_Hisp;
                SelectedOutcomes.ann_ARTCarePrEPTransitionAndSEPCost_Disc_Oth = AllOutcomes.ann_ARTCarePrEPTransitionAndSEPCost_Disc_Oth;

                OutcomeFields = fieldnames(SelectedOutcomes);

                nOutcomes = length(OutcomeFields);

                rowcount = 1;
                
                for OCount=1:nOutcomes                   
                    SelectedOutcomes.AnnualSummary(1: length(SelectedOutcomes.ann_PctAware_Overall),OCount) = SelectedOutcomes.(OutcomeFields{OCount});

                    rowcount = rowcount + 1;                   
                end
                               
        end
    end

end