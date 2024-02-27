function [Compartments, Results, TTProg] = StepStateCont(Params, TransRates, Compartments, Year, k, Results, TTProg)
%% Purpose: For the continuous version of the model: a function which 
% advances the state of the system by running the ODE and applying the 
% rates to Compartments.
% Called from: HIVEpiModel

%% 1. Declare and Initialize Variables
   
% Declare global variables that are needed from the ODE
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
    
%Initialize/clear vars
    g_TransRatesAnnual = [];
    g_RelativeInfRates = [];
    g_lambdaVonly = [];
    g_lambdaAonly = []; 
    g_lambdaVsomeAI = [];
    g_lambdaAsomeVI = [];
    g_lambdaN = [];
    g_odeStepCounter = 1;
    
    g_RelativeInfRates_PrEP_Oral_High = [];
    g_lambdaVonly_PrEP_Oral_High = [];
    g_lambdaAonly_PrEP_Oral_High = [];
    g_lambdaVsomeAI_PrEP_Oral_High = [];
    g_lambdaAsomeVI_PrEP_Oral_High = [];
    g_lambdaN_PrEP_Oral_High = [];
    
    g_RelativeInfRates_PrEP_Oral_Low = [];
    g_lambdaVonly_PrEP_Oral_Low = [];
    g_lambdaAonly_PrEP_Oral_Low = [];
    g_lambdaVsomeAI_PrEP_Oral_Low = [];
    g_lambdaAsomeVI_PrEP_Oral_Low = [];
    g_lambdaN_PrEP_Oral_Low = [];
    
    g_RelativeInfRates_PrEP_Inject_High = [];
    g_lambdaVonly_PrEP_Inject_High = [];
    g_lambdaAonly_PrEP_Inject_High = [];
    g_lambdaVsomeAI_PrEP_Inject_High = [];
    g_lambdaAsomeVI_PrEP_Inject_High = [];
    g_lambdaN_PrEP_Inject_High = [];
    
    g_RelativeInfRates_PrEP_Inject_Low = [];
    g_lambdaVonly_PrEP_Inject_Low = [];
    g_lambdaAonly_PrEP_Inject_Low = [];
    g_lambdaVsomeAI_PrEP_Inject_Low = [];
    g_lambdaAsomeVI_PrEP_Inject_Low = [];
    g_lambdaN_PrEP_Inject_Low = [];

%% 2. Prepare Settings for Differential Equation (ODE) Solver
    
    % Refine value is how granular you would like the ODE output
        % Minimum is 1 (lower is more granular)
    RefineValue = 1;
    
    % Number of times the rates function evaluates per ODE micro-step
        % Number is based on the ODE-solver and refine values selected.
    nFevalsPerTStep = 6; % Needs to be 6 if refine val = 1
    
    % Time span over which ode is solved
        % Set to equal 1 year
    tspan = [0,1];
    
    % Prepare rates function for ODE
    rateswrap = @(T,Y) RatesODE(Y, Params, TransRates, TTProg,nFevalsPerTStep,Year);
    
    % Prepare options for ODE
    options = odeset('Nonnegative', 1:Params.numComparts*Params.numStrats,'Refine',RefineValue);

%% 3. Run the ODE
    
    % Initial condition
    y0 = reshape(Compartments,Params.numComparts*Params.numStrats,1);    

    % Call ODE
    [T,Y] = ode45(rateswrap, tspan, y0, options);
    
%% 4. Reformat Output From ODE

    odeCounter = size(T,1);
    CompartmentsByODEMicroStep = reshape(Y', Params.numComparts,Params.numStrats, odeCounter);
        
    % Assign global inf rates from ODE to struct
    InfRate.RelativeInfRates = g_RelativeInfRates;
    InfRate.lambdaVonly = g_lambdaVonly;
    InfRate.lambdaAonly = g_lambdaAonly;
    InfRate.lambdaVsomeAI = g_lambdaVsomeAI;
    InfRate.lambdaAsomeVI = g_lambdaAsomeVI;
    InfRate.lambdaN = g_lambdaN;
    
    InfRate_PrEP_Oral_High.RelativeInfRates = g_RelativeInfRates_PrEP_Oral;
    InfRate_PrEP_Oral_High.lambdaVonly = g_lambdaVonly_PrEP_Oral;
    InfRate_PrEP_Oral_High.lambdaAonly = g_lambdaAonly_PrEP_Oral;
    InfRate_PrEP_Oral_High.lambdaVsomeAI = g_lambdaVsomeAI_PrEP_Oral;
    InfRate_PrEP_Oral_High.lambdaAsomeVI = g_lambdaAsomeVI_PrEP_Oral;
    InfRate_PrEP_Oral_High.lambdaN = g_lambdaN_PrEP_Oral;
    
    InfRate_PrEP_Oral_Low.RelativeInfRates = g_RelativeInfRates_PrEP_Oral;
    InfRate_PrEP_Oral_Low.lambdaVonly = g_lambdaVonly_PrEP_Oral;
    InfRate_PrEP_Oral_Low.lambdaAonly = g_lambdaAonly_PrEP_Oral;
    InfRate_PrEP_Oral_Low.lambdaVsomeAI = g_lambdaVsomeAI_PrEP_Oral;
    InfRate_PrEP_Oral_Low.lambdaAsomeVI = g_lambdaAsomeVI_PrEP_Oral;
    InfRate_PrEP_Oral_Low.lambdaN = g_lambdaN_PrEP_Oral;
    
    InfRate_PrEP_Inject_High.RelativeInfRates = g_RelativeInfRates_PrEP_Inject;
    InfRate_PrEP_Inject_High.lambdaVonly = g_lambdaVonly_PrEP_Inject;
    InfRate_PrEP_Inject_High.lambdaAonly = g_lambdaAonly_PrEP_Inject;
    InfRate_PrEP_Inject_High.lambdaVsomeAI = g_lambdaVsomeAI_PrEP_Inject;
    InfRate_PrEP_Inject_High.lambdaAsomeVI = g_lambdaAsomeVI_PrEP_Inject;
    InfRate_PrEP_Inject_High.lambdaN = g_lambdaN_PrEP_Inject;
    
    InfRate_PrEP_Inject_Low.RelativeInfRates = g_RelativeInfRates_PrEP_Inject;
    InfRate_PrEP_Inject_Low.lambdaVonly = g_lambdaVonly_PrEP_Inject;
    InfRate_PrEP_Inject_Low.lambdaAonly = g_lambdaAonly_PrEP_Inject;
    InfRate_PrEP_Inject_Low.lambdaVsomeAI = g_lambdaVsomeAI_PrEP_Inject;
    InfRate_PrEP_Inject_Low.lambdaAsomeVI = g_lambdaAsomeVI_PrEP_Inject;
    InfRate_PrEP_Inject_Low.lambdaN = g_lambdaN_PrEP_Inject;

    % Assign ODE rates to TransRates
    TransRates.AnnualRates = g_TransRatesAnnual;

    % Determine time per ODE step
    timePerODEStep = T(2:odeCounter)- T(1:odeCounter-1); 
    timePerODEStep = [timePerODEStep ; 1-sum(timePerODEStep,1)];


% Test to see if this is first run -> if not initialize
    tst = isfield(TTProg,'stepsPerYear');
    
    if tst == 0;
        TTProg.stepsPerYear = odeCounter;
        TTProg.k  = 0;
        TTProg.timePerODEStep = timePerODEStep;
    else
        TTProg.stepsPerYear = [TTProg.stepsPerYear, odeCounter];
        TTProg.timePerODEStep = [TTProg.timePerODEStep; timePerODEStep];
    end

%% 5. Collect Results

    for j = 1:odeCounter

        % Calculate the number of people who transition between each of the compartments
            % Based on the rates function
            
            for i = 1:Params.numStrats
                
                % This is a 31 x 31 matrix of the population
                C = repmat(CompartmentsByODEMicroStep(:,i,j),[1,Params.numComparts]);
                
                Transitions(:,:,i,j) = C .* (TransRates.AnnualRates(:,:,i,j) * timePerODEStep(j,1));
                
            end

        % Update step counter
        TTProg.k = 1+TTProg.k;
        

        % Calculate number of people tested during the tstep
            % Compartments at the current ODE timestep
            CompartmentsTS = CompartmentsByODEMicroStep(:,:,j);
            
            
            % Calc number of tests and test/diag rates
            TTProg = CalcNumIntnAffected(CompartmentsTS,Params,TTProg,Year,j);

        % Record Results
        Results = CollectResults( Params, CompartmentsByODEMicroStep, Year, InfRate, InfRate_PrEP_Oral_High, InfRate_PrEP_Oral_Low, InfRate_PrEP_Inject_High, InfRate_PrEP_Inject_Low, j, Transitions, Results, TTProg);
        
    end

%% 6. Update Variables
    
    % Update Compartment variable to the latest distribution of the
    % population
    Compartments = CompartmentsByODEMicroStep(:,:,odeCounter);
    
end