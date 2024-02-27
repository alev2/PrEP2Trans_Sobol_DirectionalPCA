function [x0,A,b,Aeq,beq,lb,ub,options]=GenOptComponents(Params)
%% Set up supporting variables
    
    % Set variables based on if run is targeted to YMSM only or not
    if Params.TargetIntnsToYMSM == 1
        % Establish number of DVs per allocation period
        nNumDVsPerAlloc = Params.numIntns_YMSM+Params.alloc_ConsiderARTAndCareCosts;
        % Set initial allocation variable
        initAlloc = Params.initalloc(Params.YMSMint,:);               
    else
        % Establish number of DVs per allocation period
        nNumDVsPerAlloc = Params.numIntns+Params.alloc_ConsiderARTAndCareCosts;
        % Set initial allocation variable
        initAlloc = Params.initalloc;
    end
    
    % Set up indicators for DVs by allocation period
    if Params.numAllocationPeriods == 1
        indAllocPeriod = [ones(1,nNumDVsPerAlloc)];
    elseif Params.numAllocationPeriods == 2
        indAllocPeriod = [ones(1,nNumDVsPerAlloc) zeros(1,nNumDVsPerAlloc);
             zeros(1,nNumDVsPerAlloc) ones(1,nNumDVsPerAlloc)];
    else % Params.numAllocationPeriods == 3;
        indAllocPeriod = [ones(1,nNumDVsPerAlloc) zeros(1,2*(nNumDVsPerAlloc));
             zeros(1,nNumDVsPerAlloc) ones(1,nNumDVsPerAlloc) zeros(1,nNumDVsPerAlloc);
             zeros(1,2*(nNumDVsPerAlloc)) ones(1,nNumDVsPerAlloc)];
    end

% INDICATOR: Intns by type
%     if Params.TargetIntnsToYMSM == 1
%         
%         ind_IntnsByIntnType = ...
%             [ones(1,Params.numIntns_Testing_YMSM) zeros(1,nNumDVsPerAlloc-Params.numIntns_Testing_YMSM);
%             zeros(1,Params.numIntns_Testing_YMSM) ones(1,Params.numIntns_LTCatDiag) zeros(1,nNumDVsPerAlloc-Params.numIntns_Testing_YMSM-Params.numIntns_LTCatDiag);
%             zeros(1,Params.numIntns_Testing_YMSM+Params.numIntns_LTCatDiag) ones(1,Params.numIntns_LTCafterDiag) zeros(1,nNumDVsPerAlloc-Params.numIntns_Testing_YMSM-Params.numIntns_LTCatDiag-Params.numIntns_LTCafterDiag);
%             zeros(1,Params.numIntns_Testing_YMSM+Params.numIntns_LTCatDiag+Params.numIntns_LTCafterDiag) ones(1,Params.numIntns_ARTPrescription) zeros(1,nNumDVsPerAlloc-Params.numIntns_Testing_YMSM-Params.numIntns_LTCatDiag-Params.numIntns_LTCafterDiag-Params.numIntns_ARTPrescription);
%             zeros(1,Params.numIntns_Testing_YMSM+Params.numIntns_LTCatDiag+Params.numIntns_LTCafterDiag+Params.numIntns_ARTPrescription) ones(1,Params.numIntns_ARTAdherence) zeros(1,nNumDVsPerAlloc-Params.numIntns_Testing_YMSM-Params.numIntns_LTCatDiag-Params.numIntns_LTCafterDiag-Params.numIntns_ARTPrescription-Params.numIntns_ARTAdherence);           
%             zeros(1,Params.numIntns_YMSM-Params.numIntns_PrEP_YMSM) ones(1,Params.numIntns_PrEP_YMSM) zeros(1,nNumDVsPerAlloc - Params.numIntns_YMSM)];
%         if Params.alloc_ConsiderARTAndCareCosts == 1
%             ind_IntnsByIntnType = [ind_IntnsByIntnType zeros(1)];
%         end
%         
%     else
%         
%         ind_IntnsByIntnType = ...
%             [ones(1,Params.numIntns_Testing) zeros(1,nNumDVsPerAlloc-Params.numIntns_Testing);
%             zeros(1,Params.numIntns_Testing) ones(1,Params.numIntns_LTCatDiag) zeros(1,nNumDVsPerAlloc-Params.numIntns_Testing-Params.numIntns_LTCatDiag);
%             zeros(1,Params.numIntns_Testing+Params.numIntns_LTCatDiag) ones(1,Params.numIntns_LTCafterDiag) zeros(1,nNumDVsPerAlloc-Params.numIntns_Testing-Params.numIntns_LTCatDiag-Params.numIntns_LTCafterDiag);
%             zeros(1,Params.numIntns_Testing+Params.numIntns_LTCatDiag+Params.numIntns_LTCafterDiag) ones(1,Params.numIntns_ARTPrescription) zeros(1,nNumDVsPerAlloc-Params.numIntns_Testing-Params.numIntns_LTCatDiag-Params.numIntns_LTCafterDiag-Params.numIntns_ARTPrescription);
%             zeros(1,Params.numIntns_Testing+Params.numIntns_LTCatDiag+Params.numIntns_LTCafterDiag+Params.numIntns_ARTPrescription) ones(1,Params.numIntns_ARTAdherence) zeros(1,nNumDVsPerAlloc-Params.numIntns_Testing-Params.numIntns_LTCatDiag-Params.numIntns_LTCafterDiag-Params.numIntns_ARTPrescription-Params.numIntns_ARTAdherence);
%             zeros(1,Params.numIntns-Params.numIntns_PrEP-Params.numIntns_SEP) ones(1,Params.numIntns_SEP) zeros(1,nNumDVsPerAlloc-Params.numIntns_Testing-Params.numIntns_LTCatDiag-Params.numIntns_LTCafterDiag-Params.numIntns_ARTPrescription-Params.numIntns_ARTAdherence-Params.numIntns_SEP);
%             zeros(1,Params.numIntns-Params.numIntns_PrEP) ones(1,Params.numIntns_PrEP) zeros(1,nNumDVsPerAlloc - Params.numIntns)];
%         if Params.alloc_ConsiderARTAndCareCosts == 1
%             ind_IntnsByIntnType = [ind_IntnsByIntnType zeros(1)];
%         end
%     
%     end

%% x0: initial solution; can be a scalar, vector, or matrix.

    if and(Params.optObjNum ~= 4,Params.alloc_ConsiderARTAndCareCosts == 1)
        x0 = [initAlloc(:,1); Params.alloc_TotalARTAndCareCosts];
        if Params.numAllocationPeriods == 2
            x0 = [x0; initAlloc(:,2); alloc_TotalARTAndCareCosts];
        elseif Params.numAllocationPeriods == 3
            x0 = [x0; initAlloc(:,2); alloc_TotalARTAndCareCosts;...
                initAlloc(:,3); alloc_TotalARTAndCareCosts];
        end
    else
        x0 = [initAlloc(:,1)];
        if Params.numAllocationPeriods == 2
            x0 = [x0; initAlloc(:,2)];
        elseif Params.numAllocationPeriods == 3
            x0 = [x0; initAlloc(:,2);initAlloc(:,3)];
        end
    end

%% A / b: defines linear inequalities A*x <= b. If no inequalities exist, set A = [] and b = []. If A large w/ few nonzero entries, make sparse.

    % Budget constraint: Allocate only as much funding as is available
    %       Sum of allocations (DVs) over all grantees, B & C <= Budget
    ScaleFactor = 1;
    if or(Params.optObjNum == 4, Params.alloc_ConsiderARTAndCareCosts == 1) 
        %Don't need linear budget constraints if minimizing cost to hit
        % hit incidence targets (no budget constraint) or
        % including treatment costs in budget (in which case non-linear
        %budget contraints are used)
        A = [];
        b = [];
    else
        % A starts with one row per allocation period with indicators for
        % DVs by allocation period
        A = indAllocPeriod*ScaleFactor;
        b = Params.alloc_TotalIntnBudget(1)*ScaleFactor;
        if Params.numAllocationPeriods == 2 
            b = [b; 
                Params.alloc_TotalIntnBudget(2)*ScaleFactor];
        elseif Params.numAllocationPeriods == 3
            b = [b;
                Params.alloc_TotalIntnBudget(2)*ScaleFactor;
                Params.alloc_TotalIntnBudget(3)*ScaleFactor];
        end
    end

        
%% Aeq / beq: defines linear equalities Aeq*x = beq. If no equalities exist, set Aeq = [] and beq = []
    Aeq=[];
    beq=[];
        
%% lb / ub: lower and upper bounds on the design variables in x, so that 
% the solution is always in the range lb <= x <= ub. 
% If no bounds exist, set lb = [] and/or ub = [].
% If x(i) is unbounded below, set lb(i) = -Inf, and if x(i) is unbounded above, set ub(i) = Inf.

    % Minimum and maximum annual allocation of budget to each intervention
    % for each allocation period
    if Params.TargetIntnsToYMSM == 1
        
        lb = zeros(Params.numIntns_YMSM*Params.numAllocationPeriods,1);
        ub = zeros(Params.numIntns_YMSM*Params.numAllocationPeriods,1);
        for period = 1:Params.numAllocationPeriods
           
            lb((period-1)*nNumDVsPerAlloc+1:period*nNumDVsPerAlloc,1) = ...
                Params.minalloc(Params.YMSMint,period);
            ub((period-1)*nNumDVsPerAlloc+1:period*nNumDVsPerAlloc,1) = ...
                Params.maxalloc(Params.YMSMint,period);
        end
        
    else
        
        lb = zeros(Params.numIntns*Params.numAllocationPeriods,1);
        ub = zeros(Params.numIntns*Params.numAllocationPeriods,1);
        for period = 1:Params.numAllocationPeriods
            lb((period-1)*nNumDVsPerAlloc+1:period*nNumDVsPerAlloc,1) = ...
                Params.minalloc(:,period);
            ub((period-1)*nNumDVsPerAlloc+1:period*nNumDVsPerAlloc,1) = ...
                Params.maxalloc(:,period);
        end
    
    end
    


%%	nonlcon: defines nonlinear inequalities c(x) or equalities ceq(x). 
%       fmincon optimizes such that c(x) <= 0 and ceq(x) = 0. 
%       accepts a vector x and returns the two vectors c and ceq. 
%       c is a vector that contains the nonlinear inequalities evaluated at x, 
%       and ceq is a vector that contains the nonlinear equalities evaluated at x. 
%       nonlcon should be specified as a function handle to a file or to 
%       an anonymous function, such as 
%       mycon: x = fmincon(@myfun,x0,A,b,Aeq,beq,lb,ub,@mycon), 
%       where mycon is a MATLAB function such as 
%       function [c,ceq] = mycon(x)
%        c = ...     % Compute nonlinear inequalities at x.
%        ceq = ...   % Compute nonlinear equalities at x.
% These are set in HIVEpiModel.m
                        
%% options

        % options = optimoptions(@fmincon,'Display','iter-detailed','Algorithm','interior-point','TolFun',0.01); %'interior-point','sqp','active-set'
        options = optimoptions('fmincon','Display','iter-detailed','Algorithm','interior-point','TolFun',0.001,'MaxFunEvals',15000); %'interior-point','sqp','active-set'

end

%% Trash code        
%     %CostTestPP_ByCohort (1 x numStrats?) = average cost per person for testing
%     %NEEDS TO BE CHECKED
%         numberUnawareWilling = sum(Compartments([Params.B1,Params.C1,...
%                 Params.D1,Params.E1,Params.F1,Params.F0],:))';         
%         totalWillingToBeTested = numberUnawareWilling'+Compartments(Params.A1,:);
%         probTestHIVNegByCohort = Compartments(Params.A1,:)./max(totalWillingToBeTested,1);
%         probTestHIVPosByCohort = ...
%             (totalWillingToBeTested - Compartments(Params.A1,:))./max(totalWillingToBeTested,1);        
%         CostTestPP_ByCohort = Params.costPP_TestNegRapid * probTestHIVNegByCohort ...
%                     .* Params.tt_pctRapid_r_3' + ...
%                     Params.costPP_TestPosRapid * probTestHIVPosByCohort .* ...
%                     Params.tt_pctRapid_r_3'  + ...
%                     Params.costPP_TestNegConv * probTestHIVNegByCohort .* ...
%                     (1 - Params.tt_pctRapid_r_3)' + ...
%                     Params.costPP_TestPosConv * probTestHIVPosByCohort .* ...
%                     (1 - Params.tt_pctRapid_r_3)';
