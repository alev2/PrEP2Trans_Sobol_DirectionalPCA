function [x0,A,b,Aeq,beq,lb,ub,options]=Calib_GenOptComponents(CalibParams, initValueSource, algorithm, calibStopCondition, calibTimeToRun, calibOVStop)
%% Set up supporting variables
    
    % Initialize / pre-allocate variables
    ParamFields = fieldnames(CalibParams);
    nCalibratedParameters = length(ParamFields);
    LB=1;
    UB=2;
    
%% x0: initial solution; can be a scalar, vector, or matrix.

x0(nCalibratedParameters,1)=0;
if initValueSource == 1 %base values 
    for i = 1:nCalibratedParameters
        x0(i) = CalibParams.(ParamFields{i}).baseValue;
    end    
else  %random
    for i = 1:nCalibratedParameters
        x0(i) = CalibParams.(ParamFields{i}).range(LB)+ ...
            rand(1)*(CalibParams.(ParamFields{i}).range(UB)-CalibParams.(ParamFields{i}).range(LB));
    end       
     
end
        
%% A / b: defines linear inequalities A*x <= b. If no inequalities exist, set A = [] and b = []. If A large w/ few nonzero entries, make sparse.

    % Budget constraint: Allocate only as much funding as is available
    %       Sum of allocations (DVs) over all grantees, B & C <= Budget
    A = [];
    b = [];
  
%% Aeq / beq: defines linear equalities Aeq*x = beq. If no equalities exist, set Aeq = [] and beq = []
    
    Aeq=[];
    beq=[];
        
%% lb / ub: lower and upper bounds on the design variables in x, so that the solution is always in the range lb <= x <= ub. If no bounds exist, set lb = [] and/or ub = [].If x(i) is unbounded below, set lb(i) = -Inf, and if x(i) is unbounded above, set ub(i) = Inf.

    for i = 1:nCalibratedParameters
        lb(i) = CalibParams.(ParamFields{i}).range(LB);
        ub(i) = CalibParams.(ParamFields{i}).range(UB);
    end   

%% options

%stop condition
%if stop condition is 1, then calibTimeToRun is used and the algorithms are
%ran for a certain amount of time.  If stop condition is 2, then the
%algorithms will stop absed on objective value, this is done in the nested
%function ObjValue_Calib in HIVEpiMode.m, so here the time will be set to
%infinite.

    if calibStopCondition == 2
        calibTimeToRun = Inf;
    end
    
    if algorithm == 1
        
        options = psoptimset('Display','iter','CompletePoll','Off','PollingOrder','Success','TimeLimit',calibTimeToRun);

    elseif algorithm == 2 %simulannealbnd
        
        options = saoptimset('Display','iter','TimeLimit',calibTimeToRun,'AnnealingFcn',@annealingboltz);
          
    else %fmincon
        
        options = optimoptions('fmincon','Display','iter-detailed','Algorithm','interior-point','TolFun',0.0001,'MaxFunEvals',10000); %'interior-point','sqp','active-set'
      
    end

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

    % Constraint: Only individuals in eligible populations and model states receive each intervention 
    % and Maximum intervention reach for eligible individuals
%                 yCoeff(i,j)=[1/intncost]*[includePop(i)] 
%     for i = 1:Params.numStrats 
%         yCoeff(i,1)= 1/CostTestPP_ByCohort(i);
%         for j = 1:Params.numIntns
%         end
%     end
%     
%     costPP_NotifyPosConv
% costPP_LTCFirst
% costPP_LTCAfterFirst
% costPP_ARTInitiation
% costPP_TxAdherence


%   
%     
% 
%     %includePop (1 x numStrats) = indicators that subpop i in subset to be considered in opt (for abridged version, all ones)
%     %NEEDS TO BE CHECKED
%     includePop = ones(Params.numStrats)';    

%     %	nonlcon: defines nonlinear inequalities c(x) or equalities ceq(x). 
%     %       fmincon optimizes such that c(x) <= 0 and ceq(x) = 0. 
%     %       accepts a vector x and returns the two vectors c and ceq. 
%     %       c is a vector that contains the nonlinear inequalities evaluated at x, 
%     %       and ceq is a vector that contains the nonlinear equalities evaluated at x. 
%     %       nonlcon should be specified as a function handle to a file or to 
%     %       an anonymous function, such as 
%     %       mycon: x = fmincon(@myfun,x0,A,b,Aeq,beq,lb,ub,@mycon), 
%     %       where mycon is a MATLAB function such as 
%     %       function [c,ceq] = mycon(x)
%     %        c = ...     % Compute nonlinear inequalities at x.
%     %        ceq = ...   % Compute nonlinear equalities at x.
% 
% 
%     % options = optimoptions(@fmincon,'Display','iter-detailed','Algorithm','interior-point','TolFun',0.01); %'interior-point','sqp','active-set'
% 
%     options = optimoptions('fmincon','Display','iter-detailed','Algorithm','interior-point','TolFun',0.01); %'interior-point','sqp','active-set'
