function [ObjValue, OOBObjValue, TEObjValue] = Calib_calcObjValue(CalibOutputStruct, CalibTargets, calibObjValueType)
%% Description of file

%Calculates aggregated measure of model error according to a weighted
%percent difference between model output and target datapoints
%
%Inputs:    outcomes -  Model outcome struct
%           weight - A vector containing the weights for the given outputs 
%
%Outputs:   objValue -  The calculated aggregated measure of model
%               error according to a weighted percent difference between
%               model output and target datapoints
%
%Note: To exclude an output from the error calculation, set its weight to 0 and target to 1

%% 1. Initialize variables / other set-up
  
    % User setting for which objective value to use
    nObjFxnChoice = calibObjValueType; % 1= OOBObjValue, 2 = TEObjValue
    
    % Collect names and count number of fields in CalibTargets
    TargetFields = fieldnames(CalibTargets);
    nTargetOutcomes = length(TargetFields);
    ResultsFields = fieldnames(CalibOutputStruct);
    nResults = length(ResultsFields);
    
    % Pre-allocate / initialize
    weights(nTargetOutcomes,1)=0;
    modelOutcomeValues(nResults,1)=0;
    targetValues(nTargetOutcomes,1)=0;
    target_LB(nTargetOutcomes,1)=0;
    target_UB(nTargetOutcomes,1)=0;    
    penalty(nTargetOutcomes,1)=0;    
    LB = 1;
    UB = 2;
    
    % Set value of penalty if out of bounds
    OOBPenalty = 1000;
      
%% 2. Pull relevant outcomes from the model
 
    % Pull weights and target values associated with each target outcomes
    for i = 1:nTargetOutcomes
        weights(i,1)=CalibTargets.(TargetFields{i}).weight;
        targetValues(i,1)=CalibTargets.(TargetFields{i}).targetValue;
        target_LB(i,1)=CalibTargets.(TargetFields{i}).range(LB);
        target_UB(i,1)=CalibTargets.(TargetFields{i}).range(UB);
    end
    
    % Re-scale weights
    sumOfWeights = sum(weights);
    if sumOfWeights == 0 
        sumOfWeights = 1;
    end
    normWeights = weights./sumOfWeights;
    
    
    % Pull model outcomes resulting with each target outcomes into modelOutcomeValues
    for i = 1:nResults
        modelOutcomeValues(i,1)= CalibOutputStruct.(TargetFields{i}).paramValue;
        if or(modelOutcomeValues(i,1)<target_LB(i,1),modelOutcomeValues(i,1)>target_UB(i,1))
            penalty(i,1)= OOBPenalty;
        else
            penalty(i,1)=1;            
        end
    end    
        
%% 3. Calculate aggregated error measures
    
    % OOB penalty objective measure
    errorMeasureNum = modelOutcomeValues-targetValues;
    errorMeasureDen = (modelOutcomeValues+targetValues)./2;
    errorMeasure = penalty.*(errorMeasureNum./errorMeasureDen).^2;
    OOBObjValue = transpose(errorMeasure)*normWeights;

    % Target error objective value
    absErrors = abs(modelOutcomeValues-targetValues);
    relErrors = absErrors./targetValues;
    TEObjValue = transpose(relErrors)*normWeights;
    
    switch nObjFxnChoice
        case 1 %OOBObjValue
            ObjValue = OOBObjValue;
        case 2 %TEObjValue
            ObjValue = TEObjValue;
    
    end 

end 