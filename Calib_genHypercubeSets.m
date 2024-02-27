function [CalibParams] = Calib_genHypercubeSets(nSamples, CalibParams)
%% Purpose: This code generates parameter sets using Latin Hypercube Sampling.

%% 1. Create set of random permutation of numbers between 1 and nSamples

%For each calibration parameter, we first create a set representing a 
%random permutation of numbers between 1 and the total number of samples 
%(N). The numbers in each set will represent the latin hypercube cell 
%number from which we will pull a value for the parameter in a specific 
%sample.

    field = fieldnames(CalibParams);
    N = length(field);
    for i = 1:N
        CalibParams.(field{i}).indexSet = randperm(nSamples);
    end
       
%% 2. Determine sets of parameters

% For a specific sample set and calibration parameter, we determine the 
% parameter’s sampled value as follows. First, if the minimum and maximum 
% calibration range for the parameter are equal, then, we set the parameter’s
% sampled value to their minimum/maximum range. If they are not equal, we 
% determine the size of each latin hypercube cell assuming there are N 
% equally distributed cells between the minimum and maximum calibration 
% range for the parameter. We then assume each of these cells is numbered 
% from 1 to N, increasing in value from the minimum to the maximum range. 
% Finally, based on the calibration parameter’s specified index for this 
% specific sample (determined in the code directly above), the parameter’s 
% sample value is set to a random value within the bounds of that 
% corresponding cell. We repeat this process for all calibration parameters 
% and for each sample set to be generated.

    for i = 1:nSamples
        for j = 1:N
            index = CalibParams.(field{j}).indexSet(i);
            min = CalibParams.(field{j}).range(1);
            max = CalibParams.(field{j}).range(2);
            if min == max
                CalibParams.(field{j}).paramValue(i) = min;
            else
                size = (max-min)/nSamples;
                cellmin = min + (index-1)*size;
                CalibParams.(field{j}).paramValue(i) = cellmin + rand()*size;
            end
        end
    end
end