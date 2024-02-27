

try 
    outTest=HIVEpiModel('HOPE Model V10_05_LB20230224_2_20231005.xlsm');
catch
    fprintf('Error caught!')
    outTest=HIVEpiModel('HOPE Model V10_05_LB20230224_2_20231005_Fresh.xlsm');
end