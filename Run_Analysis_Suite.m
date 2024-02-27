%This generates all of the relevant simulations for the H/E analysis and
%puts them cleanly in an excel file.

%If the excel file + sheet names are defined in such a way, this will
%automatically populate an existing 
% spreadsheet of the format "Disparities analysis results - 2023_03_01.xlsx 
%with the new results - this can make running and compiling results very
%easy. Simply make a new copy of that spreadsheet, name it accordingly with
%indicating the current run version, and the results will be automatically
%populated and tabulated according to the existing excel functions present
%in the original document.



%HOPE model version and file naming format
%baseString='.\Outcomes_LB20230224_2_4_25\inputFiles\';
baseString='.\Outcomes_COVID\inputFiles\';
baseString=strcat(baseString,'HOPE Model V09_12_2023_10_3_Calib_');
simulationType='RTI_'; 
%simulationType=''; 
calibSet='LB20230224_2';
%dateTag='_4_25';
dateTag='';
%simulationType='RTI_';
%calibSet='LB3';


%where the output file will be saved and file prefix
%directoryPath=strcat('.\Outcomes_',calibSet,dateTag,'\outcome');
directoryPath=strcat('.\Outcomes_COVID','\outcome');

%where the excel table will be saved and/or the file name
%directoryPath_Table=strcat('.\Outcomes_',calibSet,dateTag,'\Disparities analysis results - 2023_10_3_COVID.xlsx');    
directoryPath_Table=strcat('.\Outcomes_COVID\Disparities analysis results - 2023_10_5_COVID.xlsx');    
%%%directoryPath_Table=strcat('.\Outcomes_',calibSet,'\Disparities analysis results -Copy.xlsx');    




%Names of the different scenarios at the end of their corresponding excel
%files
% scenarioSet={...
%     '_BaseScenarioA_2',...
%     '_ScenarioB',...
%     '_ScenarioC',...
%     '_ScenarioD',...
%     '_ScenarioE'...
% };
scenarioSet={...
    '_BaseScenarioA_COVID',...
    '_ScenarioB_COVID',...
    '_ScenarioC_COVID',...
    '_ScenarioD_COVID',...
    '_ScenarioE_COVID'...
};

%How the different sheets in the excel sheet should be named when saving
%the results
sheetNames={...
    'A. Results',...
    'B. Results',...
    'C. Results',...
    'D. Results',...
    'E. Results'...
};

%No need to touch
simulationPathBase=strcat(baseString, simulationType, calibSet);
dotXLS='.xlsm';
%variable names for the outcomes of interest. 
load('.\variableNames.mat');
%load('.\variableNamesUndisc.mat');
discCell1='G73';
discCell2='G75';
discRate=0;
changeDiscRate=1;

for i=1:size(scenarioSet,2)

    fprintf('Simulation %s. \n',replace(scenarioSet{i},{'Scenario','Base','_'},{'','',''}));
    hopePath=strcat(simulationPathBase,scenarioSet{i},dotXLS);
    filePath=strcat(directoryPath,scenarioSet{i});

    if(changeDiscRate==1)
        fprintf('Changing the discount rates to %g.\n',discRate);
        xlswrite(hopePath,discRate,'Model Settings',discCell1);
        xlswrite(hopePath,discRate,'Model Settings',discCell2);
    end

    outcome=HIVEpiModel(hopePath);    
    save(filePath,'outcome');

    writetable(...
        generate_RelevantTable(outcome,variableNames),...
        directoryPath_Table, ...
        'Sheet',sheetNames{i},...
        'Range','A2:AK28',...
        'WriteVariableNames',1);

end


