%This generates all of the relevant simulations for the HOPE analysis for
%transmission groups.

%sheet in HOPE file where we change the mixing matrix
inputSheetName='Behaviors';

%Excel files where we output our results
tablePathFull={...
    'PrEPInterventionScenarios_MixingMatrixBaseline.xlsx',...
    'PrEPInterventionScenarios_MixingMatrix1.xlsx',...
    'PrEPInterventionScenarios_MixingMatrix2.xlsx',...
    'PrEPInterventionScenarios_MixingMatrix3.xlsx',...
    'PrEPInterventionScenarios_MixingMatrix4.xlsx'};



%function to make the excel sheet ranges easier.
rowRange=...
    @(Col1,Col2,numberRange)strcat(Col1,num2str(numberRange),':',Col2,num2str(numberRange));
%indices
outputRows_5Years=[8;4;15;26;37];
outputRows_10Years=[9;6;17;28;39];
outputRows_20Years=[10;8;19;30;41];


%cells containing entries in mixing matrix to be changed - MUST ACCOUNT NOW
%FOR HIGH/LO RISK STRUCTURES
mixCells={'H749','G750','G752','H751','I753','H753'};
yearArray=(2010:1:2042)';


%mixing matrices (for reference) - NOW MUST ACCOUNT FOR RISK STRUCTURES
%baseline:
% HETFLo: 99.9% HETM
% HETMLo: 99.9% HETF
% HETFHi: 97.6% HETM (2.1% MSM)
% HETMHi: 95.5% HETF
% MSM: 13.1% HETF (86.6% MSM)

%mix 1:
% HETFLo: 99.8% HETM
% HETMLo: 99.8% HETF 
% HETFHi: 96.6% HETM (3.1% MSM(
% HETMHi: 94.5% HETF 
% MSM: 15.1% HETF (84.6% MSM)

%Mix 2:
% HETFLo: 99.7% HETM
% HETMLo: 99.7% HETF
% HETFHi: 95.6% HETM (4.1% MSM)
% HETMHi: 93.5% HETF
% MSM: 17.1% HETF (82.6% MSM)

%Mix 3:
% HETFLo: 99.95% HETM
% HETMLo: 99.95% HETF
% HETFHi: 98.6% HETM (1.1% MSM)
% HETMHi: 96.5% HETF 
% MSM: 11.1% HETF (88.6% MSM)

%Mix 4:
% HETFLo: 99.99% HETM 
% HETMLo: 99.99% HETF
% HETFHi: 99.6% HETM (.1% MSM)
% HETMHi: 97.5% HETF
% MSM: 9.1% HETF (90.6% MSM)

%baseline mixing matrix 
baseValues=[.999  .999 .9764764333 .955 .866105747187  .1314023621 ];

%alternate mixing matrices
mixingMatFull=[...
 .999  .999 .9764764333 .955 .866105747187  .1314023621; ... baseline
 .998  .998 .9664764333 .945 .846105747187  .1514023621; %Mixing matrix 1
 .997  .997 .9564764333 .935 .826105747187  .1714023621;...%Mixing matrix 2
 .9995  .9995 .9864764333 .965 .886105747187  .1114023621;... %Mixing matrix 3
  .9999  .9999 .9964764333 .975 .906105747187  .0914023621]; %Mixing matrix 4


%names of different PrEP intervention scenario files.
%There are three - comment the one you are not currently running.

%scenario 1
%sheet on each excel file where the results are placed
%outputSheetName='Scenario 1';
%outcomeTablePath='.\Scenario1\';
%scenarioSet={...
%    'HOPE Model V09_04_LB906_2_2022_04_08.xlsm',...
%    'PrEPToMSMThreeTimes.xlsm',...
%    'PrEPHETFemalesThreeTimes.xlsm',...
%    'PrEPHETMalesThreeTimes.xlsm',...
%    'PrEPPWIDThreeTimes.xlsm',...
%};

%scenario2
%sheet on each excel file where the results are placed
% outputSheetName='Scenario 2';
% outcomeTablePath='.\Scenario2\';
% scenarioSet={...
%     'HOPE Model V09_04_LB906_2_2022_04_08.xlsm',...
%     'PrEPMSMHundredThousand.xlsm',...
%     'Copy_of_PrEPHETFemalesHundredThousand.xlsm',...
%     'PrEPHETMalesHundredThousand.xlsm',...
%     'PrEPPWIDHundredThousand.xlsm',...
% };

%scenario3
%sheet on each excel file where the results are placed
outputSheetName='ResultsTable';
outcomeTablePath='.\Scenario3\';
scenarioSet={...
   'LatestHOPE Model V10_04_LB20230224_2_20231031_wTargetMultipliers.xlsm',...
   'PrEPMSMTwoTimes.xlsm',...
   'PrEPHETFTwoTimes.xlsm',...
   'PrEPHETMTwotimes.xlsm',...
   'PrEPPWIDTwotimes.xlsm',...
};


interventionLabel={...
    'StatusQuo',...
    'HETF2x',...
    'MSM2x',...
    'HETM2x',...
    'PWID2x',...
};
load('Params.mat');

varNms={'Year','IncidenceMSM','IncidenceHETM','IncidenceHETF','IncidencePWID','IncidenceTotal','PersonYears'};


%J-loop: we iterate through the mixing matrix configurations
for j=1:size(tablePathFull,2)


    %get the output excel file and mixing matrix
    tablePath=tablePathFull{j};
    mixingMat=mixingMatFull(j,:);
    
    %i-loop: we iterate through the different PrEP intervention
    %configurations

    for i=1:size(scenarioSet,2)

        caseIdentifier=strcat( ...
            outcomeTablePath,...
               strrep(...
                   strrep(...
                       tablePath,'PrEPInterventionScenarios_',''...
                       )...
                 ,'.xlsx',''...
                 ), interventionLabel{i})

        %get the intervention scenario
        hopePath=scenarioSet{i};
        
        %redefine the mixing matrix
        xlswrite(hopePath,mixingMat(1),inputSheetName, mixCells{1} );
        xlswrite(hopePath,mixingMat(2),inputSheetName, mixCells{2} );
        xlswrite(hopePath,mixingMat(3),inputSheetName, mixCells{3} );
        xlswrite(hopePath,mixingMat(4),inputSheetName, mixCells{4} );
        xlswrite(hopePath,mixingMat(5),inputSheetName, mixCells{5} );
        xlswrite(hopePath,mixingMat(6),inputSheetName, mixCells{6} );

        %run simulation and get results
        Outcome=HIVEpiModel(hopePath);
        getOutcomes;
            
        %assemble relevant table
        outcomeTable=[IncidenceMSM IncidenceHETM IncidenceHETF IncidencePWID IncidenceTotal PersonYears];
        OutcomeTable=array2table([yearArray outcomeTable],'VariableNames',varNms);

        save(...
            strcat(...
            outcomeTablePath,...            
               strrep(...
                   strrep(...
                       tablePath,'PrEPInterventionScenarios_',''...
                       )...
                 ,'.xlsx',''...
                 ), interventionLabel{i}...
               ), 'OutcomeTable'...
             );

        %save base case values (bc in excel sheet we report infections
        %saved, so we need these at each step!
        if(i==1)
            baseCase_5=round(sum(outcomeTable(index2023:index2027,:)));
            baseCase_10=round(sum(outcomeTable(index2023:index2032,:)));
            baseCase_20=round(sum(outcomeTable(index2023:index2042,:)));        
    
            %write base case
            xlswrite(tablePath,baseCase_5(1:end-1),outputSheetName, ...
                rowRange('Q','U',outputRows_5Years(i)));
            %xlswrite(tablePath,baseCase_5(end),outputSheetName, ...
            %    rowRange('I','I',outputRows_5Years(i)));
            
            xlswrite(tablePath,baseCase_10(1:end-1),outputSheetName, ...
                rowRange('Q','U',outputRows_10Years(i)));
            %xlswrite(tablePath,baseCase_10(end),outputSheetName, ...
            %    rowRange('I','I',outputRows_10Years(i)));
    
            xlswrite(tablePath,baseCase_20(1:end-1),outputSheetName, ...
                rowRange('D','H',outputRows_20Years(i)));
            %xlswrite(tablePath,baseCase_20(end),outputSheetName, ...
            %    rowRange('I','I',outputRows_20Years(i)));
          else
            %here we are not in base case: we save infections saved!
            infSaved_5=abs(baseCase_5-round(sum(outcomeTable(index2023:index2027,:))));
            infSaved_10=abs(baseCase_10-round(sum(outcomeTable(index2023:index2032,:))));
            infSaved_20=abs(baseCase_20-round(sum(outcomeTable(index2023:index2042,:))));
    
            xlswrite(tablePath,infSaved_5(1:end-1),outputSheetName, ...
                rowRange('B','F',outputRows_5Years(i)));
            %xlswrite(tablePath,infSaved_5(end),outputSheetName, ...
            %    rowRange('I','I',outputRows_5Years(i)));
            
            xlswrite(tablePath,infSaved_10(1:end-1),outputSheetName, ...
                rowRange('B','F',outputRows_10Years(i)));
            %xlswrite(tablePath,infSaved_10(end),outputSheetName, ...
            %    rowRange('I','I',outputRows_10Years(i)));
    
            xlswrite(tablePath,infSaved_20(1:end-1),outputSheetName, ...
                rowRange('B','F',outputRows_20Years(i)));
            %xlswrite(tablePath,infSaved_20(end),outputSheetName, ...
            %    rowRange('I','I',outputRows_20Years(i)));    
        end
        
    
    
    end

    
end

%done!