function CalibOutputStruct = Calib_collectResults(Results,CalibOutputStruct, sample)

% Purpose: record key outcomes after each calibration run

% Called from: HIVEpiModel.m


% Determine last outcome year
    %LastOutcomeYr = Results.set_LastOutcomeYr;

        % 2010 outcomes-removed. Bates 6/11/19        

        % HIV Prevalence in 2016. Allaire 6/11/19. Changed to 2018 4/21/20
        % Updated to 2019. Clinkscales 6/29/2021
        CalibOutputStruct.calib_HIVPrevalence2019.paramValue(sample)= ...
            Results.calib_HIVPrevalence2019;

        %Prevalence by Race/Ethnicity in 2016 (updated to 2016; Allaire
        %6/13/2019). Updated to 2018 4/21/20. Bates
        %Updated to 2019. Clinkscales 6/29/2021
        CalibOutputStruct.calib_HIVPrevalence2019_B.paramValue(sample)= ...
            Results.calib_HIVPrevalence2019_B;
        CalibOutputStruct.calib_HIVPrevalence2019_H.paramValue(sample)= ...
            Results.calib_HIVPrevalence2019_H;
        CalibOutputStruct.calib_HIVPrevalence2019_O.paramValue(sample)= ...
            Results.calib_HIVPrevalence2019_O;
        
        %Prevalence by transmission group in 2018. Added 4/21/20
        %Updated to 2019. Clinkscales 6/29/2021
        CalibOutputStruct.calib_HIVPrevalence2019_HETM.paramValue(sample)= ...
            Results.calib_HIVPrevalence2019_HETM;
        CalibOutputStruct.calib_HIVPrevalence2019_HETF.paramValue(sample)= ...
            Results.calib_HIVPrevalence2019_HETF;
        CalibOutputStruct.calib_HIVPrevalence2019_MSM.paramValue(sample)= ...
            Results.calib_HIVPrevalence2019_MSM; 
        CalibOutputStruct.calib_HIVPrevalence2019_IDUM.paramValue(sample)= ...
            Results.calib_HIVPrevalence2019_IDUM;
        CalibOutputStruct.calib_HIVPrevalence2019_IDUF.paramValue(sample)= ...
            Results.calib_HIVPrevalence2019_IDUF;        
            
        %Prevalence by age in 2016. Bates 6/13/19
        CalibOutputStruct.calib_HIVPrevalence2016_13_24.paramValue(sample)= ...
            Results.calib_HIVPrevalence2016_13_24;
        CalibOutputStruct.calib_HIVPrevalence2016_25_34.paramValue(sample)= ...
            Results.calib_HIVPrevalence2016_25_34;
        CalibOutputStruct.calib_HIVPrevalence2016_35_44.paramValue(sample)= ...
            Results.calib_HIVPrevalence2016_35_44;
        CalibOutputStruct.calib_HIVPrevalence2016_45_54.paramValue(sample)= ...
            Results.calib_HIVPrevalence2016_45_54;
        CalibOutputStruct.calib_HIVPrevalence2016_55_64.paramValue(sample)= ...
            Results.calib_HIVPrevalence2016_55_64;
        CalibOutputStruct.calib_HIVPrevalence2016_65.paramValue(sample)= ...
            Results.calib_HIVPrevalence2016_65;
        
        %Prevalence by age and r/e in 2016. Bates 6/18/19
        CalibOutputStruct.calib_HIVPrevalence2016_B_13_24.paramValue(sample)= ...
            Results.calib_HIVPrevalence2016_B_13_24;
        CalibOutputStruct.calib_HIVPrevalence2016_H_13_24.paramValue(sample)= ...
            Results.calib_HIVPrevalence2016_H_13_24;
        CalibOutputStruct.calib_HIVPrevalence2016_O_13_24.paramValue(sample)= ...
            Results.calib_HIVPrevalence2016_O_13_24;
        
        CalibOutputStruct.calib_HIVPrevalence2016_B_25_34.paramValue(sample)= ...
            Results.calib_HIVPrevalence2016_B_25_34;
        CalibOutputStruct.calib_HIVPrevalence2016_H_25_34.paramValue(sample)= ...
            Results.calib_HIVPrevalence2016_H_25_34;
        CalibOutputStruct.calib_HIVPrevalence2016_O_25_34.paramValue(sample)= ...
            Results.calib_HIVPrevalence2016_O_25_34;
        
        % new age groups - removed 13-34 and 35-64. Bates 6/13/19. 

        %On ART, Not VLS- removed 4/7/17 Laurel Bates
        %VLS among PLWH. 
        %VLS among Diagnosed. 
       
        % Continuum of Care in 2015
        
        %Unaware. - removed. Bates 6/13/19
        
        %Diagnosed.  Updated to 2016- Bates 6/14/19. By race updated to 2019- Bates
        %11/22/22
        CalibOutputStruct.calib_TTdist2019_Diagnosed_B.paramValue(sample)  = ...
            Results.calib_TTdist2019_Diagnosed_B;
        CalibOutputStruct.calib_TTdist2019_Diagnosed_H.paramValue(sample)  = ...
            Results.calib_TTdist2019_Diagnosed_H;
        CalibOutputStruct.calib_TTdist2019_Diagnosed_O.paramValue(sample)  = ...
            Results.calib_TTdist2019_Diagnosed_O;
        %By age. Updated to 2016- Bates 6/17/19
        CalibOutputStruct.calib_TTdist2016_Diagnosed_13_24.paramValue(sample)= ...
            Results.calib_TTdist2016_Diagnosed_13_24;
        CalibOutputStruct.calib_TTdist2016_Diagnosed_25_34.paramValue(sample)= ...
            Results.calib_TTdist2016_Diagnosed_25_34;
        CalibOutputStruct.calib_TTdist2016_Diagnosed_35_44.paramValue(sample)= ...
            Results.calib_TTdist2016_Diagnosed_35_44;
        CalibOutputStruct.calib_TTdist2016_Diagnosed_45_54.paramValue(sample)= ...
            Results.calib_TTdist2016_Diagnosed_45_54;
        CalibOutputStruct.calib_TTdist2016_Diagnosed_55_64.paramValue(sample)= ...
            Results.calib_TTdist2016_Diagnosed_55_64;
        CalibOutputStruct.calib_TTdist2016_Diagnosed_65.paramValue(sample)= ...
            Results.calib_TTdist2016_Diagnosed_65;                        
          
        %By race and age (edited by JC on 12/15/2017)
        %Updated to 2016 and changed age groups to 13-24. Bates 6/17/19
        CalibOutputStruct.calib_TTdist2016_Diagnosed_MSM_B_13_24.paramValue(sample) = ...
            Results.calib_TTdist2016_Diagnosed_MSM_B_13_24;
        CalibOutputStruct.calib_TTdist2016_Diagnosed_MSM_H_13_24.paramValue(sample) = ...
            Results.calib_TTdist2016_Diagnosed_MSM_H_13_24;
        CalibOutputStruct.calib_TTdist2016_Diagnosed_MSM_O_13_24.paramValue(sample) = ...
            Results.calib_TTdist2016_Diagnosed_MSM_O_13_24;
        
        CalibOutputStruct.calib_TTdist2016_Diagnosed_MSM_B_25_34.paramValue(sample) = ...
            Results.calib_TTdist2016_Diagnosed_MSM_B_25_34;
        CalibOutputStruct.calib_TTdist2016_Diagnosed_MSM_H_25_34.paramValue(sample) = ...
            Results.calib_TTdist2016_Diagnosed_MSM_H_25_34;
        CalibOutputStruct.calib_TTdist2016_Diagnosed_MSM_O_25_34.paramValue(sample) = ...
            Results.calib_TTdist2016_Diagnosed_MSM_O_25_34;
        
        %By transmission group 6/17/19. Updated to 2018 4/21/20. Updated to
        %2019. Clinkscales 6/29/2021
        CalibOutputStruct.calib_TTdist2019_Diagnosed_HET.paramValue(sample) = ...
            Results.calib_TTdist2019_Diagnosed_HET;
        CalibOutputStruct.calib_TTdist2019_Diagnosed_MSM.paramValue(sample) = ...
            Results.calib_TTdist2019_Diagnosed_MSM;
        CalibOutputStruct.calib_TTdist2019_Diagnosed_PWID.paramValue(sample) = ...
            Results.calib_TTdist2019_Diagnosed_PWID;
        
        % total. Updated to 2016 6/17/19. Updated to 2018 4/21/20
        CalibOutputStruct.calib_TTdist2019_Diagnosed_Total.paramValue(sample)= ...
            Results.calib_TTdist2019_Diagnosed_Total;
              
        %In Care. Removed Bates 6/17/19. Added with 2019 values. Bates 11/4/2022
        CalibOutputStruct.calib_TTdist2019_InCareAmongDiag_B.paramValue(sample)  = ...
            Results.calib_TTdist2019_InCareAmongDiag_B;
        CalibOutputStruct.calib_TTdist2019_InCareAmongDiag_H.paramValue(sample)  = ...
            Results.calib_TTdist2019_InCareAmongDiag_H;
        CalibOutputStruct.calib_TTdist2019_InCareAmongDiag_O.paramValue(sample)  = ...
            Results.calib_TTdist2019_InCareAmongDiag_O;
        CalibOutputStruct.calib_TTdist2019_InCareAmongDiag_Total.paramValue(sample)  = ...
            Results.calib_TTdist2019_InCareAmongDiag_Total;


        %Aware, No Care Removed. Bates 6/17/19

        %On ART, Not VLS- removed 4/7/17 Laurel Bates
        
        %On ART Removed. Bates 6/17/19

        %VLS among PLWH. Removed Bates 6/19/19 
        
        %VLS among Diagnosed. Changed to 2015. Changed to 2017 4/21/20.
        %Changed to 2018 5/29/20. Updated to 2019 11/22/22
        CalibOutputStruct.calib_TTdist2019_VLSamongdiag_B.paramValue(sample)  = ...
            Results.calib_TTdist2019_VLSamongdiag_B;
        CalibOutputStruct.calib_TTdist2019_VLSamongdiag_H.paramValue(sample)  = ...
            Results.calib_TTdist2019_VLSamongdiag_H;
        CalibOutputStruct.calib_TTdist2019_VLSamongdiag_O.paramValue(sample)  = ...
            Results.calib_TTdist2019_VLSamongdiag_O;
        
        %VLS among transmission group. Added 4/21/20. Bates. Changed to 2018 5/29/20
        CalibOutputStruct.calib_TTdist2018_VLSamongdiag_HET.paramValue(sample)  = ...
            Results.calib_TTdist2018_VLSamongdiag_HET;
        CalibOutputStruct.calib_TTdist2018_VLSamongdiag_MSM.paramValue(sample)  = ...
            Results.calib_TTdist2018_VLSamongdiag_MSM;
        CalibOutputStruct.calib_TTdist2018_VLSamongdiag_PWID.paramValue(sample)  = ...
            Results.calib_TTdist2018_VLSamongdiag_PWID;        
        
        %by age- added 6/20/19
        CalibOutputStruct.calib_TTdist2015_VLSamongdiag_13_24.paramValue(sample)  = ...
            Results.calib_TTdist2015_VLSamongdiag_13_24;
        CalibOutputStruct.calib_TTdist2015_VLSamongdiag_25_34.paramValue(sample)  = ...
            Results.calib_TTdist2015_VLSamongdiag_25_34;
        CalibOutputStruct.calib_TTdist2015_VLSamongdiag_35_44.paramValue(sample)= ...
            Results.calib_TTdist2015_VLSamongdiag_35_44;
        CalibOutputStruct.calib_TTdist2015_VLSamongdiag_45_54.paramValue(sample)= ...
            Results.calib_TTdist2015_VLSamongdiag_45_54;
        CalibOutputStruct.calib_TTdist2015_VLSamongdiag_55_64.paramValue(sample)= ...
            Results.calib_TTdist2015_VLSamongdiag_55_64;
        CalibOutputStruct.calib_TTdist2015_VLSamongdiag_65.paramValue(sample)= ...
            Results.calib_TTdist2015_VLSamongdiag_65;
        
        %Young MSM by race - added 8/15/17. Changed to 2016 6/19/19. Bates.
        %Changed age group to 13-24 instead of 18-24. 
        CalibOutputStruct.calib_TTdist2016_VLSamongdiag_13_24_MSM_B.paramValue(sample)  = ...
            Results.calib_TTdist2016_VLSamongdiag_13_24_MSM_B;
        CalibOutputStruct.calib_TTdist2016_VLSamongdiag_13_24_MSM_H.paramValue(sample)  = ...
            Results.calib_TTdist2016_VLSamongdiag_13_24_MSM_H;
        CalibOutputStruct.calib_TTdist2016_VLSamongdiag_13_24_MSM_O.paramValue(sample)  = ...
            Results.calib_TTdist2016_VLSamongdiag_13_24_MSM_O;
        
        CalibOutputStruct.calib_TTdist2016_VLSamongdiag_25_34_MSM_B.paramValue(sample)  = ...
            Results.calib_TTdist2016_VLSamongdiag_25_34_MSM_B;
        CalibOutputStruct.calib_TTdist2016_VLSamongdiag_25_34_MSM_H.paramValue(sample)  = ...
            Results.calib_TTdist2016_VLSamongdiag_25_34_MSM_H;
        CalibOutputStruct.calib_TTdist2016_VLSamongdiag_25_34_MSM_O.paramValue(sample)  = ...
            Results.calib_TTdist2016_VLSamongdiag_25_34_MSM_O;
        
        %total. Changed to 2017 4/21/20. Changed to 2018 5/29/20. Updated
        %to 2019 11/22/22
        CalibOutputStruct.calib_TTdist2019_VLSamongdiag_Total.paramValue(sample)  = ...
            Results.calib_TTdist2019_VLSamongdiag_Total;
    
        
        % Aware PLWH Deaths. Changed to 2016. 6/14/19. Laurel
        % Updated to 2019. Clinkscales 6/29/2021
        CalibOutputStruct.calib_AwarePLWHDeaths2019.paramValue(sample)= ...
            Results.calib_AwarePLWHDeaths2019;

        CalibOutputStruct.calib_AwarePLWHDeaths2016_45_54.paramValue(sample)  = ...
            Results.calib_AwarePLWHDeaths2016_45_54;        
        CalibOutputStruct.calib_AwarePLWHDeaths2016_55_64.paramValue(sample)  = ...
            Results.calib_AwarePLWHDeaths2016_55_64;        
        CalibOutputStruct.calib_AwarePLWHDeaths2016_65.paramValue(sample)  = ...
            Results.calib_AwarePLWHDeaths2016_65;
        
        % PLWH AIDS deaths in 2019, Viguerie. Aware from surveillance, and unaware is a guess (it should be low). Added 12/9/2022.
        CalibOutputStruct.calib_AwarePLWHAIDSDeaths2019.paramValue(sample)  = ...
            Results.calib_AwarePLWHAIDSDeaths2019;
        CalibOutputStruct.calib_UnawarePLWHAIDSDeaths2019.paramValue(sample)  = ...
            Results.calib_UnawarePLWHAIDSDeaths2019;


        % New Diagnoses. Removed results by race 6/17/19. 
        %Changed to 2016. 6/19/19
        % Added two new 2019 NewDiag. outcomes: Clinskcales 04/21/2022
        CalibOutputStruct.calib_NewDiagnoses2016_HETM.paramValue(sample) = ...
            Results.calib_NewDiagnoses2016_HETM;

        CalibOutputStruct.calib_NewDiagnoses2016_HETF.paramValue(sample)  = ...
            Results.calib_NewDiagnoses2016_HETF;

        CalibOutputStruct.calib_NewDiagnoses2016_MSM.paramValue(sample)  = ...
            Results.calib_NewDiagnoses2016_MSM;

        CalibOutputStruct.calib_NewDiagnoses2016_IDUM.paramValue(sample)  = ...
            Results.calib_NewDiagnoses2016_IDUM;

        CalibOutputStruct.calib_NewDiagnoses2016_IDUF.paramValue(sample)  = ...
            Results.calib_NewDiagnoses2016_IDUF;
   
        CalibOutputStruct.calib_NewDiagnoses2016_Total.paramValue(sample)  = ...
            Results.calib_NewDiagnoses2016_Total;
        
    %Clinkscales 04/21/2022: Added NewDiagnoses2019 calibration targets below  
    
        CalibOutputStruct.calib_NewDiagnoses2019_PctB1C1.paramValue(sample)  = ...
            Results.calib_NewDiagnoses2019_PctB1C1;
    
        CalibOutputStruct.calib_NewDiagnoses2019_Total.paramValue(sample)  = ...
            Results.calib_NewDiagnoses2019_Total;


    % New Infections
        % New Infections in 2016     
        % By race changed to 2018 4/21/20. Updated to 2019. Clinkscales
        % 6/29/2021
        CalibOutputStruct.calib_NewInfections2019_B.paramValue(sample) = ...
            Results.calib_NewInfections2019_B;

        CalibOutputStruct.calib_NewInfections2019_H.paramValue(sample) = ...
            Results.calib_NewInfections2019_H;

        CalibOutputStruct.calib_NewInfections2019_O.paramValue(sample) = ...
            Results.calib_NewInfections2019_O;
        
        %By transmission group. Added 4/21/20.
        CalibOutputStruct.calib_NewInfections2019_HETM.paramValue(sample) = ...
            Results.calib_NewInfections2019_HETM;

        CalibOutputStruct.calib_NewInfections2019_HETF.paramValue(sample) = ...
            Results.calib_NewInfections2019_HETF;

        CalibOutputStruct.calib_NewInfections2019_MSM.paramValue(sample) = ...
            Results.calib_NewInfections2019_MSM;        

        CalibOutputStruct.calib_NewInfections2019_IDUM.paramValue(sample) = ...
            Results.calib_NewInfections2019_IDUM;

        CalibOutputStruct.calib_NewInfections2019_IDUF.paramValue(sample) = ...
            Results.calib_NewInfections2019_IDUF;        
        
        
        %MSM by age- Added 6/20/17. Updated to 2016 and 18-24 to 13-24. Bates
        %6/20/16
        CalibOutputStruct.calib_NewInfections_2016_MSM_13_24.paramValue(sample)  = ...
            Results.calib_NewInfections_2016_MSM_13_24;

        CalibOutputStruct.calib_NewInfections_2016_MSM_13_24_B.paramValue(sample)  = ...
            Results.calib_NewInfections_2016_MSM_13_24_B;

        CalibOutputStruct.calib_NewInfections_2016_MSM_13_24_H.paramValue(sample)  = ...
            Results.calib_NewInfections_2016_MSM_13_24_H;

        CalibOutputStruct.calib_NewInfections_2016_MSM_13_24_O.paramValue(sample)  = ...
            Results.calib_NewInfections_2016_MSM_13_24_O;

        CalibOutputStruct.calib_NewInfections_2016_MSM_25_34.paramValue(sample)  = ...
            Results.calib_NewInfections_2016_MSM_25_34;

        CalibOutputStruct.calib_NewInfections_2016_MSM_25_34_B.paramValue(sample)  = ...
            Results.calib_NewInfections_2016_MSM_25_34_B;

        CalibOutputStruct.calib_NewInfections_2016_MSM_25_34_H.paramValue(sample)  = ...
            Results.calib_NewInfections_2016_MSM_25_34_H;

        CalibOutputStruct.calib_NewInfections_2016_MSM_25_34_O.paramValue(sample)  = ...
            Results.calib_NewInfections_2016_MSM_25_34_O;       
        %Updated to 2016. Bates. 6/20/19. Updated total to 4/21/20
        CalibOutputStruct.calib_NewInfections2019_Total.paramValue(sample)  = ...
            Results.calib_NewInfections2019_Total;    
    
        CalibOutputStruct.calib_NewInfections2016_45_54.paramValue(sample)  = ...
            Results.calib_NewInfections2016_45_54;        
        CalibOutputStruct.calib_NewInfections2016_55_64.paramValue(sample)  = ...
            Results.calib_NewInfections2016_55_64;        
        CalibOutputStruct.calib_NewInfections2016_65.paramValue(sample)  = ...
            Results.calib_NewInfections2016_65;
        
        %Number on PrEP. Added by JC 2021/08/02. 
        % Modified to 2019 and 2021. Bates 02/10/23
        
        CalibOutputStruct.calib_NumOnPrEP2019_Male.paramValue(sample) = Results.calib_NumOnPrEP2019_Male;
            
        CalibOutputStruct.calib_NumOnPrEP2019_Female.paramValue(sample) = Results.calib_NumOnPrEP2019_Female;

        CalibOutputStruct.calib_NumOnPrEP2019_Black.paramValue(sample) = Results.calib_NumOnPrEP2019_Black;

        CalibOutputStruct.calib_NumOnPrEP2019_Hispanic.paramValue(sample) = Results.calib_NumOnPrEP2019_Hispanic;

        CalibOutputStruct.calib_NumOnPrEP2019_Other.paramValue(sample) = Results.calib_NumOnPrEP2019_Other;

        CalibOutputStruct.calib_NumOnPrEP2019_Total.paramValue(sample) = Results.calib_NumOnPrEP2019_Total;

        CalibOutputStruct.calib_NumOnPrEP2021_Male.paramValue(sample) = Results.calib_NumOnPrEP2021_Male;
            
        CalibOutputStruct.calib_NumOnPrEP2021_Female.paramValue(sample) = Results.calib_NumOnPrEP2021_Female;

        CalibOutputStruct.calib_NumOnPrEP2021_Black.paramValue(sample) = Results.calib_NumOnPrEP2021_Black;

        CalibOutputStruct.calib_NumOnPrEP2021_Hispanic.paramValue(sample) = Results.calib_NumOnPrEP2021_Hispanic;

        CalibOutputStruct.calib_NumOnPrEP2021_Other.paramValue(sample) = Results.calib_NumOnPrEP2021_Other;

        CalibOutputStruct.calib_NumOnPrEP2021_Total.paramValue(sample) = Results.calib_NumOnPrEP2021_Total;
        
        %Changed from 2006 to 2010 and 2015 to 2016. Bates 6/20/19. Changed
        %to 2018 4/21/20. Updated to 2019. Clinkscales 6/29/2021
        CalibOutputStruct.calib_OverallPrev2019v2010.paramValue(sample) = ...
            Results.calib_OverallPrev2019v2010;
    
        CalibOutputStruct.calib_HETPrev2019v2010_LR.paramValue(sample) = ...
            Results.calib_HETPrev2019v2010_LR;

        CalibOutputStruct.calib_HETPrev2019v2010_HR.paramValue(sample) = ...
            Results.calib_HETPrev2019v2010_HR;
        
        % 2016 population targets. Updated to 2019. Bates 11/22/22
    
        CalibOutputStruct.calib_populationMaleBlk2019.paramValue(sample) = ...
            Results.calib_populationMaleBlk2019;
        CalibOutputStruct.calib_populationMaleHisp2019.paramValue(sample) = ...
            Results.calib_populationMaleHisp2019;
        CalibOutputStruct.calib_populationMaleOth2019.paramValue(sample) = ...
            Results.calib_populationMaleOth2019;

        CalibOutputStruct.calib_populationFemaleBlk2019.paramValue(sample) = ...
            Results.calib_populationFemaleBlk2019;
        CalibOutputStruct.calib_populationFemaleHisp2019.paramValue(sample) = ...
            Results.calib_populationFemaleHisp2019;
        CalibOutputStruct.calib_populationFemaleOth2019.paramValue(sample) = ...
            Results.calib_populationFemaleOth2019;
        
        %2018 population targets. Updated to 2019 values. Clinkscales
        %6/29/2021
        CalibOutputStruct.calib_population_1334_2019.paramValue(sample) = ...
            Results.calib_population_1334_2019;
        CalibOutputStruct.calib_population_3564_2019.paramValue(sample) = ...
            Results.calib_population_3564_2019;
        CalibOutputStruct.calib_population_65_2019.paramValue(sample) = ...
            Results.calib_population_65_2019;    
    
end