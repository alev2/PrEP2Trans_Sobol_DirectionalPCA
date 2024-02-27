load('Params.mat');
PersonYears=sum(Outcome.ann_PersonYears_onPrEP,2)
ann_NumberEligForPrEP_MSM = sum(Outcome.ann_NumberEligForPrEP * Params.popIndicator(:,Params.pop_MSM),2);  
ann_NumberEligForPrEP_IDU = sum(Outcome.ann_NumberEligForPrEP * Params.popIndicator(:,Params.pop_IDU),2);
ann_NumberEligForPrEP_HRHET_M = sum(Outcome.ann_NumberEligForPrEP * Params.popSexIndicator(:,Params.popSex_HETM),2);
ann_NumberEligForPrEP_HRHET_F = sum(Outcome.ann_NumberEligForPrEP * Params.popSexIndicator(:,Params.popSex_HETF),2);
NumberOnPrEPMSM=sum(Outcome.ann_NumberOnPrEP_HighRiskMSM,2)
NumberOnPrEPPWID=sum(Outcome.ann_NumberOnPrEP_HighRiskIDUs,2)
NumberOnPrEPHETM=sum(Outcome.ann_NumberOnPrEP_HighRiskHETs_M,2)
NumberOnPrEPHETF=sum(Outcome.ann_NumberOnPrEP_HighRiskHETs_F,2)
CoverageMSM=sum(Outcome.ann_PctEligOnPrEP_HRMSM,2)
COveragePWID=sum(Outcome.ann_PctEligOnPrEP_PWID,2)
ann_PctEligOnPrEP_HRHET_M=sum(Outcome.ann_NumberOnPrEP_HighRiskHETs_M,2) ./sum(Outcome.ann_NumberEligForPrEP * Params.popSexIndicator(:,Params.popSex_HETM),2);
ann_PctEligOnPrEP_HRHET_F=sum(Outcome.ann_NumberOnPrEP_HighRiskHETs_F,2) ./sum(Outcome.ann_NumberEligForPrEP * Params.popSexIndicator(:,Params.popSex_HETF),2);
IncidenceMSM=sum(Outcome.ann_NewInfectionsMSM,2)
IncidencePWID=sum(Outcome.ann_NewInfectionsIDU,2)
IncidenceHETM=sum(Outcome.ann_NewInfectionsHETM,2)
IncidenceHETF=sum(Outcome.ann_NewInfectionsHETF,2)
IncidenceTotal=sum(Outcome.ann_TotalNewInfections,2)


OnPrEP_HighRiskMSM_Blk = sum(Outcome.ann_NumberOnPrEP_HighRiskMSM .* Params.raceIndicator(:,Params.race_B)',2);
OnPrEP_HighRiskMSM_Hisp =sum( Outcome.ann_NumberOnPrEP_HighRiskMSM .* Params.raceIndicator(:,Params.race_H)',2);
OnPrEP_HighRiskMSM_Oth = sum(Outcome.ann_NumberOnPrEP_HighRiskMSM .* Params.raceIndicator(:,Params.race_O)',2);
OnPrEP_TotalMSM=sum(Outcome.ann_NumberOnPrEP_HighRiskMSM,2);
OnPrEP_HighRiskHETs_F_Blk = sum(Outcome.ann_NumberOnPrEP_HighRiskHETs_F .* Params.raceIndicator(:,Params.race_B)',2);
OnPrEP_HighRiskHETs_F_Hisp = sum(Outcome.ann_NumberOnPrEP_HighRiskHETs_F .* Params.raceIndicator(:,Params.race_H)',2);
OnPrEP_HighRiskHETs_F_Oth = sum(Outcome.ann_NumberOnPrEP_HighRiskHETs_F .* Params.raceIndicator(:,Params.race_O)',2);
OnPrEP_TotalHETF=sum(Outcome.ann_NumberOnPrEP_HighRiskHETs_F,2);
OnPrEP_HighRiskHETs_M_Blk = sum(Outcome.ann_NumberOnPrEP_HighRiskHETs_M .* Params.raceIndicator(:,Params.race_B)',2);
OnPrEP_HighRiskHETs_M_Hisp = sum(Outcome.ann_NumberOnPrEP_HighRiskHETs_M .* Params.raceIndicator(:,Params.race_H)',2);
OnPrEP_HighRiskHETs_M_Oth = sum(Outcome.ann_NumberOnPrEP_HighRiskHETs_M .* Params.raceIndicator(:,Params.race_O)',2);
OnPrEP_TotalHETM=sum(Outcome.ann_NumberOnPrEP_HighRiskHETs_M,2);
OnPrEP_HighRiskIDUs_Blk = sum(Outcome.ann_NumberOnPrEP_HighRiskIDUs .* Params.raceIndicator(:,Params.race_B)',2);
OnPrEP_HighRiskIDUs_Hisp = sum(Outcome.ann_NumberOnPrEP_HighRiskIDUs .* Params.raceIndicator(:,Params.race_H)',2);
OnPrEP_HighRiskIDUs_Oth = sum(Outcome.ann_NumberOnPrEP_HighRiskIDUs .* Params.raceIndicator(:,Params.race_O)',2);
OnPrEP_TotalPWID=sum(Outcome.ann_NumberOnPrEP_HighRiskIDUs,2);



index2023=14;
index2027=18;
index2032=23;
index2042=33;
