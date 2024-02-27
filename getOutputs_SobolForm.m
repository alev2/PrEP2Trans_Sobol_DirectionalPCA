

hopePathSim='./HOPE Model V10_04_BaseCase.xlsm';
outTest=HIVEpiModel(hopePathSim);



ann_NewInfections_Blk=sum(outTest.ann_NewInfections_Blk,2);
ann_NewInfections_Hisp=sum(outTest.ann_NewInfections_Hisp,2);
ann_NewInfections_Oth=sum(outTest.ann_NewInfections_Oth,2);

ann_TotalNewInfections_PerPerson_Blk=ann_NewInfections_Blk./outTest.ann_popSize_Blk;
ann_TotalNewInfections_PerPerson_Hisp=ann_NewInfections_Hisp./outTest.ann_popSize_Hisp;
ann_TotalNewInfections_PerPerson_Oth=ann_NewInfections_Oth./outTest.ann_popSize_Oth;

irr_Blk=ann_TotalNewInfections_PerPerson_Blk./ann_TotalNewInfections_PerPerson_Oth;
irr_Hisp=ann_TotalNewInfections_PerPerson_Hisp./ann_TotalNewInfections_PerPerson_Oth;

spending=outTest.ann_ARTCarePrEPTransitionAndSEPCost_Disc;

outputMatrix=[ann_NewInfections_Blk ann_NewInfections_Oth ann_NewInfections_Hisp irr_Blk irr_Hisp spending];