%this defines some variables so we can do the sobol analysis quickly,
%avoiding r

% initialize paramaters that we do not vary for HIV problem, as well as
% the sparse grid library.
%basePath='C:\Users\xjm9\OneDrive - CDC\MATLAB\V10.05_HopeAndSobol\Detailed_Sobol_Example';
%modelName='HOPE Model V10_05_LB20230224_2_20231005_Fresh.xlsm';

%addpath(genpath('./sparse-grids-matlab-kit-main/'));

modelName='HopeModel_BaseCase.xlsm';
hopePath=strcat('.\',modelName);

sheetName='TTProgression';


yearRange=(2010:1:2042)';

[yearInds,~]=find(yearRange>=2023);

%%%%%TESTING (Risk ratio!)
cells1to2=...
    {
        'Q42',...black
        'Q43',...hiwp/latino
        'Q44'... white
    };

inds1to2=...
    [
        959-789; %B
        960-789; %HL
    ];

%%%%%LINKAGE
cells2to3_atDiag=...
    {
      'Q315',... black
      'Q316',... hisp/latino
      'Q317'...  white
    };

inds2to3_atDiag=...
    [
        1011-789; %B
        1012-789; %HL    
    ];

cells2to3_afterDiag=...
    {
      'Q373',... black
      'Q374',... hisp/latino
      'Q375'...  white
    };


inds2to3_afterDiag=...
    [
        1023-789;
        1024-789;
    ];

%dropping out of care
cells3to2=...
     {
      'Q412',... black
      'Q413',... hisp/latino
      'Q414'... white
     };

inds3to2=...
    [
        261; %B
        262; %HL
    ];


inds3to4=...
    [
       1188-789; %B
       1189-789; %HL
       1190-789;
    ];


%%%%%%VLS/ART
%dropping off ART
cells4to3=...
     {
      'Q444',... black
      'Q445',... hisp/latino
      'Q446'... white
     };

inds4to3=...
    [
        279;
        280;
        281;
    ];

%losing VLS
cells5to4=...
     {
      'Q479',... black
      'Q480',... hisp/latino
      'Q481'... white
     };

inds5to4=...
    [
        316;
        317;    
        318;
    ];

%becoming VLS
cells4to5=...
     {
      'Q639',... black
      'Q640',... hisp/latino
      'Q641'... white
     };

inds4to5=...
    [
        408;
        409;
    ];


%%%%PrEP
cellsPrEP={...
    'T836'...black
    'T837'...hisp/latino
    'T838'...white
    };

indsPrEPBlack=[
            1319-789;
            1322-789;
            1325-789;
            1328-789;
        ];


indsPrEPHisp=[
            1320-789;
            1323-789;
            1326-789;
            1329-789;
        ];


indsDiscRate=810-789;


indsPrepHETF2023=[
        1322-789;
        1323-789;
        1324-789;
    ];
indsPrepHETF2024=[
        1334-789;
        1335-789;
        1336-789;
    ];

indsPrepHETM2023=[
        1319-789;
        1320-789;
        1321-789;
    ];
indsPrepHETM2024=[
        1331-789;
        1332-789;
        1333-789;
    ];


indsPrepMSM2023=[
        1325-789;
        1326-789;
        1327-789;
    ];
indsPrepMSM2024=[
        1337-789;
        1338-789;
        1339-789;
    ];


indsPrepPWID2023=[
        1328-789;
        1329-789;
        1330-789;
    ];
indsPrepPWID2024=[
        1340-789;
        1341-789;
        1342-789;
    ];



ind1to2HET=956-789;
ind1to2MSM=957-789;
ind1to2PWID=958-789;

ind5to4HET=1117-789;
ind5to4MSM=1118-789;
ind5to4PWID=1119-789;

indsDiscRate=810-789;



inds_LoHETM_HETF=1762-789;
inds_LoHETM_PWID=1763-789;

inds_LoHETF_HETM=1764-789;
inds_LoHETF_MSM=1765-789;
inds_LoHETF_PWID=1766-789;

inds_HiHETM_HETF=1767-789;
inds_HiHETM_PWID=1768-789;

inds_HiHETF_HETM=1769-789;
inds_HiHETF_MSM=1770-789;
inds_HiHETF_PWID=1771-789;

inds_MSM_HETF=1772-789;
inds_MSM_MSM=1773-789;
inds_MSM_PWID=1774-789;


mixingInds=[...
    inds_LoHETM_HETF;
    inds_LoHETM_PWID;
    inds_LoHETF_HETM;
    inds_LoHETF_MSM;
    inds_LoHETF_PWID;
    inds_HiHETM_HETF;
    inds_HiHETM_PWID;
    inds_HiHETF_HETM;
    inds_HiHETF_MSM;
    inds_HiHETF_PWID;
    inds_MSM_HETF;
    inds_MSM_MSM;
    inds_MSM_PWID;   
];






