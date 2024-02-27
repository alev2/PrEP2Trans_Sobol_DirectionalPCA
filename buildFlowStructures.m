%Define progressions with topological structure
diagRates_Topology=...
    [(1-F_a)*diagRates(1) 2 4;
     (1-F_u)*diagRates(2) 3 4;
     F_a*diagRates(1) 2 5;
     F_u*diagRates(2) 3 5];

diseaseProgression_Topology=...
    [sigma_a 2 3;
     kappa   3 7;
     kappa   4 7];

treatmentProgression_Topology=...
    [xi 4 5;    
     eta 5 4;
     nu 5 6;
     omega 6 5;
     etaV 6 4;
     alpha 7 5;
    ];
    
flowStruct=struct(...
    'diagRates_Topology',diagRates_Topology,...
    'diseaseProgression_Topology',diseaseProgression_Topology,...
    'treatmentProgression_Topology',treatmentProgression_Topology,...
    'mortalityRates',mortalityRates);


mortalityTop=[mortalityRates (1:7)' [1;8*ones(6,1)]];

fullFlowTopo=[diagRates_Topology;...
              diseaseProgression_Topology;...
              treatmentProgression_Topology;...
              mortalityTop];

fullFlowTopo=[fullFlowTopo ones(size(fullFlowTopo,1))];
flowParams=fullFlowTopo(:,1);
fullFlowTopo=fullFlowTopo(:,2:end);
fullFOITopo=[1 2 1];
numStrats=1;
numFlowStages=8;
infectedStates=[(2:7)'];
absorbingStates=[8];
PsiK=Psi(2:end);
