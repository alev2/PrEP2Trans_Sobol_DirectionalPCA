%%%This runs the SIRD solution

addpath('./sparse-grids-matlab-kit-main/')
addpath('./sparse-grids-matlab-kit-main/main/')
addpath('./sparse-grids-matlab-kit-main/src/')
addpath('./sparse-grids-matlab-kit-main/tools/')
addpath('./sparse-grids-matlab-kit-main/docs-examples/')
addpath('../');
%addpath('../Data/GA/');

lw=2;
fs=15;
posMat=[-500 0 800 500];

%load state
stateName='GA';
load(strcat('../Data/',stateName,'/TotalPopulation.mat'));
load(strcat('../Data/',stateName,'/survMort.mat'));
load(strcat('../Data/',stateName,'/survPrevRange.mat'));
load(strcat('../Data/',stateName,'/survIncRange.mat'));

%system dimension
sysDim=7;

%time variables
num_Yrs=10; first_Yr=2017; year_array=(first_Yr:1:(first_Yr+num_Yrs-1))';
t_0=0; t_end=num_Yrs; num_Yrs*365+1;

%for computing annual incidence
prev_Indices=(0:num_Yrs)*365+1;

%general mortality annual rate
mu=-log(1-.004);

%PWH mortality rate multiplier (VLS)
m_v=1.0;

%PWH mortality rate multipliers (added on top of VLS, for non-VLS cases)
m_a=1.0; m_u=2.3; m_nt=2.3; m_t=1.2; m_p=8;
mortalityMultipliers=[m_a;m_u;m_nt;m_t;m_v;m_p];

%length of acute phase
sigma_a=1/(60/365);

%diagnosis rate multiplier
diagMultA=1.0; diagMultP=4.0;

%initiation care at diag.
F_a=0.86; F_u=0.86;

%chronic to AIDS
kappa=1/(10.);

%rate of obtaining VLS
nu=-log(1-.4);

%rate of dropping off treatment
eta=-log(1-.378);

%rate of losing VLS
omega=-log(1-.128);

%rate of treatment initiation/re-engagement
xi=-log(1-.44);

%rate of treatment initiation for AIDS
alpha=-log(1-.92);

%prep parameters
prepFrac=.01;
prepEfficacy=.96;

% %contact rate
% numPartners=1.8;
% actsPerPartner=54.3;
% transPerAct=.0015;
% beta=-log(1-...
%            (1-(1-transPerAct)^(numPartners*actsPerPartner))...
%            );

%FOI multipliers
psi_a=2.6; psi_u=1.0; psi_nt=.5; psi_t=.2; psi_v=0.; psi_p=2.6;
Psi=[0; psi_a; psi_u; psi_nt; psi_t; psi_v; psi_p];

%necessary to compute diagnosis and mortality rates
pctAIDSOverall=.01;
pctAIDSamongUndiag=.0;
pctAcuteAmongUndiag=.04;

%struct for mult and diag rates
diagMortParams=struct(...
    'mu',mu,...
    'm_a',m_a,'m_u',m_u,'m_nt',m_nt,'m_t',m_t,'m_v',m_v,'m_p',m_p,...
    'diagMultA',diagMultA,'diagMultP',diagMultP,...
    'pctAIDSOverall',pctAIDSOverall,'pctAIDSamongUndiag',pctAIDSamongUndiag,'pctAcuteAmongUndiag',pctAcuteAmongUndiag);

incParams=struct(...
    'Psi',Psi,...
    'pctAIDSOverall',pctAIDSOverall,'pctAIDSamongUndiag',pctAIDSamongUndiag,'pctAcuteAmongUndiag',pctAcuteAmongUndiag);


%get the mortality and diagnosis rates based on state-specific data
[mu_VLS,gamma_u,mortalityRates,diagRates,...
    continuumPct, initialCond] = obtain_MortalityDiagnosisRates(stateName,diagMortParams);

%get incidence rates based on state-specific data
[incRates, beta, totalIncRate] = obtain_IncidenceRates(stateName,incParams);

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
     alpha 7 5;
    ];
    
flowStruct=struct(...
    'diagRates_Topology',diagRates_Topology,...
    'diseaseProgression_Topology',diseaseProgression_Topology,...
    'treatmentProgression_Topology',treatmentProgression_Topology,...
    'mortalityRates',mortalityRates);

%assemble fixed flow matrix (no force-of-infection)
% [V,diagMat,hivProgMat,treatProgMat,mortalityMat]=...
%     build_Flow_Matrix(sysDim,flowStruct);

%initial condition
y0=[TotalPopulation-sum(initialCond);initialCond;0];
%y0=[TotalPopulation-sum(initialCond);initialCond];

%FOI multipliers

%put everything in parameter struct
% problem_Parameters_Alternate=struct(...
%     'sysDim',sysDim,...
%     'V',V,...
%     'prepFrac',prepFrac,...
%     'prepEfficacy',prepEfficacy,...
%     'beta',beta,...
%     'mortalityRates',mortalityRates,...
%     'Psi',Psi);

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

hiv_Model=homogenousSystem(...
    fullFlowTopo,...
    fullFOITopo,...
    numStrats,...
    numFlowStages,...
    infectedStates,...
    absorbingStates,...
    flowParams,...
    beta*(1-prepFrac*prepEfficacy),...
    PsiK);

%[t,y]=ode45(@(t,y)HIV_Model_Alternate(t,y,problem_Parameters_Alternate),tVec,y0);
[t,y]=ode45(@(t,y)solve_EpidemicProblem(t,y,hiv_Model),tVec,y0);

%collect outputs
prevalence=y*[0;1;1;1;1;1;1;0];
totalDeaths=y(:,end);
prevChange=round(diff(prevalence(prev_Indices)));
annualDeaths=round(diff(totalDeaths(prev_Indices)));
annualIncidence=prevChange+annualDeaths;
unaware=(y(:,2)+y(:,3));
notInCare=y(:,4);
inCare=y(:,5)+y(:,6);
VLS=y(:,6);
AIDS=y(:,7);

colorMat=[...
    .000 .447 .741;
    .850 .325 .098;
    .929 .694 .125;
    .494 .184 .556;
    .466 .674 .188;
    .301 .745 .933;
    .635 .078 .184];
stages=...
    {'Susceptible','Acute','Chronic, unaware','Aware, no ART','ART, no VLS','VLS','AIDS'};
%generate some plots
subplot(2,2,1)
plot(t+first_Yr, prevalence ,':r','LineWidth',lw);
hold on
plot([survPrevRange(1,1) survPrevRange(1,1)] ,[survPrevRange(1,2) survPrevRange(1,3)] ,'--*b','LineWidth',lw);
xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Num. size', 'Interpreter','latex','FontSize',50);
legend(['Model'],['Surveillance'],'Interpreter','latex','FontSize',fs,'Location','northwest',...
    'AutoUpdate','off');
for k=2:size(survPrevRange,1)
    plot([survPrevRange(k,1) survPrevRange(k,1)] ,[survPrevRange(k,2) survPrevRange(k,3)] ,'--*b','LineWidth',lw);
end
title('Prevalence', 'Interpreter','latex','FontSize',50);
fs=15;
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;

subplot(2,2,2)
plot(t+first_Yr,1-unaware./prevalence,':r','LineWidth',lw);
hold on
plot(t+first_Yr,inCare./prevalence,'-.m','LineWidth',lw);
plot(t+first_Yr,VLS./prevalence,'--k','LineWidth',lw);
xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('\% [-]', 'Interpreter','latex','FontSize',50);
legend(['\% aware'],['\% in care'],['\% VLS'], 'Interpreter','latex','FontSize',fs,'Location','northwest');
title('Continuum of Care', 'Interpreter','latex','FontSize',50);
fs=15;
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;

subplot(2,2,3)
plot(year_array,annualIncidence,'--or','LineWidth',lw);
hold on
plot(year_array,annualDeaths,':sk','LineWidth',lw);
xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Persons / Year', 'Interpreter','latex','FontSize',50);
plot([survIncRange(1,1) survIncRange(1,1)] ,[survIncRange(1,2) survIncRange(1,3)] ,'--*b','LineWidth',lw);
plot(survMort(:,1) ,survMort(:,2),'*m','LineWidth',lw);

legend(['Annual incidence, model'],['Annual PWH deaths, model'],...
    ['Annual incidence, surveillance'],['Annual PWH deaths, surveillance'],...
    'Interpreter','latex','FontSize',fs,'Location','northwest',...
    'AutoUpdate','off');

for k=2:size(survPrevRange,1)
    plot([survIncRange(k,1) survIncRange(k,1)] ,[survIncRange(k,2) survIncRange(k,3)] ,'--*b','LineWidth',lw);
end

title('Incidence + Deaths', 'Interpreter','latex','FontSize',50);
fs=15;
set(gca,'FontSize',fs);
ax=gca;
ax.YAxis.Exponent=0;

subplot(2,2,4)
plot(t+first_Yr,y(:,2)./prevalence,':','Color',colorMat(1,:), 'LineWidth',lw);
hold on
for k=3:(size(y,2)-1)
    plot(t+first_Yr,y(:,k)./prevalence,':','Color',colorMat(k,:),'LineWidth',lw);
end

xlabel('Year', 'Interpreter','latex','FontSize',50);
ylabel('Pct PWH population [\%]', 'Interpreter','latex','FontSize',50);
%plot([survIncRange(1,1) survIncRange(1,1)] ,[survIncRange(1,2) survIncRange(1,3)] ,'--*b','LineWidth',lw);
legend(stages(2:end),...
    'Interpreter','latex','FontSize',fs,'Location','northwest',...
    'AutoUpdate','off');

title('Total PWH Pop. \%', 'Interpreter','latex','FontSize',50);
fs=15;
set(gca,'FontSize',fs);
set(gcf,'Position',posMat);
ax=gca;
ax.YAxis.Exponent=0;




