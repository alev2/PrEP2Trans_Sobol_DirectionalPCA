function [mortRateVLS,diagRateChronic,mortalityRates,diagRates, continuumPct, initialCond, survIndicators] = obtain_MortalityDiagnosisRates(stateName,diagMortParams)

    %this is a helper function designed to define state-based mortality and
    %diagnosis rates. not it's not really necessary, but it can be useful.
    
    fileInput=strcat('./Data/',stateName);
    
    %read key data points
    load(strcat(fileInput,'/Diagnoses.mat'));
    load(strcat(fileInput,'/Prevalence.mat'));
    load(strcat(fileInput,'/PWHDeaths.mat'));
    load(strcat(fileInput,'/PctVLS.mat'));
    load(strcat(fileInput,'/PctLTC.mat'));
    load(strcat(fileInput,'/PctAware.mat'));
    
    %get continuum percentages
    pctAware=mean(PctAware);
    pctVLS=mean(PctVLS);
    pctLTC=mean(PctLTC);
    pctTNoVLS=pctLTC-pctVLS;
    pctAwareNoT=1-pctLTC;
    
    %diagnosis rates
    totalDiagRate=mean(Diagnoses./((1-PctAware).*Prevalence));    
    totalMortRate=mean(PWHDeaths./((PctAware).*Prevalence));

    %mortality multipliers
    mu=diagMortParams.mu;
	m_a=diagMortParams.m_a;
    m_u=diagMortParams.m_u;
	m_nt=diagMortParams.m_nt;
    m_t=diagMortParams.m_t;
    m_v=diagMortParams.m_v;
    m_p=diagMortParams.m_p;
    %arrange multipliers in vector
    mortalityMultipliers=[m_a;m_u;m_nt;m_t;m_v;m_p];

    
    %percetnage acute and AIDS, overall and among undiag
    pctAIDSOverall=diagMortParams.pctAIDSOverall;
    pctAIDSamongUndiag=diagMortParams.pctAIDSamongUndiag;
    pctAcuteAmongUndiag=diagMortParams.pctAcuteAmongUndiag;
    
    %diagnosis probability multipliers
    diagMultA=diagMortParams.diagMultA;
    diagMultP=diagMortParams.diagMultP;
    
    %full continuum of care pcts
    continuumPct=[(1-pctAware)*pctAcuteAmongUndiag;
                  (1-pctAware)*(1-pctAcuteAmongUndiag-pctAIDSamongUndiag);
                  pctAware*pctAwareNoT;
                  pctAware*pctTNoVLS;
                  pctAware*pctVLS];

	%AIDS should be the rest
    continuumPct=[continuumPct;(1-sum(continuumPct))];
    %initial conditions for PWH
    initialCond=round(mean(Prevalence)*continuumPct);

    %mortality rate from surveillance does not factor in all unaware deaths, so we
    %so we adjust here
    continuumPctMortRate=continuumPct(3:end)/sum(continuumPct(3:end));
    
    %here we get the acute/VLS mortality rate
    mortRateVLS=totalMortRate/(mortalityMultipliers(3:end)'*continuumPctMortRate);    
    %plugging back in with multipliers gives us all mortality rate
    mortalityRates=mortRateVLS*mortalityMultipliers;
    mortalityRates=[mu;mortalityRates];

    %similar procedure for diagnosis as for mortality
    diagMultipliers=[diagMultA; 1.0 ; diagMultP];
    undiagPcts=[pctAcuteAmongUndiag; 1-pctAcuteAmongUndiag-pctAIDSamongUndiag; pctAIDSamongUndiag];

    diagRateChronic=totalDiagRate/(diagMultipliers'*undiagPcts);    
    diagRates=diagRateChronic*diagMultipliers;

    
end

