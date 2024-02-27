function [incRates, chronicIncRate, totalIncRate] = obtain_IncidenceRates(stateName,incidence_Params)

    %this is a helper function designed to define state-based incidence
    %rates. not it's not really necessary.

    fileInput=strcat('./Data/',stateName);
    
    %read key data points
    load(strcat(fileInput,'/Diagnoses.mat'));
    load(strcat(fileInput,'/Prevalence.mat'));
    load(strcat(fileInput,'/PWHDeaths.mat'));
    load(strcat(fileInput,'/PctVLS.mat'));
    load(strcat(fileInput,'/PctLTC.mat'));
    load(strcat(fileInput,'/PctAware.mat'));
    load(strcat(fileInput,'/Incidence.mat'));
    
    %get continuum percentages
    pctAware=mean(PctAware);
    pctVLS=mean(PctVLS);
    pctLTC=mean(PctLTC);
    pctTNoVLS=pctLTC-pctVLS;
    pctAwareNoT=1-pctLTC;
    
    %total infectivity rate
    totalIncRate=mean(Incidence./Prevalence);
    
    %infectivity multipliers (throw out susceptible)
    Psi=incidence_Params.Psi(2:end);
        
    %percentage acute and AIDS, overall and among undiag
    pctAIDSamongUndiag=incidence_Params.pctAIDSamongUndiag;
    pctAcuteAmongUndiag=incidence_Params.pctAcuteAmongUndiag;
    
    
    %full continuum of care pcts
    continuumPct=[(1-pctAware)*pctAcuteAmongUndiag;
                  (1-pctAware)*(1-pctAcuteAmongUndiag-pctAIDSamongUndiag);
                  pctAware*pctAwareNoT;
                  pctAware*pctTNoVLS;
                  pctAware*pctVLS];

	%AIDS should be the rest
    continuumPct=[continuumPct;(1-sum(continuumPct))];
    
    %here we get the chronic infection rate 
    chronicIncRate=totalIncRate/(Psi'*continuumPct);    

    %plugging back in with multipliers gives us all mortality rate
    incRates=chronicIncRate*Psi;
    
    %augment for the susceptible pool for ease of vector operations (no infectivity)
    incRates=[0;incRates];
    
end
