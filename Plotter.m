% Purpose: Plot selected outcomes over time from Outcome object.
%[Outcome, p] = HIVEpiModel('HOPE_Model_AW_Testing');
function Plotter(Outcome)
    %% Plotter.m
    % Purpose: Plot selected outcomes over time from Outcome object created by HIVEpiModel.m
    % Written by Adam Walter awalter@rti.org
    % Last Updated: 09/16/2019
    
    % Current Set of Outputs:
    % 1. HIV Prevalence by Age Group
    % 2. HIV Prevalence by Race/Eth Group
    % 3. HIV Prevalence by Transmission Risk Group
    % 4. HIV Incidence by Age
    % 5. HIV Incidence by Race/Eth Group
    % 6. HIV Incidence by Transmission Risk Group
    % 7. Percent (of those with HIV) Diagnosed by Age Group (All)
    % 8. Percent (of those with HIV) Diagnosed by Transmission Risk Group
    % 9. Percent (of those with HIV) Diagnosed by Race/Eth Group
    % 10. Percent VLS among Diagnosed by Age Group
    % 11. Percent VLS along Diagnosed by Race/Eth Group
    % 12. Percent VLS Among Diagnosed by Transmission Risk Group
    % 13. Distribution of PLWH Between Compartments at End of Time Period 1
    % 14. Distribution of PLWH Between Compartments at End of Time Period 2
    % 15. Distribution of PLWH Between Compartments at End of Time Period 3
    
    
    %% 0.0 Initialize Years
    yrs = Outcome.set_FirstOutcomeYr:1:Outcome.set_LastOutcomeYr;
    yrs = yrs(:);
    
    %% 1. HIV PREVALENCE BY AGE GROUP
    
    % 1.1 Selected Outcomes:
    % ann_HIVPrevalence_13_17 
    % ann_HIVPrevalence_18_24 
    % ann_HIVPrevalence_25_34 
    % ann_HIVPrevalence_35_44 
    % ann_HIVPrevalence_45_54 
    % ann_HIVPrevalence_55_64 
    % ann_HIVPrevalence_65

    % 1.2 Calculate totals within each age group by sum across rows
    ann_HIVPrevalence_13_17_tot = sum(Outcome.ann_HIVPrevalence_13_17, 2);
    ann_HIVPrevalence_18_24_tot = sum(Outcome.ann_HIVPrevalence_18_24, 2);
    ann_HIVPrevalence_25_34_tot = sum(Outcome.ann_HIVPrevalence_25_34, 2);
    ann_HIVPrevalence_35_44_tot = sum(Outcome.ann_HIVPrevalence_35_44, 2);
    ann_HIVPrevalence_45_54_tot = sum(Outcome.ann_HIVPrevalence_45_54, 2);
    ann_HIVPrevalence_55_64_tot = sum(Outcome.ann_HIVPrevalence_55_64, 2);
    ann_HIVPrevalence_65_tot = sum(Outcome.ann_HIVPrevalence_65, 2);

    % 1.3 Create Subplot 1 within Figure 1
    figone = figure('visible','off');
    subplot(4,1,1);
    plot(yrs, ann_HIVPrevalence_13_17_tot,'DisplayName','13 to 17');
    hold on;
    plot(yrs, ann_HIVPrevalence_18_24_tot,'DisplayName','18 to 24');
    plot(yrs, ann_HIVPrevalence_25_34_tot,'DisplayName','25 to 34');
    plot(yrs, ann_HIVPrevalence_35_44_tot,'DisplayName','35 to 44');
    plot(yrs, ann_HIVPrevalence_45_54_tot,'DisplayName','45 to 54');
    plot(yrs, ann_HIVPrevalence_55_64_tot,'DisplayName','55 to 64');
    plot(yrs, ann_HIVPrevalence_65_tot,'DisplayName','65+');
    legend('Location', 'eastoutside', 'Orientation', 'vertical');
    title('HIV Prevalence by Age Group');
    hold off;

    %% 2. HIV PREVALENCE BY RACE/ETH GROUP
    
    % 2.1 Selected Outcomes:
    % ann_HIVPrevalence_Blk 
    % ann_HIVPrevalence_Hisp 
    % ann_HIVPrevalence_Oth

    % 2.2 Calculate totals within each race/eth group by sum across rows
    ann_HIVPrevalence_Blk_tot = sum(Outcome.ann_HIVPrevalence_Blk, 2);
    ann_HIVPrevalence_Hisp_tot = sum(Outcome.ann_HIVPrevalence_Hisp, 2);
    ann_HIVPrevalence_Oth_tot = sum(Outcome.ann_HIVPrevalence_Oth, 2);

    % 2.3 Create Subplot 2 within Figure 1
    subplot(4,1,2);
    plot(yrs, ann_HIVPrevalence_Blk_tot,'DisplayName','Black');
    hold on;
    plot(yrs, ann_HIVPrevalence_Hisp_tot,'DisplayName','Hispanic');
    plot(yrs, ann_HIVPrevalence_Oth_tot,'DisplayName','Other');
    legend('Location', 'eastoutside', 'Orientation', 'vertical');
    title('HIV Prevalence by Race/Eth Group');
    hold off;

    %% 3. HIV PREVALENCE BY TRANSMISSION RISK GROUP
    
    % 3.1 Selected Outcomes:
    % ann_HIVPrevalence_HET
    % ann_HIVPrevalence_MSM
    % ann_HIVPrevalence_IDU

    % 3.2 Calculate totals within each risk group by sum across rows
    ann_HIVPrevalence_HET_tot = sum(Outcome.ann_HIVPrevalence_HET, 2);
    ann_HIVPrevalence_MSM_tot = sum(Outcome.ann_HIVPrevalence_MSM, 2);
    ann_HIVPrevalence_IDU_tot = sum(Outcome.ann_HIVPrevalence_IDU, 2);

    % 3.3 Create Subplot 3 within Figure 1
    subplot(4,1,3);
    plot(yrs, ann_HIVPrevalence_HET_tot,'DisplayName','HET');
    hold on;
    plot(yrs, ann_HIVPrevalence_MSM_tot,'DisplayName','MSM');
    plot(yrs, ann_HIVPrevalence_IDU_tot,'DisplayName','IDU');
    legend('Location', 'eastoutside', 'Orientation', 'vertical');
    title('HIV Prevalence by Transmission Risk Group');
    hold off;

    %% 4. HIV INCIDENCE BY AGE GROUP

    % 4.1 Selected Outcomes:
    % ann_NewInfections_13_17
    % ann_NewInfections_18_24
    % ann_NewInfections_25_34
    % ann_NewInfections_35_44
    % ann_NewInfections_45_54
    % ann_NewInfections_55_64

    % 4.2 Calculate totals within each age group by sum across rows
    ann_NewInfections_13_17_tot = sum(Outcome.ann_NewInfections_13_17, 2);
    ann_NewInfections_18_24_tot = sum(Outcome.ann_NewInfections_18_24, 2);
    ann_NewInfections_25_34_tot = sum(Outcome.ann_NewInfections_25_34, 2);
    ann_NewInfections_35_44_tot = sum(Outcome.ann_NewInfections_35_44, 2);
    ann_NewInfections_45_54_tot = sum(Outcome.ann_NewInfections_45_54, 2);
    ann_NewInfections_55_64_tot = sum(Outcome.ann_NewInfections_55_64, 2);

    % 4.3 Create Subplot 4 within Figure 1 (Last subplot on first PDF / Figure 1)
    subplot(4,1,4);
    plot(yrs, ann_NewInfections_13_17_tot,'DisplayName','13 to 17');
    hold on;
    plot(yrs, ann_NewInfections_18_24_tot,'DisplayName','18 to 24');
    plot(yrs, ann_NewInfections_25_34_tot,'DisplayName','25 to 34');
    plot(yrs, ann_NewInfections_35_44_tot,'DisplayName','35 to 44');
    plot(yrs, ann_NewInfections_45_54_tot,'DisplayName','45 to 54');
    plot(yrs, ann_NewInfections_55_64_tot,'DisplayName','55 to 64');
    ax = gca;
    ax.YRuler.Exponent = 0;
    legend('Location', 'eastoutside', 'Orientation', 'vertical');
    title('HIV Incidence by Age Group');
    hold off;

    %% 5. HIV INCIDENCE BY RACE/ETH GROUP

    % 5.1 Selected Outcomes:
    % ann_NewInfections_Blk
    % ann_NewInfections_Hisp
    % ann_NewInfections_Oth

    % 5.2 Calculate totals within each race/eth group by sum across rows
    ann_NewInfections_Blk_tot = sum(Outcome.ann_NewInfections_Blk, 2);
    ann_NewInfections_Hisp_tot = sum(Outcome.ann_NewInfections_Hisp, 2);
    ann_NewInfections_Oth_tot = sum(Outcome.ann_NewInfections_Oth, 2);

    % 5.3 Create Subplot 1 within Figure 2
    figtwo = figure('visible','off');
    subplot(4,1,1);
    plot(yrs, ann_NewInfections_Blk_tot,'DisplayName','Black');
    hold on;
    plot(yrs, ann_NewInfections_Hisp_tot,'DisplayName','Hispanic');
    plot(yrs, ann_NewInfections_Oth_tot,'DisplayName','Other');
    ax = gca;
    ax.YRuler.Exponent = 0;
    legend('Location', 'eastoutside', 'Orientation', 'vertical');
    title('HIV Incidence by Race/Eth Group');
    hold off;

    %% 6. HIV INCIDENCE BY TRANSMISSION RISK GROUP

    % 6.1 Selected Outcomes:
    % ann_NewInfectionsHET
    % ann_NewInfectionsMSM
    % ann_NewInfectionsIDU

    % 6.2 Calculate totals within each risk group by sum across rows
    ann_NewInfectionsHET_tot = sum(Outcome.ann_NewInfectionsHET, 2);
    ann_NewInfectionsMSM_tot = sum(Outcome.ann_NewInfectionsMSM, 2);
    ann_NewInfectionsIDU_tot = sum(Outcome.ann_NewInfectionsIDU, 2);

    % 6.3 Create Subplot 2 within Figure 2
    subplot(4,1,2);
    plot(yrs, ann_NewInfectionsHET_tot,'DisplayName','HET');
    hold on;
    plot(yrs, ann_NewInfectionsMSM_tot,'DisplayName','MSM');
    plot(yrs, ann_NewInfectionsIDU_tot,'DisplayName','IDU');
    ax = gca;
    ax.YRuler.Exponent = 0;
    legend('Location', 'eastoutside', 'Orientation', 'vertical');
    title('HIV Incidence by Transmission Risk Group');
    hold off;

    %% 7. PERCENT DIAGNOSED BY AGE GROUP

    % 7.1 Selected Outcomes:
    % ann_PctAware_13_17
    % ann_PctAware_18_24
    % ann_PctAware_25_34
    % ann_PctAware_35_44
    % ann_PctAware_45_54
    % ann_PctAware_55_64
    % ann_PctAware_65

    % 7.2 Pull percentages within each age group, they are vectors
    ann_PctAware_13_17_tot = Outcome.ann_PctAware_13_17*100;
    ann_PctAware_18_24_tot = Outcome.ann_PctAware_18_24*100;
    ann_PctAware_25_34_tot = Outcome.ann_PctAware_25_34*100;
    ann_PctAware_35_44_tot = Outcome.ann_PctAware_35_44*100;
    ann_PctAware_45_54_tot = Outcome.ann_PctAware_45_54*100;
    ann_PctAware_55_64_tot = Outcome.ann_PctAware_55_64*100;
    ann_PctAware_65_tot = Outcome.ann_PctAware_65*100;

    % 7.3 Create Subplot 3 within Figure 2
    subplot(4,1,3);
    plot(yrs, ann_PctAware_13_17_tot,'DisplayName','13 to 17');
    hold on;
    plot(yrs, ann_PctAware_18_24_tot,'DisplayName','18 to 24');
    plot(yrs, ann_PctAware_25_34_tot,'DisplayName','25 to 34');
    plot(yrs, ann_PctAware_35_44_tot,'DisplayName','35 to 44');
    plot(yrs, ann_PctAware_45_54_tot,'DisplayName','45 to 54');
    plot(yrs, ann_PctAware_55_64_tot,'DisplayName','55 to 64');
    plot(yrs, ann_PctAware_65_tot,'DisplayName','65+');
    ytickformat('percentage');
    legend('Location', 'eastoutside', 'Orientation', 'vertical');
    title('Percent Diagnosed by Age Group');
    hold off;

    %% 8. PERCENT DIAGNOSED BY TRANSMISSION RISK GROUP

    % 8.1 Selected Outcomes:
    % ann_PctAware_HRHET
    % ann_PctAware_OHET
    % ann_PctAware_HRMSM
    % ann_PctAware_OMSM
    % ann_PctAware_PWID

    % 8.2 Pull percentages within each age group, they are vectors
    ann_PctAware_HRHET_tot = Outcome.ann_PctAware_HRHET*100;
    ann_PctAware_OHET_tot = Outcome.ann_PctAware_OHET*100;
    ann_PctAware_HRMSM_tot = Outcome.ann_PctAware_HRMSM*100;
    ann_PctAware_OMSM_tot = Outcome.ann_PctAware_OMSM*100;
    ann_PctAware_PWID_tot = Outcome.ann_PctAware_PWID*100;

    % 8.3 Create Subplot 4 within Figure 2 (Last plot on second PDF / Figure 2)
    subplot(4,1,4);
    plot(yrs, ann_PctAware_HRHET_tot,'DisplayName','HRHET');
    hold on;
    plot(yrs, ann_PctAware_OHET_tot,'DisplayName','OHET');
    plot(yrs, ann_PctAware_HRMSM_tot,'DisplayName','HRMSM');
    plot(yrs, ann_PctAware_OMSM_tot,'DisplayName','OMSM');
    plot(yrs, ann_PctAware_PWID_tot,'DisplayName','PWID');
    ytickformat('percentage');
    legend('Location', 'eastoutside', 'Orientation', 'vertical');
    title('Percent Diagnosed by Transmission Risk Group');
    hold off;

    %% 9. PERCENT DIAGNOSED BY RACE/ETH GROUP
    
    % 9.1 Selected Outcomes:
    % ann_PctAware_Blk
    % ann_PctAware_Hisp
    % ann_PctAware_Oth

    % 9.2 Pull percentages within each race/eth group, they are vectors
    ann_PctAware_Blk_tot = Outcome.ann_PctAware_Blk*100;
    ann_PctAware_Hisp_tot = Outcome.ann_PctAware_Hisp*100;
    ann_PctAware_Oth_tot = Outcome.ann_PctAware_Oth*100;
    
    % 9.3 Create subplot 1 within Figure 3 (First plot on third PDF / Figure 3)
    figthree = figure('visible','off');
    subplot(4,1,1);
    plot(yrs, ann_PctAware_Blk_tot,'DisplayName','Black');
    hold on;
    plot(yrs, ann_PctAware_Hisp_tot,'DisplayName','Hispanic');
    plot(yrs, ann_PctAware_Oth_tot,'DisplayName','Other');
    ytickformat('percentage');
    legend('Location', 'eastoutside', 'Orientation', 'vertical');
    title('Percent Diagnosed by Race/Eth Group');
    hold off;
    
    %% 10. PERCENT VLS AMONG DIAGNOSED BY AGE GROUP
    
    % 10.1 Selected Outcomes
    % ann_PctVLSamongdiag_13_17
    % ann_PctVLSamongdiag_18_24
    % ann_PctVLSamongdiag_25_34
    % ann_PctVLSamongdiag_35_44
    % ann_PctVLSamongdiag_45_54
    % ann_PctVLSamongdiag_55_64
    % ann_PctVLSamongdiag_65

    % 10.2 Pull percentages within each age group, they are vectors
    ann_PctVLSamongdiag_13_17_tot = Outcome.ann_PctVLSamongdiag_13_17*100;
    ann_PctVLSamongdiag_18_24_tot = Outcome.ann_PctVLSamongdiag_18_24*100;
    ann_PctVLSamongdiag_25_34_tot = Outcome.ann_PctVLSamongdiag_25_34*100;
    ann_PctVLSamongdiag_35_44_tot = Outcome.ann_PctVLSamongdiag_35_44*100;
    ann_PctVLSamongdiag_45_54_tot = Outcome.ann_PctVLSamongdiag_45_54*100;
    ann_PctVLSamongdiag_55_64_tot = Outcome.ann_PctVLSamongdiag_55_64*100;
    ann_PctVLSamongdiag_65_tot = Outcome.ann_PctVLSamongdiag_65*100;

    % 10.3 Create Subplot 2 within Figure 3
    subplot(4,1,2);
    plot(yrs, ann_PctVLSamongdiag_13_17_tot,'DisplayName','13 to 17');
    hold on;
    plot(yrs, ann_PctVLSamongdiag_18_24_tot,'DisplayName','18 to 24');
    plot(yrs, ann_PctVLSamongdiag_25_34_tot,'DisplayName','25 to 34');
    plot(yrs, ann_PctVLSamongdiag_35_44_tot,'DisplayName','35 to 44');
    plot(yrs, ann_PctVLSamongdiag_45_54_tot,'DisplayName','45 to 54');
    plot(yrs, ann_PctVLSamongdiag_55_64_tot,'DisplayName','55 to 64');
    plot(yrs, ann_PctVLSamongdiag_65_tot,'DisplayName','65+');
    ytickformat('percentage');
    legend('Location', 'eastoutside', 'Orientation', 'vertical');
    title('Percent VLS among Diagnosed by Age Group');
    hold off;
    
    %% 11. PERCENT VLS AMONG DIAGNOSED BY RACE/ETH GROUP
    
    % 11.1 Selected Outcomes:
    % ann_PctVLSamongdiag_Blk
    % ann_PctVLSamongdiag_Hisp
    % ann_PctVLSamongdiag_Oth

    % 11.2 Pull percentages within each race/eth group, they are vectors
    ann_PctVLSamongdiag_Blk_tot = Outcome.ann_PctVLSamongdiag_Blk*100;
    ann_PctVLSamongdiag_Hisp_tot = Outcome.ann_PctVLSamongdiag_Hisp*100;
    ann_ann_PctVLSamongdiag_Oth_tot = Outcome.ann_PctVLSamongdiag_Oth*100;
    
    % 11.3 Create Subplot 3 within Figure 3
    subplot(4,1,3);
    plot(yrs, ann_PctVLSamongdiag_Blk_tot,'DisplayName','Black');
    hold on;
    plot(yrs, ann_PctVLSamongdiag_Hisp_tot,'DisplayName','Hispanic');
    plot(yrs, ann_ann_PctVLSamongdiag_Oth_tot,'DisplayName','Other');
    ytickformat('percentage');
    legend('Location', 'eastoutside', 'Orientation', 'vertical');
    title('Percent VLS among Diagnosed by Race/Eth Group');
    hold off;
    
    %% 12. PERCENT VLS AMONG DIAGNOSED BY TRANSMISSION RISK GROUP
    
    % 12.1 Selected Outcomes:
    % ann_PctVLSamongdiag_HET
    % ann_PctVLSamongdiag_MSM
    % ann_PctVLSamongdiag_IDU

    % 12.2 Pull percentages within each age group, they are vectors
    ann_PctVLSamongdiag_HET_tot = Outcome.ann_PctVLSamongdiag_HET*100;
    ann_PctVLSamongdiag_MSM_tot = Outcome.ann_PctVLSamongdiag_MSM*100;
    ann_PctVLSamongdiag_IDU_tot = Outcome.ann_PctVLSamongdiag_IDU*100;

    % 12.3 Create Subplot 4 within Figure 3 (Last plot on third PDF / Figure 3)
    subplot(4,1,4);
    plot(yrs, ann_PctVLSamongdiag_HET_tot,'DisplayName','HET');
    hold on;
    plot(yrs, ann_PctVLSamongdiag_MSM_tot,'DisplayName','MSM');
    plot(yrs, ann_PctVLSamongdiag_IDU_tot,'DisplayName','IDU');
    ytickformat('percentage');
    legend('Location', 'eastoutside', 'Orientation', 'vertical');
    title('Percent VLS among Diagnosed by Transmission Risk Group');
    hold off;

    %% 0.1 SETUP FOR TABLES
    
    % The goal is to make the tables plot dynamically adjust to the number
    % of variables in the plot. 
    
    numplots = 0; 
    if isfield(Outcome, 'DistCompartsLastYrPeriod1')
        numplots = numplots + 1;
    end
    if isfield(Outcome, 'DistCompartsLastYrPeriod2')
        numplots = numplots + 1;
    end
    if isfield(Outcome, 'DistCompartsLastYrPeriod3')
        numplots = numplots + 1;
    end
    if isfield(Outcome, 'DistCompartsLastYrPeriod4')
        numplots = numplots + 1;
    end
    if isfield(Outcome, 'DistCompartsLastYr')
        numplots = numplots + 1;
    end
    numplots = max(2, numplots);
    figfour = figure('visible','off');
    plotposition = 1;
    %% 13. DISTRIBUTION OF PLWH BETWEEN COMPARTMENTS AT END OF TIME PERIOD 1
    
    % 13.1 Selected Outcome:
    % DistCompartsLastYrPeriod1
    
    if isfield(Outcome, 'DistCompartsLastYrPeriod1')
        % 13.2 Convert to Percentages
        eotP1 = Outcome.DistCompartsLastYrPeriod1*100;

        % 13.3 Create Heat Map of Matrix Data
        % Make Title 
        y1 = min(Outcome.set_PeriodTwoStartYear-1, Outcome.set_ModelStartYear+Outcome.set_TimeHorizon-1);
        p1Title = sprintf('Distribution of PLWH Between Compartments at End of Time Period 1 (%d)', y1);
        subplot(numplots, 1, plotposition);
        heatxvalues = {'Acute', 'CD4 > 500', 'CD4 350-500', 'CD4 200-350', 'CD4 < 200'};
        heatyvalues = {'Not Diagnosed', 'Diagnosed', 'Linked to Care', 'On ART not VLS', 'VLS'};
        h1 = heatmap(heatxvalues, heatyvalues, eotP1);
        h1.Title = p1Title;
        h1.XLabel = 'Disease Stages';
        h1.YLabel = 'Continuum Stages';
        h1.CellLabelFormat = '%.3f%%';
        plotposition = plotposition + 1;
    end
    
    %% 14. DISTRIBUTION OF PLWH BETWEEN COMPARTMENTS AT END OF TIME PERIOD 2
    
    % 14.1 Selected Outcome:
    % DistCompartsLastYrPeriod2
    
    % Wrap in exist statement because it may not exist.
    if isfield(Outcome, 'DistCompartsLastYrPeriod2')
        % 14.2 Convert to Percentages
        eotP2 = Outcome.DistCompartsLastYrPeriod2*100;

        % 14.3 Create Heat Map of Matrix Data
        % Make Title
        y2 = min(Outcome.set_PeriodThreeStartYear-1, Outcome.set_ModelStartYear+Outcome.set_TimeHorizon-1);
        p2Title = sprintf('Distribution of PLWH Between Compartments at End of Time Period 2 (%d)', y2);
        subplot(numplots, 1, plotposition);
        heatxvalues = {'Acute', 'CD4 > 500', 'CD4 350-500', 'CD4 200-350', 'CD4 < 200'};
        heatyvalues = {'Not Diagnosed', 'Diagnosed', 'Linked to Care', 'On ART not VLS', 'VLS'};
        h2 = heatmap(heatxvalues, heatyvalues, eotP2);
        h2.Title = p2Title;
        h2.XLabel = 'Disease Stages';
        h2.YLabel = 'Continuum Stages';
        h2.CellLabelFormat = '%.3f%%';
        plotposition = plotposition + 1;
    end
    
    %% 14. DISTRIBUTION OF PLWH BETWEEN COMPARTMENTS AT END OF TIME PERIOD 3
    
    % 14.1 Selected Outcome:
    % DistCompartsLastYrPeriod3
    
    % Wrap in exist statement because it may not exist.
    if isfield(Outcome, 'DistCompartsLastYrPeriod3')
        % 14.2 Convert to Percentages
        eotP3 = Outcome.DistCompartsLastYrPeriod3*100;

        % 14.3 Create Heat Map of Matrix Data
        % Make Title
        y3 = min(Outcome.set_PeriodFourStartYear-1, Outcome.set_ModelStartYear+Outcome.set_TimeHorizon-1);
        p3Title = sprintf('Distribution of PLWH Between Compartments at End of Time Period 3 (%d)', y3);
        subplot(numplots, 1, plotposition);
        heatxvalues = {'Acute', 'CD4 > 500', 'CD4 350-500', 'CD4 200-350', 'CD4 < 200'};
        heatyvalues = {'Not Diagnosed', 'Diagnosed', 'Linked to Care', 'On ART not VLS', 'VLS'};
        h3 = heatmap(heatxvalues, heatyvalues, eotP3);
        h3.Title = p3Title;
        h3.XLabel = 'Disease Stages';
        h3.YLabel = 'Continuum Stages';
        h3.CellLabelFormat = '%.3f%%';
        plotposition = plotposition + 1;
    end
    
    %% 14. DISTRIBUTION OF PLWH BETWEEN COMPARTMENTS AT END OF TIME PERIOD 4
    
    % 14.1 Selected Outcome:
    % DistCompartsLastYrPeriod4
    
    % Wrap in exist statement because it may not exist.
    if isfield(Outcome, 'DistCompartsLastYrPeriod4')
        % 14.2 Convert to Percentages
        eotP4 = Outcome.DistCompartsLastYrPeriod4*100;

        % 14.3 Create Heat Map of Matrix Data
        % Make Title
        y4 = min(Outcome.set_PeriodFiveStartYear-1, Outcome.set_ModelStartYear+Outcome.set_TimeHorizon-1);
        p4Title = sprintf('Distribution of PLWH Between Compartments at End of Time Period 4 (%d)', y4);
        subplot(numplots, 1, plotposition);
        heatxvalues = {'Acute', 'CD4 > 500', 'CD4 350-500', 'CD4 200-350', 'CD4 < 200'};
        heatyvalues = {'Not Diagnosed', 'Diagnosed', 'Linked to Care', 'On ART not VLS', 'VLS'};
        h4 = heatmap(heatxvalues, heatyvalues, eotP4);
        h4.Title = p4Title;
        h4.XLabel = 'Disease Stages';
        h4.YLabel = 'Continuum Stages';
        h4.CellLabelFormat = '%.3f%%';
        plotposition = plotposition + 1;
    end
    
    %% 15. DISTRIBUTION OF PLWH BETWEEN COMPARTMENTS IN LAST YEAR
    
    % 15.1 Selected Outcome:
    % DistCompartsLastYr
    
    % 15.2 Convert to Percentages
    if isfield(Outcome, 'DistCompartsLastYr')
        eotP5 = Outcome.DistCompartsLastYr*100;

        % 15.3 Create Heat Map of Matrix Data
        % Make Title
        p5Title = sprintf('Distribution of PLWH Between Compartments In Last Model Year (%d)', Outcome.set_LastOutcomeYr);
        subplot(numplots, 1, plotposition);
        heatxvalues = {'Acute', 'CD4 > 500', 'CD4 350-500', 'CD4 200-350', 'CD4 < 200'};
        heatyvalues = {'Not Diagnosed', 'Diagnosed', 'Linked to Care', 'On ART not VLS', 'VLS'};
        h5 = heatmap(heatxvalues, heatyvalues, eotP5);
        h5.Title = p5Title;
        h5.XLabel = 'Disease Stages';
        h5.YLabel = 'Continuum Stages';
        h5.CellLabelFormat = '%.3f%%';    
    end
    
    %% 16. SAVE FIGS AS PDF
    filename1 = sprintf('PlotsA_%s', datestr(now,'mm-dd-yy_HHMMPM'));
    filename2 = sprintf('PlotsB_%s', datestr(now,'mm-dd-yy_HHMMPM'));
    filename3 = sprintf('PlotsC_%s', datestr(now,'mm-dd-yy_HHMMPM'));
    filename4 = sprintf('PlotsD_%s', datestr(now,'mm-dd-yy_HHMMPM'));
    print(figone, filename1, '-dpdf', '-fillpage');
    print(figtwo, filename2, '-dpdf', '-fillpage');
    print(figthree, filename3, '-dpdf', '-fillpage');
    print(figfour, filename4, '-dpdf', '-fillpage');

end


































