% This code generates all the model-related non-schematic figures in the
% manuscript: "Quantitative Models of Phage-Antibiotics Combination Therapy", Rodriguez-Gonzalez et al. biorxiv 633784 (2019).

clear
clc
close all

%% Main figures

% Fig. 2 Immunophage Model
% Caption: Dynamics of the immunophage therapy model against two different bacterial inoculums.
% We simulate the phage therapy model developed by (Leung & Weitz, JTB (2017)) against two infection settings.
% In the first infection setting (a), a phage-sensitive bacterial inoculum, BP (orange solid line), 
% is challenged with phage (blue dashed line) in an immunocompetence context. In the second scenario (b),
% antibiotic-sensitive bacteria, BA (green solid line), are challenged with phage in presence of an 
% active immune response (purple dashed line). The initial bacterial density and the initial phage dose are,
% B_0 = 7.4x10^7 CFU/g and P_0 = 7.4x10^8 PFU/g, respectively. The simulation run is 96
% hours with phage being administered 2 hours after the infection.

%Fig. 2A
figure(1)
run([pwd '/code_comb_therapy/immunocompetent/HM-R_model/phage_immune/RHM_WT.m'])
saveas(gcf, [pwd '/figures/time_series/Fig2A'], 'fig')
saveas(gcf,[pwd '/figures/time_series/Fig2A'],'epsc')
clear

%Fig. 2B
figure(2)
run([pwd '/code_comb_therapy/immunocompetent/HM-R_model/phage_immune/RHM_WT_BA.m'])
saveas(gcf, [pwd '/figures/time_series/Fig2B'], 'fig')
saveas(gcf,[pwd '/figures/time_series/Fig2B'],'epsc')
clear

% Fig. 3 Phage-Antibiotic combination therapy + Immune response
% Caption: Outcomes of the phage-antibiotic combination therapy model for  
% two different infection settings. We simulate the combined effects of phage and 
% antibiotics in an immunocompetent host infected with phage-sensitive bacteria (a), 
% BP (orange solid line). In (b), the host is infected with antibiotic-sensitive bacteria,
% BA (green solid line). The dynamics of the phages (blue dashed line) and innate immunity
% (purple dashed line) are shown for each infection setting. Initial bacterial density and phage 
% dose are, B_0 = 7.4x10^7 CFU/g and P_0 = 7.4x10^8 PFU/g, respectively.
% The simulation run is 96 hours (4 days). Antibiotic and phage are administered 2 hours after the 
% beginning infection. Ciprofloxacin is maintained at a constant concentration of 0.0350 ug/ml 
% during the simulation. The carrying capacity of the bacteria is K_C = 10^10 CFU/g.


% Fig. 3A
figure(3)
run([pwd '/code_comb_therapy/immunocompetent/HM-R_model/Phage_antibio_immune/RHM_WT.m'])
saveas(gcf, [pwd '/figures/time_series/Fig3A'], 'fig')
saveas(gcf,[pwd '/figures/time_series/Fig3A'],'epsc')
clear

%Fig. 3B
figure(4)
run([pwd '/code_comb_therapy/immunocompetent/HM-R_model/Phage_antibio_immune/RHM_WT_BA.m'])
saveas(gcf, [pwd '/figures/time_series/Fig3B'], 'fig')
saveas(gcf,[pwd '/figures/time_series/Fig3B'],'epsc')
clear



% Fig. 4 Phage-Antibiotic combination therapy dynamics without immune response
% Caption: Bacterial dynamics given joint exposure to phages and antibiotic. We simulate
% bacterial growth for 96 hours in exposure to phage (blue dashed line) and antibiotic (data not shown)
% added two hours after the beginning of the inoculation. The combination of phage and antibiotic is 
% tested against two different bacterial inoculum. The first inoculum consisted of exclusively phage-sensitive 
% bacteria (a-c), BP (orange solid line). The second inoculum consisted of antibiotic-sensitive bacteria (d-f),
% BA (green solid line). Additionally, we test three different models of phage infection, heterogeneous mixing (a, d),
% phage saturation (b,e) and linear infection models (c, f). The initial bacterial density and phage dose are,
% B_0 = 7.4x10^7 CFU/g and P_0 = 7.4x10^8 PFU/g, respectively. Ciprofloxacin is maintained at a 
% constant concentration of 2.5 x MIC (i.e., 0.0350 ug/ml) during the simulations.

% Fig. 4A
figure(5)
run([pwd '/code_comb_therapy/neutropenic/HM_model/RHM_neutropenia.m'])
saveas(gcf, [pwd '/figures/time_series/Fig4A'], 'fig')
saveas(gcf,[pwd '/figures/time_series/Fig4A'],'epsc')
clear

% Fig. 4B
figure(6)
run([pwd '/code_comb_therapy/neutropenic/PS_model/RPS_neutropenia.m'])
saveas(gcf, [pwd '/figures/time_series/Fig4B'], 'fig')
saveas(gcf,[pwd '/figures/time_series/Fig4B'],'epsc')
clear

% Fig. 4C
figure(7)
run([pwd '/code_comb_therapy/neutropenic/LI_model/LIR_neutropenia.m'])
saveas(gcf, [pwd '/figures/time_series/Fig4C'], 'fig')
saveas(gcf,[pwd '/figures/time_series/Fig4C'],'epsc')
clear

% Fig. 4D
figure(8)
run([pwd '/code_comb_therapy/neutropenic/HM_model/RHM_neutropenia_BA.m'])
saveas(gcf, [pwd '/figures/time_series/Fig4D'], 'fig')
saveas(gcf,[pwd '/figures/time_series/Fig4D'],'epsc')
clear

% Fig. 4E
figure(9)
run([pwd '/code_comb_therapy/neutropenic/PS_model/RPS_neutropenia_BA.m'])
saveas(gcf, [pwd '/figures/time_series/Fig4E'], 'fig')
saveas(gcf,[pwd '/figures/time_series/Fig4E'],'epsc')
clear

% Fig. 4F
figure(10)
run([pwd '/code_comb_therapy/neutropenic/LI_model/LIR_neutropenia_BA.m'])
saveas(gcf, [pwd '/figures/time_series/Fig4F'], 'fig')
saveas(gcf,[pwd '/figures/time_series/Fig4F'],'epsc')
clear

% Fig. 5 Robustness analyis of the phage-antibiotic combination therapy model
% Caption: Outcomes of the robustness analysis for different antimicrobial strategies.
% We simulate the exposure of bacteria to different antimicrobial strategies, such as antibiotic-only (a),
% antibiotic + innate immunity (b),  phage + antibiotic (c), and phage-antibiotic combination in presence of 
% innate immunity (d). The heatmaps show the bacterial density at 96 hours post infection. Colored regions represent
% bacteria persistence (e.g., orange areas ~ 10^9 CFU/g and yellow areas ~10^10 CFU/g) while 
% the white regions represent pathogen clearance. We vary the concentration of ciprofloxacin (MIC = 0.014 ug/ml), 
% ranging from 0.1 MIC (0.0014 ug/ml) to 10 MIC (0.14 ug/ml), and the bacterial composition of the inoculum 
% ranging from 100% phage-sensitive bacteria (0% BA) to 100% antibiotic-sensitive bacteria (100% BA). Initial 
% bacterial density and phage dose (c-d) are, B0 = 7.4x10^7 CFU/g and P0 = 7.4x10^8 PFU/g, respectively. 
% Phage and antibiotic are administered two hours after the beginning of the infection.

% Fig. 5A
figure(11)
run([pwd '/code_comb_therapy/robustness_analysis/antibiotic_therapeutic.m'])
saveas(gcf, [pwd '/figures/heatmaps/Fig5A'], 'fig')
saveas(gcf,[pwd '/figures/heatmaps/Fig5A'],'epsc')
clear

% Fig. 5B
figure(13)
run([pwd '/code_comb_therapy/robustness_analysis/antibiotic_immune.m'])
saveas(gcf, [pwd '/figures/heatmaps/Fig5B'], 'fig')
saveas(gcf,[pwd '/figures/heatmaps/Fig5B'],'epsc')
clear

% Fig. 5C
figure(12)
run([pwd '/code_comb_therapy/robustness_analysis/phage_antibio_comb.m'])
saveas(gcf, [pwd '/figures/heatmaps/Fig5C'], 'fig')
saveas(gcf,[pwd '/figures/heatmaps/Fig5C'],'epsc')
clear


% Fig. 5D
figure(14)
run([pwd '/code_comb_therapy/robustness_analysis/combination_immune.m'])
saveas(gcf, [pwd '/figures/heatmaps/Fig5D'], 'fig')
saveas(gcf,[pwd '/figures/heatmaps/Fig5D'],'epsc')
clear


%% Supplementary figures


% Fig S1. Bacterial dynamics given exposure to low levels of antibiotic.
% Caption: We simulate the effects of combination therapy
% plus innate immunity on inocula with non-trivial levels of BA. First, an inoculum composed of 95% BP and 5% BA is treated
% with phage and different levels of antibiotic, 0, 0.01, and 0.001 X MIC (a, b, and c respectively). The same treatment is
% applied for an inoculum composed of 80% BP and 20% BA (d, e, and f, respectively). Initial bacterial and phage density are,
% B0 = 7.4 x 10^7 CFU/g and P0 = 7.4 x 10^8 PFU/g. Phage and antibiotic are administered two hours after infection. The
% simulation run was 96 hours (4 days).

% Fig. S1 A-B-C-D-E-F
run([pwd '/code_comb_therapy/immunocompetent/HM-R_model/Phage_antibio_immune/FigS1_lowanti/RHM_WT.m'])
letters = {'A' 'B' 'C' 'D' 'E' 'F'};
for fig_num = 1:6,
    saveas(figure(fig_num),[pwd '/figures/time_series/FigS1' letters{fig_num}], 'fig')
    saveas(figure(fig_num),[pwd '/figures/time_series/FigS1' letters{fig_num}],'epsc')
end
clear

% Fig S2. Accounting for variations in the concentration of antibiotic,
% from sub-MIC to MIC levels.
% Caption: We choose a particular inoculum composition (red boxes, 50%BP and 50%BA) from our heatmap 
% and zoom in at the dynamics level. Bacterial dynamics correspond to different antibiotic levels: 
% 0.1, 0.5, and 1 X MIC levels (a, b, and c, respectively ). The colored areas on the heatmap indicate 
% bacterial presence while the white areas indicate infection clearance after 96 hrs of treatment. 
% Phage and antibiotic are administered two hours after infection. Initial bacterial and phage density, 
% B0 = 7.4 x10^7 CFU/g and P0 = 7.4 x10^8 PFU/g, respectively.

% Fig. S2 A-B-C
run([pwd '/code_comb_therapy/immunocompetent/HM-R_model/Phage_antibio_immune/FigS2_SubtoMIC_levels/RHM_WT.m'])
letters = {'C' 'B' 'A'};
for fig_num = 1:3,
    saveas(figure(fig_num),[pwd '/figures/time_series/FigS2' letters{fig_num}], 'fig')
    saveas(figure(fig_num),[pwd '/figures/time_series/FigS2' letters{fig_num}],'epsc')
end
clear

% Fig. S2 A+P+I(Heatmap)
figure(23)
run([pwd '/code_comb_therapy/robustness_analysis/Figure_s2_tripartiteffect.m'])
saveas(gcf,[pwd '/figures/heatmaps/FigS2_heatmap'], 'fig')
saveas(gcf,[pwd '/figures/heatmaps/FigS2_heatmap'],'epsc')
clear

% Fig S3. Bacterial density after 96 h of combined treatment with intermediate innate immunity levels.
% Caption: We extend our robustness analysis of Figure 5 (bottom) to account for intermediate levels of
% innate immune activation in the context of combined therapy. We vary the levels of innate immune response 
% activation from 20% to 100% (a, b, c, d, and e). Bacterial density is calculated after 96 hr of treatment. 
% Colored regions represent bacterial presence while white regions indicate infection clearance. Phage and 
% antibiotic are administered two hours after infection. Antibiotic levels vary from 0.1 to 10 X MIC 
% (MIC of ciprofloxacin = 0.014ug/ml). Initial bacterial and phage density are, B0 = 7.4 x 10^7 CFU/g and P0 = 7.4 x 10^8 PFU/g, respectively.


% Fig. S3-Heatmaps
figure(24)
run([pwd '/code_comb_therapy/robustness_analysis/phage_anti_immunomodulation.m'])
saveas(gcf,[pwd '/figures/heatmaps/FigS3'], 'fig')
saveas(gcf,[pwd '/figures/heatmaps/FigS3'],'epsc')
clear

% Fig S4. Time delays in the application of the combined treatment.
% Caption: We extend our robustness analysis of Figure 5 (bottom) to account for time delays in the start 
% of the combined treatment. Phage and antibiotic were administered simultaneously 2, 4, 6, 8, and 10 hours 
% after the beginning of the infection in the presence (Fig. S6-top) or absence (Fig. S6-bottom) of innate 
% immunity. Colored regions on the heatmaps indicate bacterial presence while white regions indicate infection
% clearance. Antibiotic levels vary from 0.1 to 10 X MIC (MIC of ciprofloxacin = 0.014ug/ml). Initial bacterial 
% and phage density are, B0 = 7.4 x10^7 CFU/g and P0 = 7.4 x 10^8 PFU/g, respectively.

% Fig. S4-Left column (immunocompetence)
% figure(25)
run([pwd '/code_comb_therapy/robustness_analysis/combination_immune_variedtiming.m'])

save_label = {'A', 'B', 'C', 'D', 'E'};
for fig_num = 1:5,
    saveas(figure(fig_num),[pwd '/figures/heatmaps/FigS4' save_label{fig_num}],'fig')
    saveas(figure(fig_num),[pwd '/figures/heatmaps/FigS4' save_label{fig_num}],'epsc')
end
clear

% Fig. S4-Right column (neutropenic)
% figure(26)
run([pwd '/code_comb_therapy/robustness_analysis/phage_anti_variedtiming.m'])

save_label = {'F', 'G', 'H', 'I', 'J'};
for fig_num = 1:5,
    saveas(figure(fig_num),[pwd '/figures/heatmaps/FigS4' save_label{fig_num}],'fig')
    saveas(figure(fig_num),[pwd '/figures/heatmaps/FigS4' save_label{fig_num}],'epsc')
end
clear


% Fig S5. Parameter sensitivity analysis results.
% Caption: We show the distribution of the fraction of complete elimination for two therapeutic regimes
% A + P + I (blue) and A + P (red). We performed 1000 runs using perturbed parameter sets (\theta_per) and
% calculate the fraction of bacterial elimination for the two regimes. Moreover, we show the fraction of complete elimination for
% A + P + I (blue square) and A + P (red triangle) using the reference parameter set (\theta_ref ).

% NOTE: It takes ~8hrs to run the parameter sensitivity analysis for both regimes for 1000 iterations. So I'll just provide
% the data generated from the sensitivity analysis so I (and you) can generate the
% Figure S7 from the paper.

% Fig. S5
figure(27)
run([pwd '/code_comb_therapy/robustness_analysis/FigS5_param_sensitivity.m'])
saveas(gcf,[pwd '/figures/distributions/FigS5'], 'fig')
saveas(gcf,[pwd '/figures/distributions/FigS5'],'epsc')
clear

% Uncomment the following line of code if you want to run yourself the paramater sensitivity analysis for 1000 iterations
% the following line of code will also generate the figures from the
% generated data.

% run([pwd '/code_comb_therapy/robustness_analysis/param_sens_AP_API.m'])

% Fig S6. Robustness analysis of the partial resistance model
% Caption: Bacteria grew for 96 hours exposed to diferent antimicrobial strategies, antibiotic-only (a), antibiotic + innate immunity (b), phages and antibiotic (c),
% and phage-antibiotic combination plus active innate immunity (d). The heatmaps show the bacterial density at 96 hours post
% infection. Colored regions represent bacteria persistence while the clearance of the infection is represented by white regions.
% For the partial resistance model, phage can infect BA with an infection constant of \phi\delta_P , being \delta_P = 0.1. Moreover, the phagesensitive
% strain, BP , is slightly sensitive to the antibiotic with an MIC of 0.172 ug/ml (12 times higher than MIC for BA strain). We simulate different concentrations
% of ciprofloxacin (CP); using the MIC of CP for the BA strain as a reference (MIC = 0.014 ug/ml). Fold changes in MIC go
% from 0.1 MIC (0.014 ug/ml) to 10 MIC (0.14 ug/ml). Furthermore, bacterial composition of the inoculum ranges from 100%
% phage-sensitive bacteria (0% BA) to 100% antibiotic-sensitive bacteria (100% BA). Initial bacterial density and phage dose
% (c-d) are, B0 = 7.4 x 10^7 CFU/g and P0 = 7.4 x 10^8 PFU/g, respectively. Phage and antibiotic are administered two hours
% after the beginning of the infection.


% Fig. S6A
figure(19)
run([pwd '/code_comb_therapy/extended_model/neutropenic/heatmap_antibio_partres.m'])
saveas(gcf, [pwd '/figures/heatmaps/FigS6A'], 'fig')
saveas(gcf,[pwd '/figures/heatmaps/FigS6A'],'epsc')
clear

% Fig. S6B
figure(20)
run([pwd '/code_comb_therapy/extended_model/immunocompetent/heatmap_anti_immune_partres.m'])
saveas(gcf, [pwd '/figures/heatmaps/FigS6B'], 'fig')
saveas(gcf,[pwd '/figures/heatmaps/FigS6B'],'epsc')
clear

% Fig. S6C
figure(21)
run([pwd '/code_comb_therapy/extended_model/neutropenic/heatmap_anti_phage_partres.m'])
saveas(gcf, [pwd '/figures/heatmaps/FigS6C'], 'fig')
saveas(gcf,[pwd '/figures/heatmaps/FigS6C'],'epsc')
clear

% Fig. S6D
figure(22)
run([pwd '/code_comb_therapy/extended_model/immunocompetent/heatmap_anti_phage_immune_partres.m'])
saveas(gcf, [pwd '/figures/heatmaps/FigS6D'], 'fig')
saveas(gcf,[pwd '/figures/heatmaps/FigS6D'],'epsc')
clear

% Fig. S7 Bacterial dynamics of the partial resistance model
% Caption: Time series of the partial resistance model in presence and absence of host immune response. We simulate
% the combination of phage and antibiotic against a phage-sensitive (a) and an antibiotic-sensitive (b) bacterial inoculum in the
% absence of the immune response. Moreover, we simulate a within-host scenario where the combined therapy interact with the
% immune response (purple dashed line) and phage-sensitive bacteria (c) or antibiotic-sensitive bacteria (d). Here, phage (blue
% dashed line) and antibiotic are administered two hours after the infection. Initial conditions are, B0 = 7.4 x 10^7 CFU/g and
% P0 = 7.4x 10^8 PFU/g. The concentration of antibiotic (0.0350 ug/ml = 2.5 x MIC for BA strain) is maintained constant during
% the simulation, data not shown. The simulation run was 96 hours (4 days).

% Fig. S7A
figure(15)
run([pwd '/code_comb_therapy/extended_model/neutropenic/RHM_neutropenia.m'])
saveas(gcf, [pwd '/figures/time_series/FigS7A'], 'fig')
saveas(gcf,[pwd '/figures/time_series/FigS7A'],'epsc')
clear

% Fig. S7B
figure(16)
run([pwd '/code_comb_therapy/extended_model/neutropenic/RHM_neutropenia_BA.m'])
saveas(gcf, [pwd '/figures/time_series/FigS7B'], 'fig')
saveas(gcf,[pwd '/figures/time_series/FigS7B'],'epsc')
clear

% Fig. S7C
figure(17)
run([pwd '/code_comb_therapy/extended_model/immunocompetent/RHM_WT.m'])
saveas(gcf, [pwd '/figures/time_series/FigS7C'], 'fig')
saveas(gcf,[pwd '/figures/time_series/FigS7C'],'epsc')
clear

% Fig. S7D
figure(18)
run([pwd '/code_comb_therapy/extended_model/immunocompetent/RHM_WT_BA.m'])
saveas(gcf, [pwd '/figures/time_series/FigS7D'], 'fig')
saveas(gcf,[pwd '/figures/time_series/FigS7D'],'epsc')
clear
