% This code generates all the model-related non-schematic figures in the
% paper Rodriguez-Gonzalez, et al., Quantiative immunocompetent Models of
% Phage-Antibiotic Combination Therapy, bioRxiv 2019

clear
clc
close all

%% Main figures

% Fig. 2 Immunophage Model
% Caption: Dynamics of the immunophage therapy model against two different bacterial inoculums.
% We simulate the phage therapy model developed by \cite{Leung2017} against two infection settings.
% In the first infection setting (a), a phage-sensitive bacterial inoculum, BP (orange solid line), 
% is challenged with phage (blue dashed line) in an immunocompetence context. In the second scenario (b),
% antibiotic-sensitive bacteria, BA (green solid line), are challenged with phage in presence of an 
% active immune response (purple dashed line). The initial bacterial density and the initial phage dose are,
% B_0 = 7.4x10^7 CFU/g and P_0 = 7.4x10^8 PFU/g, respectively. The simulation run is 96
% hours with phage being administered 2 hours after the infection. The bacterial carrying capacity is, K_C = 10^10 CFU/g.

%Fig. 2A
figure(1)
run([pwd '/code_comb_therapy/immunocompetent/HM-R_model/Immuno_phage/RHM_WT.m'])
saveas(gcf, [pwd '/figures/time_series/Fig2A'], 'fig')
saveas(gcf,[pwd '/figures/time_series/Fig2A'],'epsc')
clear

%Fig. 2B
figure(2)
run([pwd '/code_comb_therapy/immunocompetent/HM-R_model/Immuno_phage/RHM_WT_BA.m'])
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
% bacteria persistence (e.g., blue areas ~ 10^9 CFU/g and yellow areas ~10^10 CFU/g) while 
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
figure(12)
run([pwd '/code_comb_therapy/robustness_analysis/phage_antibio_comb.m'])
saveas(gcf, [pwd '/figures/heatmaps/Fig5B'], 'fig')
saveas(gcf,[pwd '/figures/heatmaps/Fig5B'],'epsc')
clear

% Fig. 5C
figure(13)
run([pwd '/code_comb_therapy/robustness_analysis/antibiotic_immune.m'])
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

% Fig. S1 Bacterial dynamics of the partial resistance model
% Caption: Time series of the partial resistance model in presence and absence of host immune response. We simulate
% the combination of phage and antibiotic against a phage-sensitive (a) and an antibiotic-sensitive (b) bacterial inoculum in the
% absence of the immune response. Moreover, we simulate a within-host scenario where the combined therapy interact with the
% immune response (purple dashed line) and phage-sensitive bacteria (c) or antibiotic-sensitive bacteria (d). Here, phage (blue
% dashed line) and antibiotic are administered two hours after the infection. Initial conditions are, B0 = 7.4 x 10^7 CFU/g and
% P0 = 7.4x 10^8 PFU/g. The concentration of antibiotic (0.0350 ug/ml = 2.5 x MIC for BA strain) is maintained constant during
% the simulation, data not shown. The simulation run was 96 hours (4 days).

% Fig. S1A
figure(15)
run([pwd '/code_comb_therapy/extended_model/neutropenic/RHM_neutropenia.m'])
saveas(gcf, [pwd '/figures/time_series/FigS1A'], 'fig')
saveas(gcf,[pwd '/figures/time_series/FigS1A'],'epsc')
clear

% Fig. S1B
figure(16)
run([pwd '/code_comb_therapy/extended_model/neutropenic/RHM_neutropenia_BA.m'])
saveas(gcf, [pwd '/figures/time_series/FigS1B'], 'fig')
saveas(gcf,[pwd '/figures/time_series/FigS1B'],'epsc')
clear

% Fig. S1C
figure(17)
run([pwd '/code_comb_therapy/extended_model/immunocompetent/RHM_WT.m'])
saveas(gcf, [pwd '/figures/time_series/FigS1C'], 'fig')
saveas(gcf,[pwd '/figures/time_series/FigS1C'],'epsc')
clear

% Fig. S1D
figure(18)
run([pwd '/code_comb_therapy/extended_model/immunocompetent/RHM_WT_BA.m'])
saveas(gcf, [pwd '/figures/time_series/FigS1D'], 'fig')
saveas(gcf,[pwd '/figures/time_series/FigS1D'],'epsc')
clear

% Fig S2. Robustness analysis of the partial resistance model
% Caption: Bacteria grew for 96 hours exposed to diferent antimicrobial strategies, antibiotic-only (a), antibiotic + innate immunity (b), phages and antibiotic (c),
% and phage-antibiotic combination plus active innate immunity (d). The heatmaps show the bacterial density at 96 hours post
% infection. Colored regions represent bacteria persistence while the clearance of the infection is represented by white regions.
% For the partial resistance model, phage can infect BA with an infection constant of \phi\delta_P , being \delta_P = 0.1. Moreover, the phagesensitive
% strain, BP , is slightly sensitive to the antibiotic with an MIC of 0.172 ug/ml. We simulate different concentrations
% of ciprofloxacin (CP); using the MIC of CP for the BA strain as a reference (MIC = 0.014 ug/ml). Fold changes in MIC go
% from 0.1 MIC (0.014 ug/ml) to 10 MIC (0.14 ug/ml). Furthermore, bacterial composition of the inoculum ranges from 100%
% phage-sensitive bacteria (0% BA) to 100% antibiotic-sensitive bacteria (100% BA). Initial bacterial density and phage dose
% (c-d) are, B0 = 7.4 x 10^7 CFU/g and P0 = 7.4 x 10^8 PFU/g, respectively. Phage and antibiotic are administered two hours
% after the beginning of the infection.


% Fig. S2A
figure(19)
run([pwd '/code_comb_therapy/extended_model/neutropenic/heatmap_antibio_partres.m'])
saveas(gcf, [pwd '/figures/heatmaps/FigS2A'], 'fig')
saveas(gcf,[pwd '/figures/heatmaps/FigS2A'],'epsc')
clear

% Fig. S2B
figure(20)
run([pwd '/code_comb_therapy/extended_model/immunocompetent/heatmap_anti_immune_partres.m'])
saveas(gcf, [pwd '/figures/heatmaps/FigS2B'], 'fig')
saveas(gcf,[pwd '/figures/heatmaps/FigS2B'],'epsc')
clear

% Fig. S2C
figure(21)
run([pwd '/code_comb_therapy/extended_model/neutropenic/heatmap_anti_phage_partres.m'])
saveas(gcf, [pwd '/figures/heatmaps/FigS2C'], 'fig')
saveas(gcf,[pwd '/figures/heatmaps/FigS2C'],'epsc')
clear

% Fig. S2D
figure(22)
run([pwd '/code_comb_therapy/extended_model/immunocompetent/heatmap_anti_phage_immune_partres.m'])
saveas(gcf, [pwd '/figures/heatmaps/FigS2D'], 'fig')
saveas(gcf,[pwd '/figures/heatmaps/FigS2D'],'epsc')
clear