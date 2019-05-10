%% Code for the robustness analysis of the partial resistance model: Antibiotic + Phage + Innate immunity therapeutic
% How variations in initial conditions (i.e., antibiotic concentration and
% bacterial composition of the inoculum) affect model outcomes
% Immunocompetent host
% Dependencies: (1) rhmODE.m (2)simRHM_WT.m (3) myEventsFcn.m


clear
clc
close all

% Immunocompetence parameters:
Ki = 2.4e7; % maximum capacity of immune response
Io = 2.7e6; % initial immune response
P = 7.4e8; % phage treatment
%P = 0; % no phage treatment

% Antibiotic parameters for Ciprofloxacin
anti_name = 'CP';
MIC = 0.014;  % ug/ml, MIC of ciprofloxacin for BA strain
up_MIC_10 = [10:-0.5:1.5]; % from 1.5 MIC to 10 MIC
sub_MIC_0= [1:-.05:0.1]; % from 0.1 MIC to 1 MIC
anti_conc = [up_MIC_10 sub_MIC_0];

% Inoculum proportions
B_total = 7.4e7;
Bp_proportions = [1:-0.05:0];

                     
% delta parameters (modulate phage resistance of BA strain)
delta_P = 1/10;

% Simulate the combination therapy + immune response model for varied
% antibiotic concentrations and mixed bacterial inoculums

matrix = zeros(length(anti_conc), length(Bp_proportions));
matrix_threshold = zeros(length(anti_conc), length(Bp_proportions));
for i = 1:length(anti_conc)
        anti_dose = anti_conc(i)*MIC;
        
        for j = 1:length(Bp_proportions)
            proportion = Bp_proportions(j);
            
            Bp = B_total * proportion;
            Ba = B_total - Bp;
            [y, TB, time] = simRHM_WT(Ki, Io, Bp, Ba, P, anti_dose, anti_name, delta_P);
            ind_96h = find(time >= 96 & time < 97);
            if length(ind_96h) > 1
                ind_96h = ind_96h(1);
            end
            Bp_steady = y(ind_96h, 1);
            Ba_steady = y(ind_96h, 2);
            Btotal_steady = Bp_steady + Ba_steady;
            matrix(i,j) = Btotal_steady;
            if Ba_steady > Bp_steady
                ba = 1;
            else
                ba = 0;
            end
            matrix_threshold(i,j) = ba;
            
        end
end


% Heatmap configuration settings for log10 scale data
%fig = figure(1);
%cmap = colormap(parula(50));
%cmap(1,:) = [1 1 1]; % map white color to 0 Bacterial density
%mat = imagesc(matrix);
%c = log10(matrix(find(mat.CData ~=0)));
%matrix(find(mat.CData ~=0)) = c;
%mat = imagesc(matrix);
%hold on
%colormap(cmap);
%h = colorbar;
%caxis([0 10]) % white color mapped to 0, yellow (last color of colormap) mapped to 10

% Heatmap configuration settings for raw values
fig = figure(1);
cmap = colormap(parula(100));
cmap = [1,1,1; cmap]; % map white color to 0 Bacterial density
mat = imagesc(matrix);
hold on
colormap(cmap);
h = colorbar;
caxis([0 1e10]) % white color mapped to 0, yellow (last color of colormap) mapped to 1e10

ylabel(h, 'Bacterial density (g^{-1})', 'FontSize', 16, 'fontweight', 'bold')
xlim = get(gca, 'Xlim');
%title(h, 'log(B)', 'FontSize', 16, 'fontweight', 'bold')
plot(xlim, [19 19], '--k', 'Linewidth', 2 )
hold off

xlabel('B_{A} inoculum proportion','FontSize',20)
ylabel('Fold Change in MIC','FontSize',20, 'Position', [-1.5 19])
set(gca,'XTick',[1:2:21])
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel',[0:0.1:1], 'FontSize', 16)


concentrations = sort(anti_conc, 'descend');
set(gca,'YTick', [1:18:37], 'yticklabel',[],'ytickmode','manual');
yt = get(gca, 'YTick');
for j=1:length(yt)
  h=text(-.5 , yt(j), ['10^{' num2str(log10( concentrations(yt(j)))) '}']);
  set(h,'HorizontalAlignment','center','fontsize',get(gca,'fontsize'), 'fontweight', 'bold')
end
set(gca, 'fontweight', 'bold')
title({"Antibiotic + Phages + Immune response"; "MIC_{CP} for B_{P} = 0.172 \mug/ml, MIC_{CP} for B_{A} = 0.014 \mug/ml"}, 'FontSize', 20)

s = size(matrix_threshold);
ncol = s(2);
threshold = [];
for i = 1:ncol
    ind = find(matrix_threshold(:,i) == 1);
    if isempty(ind)
        ind = [38];
    end
    threshold = [threshold ind(1)];
end
hold on
plot([1:22], [threshold threshold(end)], '-k', 'Linewidth', 3 )
hold off
%set(gcf, 'position', [1000 827 637 511], 'units', 'pixels')