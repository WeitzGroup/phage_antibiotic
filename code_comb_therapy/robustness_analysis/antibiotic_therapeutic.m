% Robustness analysis of the combination therapy model: 
% This code simulate the effect of antibiotic against mixed bacterial
% inoculums and generates a heatmap.
% Neutropenic host
% Dependencies: (1) rhmODE.m (2)simRHM.m (3) myEventsFcn.m

clear
clc
close all

% Neutropenia parameters:
Ki = 2.4e7; 
Io = 0; % Initial immune response
P = 0; %  no phage treatment
%P = 7.4e8; % phage treatment


% Antibiotic parameters for Ciprofloxacin
anti_name = 'CP';
MIC = 0.014;  % ug/ml, MIC of ciprofloxacin for BA strain

up_MIC_10 = [10:-0.5:1.5]; % from 1.5 MIC to 10 MIC
sub_MIC_0= [1:-.05:0.1]; %  from 0.1 MIC to 1 MIC
anti_conc = [up_MIC_10 sub_MIC_0];

% Inoculum proportions
B_total = 7.4e7; % CFU/g,  initial bacterial density
Bp_proportions = [1:-0.05:0];

                     
% Antibiotic parameters for Ceftazidime

%suggested_dosage = [0:10:600]*0.025;
%dose = suggested_dosage(7);
%dosing_interval = 1;
%anti_name = 'CAZ';

% Simulate the antibiotic therapeutic model for varied
% antibiotic concentrations and mixed bacterial inoculums

matrix = zeros(length(anti_conc), length(Bp_proportions));
matrix_threshold = zeros(length(anti_conc), length(Bp_proportions));
time_v = [];
for i = 1:length(anti_conc)
        anti_dose = anti_conc(i)*MIC;
        
        for j = 1:length(Bp_proportions)
            proportion = Bp_proportions(j);
            
            Bp = B_total * proportion;
            Ba = B_total - Bp;
            [y, TB, time] = simRHM(Ki, Io, Bp, Ba, P, anti_dose, anti_name);
            time_v = [time_v time(end)];
            Bp_steady = y(end, 1); % Take the bacterial density at 96 hours
            Ba_steady = y(end, 2);
            Btotal_steady = Bp_steady + Ba_steady;
            matrix(i,j) = Btotal_steady;
            if Bp_steady > Ba_steady
                bp = 1;
            else
                bp = 0;
            end
            matrix_threshold(i,j) = bp;
            
        end
end

m = matrix;
%---- Creat heatmap -----%

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

% Axis settings for the heatmaps
ylabel(h, 'Bacterial density (g^{-1})', 'FontSize', 16, 'fontweight', 'bold')
%title(h, 'log(B)', 'FontSize', 16, 'fontweight', 'bold')
xlim = get(gca, 'Xlim');
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
  set(h,'HorizontalAlignment','center','fontsize',get(gca,'fontsize'), 'fontweight','bold')
end
set(gca, 'fontweight', 'bold')
title('Antibiotic', 'FontSize', 20)

s = size(matrix_threshold);
ncol = s(2);
threshold = [];
for i = 1:ncol
    ind = find(matrix_threshold(:,i) == 1);
    if ind(end) == 37
       ind(end) = 38;
    end
    threshold = [threshold ind(end)];
end
hold on
plot([0:22], [threshold(1) threshold threshold(end)], '-k', 'Linewidth', 3 )
hold off