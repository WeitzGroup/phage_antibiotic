%% Code for the robustness analysis of the partial resistance model: Antibiotic therapeutic
% How variations in initial conditions (i.e., antibiotic concentration and
% bacterial composition of the inoculum) affect model outcomes
% Neutropenic host
% Dependencies: (1) rhmODE.m (2)simRHM.m (3) myEventsFcn.m


clear
clc
close all

% Neutropenia parameters:
Ki = 2.4e7;
Io = 0;
%P = 7.4e8; % phage treatment
P = 0; % no phage treatment

% Antibiotic parameters for Ciprofloxacin
anti_name = 'CP';
MIC = 0.014;  % ug/ml, MIC of ciprofloxacin for BA strain

% Antibiotic concentrations, from 0.1 to 10 MIC (of cipro for BA PAPS
% strain)
anti_conc = linspace(log10(10), log10(.1), 41);
anti_conc = 10.^anti_conc;

% Inoculum proportions
B_total = 7.4e7;
Bp_proportions = [1:-0.05:0];

                     
% delta parameter, modulate phage resistance for BA strain
delta_P = 1/10;

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
            [y, TB, time] = simRHM(Ki, Io, Bp, Ba, P, anti_dose, anti_name, delta_P);
            ind_96h = find(time >= 96 & time < 97);
            if length(ind_96h) > 1
                ind_96h = ind_96h(1);
            end
            time_v = [time_v time(ind_96h)];
            Bp_steady = y(ind_96h, 1);
            Ba_steady = y(ind_96h, 2);
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
%mat = imagesc(matrix);
matrix(find(matrix <= 1)) = 1;
matrix = log10(matrix);
mat = imagesc(matrix);
hold on
colormap(cmap);
h = colorbar;
caxis([0 10]) % white color mapped to 0, yellow (last color of colormap) mapped to 1e10

% Axis settings for the heatmaps
ylabel(h, 'log_{10} Bacterial density', 'FontSize', 16, 'fontweight', 'bold')
%title(h, 'log10 (B)', 'FontSize', 16, 'fontweight', 'bold')
xlim = get(gca, 'Xlim');
plot(xlim, [21 21], '--k', 'Linewidth', 2 )
hold off

xlabel('B_{A} inoculum proportion','FontSize',20)
ylabel('Fold Change in MIC','FontSize',20, 'Position', [-2.2 20])
set(gca,'XTick',[1:5:21])
xt = get(gca, 'XTick');
set(gca, 'XTick', xt, 'XTickLabel',[0:0.25:1], 'FontSize', 16)


concentrations = sort(anti_conc, 'descend');
set(gca,'YTick', [1:10:41], 'yticklabel',[],'ytickmode','manual');
yt = get(gca, 'YTick');
for j=1:length(yt)
  h=text(-1 , yt(j), ['10^{' num2str(log10( concentrations(yt(j)))) '}']);
  set(h,'HorizontalAlignment','center','fontsize',get(gca,'fontsize'), 'fontweight','bold')
end
set(gca, 'fontweight', 'bold')

title({"Antibiotic"; "MIC_{CP} for B_{P} = 0.172 \mug/ml, MIC_{CP} for B_{A} = 0.014 \mug/ml"}, 'FontSize', 20)

s = size(matrix_threshold);
ncol = s(2);
threshold = [];
for i = 1:ncol
    ind = find(matrix_threshold(:,i) == 1);
    if ind(end) == 41
       ind(end) = 42;
    end
    threshold = [threshold ind(end)];
end
hold on
plot([0:22], [threshold(1) threshold threshold(end)], '-k', 'Linewidth', 3 )
hold off
%text(20.1, 40, 'B_{A}', 'FontSize',16,'fontweight', 'bold')
ta = annotation('textarrow',[0.66 0.74], [0.25 0.17], 'units', 'normalized', 'String','B_{A}','FontSize',16, 'fontweight', 'bold', 'Linewidth',2);
text(10, 15, 'B_{P}', 'FontSize',16,'fontweight', 'bold')
text(0.02, 0.95, 'a)', 'units', 'normalized', 'FontSize',16,'fontweight', 'bold')
axis square