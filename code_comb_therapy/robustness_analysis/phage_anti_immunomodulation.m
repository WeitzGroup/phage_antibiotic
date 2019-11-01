% Robustness analysis of the combination therapy model.
% This code simulates the phage-antibiotic combination therapy + immune
% modulated response agains mixed-strain infections.
% A heatmap is generated.
% Immunocompetent host
% Dependencies: (1) rhmODE.m (2)simRHM_WT.m (3) myEventsFcn.m


clear
clc
close all

percents = [0.2 0.4 0.6 0.8 1]; % percent of immune activation
fig_label = {'a)', 'b)', 'c)', 'd)', 'e)'};
tic
for perc = 1:length(percents),
    
    % Immunocompetence parameters:
    Ki = 2.4e7*percents(perc); % maximum capacity of immune response
    Io = 2.7e6*percents(perc); % Initial immune response
    P = 7.4e8; % phage treatment
    %P = 0; % no phage treatment

    % Antibiotic parameters for Ciprofloxacin
    anti_name = 'CP';
    MIC = 0.014;  % ug/ml, MIC of ciprofloxacin for antibiotic-sensitive strain (BA)
    
    % Antibiotic concentrations, from 0.1 to 10 MIC (of cipro for BA PAPS
    % strain)
    anti_conc = linspace(log10(10), log10(.1), 41);
    anti_conc = 10.^anti_conc;

    % Inoculum proportions
    B_total = 7.4e7;
    Bp_proportions = [1:-0.05:0];


    % Simulate the phage-antibiotic combination + immune response model for varied
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
                [y, TB, time] = simRHM_WT(Ki, Io, Bp, Ba, P, anti_dose, anti_name);
                time_v = [time_v time(end)];
                
                Bp_steady = y(end, 1);
                Ba_steady = y(end, 2);
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


    m = matrix;
    %%%--- Create the heatmap ----%%%

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

    %fig = figure(1);
    subplot(2, 3, perc);
    cmap = colormap(parula(100));
    cmap = [1, 1, 1; cmap];
    matrix(find(matrix <= 1)) = 1;
    matrix = log10(matrix);
    imagesc(matrix);
    colormap(cmap);
    hold on
    h = colorbar;
    caxis([0 10])
    
    
    ylabel(h, 'Bacterial density (g^{-1})', 'FontSize', 16, 'fontweight', 'bold')
    title(h, 'log10 (B)', 'FontSize', 16, 'fontweight', 'bold')
    xlim = get(gca, 'Xlim');
    plot(xlim, [21 21], '--k', 'Linewidth', 2 )
    hold off

    xlabel('B_{A} inoculum proportion','FontSize',20)
    ylabel('Fold Change in MIC','FontSize',20, 'Position', [-2.7 20])
    set(gca,'XTick',[1:5:21])
    xt = get(gca, 'XTick');
    set(gca, 'XTick', xt, 'XTickLabel',[0:0.25:1], 'FontSize', 16)


    concentrations = sort(anti_conc, 'descend');
    set(gca,'YTick', [1:10:41], 'yticklabel',[],'ytickmode','manual');
    yt = get(gca, 'YTick');
    for j=1:length(yt)
      h=text(-1.1 , yt(j), ['10^{' num2str(log10( concentrations(yt(j)))) '}']);
      set(h,'HorizontalAlignment','center','fontsize',get(gca,'fontsize'), 'fontweight','bold')
    end
    set(gca, 'fontweight', 'bold')
    
    title(['A + P + ' int2str(percents(perc)*100) '% I'], 'FontSize', 20)

    s = size(matrix_threshold);
    ncol = s(2);
    threshold = [];
    for i = 1:ncol
        ind = find(matrix_threshold(:,i) == 1);
        if isempty(ind)
            ind = [42];
        end
        threshold = [threshold ind(1)];
    end
    hold on
    plot([0:22], [threshold(1) threshold threshold(end)], '-k', 'Linewidth', 3 )
    hold off
    
    if percents(perc) < 0.8
        bp_text=text(10 , 15, 'B_P', 'fontsize',  16,'fontweight', 'bold');
    end
    ba_text=text(10 , 35, 'B_A','fontsize', 16,'fontweight', 'bold' );
    axis square
    text(0.02, 0.95, fig_label{perc}, 'units', 'normalized', 'FontSize',16,'fontweight', 'bold')

end
set(gcf, 'position', [299 434 1446 791])
toc
