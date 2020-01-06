% Code to simulate Phage-Antibiotic combination therapy on heterogeneous mixing model (ODE)
% Inoculum: Antibiotic-sensitive bacteria (BP)
% Phage and Antibiotic added two hours after infection
% Dependencies: (1) rhmODE.m (2)simRHM_WT.m (3) myEventsFcn.m

clear
clc
close all

% Immunocompetence parameters:
Ki = 2.4e7; % Carrying capacity of the immune response
Io = 2.7e6; % Initial immune response
B = 7.4e7;  % Initial bacterial inoculum
P = 7.4e8; % phage treatment
%P = 0; % no phage treatment

% Antibiotic parameters for Ciprofloxacin
dose = 0.014*2.5; % ug/ml
anti_name = 'CP'; 



% Simulate phage-antibiotic combination therapy against an
% antibiotic-sensitive inoculum
[y, TB, time] = simRHM_WT(Ki, Io, 0, B, P, dose, anti_name);


%----------------------------------------
% plotting
%----------------------------------------

% Do you want to show long run simulation (300 h)?
long_run = 0; % 1 for long run sim.

% Plot default values
set(0,'DefaultAxesLinewidth',2)
set(0, 'DefaultAxesFontName', 'Arial')

% Define color vectors
Bvector = [235 140 87]./255;
Rvector = [204 204 0]./255;
Pvector = [0 171 169]./255;
Ivector = [159 0 197]./255;
Avector = [0 0 0];
Kvector = [1 0 0];

figure(1)
semilogy(time,y(:,1),'Color', Bvector, 'Linewidth',3.5);
hold on;
semilogy(time,y(:,2),'Color', Rvector,'Linewidth',3.5)
semilogy(time,y(:,3),'LineStyle','--','Color', Pvector,'Linewidth',3.5)
semilogy(time,y(:,4),'LineStyle','--','Color', Ivector,'Linewidth',3.5)
%semilogy(time,y(:,5),'LineStyle','-','Color', Avector,'Linewidth',2)


%-----------------------------------------------------
xlabel('Hours post infection', 'FontSize', 16,'fontweight','bold')
if long_run
    set(gca,'XTick',[0:30:300])
else
    set(gca,'XTick',[0:12:96])
end
ylabel('Density (g^{-1})', 'FontSize', 16,'fontweight','bold')
h = findobj('Color',Bvector);
j = findobj('Color',Rvector);
g = findobj('Color',Pvector);
i = findobj('Color',Ivector);
d = findobj('Color', Avector);
k = findobj('Color', Kvector);

if P ~= 0
    % Legend for condition with phage treatment
    v = [g(1) i(1) h(1) j(1)];
    h_leg = legend(v, 'phage','host immunity','BP','BA', 'Location', 'southeast');
else
    % Legend for condition with NO phage treatment
    v = [i(1) h(1) j(1)];
    h_leg = legend(v,'host immunity','BP','BA', 'Location','southeast');
end

if long_run
    end_sim = 300;
else
    end_sim = 96;
end

if dose > 0 && dose < 1
    axis([0,end_sim,1,1e13])
else
    axis([0,end_sim,1,1e13])
end
legend boxoff
set(gca,'FontSize',20,'fontweight','bold')
set(h_leg, 'FontSize',20,'fontweight','normal')
set(gcf,'PaperPositionMode','manual','PaperPosition',[0.25 2.5 8 6],'PaperUnits','inches')
text(0.02, 0.95, 'b)', 'units', 'normalized', 'FontSize',16,'fontweight', 'bold')
title(["Combined therapy + Immune response against B_{A} inoculum"; "Heterogeneous mixing model"], 'FontSize', 20, 'fontweight', 'bold')
set(gcf, 'position', [440   369   602   429])
text(0.02, 0.95, 'd)', 'units', 'normalized', 'FontSize',16,'fontweight', 'bold')
