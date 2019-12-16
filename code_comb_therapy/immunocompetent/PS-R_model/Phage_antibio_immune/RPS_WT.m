% Code to simulate Phage-Antibiotic combination therapy, Phage saturation model (ODE)
% Immunocompetent host
% Inoculum: Phage-sensitive bacteria (BP)
% Phage and Antibiotic added two hours after infection
% Dependencies: (1) rpsODE.m (2)simRPS_WT.m (3) myEventsFcn.m

clear
clc
close all

% Immunocompetence parameters:
Ki = 2.4e7; % Maximum carrying capacity of the immune response
Io = 2.7e6; % Initial Immune response
B = 7.4e7;  % Bacterial inoculum 
P = 7.4e8; % phage treatment
%P = 0; % no phage treatment

% Antibiotic parameters for Ciprofloxacin

dose = 0.014*2.5; % ug/ml
anti_name = 'CP'; % ciprofloxacin

% Antibiotic parameters for Ceftazidime
%suggested_dosage = [0:10:600]*0.025;
%dose = suggested_dosage(1);
%dosing_interval = 0;
%anti_name = 'CAZ';

% Simulate combination therapy against phage-sensitive inoculum under a
% phage-saturation context
[y, TB, time] = simRPS_WT(Ki, Io, B, 0, P, dose, anti_name);

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
    h_leg = legend(v, 'phages','host immunity','BP','BA', 'Location', 'northeast');
else
    % Legend for condition with NO phage treatment
    v = [i(1) h(1) j(1)];
    h_leg = legend(v,'host immunity','BP','BA', 'Location','northeast');
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
%set(gca, 'Units','inches','Position',[1 1 3 2.5])
set(h_leg, 'FontSize',20,'fontweight','normal')
set(gcf,'PaperPositionMode','manual','PaperPosition',[0.25 2.5 8 6],'PaperUnits','inches')
set(gcf, 'position', [440   369   602   429])
title(["Combined therapy + Immune response against B_{P} inoculum"; "Phage saturation model"], 'FontSize', 20, 'fontweight', 'bold')
text(0.02, 0.95, 'b)', 'units', 'normalized', 'FontSize',16,'fontweight', 'bold')