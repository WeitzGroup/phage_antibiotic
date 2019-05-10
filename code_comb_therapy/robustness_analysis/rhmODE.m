% Phage resistant heterogeneous mixing model ODE

function [ dx ] = rhmODE( t,x,p )

dx = [0;0;0;0;0];

%-----------Parameters------------
% phage-susceptible bacteria growth rate
r = p.r;
% antibiotic-sensitive bacteria growth rate
rp = p.rp;
% total bacteria carrying capacity
Kc = p.Kc;
% nonlinear adsorption rate of phage:
phi = p.phi;
% power law exponent in phage infection:
g = p.g;
% immune response killing rate parameter:
ep = p.ep;
% bacterial conc. at which immune response is half as effective:
Kd = p.Kd;
% burst size of phage:
beta = p.beta;
% decay rate of phage:
w = p.w;
% maximum growth rate of immune response:
a = p.a;
% max capacity of immune response:
Ki = p.Ki;
% conc. of bacteria at which imm resp growth rate is half its maximum:
Kn = p.Kn;
% probability of emergence of phage-resistant mutation per cell division
m = p.m;
% probability of reversible mutation probability (Phage resistant -> Phage sensitive)
m2= p.m2;


% Maximum antibiotic-mediated kill rate
kkill = p.kkill;
% EC_50, antibiotic concentration when kkill/2
ec = p.ec;
% Steepness of the sigmoid relationship between Antibiotic and death rate
H = p.H;
% Elimination rate of antibiotic from serum
theta = p.theta;
% Antibiotic constant pulse, such that A is at equilibrium
pulse = p.pulse;


B = x(1);
R = x(2);
P = x(3);
I = x(4);
A = x(5);

% Change in susceptible bacterial population
%dx(1) = (r*B*(1-((B+R)/Kc))*(1-m))-(phi*B*(P^g))-(ep*I*B/(1+((B+R)/Kd)));
dx(1) = (r*B*(1-((B+R)/Kc))*(1-m))-(phi*B*(P^g))-(ep*I*B/(1+((B+R)/Kd)))+ rp*R*(1-((B+R)/Kc))*m2; % with reversible mutation
% Change in resistant bacterial population
%dx(2) = ((rp*R*(1-((B+R)/Kc)))+r*B*(1-((B+R)/Kc))*m)-(ep*I*R/(1+((B+R)/Kd))) -((kkill*(A^H)*R)/(ec^H + A^H));
dx(2) = ((rp*R*(1-((B+R)/Kc)))*(1-m2)+r*B*(1-((B+R)/Kc))*m)-(ep*I*R/(1+((B+R)/Kd))) -((kkill*(A^H)*R)/(ec^H + A^H)); % with reversible mutation
% Change in phage population
dx(3) = (beta*phi*B*(P^g))-(w*P);
% Change immune response
dx(4) = a*I*((B+R)/(B+R+Kn))*(1-(I/Ki));
% Change in antibiotic
dx(5) = pulse - (theta*A);

end

