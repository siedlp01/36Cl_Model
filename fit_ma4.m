% Automated model fitting for site MA4
data = load('data/datarockMA4.txt');
coll = load('data/datacolluviumMA4.txt');
 
EL = load('mag_data/datamagfieldMA4Stone2000.txt');
 
% Schlagenhauf et al. (2011) ages and slip for MA4
% age=[8200 4500 1000];
% slip=[60 450 250];
% preexp=13000;

slip=[60 450 250];
preexp=13000; 

% fixed parameters
epsilon=0; % erosion rate (mm/yr)

alpha=30; % colluvium wedge dip (degrees) 
beta=42;  % preserved scarp dip (degrees)
gamma=30; % upper eroded scarp dip (degrees)
Hfinal=770; % final height (present height) of the fault scarp of dip \beta (cm)

rho_coll = 1.5 ; % MA3 : 1.5
rho_rock = 2.64 ; % MA3 : 2.7

max_age = 8500;
age_step = 100;
%{
fit_scarp(max_age, age_step, data, coll, slip, preexp, EL, epsilon, alpha, beta, gamma, Hfinal, rho_coll, rho_rock, true, true);
%}

fit_pinned_scarp(data, coll, max_age, age_step, slip, preexp, EL, epsilon, alpha, beta, gamma, Hfinal, rho_coll, rho_rock, true)
