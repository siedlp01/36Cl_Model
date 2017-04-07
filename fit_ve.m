% Automated model fitting for site VE
data = load('data/datarockVE.txt');
data = flipud(data); % fit_scarp requires reverse order of samples: high to low
coll = load('data/datacolluviumVE.txt');
 
EL = load('mag_data/datamagfieldVEStone2000.txt');
 
% Schlagenhauf et al. (2011) ages and slip for VE
% age=[7500 4700 4400 4000 1300];
% slip=[130 110 170 250 130];
% preexp=6000; 

slip = [130 110 170 250 130];
preexp = 6000; 

% fixed parameters
epsilon=0; % erosion rate (mm/yr)

alpha=30; % colluvium wedge dip (degrees) 
beta=40;  % preserved scarp dip (degrees)
gamma=35; % upper eroded scarp dip (degrees)
Hfinal=950; % final height (present height) of the fault scarp of dip \beta (cm)

rho_coll = 1.5 ; % MA3 : 1.5
rho_rock = 2.71 ; % MA3 : 2.7

max_age = 8200;
age_step = 100;

%{
fit_scarp(max_age, age_step, data, coll, slip, preexp, EL, epsilon, alpha, beta, gamma, Hfinal, rho_coll, rho_rock, true, true);
%}

fit_pinned_scarp(data, coll, max_age, age_step, slip, preexp, EL, epsilon, alpha, beta, gamma, Hfinal, rho_coll, rho_rock, true)

