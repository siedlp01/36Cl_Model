% Automated model fitting for site VE
data = load('data/datarockVE.txt');
data = flipud(data); % fit_scarp requires reverse order of samples: high to low
coll = load('data/datacolluviumVE.txt');
 
EL = load('mag_data/datamagfieldVEStone2000.txt');
 
% Schlagenhauf et al. (2011) ages and slip for VE
% age=[7500 4700 4400 4000 1300];
% slip=[130 110 170 250 130];
% preexp=6000; 

age = [12900   3400   3400   3400   2200];
slip = [130 110 170 250 130];
preexp = 0; 

% fixed parameters
epsilon=0; % erosion rate (mm/yr)

alpha=30; % colluvium wedge dip (degrees) 
beta=40;  % preserved scarp dip (degrees)
gamma=35; % upper eroded scarp dip (degrees)
Hfinal=950; % final height (present height) of the fault scarp of dip \beta (cm)

rho_coll = 1.5 ; % MA3 : 1.5
rho_rock = 2.71 ; % MA3 : 2.7

modelscarp(data, coll, age, slip, preexp, EL, epsilon, alpha, beta, gamma, Hfinal, rho_coll, rho_rock)