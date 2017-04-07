% Automated model fitting for site MA2
data = load('data/datarockMA2.txt');
data = flipud(data); % fit_scarp requires reverse order of samples: high to low
coll = load('data/datacolluviumMA2.txt');
 
EL = load('mag_data/datamagfieldMA2Stone2000.txt');
 
% Schlagenhauf et al. (2011) ages and slip for MA2
% age=[14700 11000 10100 9700 6000 1600];
% slip=[85 125 130 130 175 210];
% preexp=9800;

slip=[85 125 130 130 175 210]; %[85 125 260 175 210];
preexp=9800; 

% fixed parameters
epsilon=0; % erosion rate (mm/yr)

alpha=35; % colluvium wedge dip (degrees) 
beta=50;  % preserved scarp dip (degrees)
gamma=35; % upper eroded scarp dip (degrees)
Hfinal=860; % final height (present height) of the fault scarp of dip \beta (cm)

rho_coll = 1.5 ; % MA3 : 1.5
rho_rock = 2.67 ; % MA3 : 2.7

max_age = 15000;
age_step = 100;
%{
fit_scarp(max_age, age_step, data, coll, slip, preexp, EL, epsilon, alpha, beta, gamma, Hfinal, rho_coll, rho_rock, true, true);
%}

fit_pinned_scarp(data, coll, max_age, age_step, slip, preexp, EL, epsilon, alpha, beta, gamma, Hfinal, rho_coll, rho_rock, true)
