% Automated model fitting for site MA3
data = load('data/datarockMA3_trunc.txt');
coll = load('data/datacolluviumMA3.txt'); 

EL = load('mag_data/datamagfieldMA3Stone2000.txt');
 
% Schlagenhauf et al. (2011) ages and slip for MA3
% age=[7900 4800 4400 4000 1300];
% slip=[190 205 160 360 200];
% preexp=2500;

slip=[190 205 160 360 200];
preexp=2500; 

% fixed parameters
epsilon=0; % erosion rate (mm/yr)

alpha=30; % colluvium wedge dip (degrees) 
beta=45;  % preserved scarp dip (degrees)
gamma=30; % upper eroded scarp dip (degrees)
Hfinal=2000; % final height (present height) of the fault scarp of dip \beta (cm)

rho_coll = 1.5 ; % MA3 : 1.5
rho_rock = 2.7 ; % MA3 : 2.7

max_age = 8000;
age_step = 100;

%{
fit_scarp(max_age, age_step, data, coll, slip, preexp, EL, epsilon, alpha, beta, gamma, Hfinal, rho_coll, rho_rock, true, true);
%}

fit_pinned_scarp(data, coll, max_age, age_step, slip, preexp, EL, epsilon, alpha, beta, gamma, Hfinal, rho_coll, rho_rock, true)
