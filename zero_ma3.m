% Automated model fitting for site MA3
data = load('data/datarockMA3_trunc.txt');
coll = load('data/datacolluviumMA3.txt'); 

EL = load('mag_data/datamagfieldMA3Stone2000.txt');
 
% Schlagenhauf et al. (2011) ages and slip for MA3
% age=[7900 4800 4400 4000 1300];
% slip=[190 205 160 360 200];
% preexp=2500;

age = [9800  4300  4200  3800  1800];
slip=[190 205 160 360 200];
preexp=0; 

% fixed parameters
epsilon=0; % erosion rate (mm/yr)

alpha=30; % colluvium wedge dip (degrees) 
beta=45;  % preserved scarp dip (degrees)
gamma=30; % upper eroded scarp dip (degrees)
Hfinal=2000; % final height (present height) of the fault scarp of dip \beta (cm)

rho_coll = 1.5 ; % MA3 : 1.5
rho_rock = 2.7 ; % MA3 : 2.7


modelscarp(data, coll, age, slip, preexp, EL, epsilon, alpha, beta, gamma, Hfinal, rho_coll, rho_rock)