% Automated model fitting for site MA1
data = load('data/datarockMA1.txt');
coll = load('data/datacolluviumMA1.txt');
 
EL = load('mag_data/datamagfieldMA1Stone2000.txt');
 
% Schlagenhauf et al. (2011) ages and slip
% age=[8000 4800 4500 4000 1300];
% slip=[40 175 250 280 60];
% preexp=4200;

age = [12000   4500   4400   3700   2400];
slip=[40 175 250 280 60];
preexp=0; % pre-exposure duration (years)

% fixed parameters
epsilon=0; % erosion rate (mm/yr)

alpha=25; % colluvium wedge dip (degrees) 
beta=40;  % preserved scarp dip (degrees)
gamma=35; % upper eroded scarp dip (degrees)
Hfinal=1500; % final height (present height) of the fault scarp of dip \beta (cm)

rho_coll = 1.5 ; % MA3 : 1.5
rho_rock = 2.67 ; % MA3 : 2.7


modelscarp(data, coll, age, slip, preexp, EL, epsilon, alpha, beta, gamma, Hfinal, rho_coll, rho_rock)