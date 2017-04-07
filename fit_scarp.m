function Nf = fit_scarp(maxAge, ageStep, data,coll,slip,preexp,EL,epsilon,alpha,beta,gamma,Hfinal,rho_coll,rho_rock, debug, graphit)
% Given the slip for each event, calculate the best fit for the ages of the
% earthquakes by looking for the best fit (least root mean square) for each
% age/event

%Psi_Cl36_Ca_0 = 48.8 ;% (at of Cl36 /g of Ca per yr)
Psi_Cl36_Ca_0 = 42.2 ;% (at of Cl36 /g of Ca per yr)
%
%--------------------------------------------------------------------------

%---------------------CONSTANTS--------------------------------------------
% Radioactive decay constant for 36Cl (a-1)
lambda36 = 2.303e-6 ;
%
% True attenuation length for fast neutron (g.cm-2)
Lambda = 208 ;
%--------------------------------------------------------------------------

%----------------EARTH MAG FIELD LOADING-----------------------------------
% Loading of Earth magnetic field variations from file 'EL'

ti = EL(:,1) ; % time period (years)
it = EL(:,2) ; % time steps (years) - should be 100 yrs
EL_f = EL(:,3) ; % scaling factor for neutrons (S_el,f)
EL_mu = EL(:,4) ; % scaling factor for muons (S_el,mu)
%--------------------------------------------------------------------------

%---------------VARIABLES INITIALIZATION-----------------------------------
[m,n] = size(data) ; if n ~= 66, error('File data must have 66 columns'), end
nc = size(coll,2) ; if nc ~= 62, error('File coll must have 62 columns'), end
nel = size(EL,2) ; if nel ~= 4, error('File EL must have 4 columns'), end

%
N_eq = length(slip) ; % number of earthquakes
age = zeros(size(slip));

%
R = sum(slip) ; % total cumulative slip
Rc = cumsum(slip) ;  Rc = [0 Rc] ;% slip added up after each earthquake
%
if Hfinal < R , error('Hfinal cannot be lower than cumulative slip R') , end
%
Hinit = Hfinal - R ; % initial height of the scarp during pre-exposure

%
%--------------------------------------------------------------------------

%--------------------SURFACE SCALING---------------------------------------
%               using scsurf.m for z>=0
% Calculates a scaling factor S_S(z>=0) every cm used for the samples at  
% surface which is normalized by S_S(z=0) after in the calculation of production 
% at surface (Parts B and C). This allows to take into account for the 
% presence of upper part of dip gamma.
%
Zs = 0:1:R ; % initialization of Zs ; one calculation point every cm  
S_S = zeros(size(Zs)) ; % initialization of S_S (Surface Scaling)
%
for i = 1:length(Zs)    % loop on Zs
    a = scsurf(Zs(i),Hfinal,Lambda,beta,gamma,rho_rock) ;
    S_S(i) = a ;
end
%--------------------------------------------------------------------------

%----------------------DEPTH SCALING FOR NEUTRONS--------------------------
%       using scdepth.m, function of Hiseg (earthquakes) for z<=0
% Calculates a scaling factor S_D(z<=0) every 10 cm fitted by fitexp.m
% (S_D=so_f.exp(-z/Lambda_f). The derived so_f and Lambda_f depend
% on the height of the scarp of dip beta which grows after each earthquake
% (Hiseg = Hinitial + Rc(i) with earthquake i), so that so_f_d_iseg and 
% Lambda_f_d_iseg are calculated iteratively.
% They are used later (parts B and C) to calculate the productions at depth
% which are then scaled to production at z=0 to derive a scaling factor
% function of z<=0.
%
so_f_diseg = zeros(1,N_eq) ; % initialization of so_f_d_iseg
Lambda_f_diseg = zeros(1,N_eq) ; % initialization of Lambda_f_d_iseg
%
for is = 1:N_eq     % earthquake loop
    Hiseg = Hinit + Rc(is) ; % Height of exhumed scarp after each earthquake
    Ziseg = 0:10:R ; % one calculation point every 10 cm is sufficient  
    Ziseg = -Ziseg ; % negative because at depth
    S_D_iseg = zeros(size(Ziseg)) ; % initialization of S_D_iseg
    for i = 1:length(Ziseg)     % loop on z
        a = scdepth(Ziseg(i),Hiseg,Lambda,alpha,beta,gamma,rho_rock,rho_coll) ;
        S_D_iseg(i) = a ;
    end
    [dd,ee] = fitexp(-Ziseg*rho_coll,S_D_iseg,Lambda) ; % fit by fitexp.m
    so_f_diseg(is) = dd ; % constant so
    Lambda_f_diseg(is) = ee ; % attenuation length for neutron in direction z 
end
%
% attenuation length perpendicular to colluvium surface after each
% earthquake (with H increasing after each earthquake):
Lambda_f_diseg = Lambda_f_diseg*sind(beta - alpha) ; 

%---------------------------
% For beta infinite plane (used in B2 and C6):
Zbeta_inf = 0:10:1000; Zbeta_inf = -Zbeta_inf ; % initialization
S_D_beta_inf = zeros(size(Zbeta_inf));
%
for i = 1:length(Zbeta_inf)     % loop on z
        a = scdepth(Zbeta_inf(i),2000,Lambda,alpha,beta,gamma,rho_rock,rho_coll) ;
        S_D_beta_inf(i) = a ;
end
[so_f_beta_inf,Lambda_f_beta_inf] = fitexp(-Zbeta_inf*rho_coll,S_D_beta_inf,Lambda) ; % fit by fitexp.m
Lambda_f_beta_inf = Lambda_f_beta_inf*sind(beta - alpha) ; % attenuation perp. to colluvium surface
%--------------------------------------------------------------------------

%-------------------ROCK SCALING FOR NEUTRONS------------------------------
%        using scrock.m (attenuation in the direction of 'e')
%
e = 0:1:100 ; % e is in cm and perpendicular to scarp surface
Se = zeros(size(e)) ; % initialization of scaling Se
for i = 1:length(e)         % Loop on e
	Se(i) = scrock(e(i),Lambda,beta,rho_rock) ;
end
%
[so_f_e,Lambda_f_e] = fitexp(e*rho_rock,Se,Lambda) ; % exponential fit
%--------------------------------------------------------------------------

%---------------VARIABLES INITIALIZATION-----------------------------------
%
% h must be in cm and integers
h = data(:,n-3) ;  % initial positions of the samples at surface (cm)- integer
Z = (R - data(:,n-3))*rho_coll ;

d = data ; % substitution of matrix data by matrix d
d(:,n-3) = Z ; % samples position along z
d(:,n-2) = data(:,n-2)*rho_rock ; % thickness converted in g.cm-2
cl36AMS = d(:,n-1) ; % sample concentration in [36Cl] measured by AMS
sig_cl36AMS = d(:,n) ; % uncertainty on [36Cl] AMS measurements


slip_gcm2 = slip*rho_coll ; % coseismic slip in g.cm-2
sc = cumsum(slip_gcm2) ; % cumulative slip after each earthquake (g.cm-2)
sc0 = [0 sc] ;

% Positions along e initially (eo)
thick = data(:,n-2) ;
th2 = (thick/2)*rho_rock ; % 1/2 thickness converted in g.cm-2
eo = zeros(size(Z)) ;
%for iseg = 1:N_eq
%    eo(Z > sc0(iseg) & Z <= sc0(iseg + 1)) = epsilon*age(iseg)*0.1*rho_rock ; % in g.cm-2
%end
%eo(length(Z)) = epsilon*age(1)*0.1*rho_rock ;
%eo = eo + th2 ; % we add the 1/2 thickness : sample position along e is given at the sample center
%--------------------------------------------------------------------------

%----- B ------------------------------------------------------------------
% comment pre-exposure part if pre-exp = 0, and uncomment the line below:
% N_in = zeros(size(Z)) ; Ni = zeros(size(Z)) ; Nf = zeros(size(Z)) ;
%-----------------------------PRE-EXPOSURE PROFILE-------------------------
%
% Calculation of [36Cl] concentration profile at the end of pre-exposure.
%
% initialization at 0

 % No : initial concentration (here = zero) before pre-exposure
Ni = zeros(size(Z)); % Ni :  
Nf = zeros(size(Z)) ; % Nf : final 36Cl concentration 

%--------------------------------------------------------------------------


%----- C ------------------------------------------------------------------
%-----------------------------SEISMIC PHASE--------------------------------
%
% -the term 'segment' is used for the samples exhumed by an earthquake-
% Calculation of [36Cl] profiles during seismic cycle.
% Separated in two stages : 
%   * when samples are at depth and progressively rising because of earthquakes
%   (moving in the direction z with their position in direction e fixed)
%   * and when samples are brought to surface and only sustaining erosion
%   (moving along the direction e)
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% FIRST EXHUMED SEGMENT is treated alone.
%
% variables initialization: 
j1 = find(Z >= sc0(1) & Z <= sc0(2)) ; % samples from first exhumed segment 

jmax = j1(length(j1));

seg1cl36AMS = d(j1, n-1);
seg1sig_cl36AMS = d(j1, n);

%
% C1 - Loop - iteration on samples (k) from first exhumed segment
minRMSw = 0;
bestAge = age(1);
for current = maxAge:-ageStep:0
    if debug == true
        display(['Calculating segment 1: using age ', num2str(current)]);
    end
    
    if preexp > 0 
        Ni = zeros(size(Z)) ;
        No = zeros(size(Z)) ;
        tt = find(ti <= (current + preexp) & ti > current) ; % epoch index corresponding to pre-exposure
        ip = it(tt) ; % corresponding intervals

        % B1 - Loop - iteration on every samples
        for j = 1:m 

            dpj = d(j,:) ;
            dpj(n-3) = dpj(n-3)*sin((beta - alpha)*pi/180) ; % in the direction perpendicular to colluvium surface
            d0 = dpj ;
            d0(n-3) = 0 ;

        % scaling is separated in two expressions: exp(-z/Lambda_f)*exp(-e/Lambda_f_e).

        % in z direction : samples are at Z(j) = R - h(j)
        % in e direction : samples are at position eo(j) = epsilon*T(j), where T(j) is the age of the earthquake which brings the sample j at surface.

            N_in = No(j) ; % initial concentration (here = zero)
            N_out=N_in;
            % B2 - LOOP - iteration on time (ii) during pre-exposure
            for ii = 1:length(tt)
                [P_cosmo,P_rad] = clrock(d(j,:),eo(j),Lambda_f_e,so_f_e,EL_f(tt(ii)),EL_mu(tt(ii)),Psi_Cl36_Ca_0,rho_rock) ;

                % scaling at depth due to the presence of the colluvium: scoll=Pcoll(j)/Pcoll(z=0)
                P_coll = clcoll(coll,dpj,Lambda_f_diseg(1),so_f_diseg(1),EL_f(tt(ii)),EL_mu(tt(ii)),Psi_Cl36_Ca_0,rho_rock,rho_coll) ;
                P_zero = clcoll(coll,d0,Lambda_f_beta_inf,so_f_beta_inf,EL_f(tt(ii)),EL_mu(tt(ii)),Psi_Cl36_Ca_0,rho_rock,rho_coll) ;
                scoll = P_coll/P_zero ; 

                P_tot = P_rad + P_cosmo*scoll; % only P (Pcosmogenic) is scalled by scoll
                N_out = N_in + (P_tot - lambda36*N_in)*ip(ii) ; % minus radioactive decrease during same time step
                N_in = N_out ;
            end    

            Ni(j) = N_out ;

        end

    end
    
    N1 = zeros(size(Z(j1))) ;
    tt = find(ti <= current) ; % epoch index more recent than first earthquake
    ip = it(tt) ; % time intervals corresponding

    eo = zeros(size(Z)) ;
    eo(Z > sc0(1) & Z <= sc0(1 + 1)) = epsilon*current*0.1*rho_rock ; % in g.cm-2
    eo(length(Z)) = epsilon*current*0.1*rho_rock ;
    eo = eo + th2 ;

    for k = 1:length(j1)
        djk = d(j1(k),:) ;  
        hjk = h(j1(k)) ;    % position of sample k (cm)
        N_in = Ni(j1(k)) ;  % initial concentration is Ni, obtained after pre-exposure
        ejk = eo(j1(k)) ;   % initial position along e is eo(j1(k)) 
    
        % C2 - Loop - iteration on  time steps ii from t1 (= age eq1) to present
        for ii = 1:length(tt)
            [P_cosmo,P_rad] = clrock(djk,ejk,Lambda_f_e,so_f_e,EL_f(tt(ii)),EL_mu(tt(ii)),Psi_Cl36_Ca_0,rho_rock) ;
        
            scorr = S_S(1+hjk)/S_S(1) ;     % surface scaling factor (scorr)
            P_tot = P_rad + P_cosmo*scorr ;           % only Pcosmogenic is scaled with scorr
            N_out = N_in + (P_tot - lambda36*N_in)*ip(ii) ; % minus radioactive decrease during same time step
        
            ejk = ejk - epsilon*ip(ii)*0.1*rho_rock ; % new position along e at each time step (g.cm-2)
            N_in = N_out ; 
        end
    
        N1(k) = N_out ;
    
    end

    seg1rmsw = ((seg1cl36AMS - N1)./seg1sig_cl36AMS).^2 ;
    seg1rmsw = sum(seg1rmsw) ; % rmsw = sum(rmsw)/m ; 
    seg1rmsw = sqrt(seg1rmsw);
    
    if debug == true
        display(['RMSw for segment 1: using age ', num2str(current), ' is ', sprintf('%0.15f', seg1rmsw)]);
    end 
    
    if minRMSw == 0 || seg1rmsw < minRMSw 
        if debug == true
            disp(['Age ', num2str(current), ' is a better fit: ']);
        end
        bestAge = current;
        minRMSw = seg1rmsw;
        Nf(j1) = N1 ;
    else
        age(1) = bestAge;
        if debug == true
            disp(['Age ', num2str(age(1)), ' is the best fit for first event']);
        end
        break;
    end
    
end

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% ITERATION ON SEGMENTS 2 to N_eq
%
% C3 - Loop - iteration on each segment (from segment 2 to N_eq=number of eq)

for iseg = 2:N_eq
    if debug == true
        disp(['Analysing segment ', num2str(iseg)]);
    end
    j = find(Z > sc0(iseg) & Z <= sc0(iseg+1)) ; % index of samples from segment iseg
    
    if isempty(j)
        if debug == true
            disp(['Empty find using previous age ']);
        end
        age(iseg) = age(iseg - 1);
        continue;
    end
    jmax = j(length(j));
    segncl36AMS = d(1:jmax, n-1);
    segnsig_cl36AMS = d(1:jmax, n);
    
    minRMSwn = -1;
    bestAgeN = age(iseg - 1);
   
    z_j = Z(j) ; % initial depth along z of these samples (g.cm-2)
    upper = age(iseg - 1);
    %if iseg == N_eq 
    %    upper = 0;
    %end
    for current = upper:-ageStep:0
        if debug == true
            disp(['Calculating segment ', num2str(iseg), ': using age ', num2str(current), ' j ', mat2str(j), ' ', num2str(jmax)]);
        end
        N_new = zeros(size(z_j)) ;
        eo = zeros(size(Z)) ;
        for s = 1:iseg
            a = age(s);
            if s == iseg 
                a = current;
            end
            eo(Z > sc0(s) & Z <= sc0(s + 1)) = epsilon*a*0.1*rho_rock ; % in g.cm-2
        end
        eo(length(Z)) = epsilon*age(1)*0.1*rho_rock ;
        eo = eo + th2 ;
        % C4 - Loop - iteration each sample from segment iseg     
        for k = 1:length(j)                                                  
        
            ejk = eo(j(k)) ; % initial position along e is stil eo.
            djk = d(j(k),:) ;
            djk(n-3) = djk(n-3)*sind(beta - alpha) ;
        
            N_in = Ni(j(k)) ; %  initial concentration is Ni
        
            % C5 - Loop - iteration on previous earthquakes
            for l = 1:iseg-1   
                ttt_limit = age(l+1);
                if l == iseg -1
                    ttt_limit = current;
                end
                ttt = find(ti <= age(l) & ti > ttt_limit) ; % epoch index 
                ipp = it(ttt) ; % time intervals corresponding
            
                % depth (along z) are modified after each earthquake
                djk(n-3) = djk(n-3) - slip(l)*rho_coll*sind(beta - alpha) ;
                d0 = djk ;
                d0(n-3) = 0 ;
            
%------------------------------            
            % C6 - DEPTH LOOP - iteration during BURIED PERIOD (T1 -> T(iseg-1))
%------------------------------ 
                for iii = 1:length(ttt)
                    [P_cosmo,P_rad] = clrock(djk,ejk,Lambda_f_e,so_f_e,EL_f(ttt(iii)),EL_mu(ttt(iii)),Psi_Cl36_Ca_0,rho_rock) ;
                
                    % scaling at depth due to the presence of the colluvium: scoll=Pcoll(j)/Pcoll(z=0)
                    P_coll = clcoll(coll,djk,Lambda_f_diseg(l+1),so_f_diseg(l+1),EL_f(ttt(iii)),EL_mu(ttt(iii)),Psi_Cl36_Ca_0,rho_rock,rho_coll) ;
                    P_zero = clcoll(coll,d0,Lambda_f_beta_inf,so_f_beta_inf,EL_f(ttt(iii)),EL_mu(ttt(iii)),Psi_Cl36_Ca_0,rho_rock,rho_coll) ;
                    scoll = P_coll/P_zero ; 
                
                    P_tot = P_rad + P_cosmo*scoll; % only P (Pcosmogenic) is scalled by scoll
                    N_out = N_in + (P_tot - lambda36*N_in)*ipp(iii) ; % minus radioactive decrease during same time step
                    N_in = N_out ;
                end
            
                N_in = N_out ; 
            
            end

            N_in = N_out ;
        
            tt = find(ti <= current) ; % epoch index more recent than earthquake iseg
            ip = it(tt) ; % time intervals corresponding
            djk = d(j(k),:) ;
            hjk = h(j(k)) ;
        
%------------------------------         
            % C7 - SURFACE LOOP - iteration during EXHUMED PERIOD 
%------------------------------ 
            for ii = 1:length(tt)
                [P_cosmo,P_rad] = clrock(djk,ejk,Lambda_f_e,so_f_e,EL_f(tt(ii)),EL_mu(tt(ii)),Psi_Cl36_Ca_0,rho_rock) ;
                
                scorr = S_S(1+hjk)/S_S(1) ; % surface scaling factor (scorr)
                P_tot = P_rad + P_cosmo*scorr ; % only Pcosmogenic is scaled with scorr
                N_out = N_in + (P_tot - lambda36*N_in)*ip(ii) ; % minus radioactive decrease during same time step
              
                ejk = ejk - epsilon*ip(ii)*0.1*rho_rock ; % new position along e at each time step (g.cm-2)
                N_in = N_out ;
            end
        
            N_new(k) = N_out ;

        end
    
        N_test = Nf;
        N_test(j) = N_new;
        segnrmsw = ((segncl36AMS - N_test(1:jmax))./segnsig_cl36AMS).^2 ;
        segnrmsw = sum(segnrmsw) ; % rmsw = sum(rmsw)/m ; 
        segnrmsw = sqrt(segnrmsw);
        
        if debug == true
            display(['RMSw for segment ', num2str(iseg), ': using age ', num2str(current), ' is ', sprintf('%0.15f', segnrmsw)]);
        end
        
        if (minRMSwn == -1 || segnrmsw == minRMSwn || segnrmsw < minRMSwn) 
            if debug == true
                disp(['Age ', num2str(current), ' is a better fit.']);
            end
            bestAgeN = current;
            age(iseg) = current;
            minRMSwn = segnrmsw;
            Nf(j) = N_new ;
        else
             age(iseg) = bestAgeN;
             if debug == true
                 disp(['Age ', num2str(age(iseg)), ' is the best fit for event ', num2str(iseg)]);
             end
             break;
        end
    end    
end

% RMSw (weighted least square) :
rmsw = ((cl36AMS - Nf)./sig_cl36AMS).^2 ;
rmsw = sum(rmsw) ; % rmsw = sum(rmsw)/m ; 
rmsw = sqrt(rmsw) ;

% AICC (Akaike Information Criterion):
% if file data contains samples from the buried part of the scarp
% then, nb_param = 2*N_eq + 3 "-2" ; 
% (and we add an earthquake at time = 0 and of slip = height of buried samples

nb_param = 2*N_eq + 3 ; % 2*N_eq + pre-exp + variance + erosion (epsilon)
hauteur = h;
if any(age == 0)
    nb_param = nb_param - 2 ;
    hauteur = h - slip(end) ;
end

aicc = ak(cl36AMS,Nf,nb_param) ;

% Chi_square
chi_square = ((cl36AMS - Nf)./sig_cl36AMS).^2 ;
chi_square = sum(chi_square) ;
chi_square = (1/(m - nb_param - 1))*chi_square;

if graphit == true
    % error bar coordinates
    X = [cl36AMS(:)-sig_cl36AMS(:) cl36AMS(:)+sig_cl36AMS(:)] ;
    X = X' ;
    Y = [hauteur(:) hauteur(:)] ;
    Y = Y' ;

% plot of calculated and measured values of 36Cl with associated error bars
    hold on;
    plot(Nf,hauteur/100,'bo',cl36AMS,hauteur/100,'k.',X,Y/100,'k-') 
    
    hp = ones(N_eq,1)*get(gca,'XLim') ;
    hp = hp' ;
    vp = cumsum(fliplr(slip)) ;
    vp = [vp' vp'] ;

    if any(age == 0) % in case of samples coming from the buried part of the scarp
        vp = vp' - slip(end) ;
    else
        vp = vp' ;
    end

% plot of segments limits in z
    plot(hp,vp/100,'b-') 
    
%hold off;

% proportion of points where |calc - ams| <= sigma
    i = find(abs(cl36AMS - Nf) <= sig_cl36AMS) ;
    i = fix(100*length(i)/m) ;

% titles and legends
    str_age = ['T (ka) = ',num2str(age/1000,'%5.1f')] ;
    str_slip = ['; S (m) = ',num2str(slip/100,'%5.2f'),'; P (ka) = ',int2str(preexp/1000)] ;
    str_geom = ['; a = ',int2str(alpha),'; b = ',int2str(beta),'; c = ',int2str(gamma)] ;
    str_rho = ['; rrock = ',num2str(rho_rock),'; rcoll = ',num2str(rho_coll)] ;
    str_rms = ['; RMSw = ',int2str(fix(rmsw)),'; Chi2 = ',int2str(chi_square),'; i = ',int2str(i),'; AICc = ',int2str(fix(aicc))] ;

    title([str_age str_slip str_geom str_rho str_rms])
    xlabel('36Cl (at.g-1)')
    ylabel('Height on fault scarp (m)')
% axis([0.5e5 5e5 -4 12]) % axis for MA3 site (Magnola fault, Italy)
end

disp([num2str(preexp), ' ', num2str(age), ' ', sprintf('%0.15f', rmsw)]);

%----------- AICC FUNCTION ---------
function aicc = ak(measurements,calculations,K)

n = length(measurements) ;
aicc = sum((measurements - calculations).^2) ;
aicc = n*log(aicc/n) + (2*n*K)/(n - K - 1) ;
%-----------------------------------


