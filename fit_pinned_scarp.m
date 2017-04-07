function age = fit_pinned_scarp(data,coll,max_age, age_step, slip,preexp,EL,epsilon,alpha,beta,gamma,Hfinal,rho_coll,rho_rock, graphit)
% Another method for determining the ages of the events. The first age is
% pinned using fit first
% be supplied to pin the model so that it fits for all events. 

[first_age, minRMSw ] = fit_first(true, max_age, age_step, data, coll, slip, preexp, EL, epsilon, alpha, beta, gamma, Hfinal, rho_coll, rho_rock);

N_eq = length(slip);
age = zeros(size(slip));

age(1) = first_age;
for quake = 2:N_eq

    quake_slip = [ slip(1:quake - 1) sum(slip(quake:N_eq)) ];
    quake_age = [age(1:quake - 1) age(quake - 1)];
    start_age = age(quake - 1);
    
    min_rmsw = 0;
    best_age = start_age;
    for current = start_age : -age_step : 0
        
        quake_age(quake) = current;
        
        [rmsw, chi_square, aicc, i] = correlate_scarp(data,coll,quake_age,quake_slip,preexp,EL,epsilon,alpha,beta,gamma,Hfinal,rho_coll,rho_rock);    
        
        disp(['quake ', num2str(quake), ' age: ', num2str(current), ' rmsw: ', num2str(rmsw), ' chi: ', num2str(chi_square), ' aicc: ', num2str(aicc), ' i: ', num2str(i)]);
        
        if min_rmsw == 0 || rmsw < min_rmsw
            min_rmsw = rmsw;
            best_age = current;
        else
            age(quake) = best_age;
            break;
        end
    end
 
end

if graphit
   
    modelscarp(data, coll, age, slip, preexp, EL, epsilon, alpha, beta, gamma, Hfinal, rho_coll, rho_rock)
    
end


