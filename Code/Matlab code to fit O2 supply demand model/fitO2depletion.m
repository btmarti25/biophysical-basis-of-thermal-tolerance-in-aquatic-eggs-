function [ssq] = fitO2depletionFINAL(x)
global tab j_list

dt=1/300; %hr
V_i = .2121748; %ml
V_o = 1.7 - V_i; %ml
Init_C_o = 18; % ug/ml
k=x(1); % half-saturation coefficient
k_e= x(2); % mass transfer coefficient ugO2/(h ugO2/ml)
ssq=0;
sse=0;
treatments = unique(tab.devDay);
for j = 1:length(treatments)
    stab = tab(tab.devDay==treatments(j_list(j)),:);
    mass = exp(-12.630 + 2.692 *log(stab.devDay(1))+ 2.838* log(12));
    mo2=  x(3)*mass^1*exp(x(4)*stab.temp(1));
    DO_points = [10,20,30,40,50,60,70,80,90,100];
    % calc initital internal DO concentration
    Init_C_i = (-k * k_e - mo2 + k_e *Init_C_o + sqrt(k^2 * k_e^2 + (mo2 - k_e *Init_C_o)^2 + 2* k* k_e* (mo2 + k_e *Init_C_o)))/(2 *k_e);
    C_i=Init_C_i;
    C_o = Init_C_o;
    t=1;
    while C_o(t) > .8
        f= C_i(t)^(1)/(k + C_i(t)^(1)); % Michaelis-Menton function
        demand(t) =  mo2 * f; % oxygen condumption ug/hr
        in_flux = k_e * (C_o(t) - C_i(t)); % ug/hr
        dC_i = (in_flux - demand(t))/V_i;
        dC_o(t) = -in_flux/V_o; % dC_o is in units of  ug/(hr ml)
        C_i(t+1)=C_i(t)+dC_i*dt;
        C_o(t+1)=C_o(t)+dC_o(t)*dt;
        t = t+1;
    end
        
    for i = 1:length(DO_points)
        % find time with the modelled C_0 closest to the comaprison C_0 levels
        [minValue,closestIndex] =  min(abs(C_o(1:(end-1))-DO_points(i)* stab.DO_satconc(1)/100))  ;
        pred_mo2 =  V_o * -dC_o(closestIndex); % convert dC_o to ug/hr by multiplying by the well volume
        obs_mo2 = stab.mo2(stab.DOsat==DO_points(i),:); % find the observed mo2s for this DO level
        if isempty(obs_mo2) == 0
            temp = sum(((pred_mo2-obs_mo2)).^2); % calc sum of squares for these observations
        end
        ssq = ssq + temp; % add ssq
    end
end
ssq;
end