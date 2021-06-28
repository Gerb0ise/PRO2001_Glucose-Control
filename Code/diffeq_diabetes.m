function dxdt = diffeq_diabetes(t,x,p_internal,c,input)
%diffeq_diabetes - Differential equations of the E-DES model.
%
%   dxdt = diffeq_diabetes(t,x,p_internal,c,input)
%
%Inputs
% t:    current time point
% x:    vector with state variables (at current time t)
% p_internal: vector of model parameters
% c:    struct of model constants
%input: struct with model inputs (time of the meal, carbohydrate content of
%       the meal, timing and dose of exogenous insulin 
%
%Output
% dxdt: vector with derivatives of state variables


% -- state variables
Mg    = x(1); %Mg: glucose mass in gut [mg]
Gpl   = x(2); %Gpl: plasma glucose concentration [mmol/L]
Ipl   = x(3); %Ipl: plasma insulin concentration [mU/L]
Gint  = x(4); %Gint: integrated plasma glucose increase (int (Gpl)) [mmol/L]

% -- model constants
f      = c.f;       %f [mmol/mg], conversion factor from mg/dl to mmol/l
vg     = c.vg;      %vg [L/kg], glucose distribution volume
gbliv  = c.gbliv;   %gbliv [mmol/L], basal liver glucose output
beta   = c.beta;    %beta [(mmol/L)/(microU/mL)], conversion from microU/ml to mmol/l
taui   = c.taui;    %taui [min], integration time constant
taud   = c.taud;    %taud [min], differentation time constant
Gthpl  = c.Gthpl;   %Gthpl [mmol/L], threshold for renal glucose removal
c1     = c.c1;      %c1 [1/min](previously k8); constant term in gren
c2     = c.c2;      %c2 [mmol/L/min] = gbliv * (Km+Gb)/Gb - k5*beta*Ibpl; constant term in gnon-it
c3     = c.c3;      %c3 [1/min] = k7 * Gb/(beta*tau_i*Ibpl).*t_integralwindow ; constant term in iliv
t_integralwindow = c.t_integralwindow; %Lower bound of moving time window of Gint

% -- model input
D       = input.D; %D (total amount of carbohydrates ingested) [mg] 

t_meal_start = input.t_meal_start;% starting time of the meal [min, counted from 0 = 0:00]
Mb      = input.Mb; %Mb (body weight) [kg]

% -- model parameters
k1    = p_internal(1); %[1/min]
k2    = p_internal(2); %[1/min]
k3    = p_internal(3); %[1/min]
k4    = p_internal(4); %[1/min]
k5    = p_internal(5); %[1/min]
k6    = p_internal(6); %[1/min]
k7    = p_internal(7); %[1/min]
k8    = p_internal(8); %[1/min]
k9   = p_internal(9);  %[1/min]
sigma = p_internal(11); %[-]
KM    = p_internal(12); %[mmol/l]
Gb    = p_internal(13); %[mmol/l]
Iplb  = p_internal(14); %[microU/ml]

%% ------------ ODEs --------------------------------------------------------
%% dMg/dt -- glucose in gut
if t>t_meal_start
    mgmeal  = sigma .* k1.^sigma .* (t-t_meal_start).^(sigma-1) .* exp(-(k1.*(t-t_meal_start)).^sigma) .* D;
else
    mgmeal = 0;
end
mgpl    = k2.*Mg;
dMg_dt   = mgmeal - mgpl;


%% dGpl/dt -- glucose in plasma
gliv    = gbliv - k3.*(Gpl-Gb) - k4.*beta.*(Ipl-Iplb);
ggut    = k2 .*(f/(vg*Mb)) .* Mg;
gnonit  = c2 .* (Gpl./(KM+Gpl));
git     = k5.*beta.*Ipl .* (Gpl./(KM+Gpl));
if Gpl > Gthpl
    gren  = c1./(vg*Mb) .* (Gpl - Gthpl);
else
    gren = 0;
end
dGpl_dt = gliv + ggut - gnonit - git - gren;


%% dIpl/dt -- insulin in plasma
% ipnc
global t_saved Gpl_saved
       Gpl_saved = 4.9;
       t_saved = 0;
       
t_lowerbound = t - t_integralwindow;
if (t > t_integralwindow) && (length(t_saved)>1) && (length(t_saved) == length(Gpl_saved))
    Gpl_lowerbound = interp1(t_saved,Gpl_saved,t_lowerbound, 'spline');
else
    Gpl_lowerbound = Gpl_saved(1);  % is called when t < t_integralwindow, or if there is no saved step yet (steps are only saved at pre-defined time points)
end
%Gpl_lowerbound = interp1(t_saved,Gpl_saved,t_lowerbound, 'spline');
dGint_dt = (Gpl-Gb) - (Gpl_lowerbound-Gb);
ipnc    = (beta.^-1).*(k6.*(Gpl-Gb) + (k7/taui).*Gint + (k7/taui).*Gb.*t_integralwindow + (k8.*taud).*dGpl_dt);

%iliv & iif
iliv    = c3.*Ipl;
iif     = k9.*(Ipl-Iplb);

dIpl_dt = ipnc - iliv - iif;


% -- catch an error where the timestep of the integration becomes too small
% MINSTEP = 1e-10; %Minimum step
% 
% persistent tprev
% 
% if isempty(tprev)
%     tprev = -inf;
% end
% timestep = t - tprev;
% tprev = t;
% 
% if (timestep > 0) && (timestep < MINSTEP)
%    error(['Stopped. Time step is too small: ' num2str(timestep)])
% end


%% -- output differential equations
dxdt = [dMg_dt
        dGpl_dt
        dIpl_dt
        dGint_dt];