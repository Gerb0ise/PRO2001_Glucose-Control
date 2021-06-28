%% Generating parameters k6 k7 k8 

clearvars; clearvars -global; close all; clc;

p0= [5.82e-1 2.20e-2 4.71];         % initial parameter values for k6 k7 k8
multiplier = 100;                   % bounds multiplier 
lbP = p0/multiplier;                % lower bound of parameter set
hbP = p0*multiplier;                % upper bound of the parameter set
total_k = 100;                      % amount of points in each parameter vector
    
k6 = exp(linspace(log(lbP(1)) , log(hbP(1)) , total_k));
k7 = exp(linspace(log(lbP(2)) , log(hbP(2)) , total_k));
k8 = exp(linspace(log(lbP(3)) , log(hbP(3)) , total_k));


%% Necessary variables for the PID model

timespan = 0:180;                   % timespan for the ODE in minutes

p = [1.35e-2                        % k1
    6.33e-1                         % k2
    5.00e-5                         % k3
    1.00e-3                         % k4
    3.80e-3                         % k5
    5.82e-1                         % k6
    2.20e-2                         % k7
    4.71                            % k8
    1.08e-2                         % k9
    2.60                            % k10
    1.35                            % sigma /index 11
    0.63                            % Km /index 12
    4.9                             % Gb /index 13
    6.7];                           % Iplb /index 14  

c = struct();                       % Constants
    c.f      = 0.005551;            % f [mmol/mg], must be equal to mgdL_to_mmolL /10 in load_data.m
    c.vg     = 17/70;               % vg [L/kg]
    c.gbliv  = 0.043;               % gbliv [mmol/L]
    c.beta   = 1;                   % beta [(mmol/L)/(microU/mL)]
    c.taui   = 31;                  % taui [min]
    c.taud   = 3;                   % taud [min]
    c.vi     = 13/70;               % distribution volume of insulin per kg bodymass, vi [L/kg]
    c.Gthpl  = 9;                   % Gthpl [mmol/L]
    c.t_integralwindow = 30;        % Lower bound of moving time window of Gint
    c.c1     = 0.1;                 % c1 [1/min](previously k8)
    c.c2     = c.gbliv.*(p(12)+p(13))/p(13) - p(5)*c.beta*p(14);          % c2 [mmol/L/min] = gbliv * (Km+Gb)/Gb - k5*beta*Ibpl; constant term in gnon-it
    c.c3     = p(7).*p(13)/(c.beta*c.taui.*p(14)).*c.t_integralwindow;    % c3 [1/min] = k7 * Gb/(beta*tau_i*Ibpl).*t_integralwindow ; constant term in iliv

input = struct();
    input.D  = 75000;               % D (total amount of carbohydrates ingested) [mg]
    input.t_meal_start = 0;         % starting time of the meal [min, counted from 0 = 0:00]
    input.Mb = 70;                  % Mb (body weight) [kg]

x0 = [0                             % Mg: glucose mass in gut [mg]
      p(13)                         % Gpl: plasma glucose concentration [mmol/L]
      p(14)                         % Ipl: plasma insulin concentration [mU/L]
      0];

  
%% Calculating peaks for a set of k6 k7 and k8 parameters

ODE_optionsG = odeset('RelTol',1e-5);        
peaks = zeros(total_k^3,1);

[xx,yy,zz] = meshgrid(k6,k7,k8);    % creates a mesh for every possible combination of k6 k7 and k8 values
x = reshape(xx,[],1);
y = reshape(yy,[],1);
z = reshape(zz,[],1);

parfor i = 1:total_k^3
    tmpp = p;
    tmpp(6) = x(i);
    tmpp(7) = y(i);
    tmpp(8) = z(i);

	[tmp,xmodelGID] = ode15s(@diffeq_diabetes,timespan,x0,ODE_optionsG,tmpp,c,input);
	peaks(i) = numel(findpeaks(xmodelGID(:,3)'));               % counts the peaks for insulin
      
    if mod(i, 1000) == 0
        disp(strcat(num2str(i/(total_k^3)*100,'%.1f'),' %'))    % displays advancement percentage
    end

end

fpath = 'C:\Users\Gerbie\Google Drive\Current courses\PRO2001 - Glucose level\edes_matlab\data';
save(fullfile(fpath,"data_n" + total_k + "_m"+ multiplier));


%% Calculating peaks for a set of k6 and k8 parameters

ODE_optionsG = odeset('RelTol',1e-5);

[xx,yy] = meshgrid(k6,k8);          % creates a set of coordinates for every possible combination of k6 and k8 values
x = reshape(xx,[],1);
y = reshape(yy,[],1);
peaks = zeros(total_k^2,1);

parfor i = 1:total_k^2
    tmpp = p;
    tmpp(6) = x(i);
    tmpp(8) = y(i);

	[tmp,xmodelGID] = ode15s(@diffeq_diabetes,timespan,x0,ODE_optionsG,tmpp,c,input);
	peaks(i) = numel(findpeaks(xmodelGID(:,3)'));               % counts the peaks for insulin
      
    if mod(i, 1000) == 0
        disp(strcat(num2str(i/(total_k^2)*100,'%.1f'),' %'))    % displays advancement percentage
    end

end

fpath = 'C:\Users\Gerbie\Google Drive\Current courses\PRO2001 - Glucose level\edes_matlab\data';
save(fullfile(fpath,"data_n" + total_k + "_m" + multiplier + "_no_k8"));


%% Plotting 3D mesh figure of peaks

fig1 = figure('units','centimeters','position',[0 0 24 18],'renderer','painters');
colormap(fig1,jet(max(peaks)-2));

peaks1 = peaks;
    peaks1(peaks1==0) = NaN;        % removes all datapoints with 0 peaks
    peaks1(peaks1==1) = NaN;        % removes all datapoints with 1 peaks
    peaks1(peaks1==2) = NaN;        % removes all datapoints with 2 peaks
    
s = scatter3(x,y,z,6,peaks1,'filled');
set(gca,'XScale','log','YScale','log','ZScale','log');

fs = 12;

title('Number of peaks as a function of k6 k7 and k8','FontSize',fs);
xlabel("k6",'FontSize',fs);
ylabel("k7",'FontSize',fs);
zlabel("k8",'FontSize',fs);

cb = colorbar;
cb.Label.String = 'Number of Peaks';
cb.Label.FontSize = fs;
cb.Ticks = linspace(3.5,14.5,13);
cb.TickLabels = num2cell(3:15);

[caz,cel] = view;
view(caz+90,cel)                   % rotates the plot 90°

fpath = 'C:\Users\Gerbie\Google Drive\Current courses\PRO2001 - Glucose level\edes_matlab\Figures';
print(fig1,fullfile(fpath,"Simu_peaks_rotated"),'-r600','-dpng')
print(fig1,fullfile(fpath,"Simu_peaks_rotated"),'-r600','-dsvg')


%% Animating the 3D peaks figure

[caz,cel] = view;
frames = 720;
set(gca,'projection','perspective')
axis vis3d
set(ax,'xlimmode','manual','ylimmode','manual','zlimmode','manual','xtickmode','manual','ytickmode','manual','ztickmode','manual','xticklabelmode','manual','yticklabelmode','manual','zticklabelmode','manual')
      
for i=1:frames
    view(caz+i*360/frames,cel)
    F(i) = getframe(gcf) ;
    drawnow
end

writerObj = VideoWriter('C:\Users\Gerbie\Google Drive\Current courses\PRO2001 - Glucose level\edes_matlab\Figures\Final\Peaks_animation2','MPEG-4');
writerObj.FrameRate = 30;

open(writerObj);

for i=1:length(F)
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end

close(writerObj);  
    

%% Plotting contour of peaks for k6 and k8 values only

fig1 = figure('units','centimeters','position',[0 0 24 18],'renderer','painters');
set(gcf,'renderer','painters');
colormap(fig1,flip(bone(max(peaks)+1)));

contourf(k6,k8,reshape(peaks,total_k,total_k),15)
set(gca,'XScale','log','YScale','log','ZScale','log');

fs = 12;

title('Number of peaks as a function of k6 and k8','FontSize',fs);
xlabel("k6",'FontSize',fs);
ylabel("k8",'FontSize',fs);
zlabel("Peaks",'FontSize',fs);

cb = colorbar;
cb.Label.String = 'Number of Peaks';
cb.Label.FontSize = fs;
cb.Ticks = linspace(0.5,13.125-0.5,15);
cb.TickLabels = num2cell(0:14);

fpath = 'C:\Users\Gerbie\Google Drive\Current courses\PRO2001 - Glucose level\edes_matlab\Figures';
print(fig1,fullfile(fpath,"Simu_peaks_noK8"),'-r600','-dpng')
print(fig1,fullfile(fpath,"Simu_peaks_noK8"),'-r600','-dsvg')


%% Plotting ODE's for glucose and insulin by varying each parameter separately

p0 = [5.82e-1 2.20e-2 4.71];            % initial parameter values for k6 k7 k8
total_k = 1000;                         % amount of points in each parameter vector

k6 = exp(linspace(log(0.001) , log(5) , total_k));
k7 = exp(linspace(log(0.001) , log(2) , total_k));
k8 = exp(linspace(log(0.001) , log(5) , total_k));

gradient = parula(total_k);
ODE_optionsG = odeset('RelTol',1e-5);
k = [k6' k7' k8'];

for j = 1:3  
    
    fig1 = figure('units','centimeters','position',[0 0 22 18],'renderer','painters');
    
    if j==1, name = 'k6';
    elseif j==2, name = 'k7';
    else, name = 'k8';
    end
    
    p(6) = p0(1);
    p(7) = p0(2);
    p(8) = p0(3);

    for i = 1:total_k

        p(j+5) = k(i,j);

        [~,xmodelGID] = ode15s(@diffeq_diabetes,timespan,x0,ODE_optionsG,p,c,input);
        
        subplot(2,1,1)
        plot(timespan,xmodelGID(:,2),'color',gradient(i,:))         % glucose
        hold on
        
        subplot(2,1,2)
        plot(timespan,xmodelGID(:,3),'color',gradient(i,:))         % insulin
        hold on

    end

    fs = 12;
    
    subplot(2,1,1)
        title("Glucose level for a set of " + name,'FontSize',fs) 
        xlabel("time, min",'FontSize',fs);
        ylabel("Plasma glucose, mmol/L",'FontSize',fs);

        ax = gca;
        ax.FontSize = fs;

        cb = colorbar;
        cb.Label.String = name;
        cb.Label.FontSize = fs;
        set(gca,'colorscale','log')
        cb.TickLabels = [min(k(:,j)) max(k(:,j))];
    
    subplot(2,1,2)
        title("Insulin level for a set of " + name,'FontSize',fs) 
        xlabel("time, min",'FontSize',fs);
        ylabel("Plasma insumin, mU/L",'FontSize',fs);

        ax = gca;
        ax.FontSize = fs;

        cb = colorbar;
        cb.Label.String = name;
        cb.Label.FontSize = fs;
        set(gca,'colorscale','log')
        cb.TickLabels = [min(k(:,j)) max(k(:,j))];
        
    fpath = 'C:\Users\Gerbie\Google Drive\Current courses\PRO2001 - Glucose level\edes_matlab\Figures';
    print(fig1,fullfile(fpath,"Simu_varying_"+ name),'-r600','-dpng')
    print(fig1,fullfile(fpath,"Simu_varying_"+ name),'-r600','-dsvg')
    
end


%% Plotting ODE for glucose and insulin for one specified k6 k7 and k8 set

fig1 = figure('units','centimeters','position',[0 0 26 18],'renderer','painters');
ODE_optionsG = odeset('RelTol',1e-5);

k6a = 11.9689;
k7a = 0.000264991;
k8a = 0.580671;

p(6) = k6a;                 % k6 value
p(7) = k7a;                 % k7 value
p(8) = k8a;                 % k8 value

[~,xmodelGID1] = ode15s(@diffeq_diabetes,timespan,x0,ODE_optionsG,p,c,input);

k6b = 0.3;
k7b = 0.000264991;
k8b = 0.580671;

p(6) = k6b;                 % k6 value
p(7) = k7b;                 % k7 value
p(8) = k8b;                 % k8 value

[~,xmodelGID2] = ode15s(@diffeq_diabetes,timespan,x0,ODE_optionsG,p,c,input);

        
subplot(2,1,1)
    plot(timespan,xmodelGID1(:,2),'linewidth',2)         % glucose1
    hold on
    plot(timespan,xmodelGID2(:,2),'linewidth',2)         % glucose2
    title("Glucose level",'FontSize',fs)
    xlabel("time, min",'FontSize',fs);
    ylabel("Plasma glucose, mmol/L",'FontSize',fs);
    ax = gca;
    ax.FontSize = fs;
    
    
subplot(2,1,2)
    plot(timespan,xmodelGID1(:,3),'linewidth',2)         % insulin1
    hold on
	plot(timespan,xmodelGID2(:,3),'linewidth',2)         % insulin2
    title("Insulin level",'FontSize',fs)
    xlabel("time, min",'FontSize',fs);
    ylabel("Plasma insumin, mU/L",'FontSize',fs);
    ax = gca;
    ax.FontSize = fs;

leg = legend({"[k6 k7 k8] = ["+k6a+" "+k7a+" "+k8a+"]","[k6 k7 k8] = [  "+k6b+"       "+k7b+" "+k8b+"]"});
leg.Position = [0.676694750271788,0.458578433037973,0.228640187554124,0.0624999983345761];

fpath = 'C:\Users\Gerbie\Google Drive\Current courses\PRO2001 - Glucose level\edes_matlab\Figures';
print(fig1,fullfile(fpath,"Simu_for_data_pts"),'-r600','-dpng')
print(fig1,fullfile(fpath,"Simu_for_data_pts"),'-r600','-dsvg')


    