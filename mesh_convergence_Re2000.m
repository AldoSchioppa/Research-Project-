close all; clc; clear variables; 

% Reference temperature
T_0=273.15; %K

% Initial global cell size
dx=0.5; %m

% Min cell size
min_dx=0.00125;

% Target CFL for the finest grid
Cfl=35;

%% Input 

% Reference length
c=1; %m

% Freestream temperature
T_inf=275.15; %K

% Reynolds number 
Re=2000;

% Mach number
M=0.85;

% Ideal gas model
gamma=1.4;
R=287;

% Dominant Str
St=[1, 0.07];

%% Outputs 

% Domain length 
L=20*2*c;

% Freestream viscosity
mu_inf=1.716e-05*(T_inf/T_0)*((T_0+110.4)/(T_inf+110.4));

% Speed of sound 
a_inf=(gamma*R*T_inf)^0.5;

% Freestream velocity
V_inf=M*a_inf;

% Freestream density 
rho_inf=Re*mu_inf/(V_inf*c);

% Freestream pressure
p_inf=rho_inf*R*T_inf;

% Vector of freestream condition
initial_conditions=[rho_inf,p_inf,mu_inf,T_inf,V_inf,a_inf];

% Convective time
conv_dt=c/V_inf;

% Dominant frequency
F=St*V_inf/c;

% Vector of dominant time step and the convective time added in
Dt=1./F;
Dt(end+1)=conv_dt;

% Division by 50 for first time step approximation
Dt=Dt/50;

% Most restrictive time step
dt=min(Dt);

% Max CFL (for the finest mesh)
CFL=dt*(V_inf+a_inf)/min_dx;

% New dt after CFL correction
dt=Cfl*min_dx/(V_inf+a_inf);

% Minimum transient time
tr=L/V_inf;

%%

% Time step vector for time convergence. Progression of 1/5.
Dt=[dt,dt/5,dt/5/5];

% Global cell size for mesh convergence. Progression of 1/2.
Dx=[dx,dx/2,dx/2/2];

% Number of time step for transient time
n_tr=tr./Dt;

%% Folders containing data
STRINGS={'dx1','dx2','dx3'; ...
         'dt1','dt2','dt3'
        };

FOLDERS={'Frequency_Analysis', ...
        'Pressure_Distribution', ...
        'Skin_Friction_Coefficient', ...
        'Von_Kármán_Shedding'};
%% Time convergence - Benchmark - Kojima et al.
close all

% Store-vectors for checking on convergence in time and space
STR={};
CL={};

for ti=1:1:length(Dt)

% LIFT
filename=strcat('D:\ResearchProject_Major\Data\Re2000\',FOLDERS{1},'\M0.85_',STRINGS{1,1},'_',STRINGS{2,ti},'_L.csv');
% Data extraction
data=readmatrix(filename);
L=data(:,2); 
% Transient estimation
transient=n_tr(ti); % To be modified.
% Erase transient data
if transient~=0
    Cl=L(transient:end)/(0.5*rho_inf*V_inf^2*c);
end

% The FFT of a time-dependent signal is a particular case of a 1D input to
% the SPOD function. In particular, uniquely passing the dt the y(t), onto
% the function as inputs retrieves the frequency spectrum.
[L_Cl,~,f_Cl] = spod(Cl,[],[],[],Dt(ti));  
% Plot
figure 
loglog(f_Cl*c/V_inf,L_Cl(:,1),'k','LineWidth',2); hold on; grid on; xlim([min(f_Cl*c/V_inf) max(f_Cl*c/V_inf)]); 
xlab=xlabel("$St [-]$"); ylab=ylabel("\textbf{FFT of $C_l$}"); 
set(xlab,'Interpreter', 'latex','FontSize',16); set(ylab,'Interpreter','latex','Rotation',90,'FontSize',14); 
set(gca,'FontSize',11);
title('Lift Frequency Spectrum')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAG 
filename=strcat('D:\ResearchProject_Major\Data\Re2000\',FOLDERS{1},'\M0.85_',STRINGS{1,1},'_',STRINGS{2,ti},'_D.csv');
% Data extraction
data=readmatrix(filename);
D=data(:,2);
% Transient estimation
transient=n_tr(ti); % To be modified.
% Erase transient data
if transient~=0
    Cd=D(transient:end)/(0.5*rho_inf*V_inf^2*c);
end

% Sampling time-step (sampling rate).
% dt_s=mean(diff(t)); 

% Frequency of acqusition
% f_s=1/dt_s;      

% Compute the FFT
[L_Cd,~,f_Cd] = spod(Cd,[],[],[],Dt(ti));  
% Plot
figure 
loglog(f_Cd*c/V_inf,L_Cd(:,1),'k','LineWidth',2); hold on; grid on; grid minor; xlim([10^-1 max(f_Cd*c/V_inf)-1]);
xlab=xlabel("$St [-]$"); ylab=ylabel("\textbf{FFT of $C_d$}"); 
set(xlab,'Interpreter', 'latex','FontSize',16); set(ylab,'Interpreter','latex','Rotation',90,'FontSize',14); 
set(gca,'FontSize',11);
title('Drag Frequency Spectrum')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VON KÁRMÁN
filename=strcat('D:\ResearchProject_Major\Data\Re2000\',FOLDERS{4},'\M0.85_',STRINGS{1,1},'_',STRINGS{2,ti},'_V.csv');
% Data extraction
data=readmatrix(filename);
V_W=data(:,2); 
% Transient estimation
transient=n_tr(ti); % To modified
% Erase transient data
if transient~=0
    M_W=V_W(transient:end)/a_inf;
end
%Compute the FFT
[L_W,~,f_W] = spod(M_W,[],[],[],Dt(ti));  
% Plot
hFig=figure;
movegui(hFig,"north")
loglog(f_W*c/V_inf,L_W,'LineWidth',2); hold on; grid on; xlim([min(f_W*c/V_inf) max(f_W*c/V_inf)-1]); axis square;
xlab=xlabel("$St [-]$"); ylab=ylabel("\textbf{FFT of velocity}"'); 
set(xlab,'Interpreter', 'latex','fontweight','bold','FontSize',16); set(ylab,'Interpreter','latex','fontweight','bold','Rotation',90,'FontSize',14); 
set(gca,'FontSize',11);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pressure_Distrubution
% filename=strcat('D:\ResearchProject_Major\Data\Re2000\',FOLDERS{2},'\M0.85_',STRINGS{1,1},'_',STRINGS{2,ti},'_p_mean.csv');
% hFig=figure; % Defining image-object
% % hFig.Name=strcat('Pressure distribution. BASE CELL SIZE=',num2str(wake_cell_size(3)),' m. TIME STEP=',num2str(time_step(3)),' s.'); % Image header
% hFig.NumberTitle = 'off'; % Number of figure off.
% movegui(hFig,"east"); % Position on the screen.
% cpplot(filename,p_inf,c,rho_inf,V_inf,true,'both'); hold on; axis square;
% %ZoomPanCallback(hFig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Skin_Friction_Coefficient
% filename=strcat('D:\ResearchProject_Major\Data\Re2000\',FOLDERS{3},'\M0.85_',STRINGS{1,1},'_',STRINGS{2,ti},'_f_mean.csv');
% hFig=figure;
% hFig.NumberTitle='off';
% movegui(hFig,"west");
% cfplot(filename,c,rho_inf,V_inf,true); axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistical convergence
% hFig=figure;
% hFig.NumberTitle='off';
% movegui(hFig,"center");
% filename=strcat('D:\ResearchProject_Major\Data\Re2000\',FOLDERS{1},'\M0.85_',STRINGS{1,1},'_',STRINGS{2,ti},'_L.csv');
% % Transient estimation
% transient=8*n_tr(ti);
% statconvergenceplot(filename,50,transient,V_inf,rho_inf,c); axis square;
% ZoomPanCallback(hFig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CL{ti}=mean(Cl);
end

% Strouhal estraction 
format long;
% Vector of Strouhal associated to the maxima
STR{1}=0.99593;
STR{2}=1.002;
STR{3}=0.9868;

%% Mesh convergence - Benchmark - Kojima et al.

close all

% Store-vectors for checking on convergence in time and space
STR={};
CL={};

for ix=1:1:length(Dx)

% LIFT
filename=strcat('D:\ResearchProject_Major\Data\Re2000\',FOLDERS{1},'\M0.85_',STRINGS{1,ix},'_',STRINGS{2,1},'_L.csv');
% Data extraction
data=readmatrix(filename);
L=data(:,2); 
% Transient estimation
transient=n_tr(1); % To be modified.
% Erase transient data
if transient~=0
    Cl=L(transient:end)/(0.5*rho_inf*V_inf^2*c);
end

% The FFT of a time-dependent signal is a particular case of a 1D input to
% the SPOD function. In particular, uniquely passing the dt the y(t), onto
% the function as inputs retrieves the frequency spectrum.
[L_Cl,~,f_Cl] = spod(Cl,[],[],[],Dt(1));  
% Plot
figure 
loglog(f_Cl*c/V_inf,L_Cl(:,1),'k','LineWidth',2); hold on; grid on; xlim([min(f_Cl*c/V_inf) max(f_Cl*c/V_inf)]); 
xlab=xlabel("$St [-]$"); ylab=ylabel("\textbf{FFT of $C_l$}"); 
set(xlab,'Interpreter', 'latex','FontSize',16); set(ylab,'Interpreter','latex','Rotation',90,'FontSize',14); 
set(gca,'FontSize',11);
title('Lift Frequency Spectrum')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DRAG 
filename=strcat('D:\ResearchProject_Major\Data\Re2000\',FOLDERS{1},'\M0.85_',STRINGS{1,ix},'_',STRINGS{2,1},'_D.csv');
% Data extraction
data=readmatrix(filename);
D=data(:,2);
% Transient estimation
transient=n_tr(1); % To be modified.
% Erase transient data
if transient~=0
    Cd=D(transient:end)/(0.5*rho_inf*V_inf^2*c);
end

% Sampling time-step (sampling rate).
% dt_s=mean(diff(t)); 

% Frequency of acqusition
% f_s=1/dt_s;      

% Compute the FFT
[L_Cd,~,f_Cd] = spod(Cd,[],[],[],Dt(1));  
% Plot
figure 
loglog(f_Cd*c/V_inf,L_Cd(:,1),'k','LineWidth',2); hold on; grid on; grid minor; xlim([10^-1 max(f_Cd*c/V_inf)-1]);
xlab=xlabel("$St [-]$"); ylab=ylabel("\textbf{FFT of $C_d$}"); 
set(xlab,'Interpreter', 'latex','FontSize',16); set(ylab,'Interpreter','latex','Rotation',90,'FontSize',14); 
set(gca,'FontSize',11);
title('Drag Frequency Spectrum')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VON KÁRMÁN
filename=strcat('D:\ResearchProject_Major\Data\Re2000\',FOLDERS{4},'\M0.85_',STRINGS{1,ix},'_',STRINGS{2,1},'_V.csv');
% Data extraction
data=readmatrix(filename);
V_W=data(:,2); 
% Transient estimation
transient=n_tr(1); % To modified
% Erase transient data
if transient~=0
    M_W=V_W(transient:end)/a_inf;
end
%Compute the FFT
[L_W,~,f_W] = spod(M_W,[],[],[],Dt(1));  
% Plot
hFig=figure;
movegui(hFig,"north")
loglog(f_W*c/V_inf,L_W,'LineWidth',2); hold on; grid on; xlim([min(f_W*c/V_inf) max(f_W*c/V_inf)-1]); axis square;
xlab=xlabel("$St [-]$"); ylab=ylabel("\textbf{FFT of velocity}"'); 
set(xlab,'Interpreter', 'latex','fontweight','bold','FontSize',16); set(ylab,'Interpreter','latex','fontweight','bold','Rotation',90,'FontSize',14); 
set(gca,'FontSize',11);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pressure_Distrubution
% filename=strcat('D:\ResearchProject_Major\Data\Re2000\',FOLDERS{2},'\M0.85_',STRINGS{1,1},'_',STRINGS{2,ti},'_p_mean.csv');
% hFig=figure; % Defining image-object
% % hFig.Name=strcat('Pressure distribution. BASE CELL SIZE=',num2str(wake_cell_size(3)),' m. TIME STEP=',num2str(time_step(3)),' s.'); % Image header
% hFig.NumberTitle = 'off'; % Number of figure off.
% movegui(hFig,"east"); % Position on the screen.
% cpplot(filename,p_inf,c,rho_inf,V_inf,true,'both'); hold on; axis square;
% %ZoomPanCallback(hFig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Skin_Friction_Coefficient
% filename=strcat('D:\ResearchProject_Major\Data\Re2000\',FOLDERS{3},'\M0.85_',STRINGS{1,1},'_',STRINGS{2,ti},'_f_mean.csv');
% hFig=figure;
% hFig.NumberTitle='off';
% movegui(hFig,"west");
% cfplot(filename,c,rho_inf,V_inf,true); axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistical convergence
% hFig=figure;
% hFig.NumberTitle='off';
% movegui(hFig,"center");
% filename=strcat('D:\ResearchProject_Major\Data\Re2000\',FOLDERS{1},'\M0.85_',STRINGS{1,1},'_',STRINGS{2,ti},'_L.csv');
% % Transient estimation
% transient=8*n_tr(ti);
% statconvergenceplot(filename,50,transient,V_inf,rho_inf,c); axis square;
% ZoomPanCallback(hFig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CL{ix}=mean(Cl);
end

% Strouhal estraction 
format long;
% Vector of Strouhal associated to the maxima
STR{1}=0;
STR{2}=0;
STR{3}=0;


%%

% SHOCK Pressure
filename=strcat('D:\ResearchProject_Major\Data\Re2000\Pressure_through_Shock\M0.89_Unsteady_Inviscuous_AoA6_p.csv');
% Data extraction
data=readmatrix(filename);
p=data(:,2); 
% Transient estimation
transient=5*n_tr(1); % To modified
% Erase transient data
if transient~=0
    p=p(transient:end)/(0.5*rho_inf*V_inf^2*c);
end
%Compute the FFT
[L_p,~,f_p] = spod(p-mean(p),[],[],[],Dt(1));  
% Plot
hFig=figure;
movegui(hFig,"north")
loglog(f_p*c/V_inf,L_p,'LineWidth',2); hold on; grid on; xlim([min(f_p*c/V_inf) max(f_p*c/V_inf)-1]); axis square;
xlab=xlabel("$St [-]$"); ylab=ylabel("\textbf{FFT of velocity}"'); 
set(xlab,'Interpreter', 'latex','fontweight','bold','FontSize',16); set(ylab,'Interpreter','latex','fontweight','bold','Rotation',90,'FontSize',14); 
set(gca,'FontSize',11);

