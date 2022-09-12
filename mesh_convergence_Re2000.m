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

for ti=1:1:length(Dt)-1

% LIFT
filename=strcat('D:\ResearchProject_Major\Data\Re2000\',FOLDERS{1},'\M0.85_',STRINGS{1,1},'_',STRINGS{2,ti},'_L.csv');
% Data extraction
data=readmatrix(filename);
L=data(:,2); 
% Transient estimation
transient=8*n_tr(ti);
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
transient=8*n_tr(ti);
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
transient=6*n_tr(1);
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
STR{3}=0;

return 
%% Time convergence - dx(1) - dt() 
close all

jj=2; % FFT length index. Defining FFT resoultion. Often times equal to the number of samples.

% Transient estimation (1 or 2 convective time w.r.t the computational domain). Given as a number of time steps.
conv_time=domain_length/u_inf; % Convective time
transient=conv_time/time_step(2); % Number of time steps

% If necessary, cut the signal from right for graphics purpose
stop=running_time/time_step(2);

% Lift FFT
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Frequency Analysis\Re2000_dx1_Dt2_L';
[Cl,FFT_Cl,nfft,f_s]=FFT(filename,f_max,3*transient,0,c,u_inf,rho_inf,a_inf,p_inf,'$C_{l}$',1);

% Manual dominant frequency estraction
format long;
% Estraction to the dominant frequency
F=FFT_Cl{jj}{1};
A=FFT_Cl{jj}{2};
f_d=F(5); % The index is to modify 
max_A=A(5);

% FFT plot
fftplot(filename,3*transient,0,nfft(jj),f_d,max_A,f_s,c,u_inf); 

% Wake velocity fluctuations FFT
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Von Kármán Shedding\dx1_Dt2_V_Fluctuations_Wake';
[~,FFT_V,nfft,f_s]=FFT(filename,f_max,transient,0,c,u_inf,rho_inf,a_inf,p_inf,'$\|{\mathbf{M}_y}\|$',3); % To modify after the first call 

% Manual dominant frequency estraction
format long;
% Estraction to the dominant frequency
F=FFT_V{jj}{1};
A=FFT_V{jj}{2};
f_d=F(43); % The index is to modify 
max_A=A(43);

% FFT plot
fftplot(filename,transient,0,nfft(jj),f_d,max_A,f_s,c,u_inf);

% Cp plot
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Pressure Distribution\dx1_Dt2_Mean_p';
hFig=figure; % Defining image-object
hFig.Name=strcat('Pressure distribution. BASE CELL SIZE=',num2str(wake_cell_size(3)),' m. TIME STEP=',num2str(time_step(2)),' s.'); % Image header
hFig.NumberTitle = 'off'; % Number of figure off.
movegui(hFig,"east"); % Position on the screen.
cpplot(filename,p_inf,c,rho_inf,u_inf,true,'both'); hold on; 
ZoomPanCallback(hFig);

% Cf plot
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Skin Friction Coefficient\dx1_Dt2_Mean_Cf';
hFig=figure;
hFig.Name=strcat('Friction Coefficient. BASE CELL SIZE=',num2str(wake_cell_size(3)),' m. TIME STEP=',num2str(time_step(2)),' s.');
hFig.NumberTitle='off';
movegui(hFig,"west");
cfplot(filename,c,rho_inf,u_inf,true); hold on; 

% Statistical convergence
hFig=figure;
hFig.Name=strcat('Free stream. BASE CELL SIZE=',num2str(wake_cell_size(3)),' m. TIME STEP=',num2str(time_step(2)),' s.');
hFig.NumberTitle='off';
movegui(hFig,"west");
statconvergenceplot(Cl,50);
ZoomPanCallback(hFig);

% Strouhal estraction 
format long;
% Vector of Strouhal associated to the maxima
STR=FFT_Cl{jj}{3};
% Store Str to check convergence in time
STR_t(2)=STR(5);
CL_t(2)=mean(Cl(:,1));
%% Re2000_60k_Dt3_FreeStream
close all

jj=2; % FFT length index. Defining FFT resoultion. Often times equal to the number of samples.

% Transient estimation (1 or 2 convective time w.r.t the computational domain). Given as a number of time steps.
conv_time=domain_length/u_inf; % Convective time
transient=conv_time/time_step(1); % Number of time steps

% If necessary, cut the signal from right for graphics purpose
% stop=running_time/time_step(1);

% Lift FFT
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Frequency Analysis\Re2000_dx1_Dt3_L';
[Cl,FFT_Cl,nfft,f_s]=FFT(filename,f_max,3*transient,0,c,u_inf,rho_inf,a_inf,p_inf,'$C_{l}$',1);

% Manual dominant frequency estraction
format long;
% Estraction to the dominant frequency
F=FFT_Cl{jj}{1};
A=FFT_Cl{jj}{2};
f_d=F(4); % The index is to modify 
max_A=A(4);

% FFT plot
fftplot(filename,3*transient,0,nfft(jj),f_d,max_A,f_s,c,u_inf); 

% Wake velocity fluctuations FFT
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Von Kármán Shedding\dx1_Dt3_V_Fluctuations_Wake';
[~,FFT_V,nfft,f_s]=FFT(filename,f_max,transient,0,c,u_inf,rho_inf,a_inf,p_inf,'$\|{\mathbf{M}_y}\|$',3); % To modify after the first call 

% Manual dominant frequency estraction
format long;
% Estraction to the dominant frequency
F=FFT_V{jj}{1};
A=FFT_V{jj}{2};
f_d=F(41); % The index is to modify 
max_A=A(41);

% FFT plot
fftplot(filename,transient,0,nfft(jj),f_d,max_A,f_s,c,u_inf);

% Cp plot
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Pressure Distribution\dx1_Dt3_Mean_p';
hFig=figure; % Defining image-object
hFig.Name=strcat('Pressure distribution. BASE CELL SIZE=',num2str(wake_cell_size(3)),' m. TIME STEP=',num2str(time_step(1)),' s.'); % Image header
hFig.NumberTitle = 'off'; % Number of figure off.
movegui(hFig,"east"); % Position on the screen.
cpplot(filename,p_inf,c,rho_inf,u_inf,true,'both'); hold on; 
ZoomPanCallback(hFig);

% Cf plot
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Skin Friction Coefficient\dx1_Dt3_Mean_Cf';
hFig=figure;
hFig.Name=strcat('Friction Coefficient. BASE CELL SIZE=',num2str(wake_cell_size(3)),' m. TIME STEP=',num2str(time_step(1)),' s.');
hFig.NumberTitle='off';
movegui(hFig,"west");
cfplot(filename,c,rho_inf,u_inf,true); hold on; 

% Statistical convergence
hFig=figure;
hFig.Name=strcat('Free stream. BASE CELL SIZE=',num2str(wake_cell_size(3)),' m. TIME STEP=',num2str(time_step(1)),' s.');
hFig.NumberTitle='off';
movegui(hFig,"west");
statconvergenceplot(Cl,50);
ZoomPanCallback(hFig);

% Strouhal estraction 
format long;
% Vector of Strouhal associated to the maxima
STR=FFT_Cl{jj}{3};
% Store Str to check convergence in time
STR_t(3)=STR(4);
CL_t(3)=mean(Cl(:,1));
%% Strouhal vs Time-step
close all;

roc=zeros(1,length(time_step));
for jj=1:length(time_step)-1
    roc(jj+1)=(abs((STR_t(jj+1)-STR_t(jj)))/STR_t(jj))*100;
end    

figure('Name', 'Time convergence.', 'NumberTitle', 'off'); hold on;

yyaxis left
plot(fliplr(time_step)*u_inf/c, STR_t,'.-','markersize',16); grid on; grid minor; axis equal;
xlab_left=xlabel("$ \Delta t \frac{u_{\infty}}{c}$"); ylab_left=ylabel("$St$",'HorizontalAlignment','right'); 
set(gca,'FontSize',10);
set(ylab_left,'Interpreter','Latex','Rotation',90,'FontSize', 13); 
set(xlab_left,'Interpreter','Latex','FontSize', 16); 
%title('Time convergence');
yyaxis right
plot(fliplr(time_step)*u_inf/c, roc,'.-','markersize',16); grid on; grid minor; axis equal;
ylab=ylabel("$ROC [\%]$",'HorizontalAlignment','left'); 
set(ylab,'Interpreter','Latex','Rotation',90,'FontSize', 13); 
% Change the digits size 
ax=gca;
ax.XLabel.FontSize=15;
ax.YLabel.FontSize=15;
%title('Strouhal number Vs Time step');
%% Cl vs time
close all;

roc=zeros(1,length(time_step));
for jj=1:length(time_step)-1
    roc(jj+1)=(abs((CL_t(jj+1)-CL_t(jj)))/CL_t(jj))*100;
end    

figure('Name', 'Time convergence.', 'NumberTitle', 'off'); hold on;

yyaxis left
plot(fliplr(time_step)*u_inf/c, CL_t,'.-','markersize',16); grid on; grid minor; axis equal;
xlab_left=xlabel("$ \Delta t \frac{u_{\infty}}{c}$"); ylab_left=ylabel("$C_{l}$",'HorizontalAlignment','right'); 
set(gca,'FontSize',10);
set(ylab_left,'Interpreter','Latex','Rotation',0,'FontSize', 13); 
set(xlab_left,'Interpreter','Latex','FontSize', 16); 
%title('Time convergence');
yyaxis right
plot(fliplr(time_step)*u_inf/c, roc,'.-','markersize',16); grid on; grid minor; axis equal;
ylab=ylabel("$ROC [\%]$",'HorizontalAlignment','left'); 
set(ylab,'Interpreter','Latex','Rotation',90,'FontSize', 13); 
% Change the digits size 
ax=gca;
ax.XLabel.FontSize=17;
ax.YLabel.FontSize=17;
%title('Strouhal number Vs Time step');
%% Re2000_136k_Dt2_FreeStream
close all

jj=2; % FFT length index. Defining FFT resoultion. Often times equal to the number of samples.

% Transient estimation (1 or 2 convective time w.r.t the computational domain). Given as a number of time steps.
conv_time=domain_length/u_inf; % Convective time
transient=conv_time/time_step(3); % Number of time steps

% If necessary, cut the signal from right for graphics purpose
% stop=running_time/time_step(1);

% Lift FFT
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Frequency Analysis\Re2000_dx2_Dt1_L';
[Cl,FFT_Cl,nfft,f_s]=FFT(filename,f_max,3*transient,0,c,u_inf,rho_inf,a_inf,p_inf,'$C_{l}$',1);

% Manual dominant frequency estraction
format long;
% Estraction to the dominant frequency
F=FFT_Cl{jj}{1};
A=FFT_Cl{jj}{2};
f_d=F(3); % The index is to modify 
max_A=A(3);

% FFT plot
fftplot(filename,3*transient,0,nfft(jj),f_d,max_A,f_s,c,u_inf); 

% Wake velocity fluctuations FFT
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Von Kármán Shedding\dx2_Dt1_V_Fluctuations_Wake';
[~,FFT_V,nfft,f_s]=FFT(filename,f_max,transient,0,c,u_inf,rho_inf,a_inf,p_inf,'$\|{\mathbf{M}_y}\|$',3); % To modify after the first call 

% Manual dominant frequency estraction
format long;
% Estraction to the dominant frequency
F=FFT_V{jj}{1};
A=FFT_V{jj}{2};
f_d=F(17); % The index is to modify 
max_A=A(17);

% FFT plot
fftplot(filename,transient,0,nfft(jj),f_d,max_A,f_s,c,u_inf);

% Cp plot
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Pressure Distribution\dx2_Dt1_Mean_p';
hFig=figure; % Defining image-object
hFig.Name=strcat('Pressure distribution. BASE CELL SIZE=',num2str(wake_cell_size(2)),' m. TIME STEP=',num2str(time_step(3)),' s.'); % Image header
hFig.NumberTitle='off'; % Number of figure off.
movegui(hFig,"east"); % Position on the screen.
cpplot(filename,p_inf,c,rho_inf,u_inf,false,'both'); hold on; 
% ZoomPanCallback(hFig);

% Cf plot
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Skin Friction Coefficient\dx2_Dt1_Mean_Cf';
hFig=figure;
hFig.Name=strcat('Friction Coefficient. BASE CELL SIZE=',num2str(wake_cell_size(2)),' m. TIME STEP=',num2str(time_step(3)),' s.');
hFig.NumberTitle='off';
movegui(hFig,"west");
cfplot(filename,c,rho_inf,u_inf,true); hold on; 

% Statistical convergence
hFig=figure;
hFig.Name=strcat('Free stream. BASE CELL SIZE=',num2str(wake_cell_size(2)),' m. TIME STEP=',num2str(time_step(3)),' s.');
hFig.NumberTitle='off';
movegui(hFig,"west");
statconvergenceplot(Cl,50);
ZoomPanCallback(hFig);

% Strouhal estraction 
format long;
% Vector of Strouhal associated to the maxima
STR=FFT_Cl{jj}{3};
STR_w=FFT_V{jj}{3};
% Store Str to check convergence in time
STR_s(2)=STR(3);
STR_s_w(2)=STR_w(17);
CL_s(2)=mean(Cl(:,1));
%% Re2000_374k_Dt2_FreeStream
close all

jj=2; % FFT length index. Defining FFT resoultion. Often times equal to the number of samples.

% Transient estimation (1 or 2 convective time w.r.t the computational domain). Given as a number of time steps.
conv_time=domain_length/u_inf; % Convective time
transient=conv_time/time_step(2); % Number of time steps

% If necessary, cut the signal from right for graphics purpose
% stop=running_time/time_step(1);

% Lift FFT
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Frequency Analysis\Re2000_dx3_Dt1_L';
[Cl,FFT_Cl,nfft,f_s]=FFT(filename,f_max,3*transient,0,c,u_inf,rho_inf,a_inf,p_inf,'$C_{l}$',1);

% Manual dominant frequency estraction
format long;
% Estraction to the dominant frequency
F=FFT_Cl{jj}{1};
A=FFT_Cl{jj}{2};
f_d=F(2); % The index is to modify 
max_A=A(2);

% FFT plot
fftplot(filename,3*transient,0,nfft(jj),f_d,max_A,f_s,c,u_inf); 

% Wake velocity fluctuations FFT
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Von Kármán Shedding\dx3_Dt1_V_Fluctuations_Wake';
[~,FFT_V,nfft,f_s]=FFT(filename,f_max,transient,4000,c,u_inf,rho_inf,a_inf,p_inf,'$\|{\mathbf{M}_y}\|$',3); % To modify after the first call 

% Manual dominant frequency estraction
format long;
% Estraction to the dominant frequency
F=FFT_V{jj}{1};
A=FFT_V{jj}{2};
f_d=F(5); % The index is to modify 
max_A=A(5);

% FFT plot
fftplot(filename,transient,4000,nfft(jj),f_d,max_A,f_s,c,u_inf);

% Cp plot
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Pressure Distribution\dx3_Dt1_Mean_p';
hFig=figure; % Defining image-object
hFig.Name=strcat('Pressure distribution. BASE CELL SIZE=',num2str(wake_cell_size(1)),' m. TIME STEP=',num2str(time_step(3)),' s.'); % Image header
hFig.NumberTitle='off'; % Number of figure off.
movegui(hFig,"east"); % Position on the screen.
cpplot(filename,p_inf,c,rho_inf,u_inf,true,'both'); hold on; 
ZoomPanCallback(hFig);

% Cf plot
filename='C:\Users\Asus\Documents\ResearchProject_Major\Data\Re2000_Definitive\Skin Friction Coefficient\dx3_Dt1_Mean_Cf';
hFig=figure;
hFig.Name=strcat('Friction Coefficient. BASE CELL SIZE=',num2str(wake_cell_size(1)),' m. TIME STEP=',num2str(time_step(3)),' s.');
hFig.NumberTitle='off';
movegui(hFig,"west");
cfplot(filename,c,rho_inf,u_inf,true); hold on; 

% Statistical convergence
hFig=figure;
hFig.Name=strcat('Free stream. BASE CELL SIZE=',num2str(wake_cell_size(1)),' m. TIME STEP=',num2str(time_step(3)),' s.');
hFig.NumberTitle='off';
movegui(hFig,"west");
statconvergenceplot(Cl,50);
ZoomPanCallback(hFig);

% Strouhal estraction 
format long;
% Vector of Strouhal associated to the maxima
STR=FFT_Cl{jj}{3};
STR_w=FFT_V{jj}{3};
% Store Str to check convergence in time
STR_s(3)=STR(2);
STR_s_w(3)=STR_w(5);
CL_s(3)=mean(Cl(:,1));
%% Mesh convergence - Free-stream - Shock wave
close all;

roc=zeros(1,length(wake_cell_size));
for jj=1:length(wake_cell_size)-1
    roc(jj+1)=(abs((STR_s(jj+1)-STR_s(jj)))/STR_s(jj))*100;
end    

figure('Name','Grid convergence.','NumberTitle','off'); hold on;

ax=gca;
yyaxis left
ax.YColor='#4DBEEE';
plot(fliplr(wake_cell_size)/delta,STR_s,'.-','markersize',16,'Color','#4DBEEE'); grid on; grid minor; 
xlab_left=xlabel("$\frac{\delta_{x}}{\delta}$"); ylab_left=ylabel("Shock wave Strouhal"); 
set(gca,'FontSize',8.5);
set(ylab_left,'Interpreter','Latex','Rotation',90,'FontSize', 13); 
set(xlab_left,'Interpreter','Latex','FontSize', 16); 
%title('Time convergence');
yyaxis right
ax.YColor='#A2142F';
plot(fliplr(wake_cell_size)/delta,roc,'.-','markersize',16,'Color','#A2142F'); grid on; grid minor; 
ylab=ylabel("$ROC [\%]$",'HorizontalAlignment','left'); 
set(ylab,'Interpreter','Latex','Rotation',90,'FontSize', 13); 
% Change the digits size 
ax=gca;
ax.XLabel.FontSize=15;
ax.YLabel.FontSize=15;
%title('Strouhal number Vs Time step');
%% Mesh convergence - Free-stream - Von Kármán street
close all;

roc=zeros(1,length(wake_cell_size));
for jj=1:length(wake_cell_size)-1
    roc(jj+1)=(abs((STR_s_w(jj+1)-STR_s_w(jj)))/STR_s_w(jj))*100;
end    

figure('Name','Grid convergence.','NumberTitle','off'); hold on;

ax=gca;
yyaxis left
ax.YColor='#4DBEEE';
plot(fliplr(wake_cell_size)/c,STR_s_w,'.-','markersize',16,'Color','#4DBEEE'); grid on; grid minor; 
xlab_left=xlabel("$\frac{\delta_{x}}{c}$"); ylab_left=ylabel("Von Kármán Strouhal"); 
set(gca,'FontSize',8.5);
set(ylab_left,'Interpreter','Latex','Rotation',90,'FontSize', 13); 
set(xlab_left,'Interpreter','Latex','FontSize', 16); 
%title('Time convergence');
yyaxis right
ax.YColor='#A2142F';
plot(fliplr(wake_cell_size)/c,roc,'.-','markersize',16,'Color','#A2142F'); grid on; grid minor; 
ylab=ylabel("$ROC [\%]$",'HorizontalAlignment','left'); 
set(ylab,'Interpreter','Latex','Rotation',90,'FontSize', 13); 
% Change the digits size 
ax=gca;
ax.XLabel.FontSize=15;
ax.YLabel.FontSize=15;
%title('Strouhal number Vs Time step');
%% Mesh convergence - Free-stream - Cl
close all;

roc=zeros(1,length(wake_cell_size));
for jj=1:length(wake_cell_size)-1
    roc(jj+1)=(abs((CL_s(jj+1)-CL_s(jj)))/CL_s(jj))*100;
end    

figure('Name','Grid convergence.','NumberTitle','off'); hold on;

ax=gca;
yyaxis left
ax.YColor='#4DBEEE';
plot(fliplr(wake_cell_size)/c,CL_s,'.-','markersize',16,'Color','#4DBEEE'); grid on; grid minor; 
xlab_left=xlabel("$\frac{\delta_{x}}{c}$"); ylab_left=ylabel("$C_{l}$"); 
set(gca,'FontSize',10);
set(ylab_left,'Interpreter','Latex','Rotation',0,'FontSize', 13); 
set(xlab_left,'Interpreter','Latex','FontSize', 16); 
%title('Time convergence');
yyaxis right
ax.YColor='#A2142F';
plot(fliplr(wake_cell_size)/c,roc,'.-','markersize',16,'Color','#A2142F'); grid on; grid minor; 
ylab=ylabel("$ROC [\%]$",'HorizontalAlignment','left'); 
set(ylab,'Interpreter','Latex','Rotation',90,'FontSize', 13); 
% Change the digits size 
ax=gca;
ax.XLabel.FontSize=16;
ax.YLabel.FontSize=17;
%title('Strouhal number Vs Time step');






