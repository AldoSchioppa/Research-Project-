
close all; clc; clear variables; 

% Reference temperature
T_0=273.15; %K

% Initial global cell size
dx=0.6; %m

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

% Vector of Mach 
M=[0.2:0.1:0.8,0.81:0.01:1];

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
rho_inf=Re*mu_inf./(V_inf*c);

% Freestream pressure
p_inf=rho_inf*R*T_inf;

% Vector of freestream condition
FREE_STREAM=[rho_inf,p_inf,mu_inf,T_inf,V_inf,a_inf];

%tr=[8000,3000,8000,8000,8000,8000,8000];

str={'M02','M03','M04','M05','M06','M07','M08',...
    'M081','M082','M083','M084','M085','M086','M087','M088','M089','M09',...
    'M091','M092','M093','M094','M095','M096','M097','M098','M099','M1',};

%% Re 2000
Cd=zeros(1,length(str));

figure(1)

for jj=1:length(str)
    filename=strcat('D:\ResearchProject_Major\Data\Re2000\Drag_Divergence\',str{jj},'_D');
    data=readmatrix(filename);
    % plot(data(:,1),data(:,2)); hold on;
    Cd(jj)=mean(data(20000:end,2))/(0.5*rho_inf(jj)*V_inf(jj)^2*c);
end

subplot(1,2,2)
semilogy(M,Cd/Cd(1),'d'); hold on;


subplot(1,2,1)
semilogy(M,Cd,'d'); hold on;

%%
Cd=zeros(1,length(str));


for jj=1:length(str)
    filename=strcat('D:\ResearchProject_Major\Data\ReINF\Drag_Divergence\',str{jj},'_D');
    data=readmatrix(filename);
    %plot(data(:,1),data(:,2)); hold on;
    Cd(jj)=mean(data(9000:end,2))/(0.5*rho_inf(jj)*V_inf(jj)^2*c);
end

subplot(1,2,2)
semilogy(M,Cd/Cd(1),'o'); grid on; grid minor; axis([min(M) max(M) 1 10^3]); axis square;
xlab=xlabel('$M$'); ylab=ylabel('$ \frac{Cd}{Cd_{M=0.2}} \;\;\;$');
set(xlab,'Interpreter','latex','FontSize',16); set(ylab,'Interpreter','latex','Rotation',0,'FontSize',18)
leg=legend({'$Re_{c}=2000$','Inviscid'},'Orientation','horizontal');
set(leg,'Interpreter','latex','FontSize',11);
title('(a)');

subplot(1,2,1)
semilogy(M,Cd,'o'); grid on; grid minor; axis([min(M) max(M) 10^-3 1]); axis square;
xlab=xlabel('$M$'); ylab=ylabel('$C_d$');
set(xlab,'Interpreter','latex','FontSize',16); set(ylab,'Interpreter','latex','Rotation',0,'FontSize',16)
leg=legend({'$Re_{c}=2000$','Inviscid'},'Orientation','horizontal');
set(leg,'Interpreter','latex','FontSize',11);
title('(b)');

%% Pressure distribution for several Mach 

symbol={'diamond','square','^','x','*','o'};
file_str={'M0.8','M0.85','M0.9','M0.95','M1'};
figure(2)

lol=subplot(1,2,1);
for jj=0:length(file_str)-1
    cpplot(strcat('D:\ResearchProject_Major\Data\Re2000\Pressure_Distribution\',file_str{jj+1},'_p_mean.csv'),p_inf(7+5*jj),c,rho_inf(7+5*jj), ...
        V_inf(7+5*jj),true,'upper',symbol{jj+1}); hold on; axis square;
end
cpplot(strcat('D:\ResearchProject_Major\Data\Re2000\Pressure_Distribution\M0.93_p_mean.csv'),p_inf(20),c,rho_inf(20), ...
        V_inf(20),true,'upper',symbol{6}); hold on; axis square;

title('(a)');
leg=legend('M=0.8','M=0.85','M=0.9','M=0.95','M=1','M=0.93','Orientation','horizontal');
set(leg,'Interpreter','latex','FontSize',11);

subplot(1,2,2)
for jj=0:length(file_str)-1
    cfplot(strcat('D:\ResearchProject_Major\Data\Re2000\Skin_Friction_Coefficient\',file_str{jj+1},'_f_mean.csv'),c,rho_inf(7+5*jj),V_inf(7+5*jj) ...
        ,true,'upper'); hold on; axis square;
end
cfplot(strcat('D:\ResearchProject_Major\Data\Re2000\Skin_Friction_Coefficient\M0.93_f_mean.csv'),c,rho_inf(20),V_inf(20) ...
        ,true,'upper'); axis square; 

title('(b)');
leg=legend('M=0.8','M=0.85','M=0.9','M=0.95','M=1','M=0.93','Orientation','horizontal');
set(leg,'Interpreter','latex','FontSize',11);

pos_shock=[0.2,0.212,0.25,0.27,0.3,0.325,0.35,0.38,0.4,0.45,0.4,0.425,0.45,0.55,0.7,0.725,0.75,0.75,0.75,0.75,0.75];

%subplot(2,2,[3,4])
figure
plot(M(7:end)',pos_shock,'diamond','color','k'); axis square; grid on; grid minor; 
xlab=xlabel('$M [-]$'); ylab=ylabel('$\frac{x_s}{c} [-] \;\;\;\;\;\;\;\;$');
set(xlab,'Interpreter','latex','FontSize',16); set(ylab,'Interpreter','latex','Rotation',0,'FontSize',18)
title('(c)');


%%

