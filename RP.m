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
FREE_STREAM=[rho_inf,p_inf,mu_inf,T_inf,V_inf,a_inf];

% Convective time
conv_dt=c/V_inf;

% Dominant frequency
F=St*V_inf/c;

% Vector of dominant time step
Dt=1./F;
Dt(end+1)=conv_dt;

% Division by 50 for first time step
Dt=Dt/50;

% Most restrictive time step, i.e. the smallest 
dt=min(Dt);

% Max CFL (for the finest mesh)
CFL=dt*(V_inf+a_inf)/min_dx;

% New dt after CFL correction
dt=Cfl*min_dx/(V_inf+a_inf);

% Minimum transient time
tr=L/V_inf;
%%
filename='D:\ResearchProject_Major\Data\Re2000\Snapshots_SHOCK\SHOCK_64184.csv';
data=readtable(filename);

display(data.Properties.VariableNames);

idx_t=1;
idx_F=3;
idx_x=4;
idx_y=5;
idx_z=6;

[x, IDX]=sort(data.(idx_x));
%x=data.(idx_x);
y=data.(idx_y); 
% y=y(IDX);

[X, Y]=meshgrid(unique(x),unique(y));
Y=flip(Y);

%%

% Specify the folder where the files live.
myFolder='D:\ResearchProject_Major\Data\Re2000\Snapshots_SHOCK';
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
    errorMessage=sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder=uigetdir(); % Ask for a new one.
    if myFolder==0
         % User clicked Cancel
         return;
    end
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.csv'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
Field={};
for jj=1:length(theFiles) % Change 
    baseFileName=theFiles(jj).name;
    fullFileName=fullfile(theFiles(jj).folder, baseFileName);
    fprintf(1,'Now reading %s\n',fullFileName);
    data = readmatrix(fullFileName);
    data = sortrows(data,[idx_x,idx_y,idx_z]);
    %data = sortrows(data);
    Field{jj} = data(:,idx_F);
    field_Snapshots(jj,:,:)=reshape(Field{jj},size(X,1),[]);
end

%% Apply POD

[U_POD, S_POD, V_POD] = pod(field_Snapshots);

% nFFT  	= 128;
% nOvlp 	= floor(nFFT/2);
% window	= hann(nFFT);

% The temporal mean of each block ends up in the zero-frequency bin when you take the Fourier transform, 
% so the mean has already been removed from the non-zero frequency components 
% and there is no need to subtract it. 
% But you can set OPTS.mean='blockwise' to remove the mean of each block if your data has a long-time trend 
% that leaks into your low frequencies (for example if your data isn't truly stationary).
OPT.mean='blockwise';

[L,P,f] = spod(field_Snapshots,[],[],[],Dt(1),OPT);

%% Plots the mode energies VS the modes' number
ModeEnergies=S_POD.^2;
ModeEnergyFraction=ModeEnergies/sum(ModeEnergies);
figure('Color','w');
bar(1:length(ModeEnergies),ModeEnergyFraction,'k'); xlim([0 100]);
title('Mode Energies');

%% Plot modes' energy VS frequency 
figure
loglog(f*c/V_inf,L); grid on; grid minor; 
xlabel('St'), ylabel('|A(St)|')

%% Animate the first n snapshots of the pressure field.
n=100;
figure('name','')
for ti=1:n
    pcolor(X,Y,squeeze(field_Snapshots(ti,:,:))); colorbar
    axis equal tight, shading interp, 
    xlabel('x'), ylabel('y')
    pause(0.05)
    drawnow
end
%% Visualize the 1st and 2nd SPOD modes at three frequencies.
% Note that the first mode has a high energy content when St~1.
figure
count = 1;
n_f = 2;
n_m = 2;
for mi = [1 2] % modes counter 
    for fi = [2 21] % frequencies counter 
        subplot(n_m,n_f,count)
        contourf(X,Y,real(squeeze(P(fi,:,:,mi))),11,'edgecolor','none'), axis equal tight, caxis(max(abs(caxis))*[-1 1]); shading interp; colorbar
        xlabel('x'), ylabel('y'), title(['St=' num2str(f(fi)*c/V_inf,'%.2f') ', mode ' num2str(mi) ', \lambda=' num2str(L(fi,mi),'%.2g')])
        xlim([min(x) max(x)]); ylim([min(y) max(y)])
        count=count+1;
    end
   % P(fi,:,:,mi), the first index are the frequencies, while the last one
   % is mode number ranked in descending order by modal energy. Lambda is
   % the energy content associated to the mi-th modes (at fi-th frequency).
end

%% Animate the same modes.
%   Note how all wavepackets travel at approximately the same phase
%   speed c_ph. The reason is that their streamwise wavenumber k_x changes 
%   with frequency such that c_ph = omega/k_x is approximately constant.
figure
nt=30;
T=1/f(21);              % period of the 21th frequency
time=linspace(0,T,nt);     % animate over one period
count=1;
n_f = 2;
n_m = 2;
for ti=1:nt
    for mi = [1 2 3 4]
        for fi = [2 21]
            subplot(n_m,n_f,count)
            pcolor(X,Y,real(squeeze(P(fi,:,:,mi)*exp(2i*pi*f(fi)*time(ti))))), shading interp, axis equal tight, caxis(max(abs(caxis))*[-1 1]); colorbar
            xlabel('x'), ylabel('y'), title(['St=' num2str(f(fi)*c/V_inf,'%.2f') ', mode ' num2str(mi) ', \lambda=' num2str(L(fi,mi),'%.2g')]);
            xlim([min(x) max(x)]); ylim([min(y) max(y)])
            count=count+1;
            hold on
        end
    end
    drawnow
    hold off
    count=1;
end

%% Frequency-time analysis.
%   Computing the mode expansion coefficients of the first nModes modes
%   using a windowing function and weights consistent with the SPOD.

nFFT=1024; % To change coherently with the number of snapshots
window=hann(nFFT);

nModes=2;
a=tcoeffs(field_Snapshots,P,hann(nFFT),[],nModes);

%% Visualize the results. Frequency-time plot (dynamic spectrum).
nt  = size(field_Snapshots,1);
t   = (1:nt)*dt;

figure

% First mode
% subplot(2,1,1)
pcolor(t,f,abs(squeeze(a(:,:,1)))); shading interp; colorbar
%daspect([100 1 1])
title('Frequency-time diagram'); % (first mode, ' num2str(sum(L(:,1))/sum(L(:))*100,'%3.1f') '% of energy)'
xlabel('t [s]'), ylabel('f [Hz]'), caxis([-0.75 0.75].*caxis)

% nModes modes
% subplot(2,1,2)
% pcolor(t,f,squeeze(abs(sum(a,nModes)))); shading interp; colorbar
% %daspect([100 1 1])
% title(['frequency-time diagram (sum of first ' num2str(nModes) ' modes, ' num2str(sum(sum(L(:,1:nModes)))/sum(L(:))*100,'%3.1f') '% of energy)'])
% xlabel('time'), ylabel('frequency'), caxis([0 0.75].*caxis)

