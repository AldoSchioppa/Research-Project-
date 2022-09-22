close all; clc; clear variables; 

% Reference temperature.
T_0=273.15; %K
% Initial global cell size.
dx=0.5; %m
% Min cell size across the domain.
min_dx=0.00125; %\m
% Target CFL for the finest grid.
Cfl=35;

%% Input 

% Reference length
c=1; %m
% Maximum airfoil thickness
delta=0.12*c;
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
St=[1,0.07];

%% Outputs 

% Domain length 
L=20*2*c; %m
% Freestream viscosity.
mu_inf=1.716e-05*(T_inf/T_0)*((T_0+110.4)/(T_inf+110.4)); %(m^2)/s
% Speed of sound.
a_inf=(gamma*R*T_inf)^0.5; %m/s
% Freestream velocity.
V_inf=M*a_inf; % m/s
% Freestream density. 
rho_inf=Re*mu_inf/(V_inf*c); %Kg/m^3
% Freestream pressure.
p_inf=rho_inf*R*T_inf; %Pa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Convective time.
conv_dt=c/V_inf; %s
% Dominant frequency.
F=St*V_inf/c; %Hz
% Vector of dominant time step.
Dt=1./F; %s
Dt(end+1)=conv_dt;
% Division by 50 for first time step.
Dt=Dt/50; %s
% Most restrictive time step.
dt=min(Dt); %s
% Max CFL (for the finest mesh).
CFL=dt*(V_inf+a_inf)/min_dx;
% New dt after CFL correction.
dt=Cfl*min_dx/(V_inf+a_inf);
% Time step vector for time convergence. Progression of 1/5.
Dt=[dt,dt/5,dt/5/5]; %s
%% Data unpacking

filename='D:\ResearchProject_Major\Data\Re2000\Snapshots_Shock_Refinement\SHOCK_2.000010e+00.csv';
d=readtable(filename,'PreserveVariableNames',true);
display(d.Properties.VariableNames);

% idx_t=1; if existent
idx_F=3; % Defining the index for the field the user shall apply the POD to.

% Grid coordinates indexes.
idx_x=4;
idx_y=5;
idx_z=6; % Not noteworthy for 2D fields.

% Coordinates extraction.
x=d.(idx_x)/c; % non-dimensionalisation with respect to the chord.
y=d.(idx_y)/c; 

% Grid generation.
[X,Y]=meshgrid(unique(x),unique(y));

% chord-wise and chord-transversal coordinates vectors.
x=linspace(min(x),max(x),size(X,1));
y=linspace(min(y),max(y),size(Y,2));

%% Data extraction

% Specify the folder where the files resire.
myFolder='D:\ResearchProject_Major\Data\Re2000\Snapshots_Shock_inv'; % To modify 
% Check to make sure that folder actually exists.  Warn user if it does not.
if ~isfolder(myFolder)
    errorMessage=sprintf(['Error: The following folder does not exist:\n%s\nPlease' ...
        'specify a new folder.'], myFolder);
    uiwait(warndlg(errorMessage));
    myFolder=uigetdir(); % Ask for a new one.
    if myFolder==0
         return;
    end
end

% Get a list of all files in the folder with the desired file name pattern.
filePattern=fullfile(myFolder, '*.csv'); % Change to whatever pattern the user needs.
theFiles=dir(filePattern);
disp('Reading the snaphots files might take a while...')
for jj=1:length(theFiles) 
    % Access the file name.
    baseFileName=theFiles(jj).name;
    % Build the path.
    fullFileName=fullfile(theFiles(jj).folder, baseFileName);
    if rem(jj,2000)==0 % displying the status each 2000 snapshots.
    fprintf(1,'Now reading %s\n',fullFileName);
    disp('------------------------------------')
    end
    % Read data.
    data=readmatrix(fullFileName);
    % Crucial step is sorting data with respect to the grid coordinates.
    data=sortrows(data,[idx_x,idx_y,idx_z]);
    % Store sorted data vector in a cell array.
    Field{jj}=data(:,idx_F);
    % Finally, reshape according to the grid discretisation.
    field_Snapshots(jj,:,:)=reshape(Field{jj},size(X,1),[]);
end

% Time vector
t=(1:length(theFiles))*Dt(1);

%% Apply POD

%load Shock_snapshots.mat

disp('It migth take a while...');

% Proper Orthogonal Decomposition
[U_POD, S_POD, V_POD] = pod(field_Snapshots);

% The temporal mean of each block ends up in the zero-frequency bin when you take the Fourier transform, 
% so the mean has already been removed from the non-zero frequency components 
% and there is no need to subtract it. 
% But you can set OPTS.mean='blockwise' to remove the mean of each block if your data has a long-time trend 
% that leaks into your low frequencies (for example if your data isn't truly stationary).
OPT.mean='blockwise';
% Spectral Proper Orthogonal Decomposition
[L,P,f] = spod(field_Snapshots,[],[],[],Dt(1),OPT);

%% Modes energy VS rank

% The pod.m function outputs a diagonal matrix S_POD whose coefficients are
% the square root of the modes energy in descending order.
ModeEnergies=S_POD.^2;
% Normalisation w.r.t. the total energy.
ModeEnergyFraction=ModeEnergies/sum(ModeEnergies);
% Plot
figure('Color','w');
bar(1:length(ModeEnergies),ModeEnergyFraction,'k'); 
xlim([0 100]);
title('Modes Energy');

%% Frequency spectrum

figure
loglog(f*c/V_inf,L); grid on; grid minor; 
xlab=xlabel("$St [-]$"); ylab=ylabel('\textbf{FFT of} pressure'); 
set(xlab,'Interpreter', 'latex','FontSize',16); 
set(ylab,'Interpreter','latex','Rotation',90,'FontSize',14); 
set(gca,'FontSize',11);

%% Animate the first n snapshots of the pressure field.

n=100; % Number of snapshost to visualise
figure
for ti=1:n
    pcolor(X,flip(Y),squeeze(field_Snapshots(ti,:,:))); colorbar;
    axis tight equal; shading interp;
    set(gca,'FontSize',11);
    xlab=xlabel("$\frac{x}{c}$"); ylab=ylabel("$\frac{y}{c}$"); 
    set(xlab,'Interpreter', 'latex','FontSize',16); 
    set(ylab,'Interpreter','latex','Rotation',0,'FontSize',16); 
    pause(0.05)
    drawnow
end

%% Visualize the SPOD modes at multiple frequencies.

% P(fi,:,:,mi), the first index loops over the frequencies, while the last one
% is mode number ranked in descending order by modal energy. Lambda is
% the energy content associated to the mi-th modes at fi-th frequency.

figure
count=1;
n_f=2; % Number of frequencies.
n_m=2; % Number of modes.
for mi = [1 2] % modes counter.
    for fi = [2 21] % frequencies counter.
        subplot(n_m,n_f,count)
        contourf(X,Y,real(squeeze(P(fi,:,:,mi))),11,'edgecolor','none'); axis equal tight; shading interp;
        caxis(max(abs(caxis))*[-1 1]); 
        colorbar;
        colormap(gray(8192));
        %brighten(.5);
        set(gca,'FontSize',11);
        xlab=xlabel("$\frac{x}{c}$"); ylab=ylabel("$\frac{y}{c} \;\;\;\;$"); 
        set(xlab,'Interpreter', 'latex','FontSize',16); 
        set(ylab,'Interpreter','latex','Rotation',0,'FontSize',16); 
        title(['St=' num2str(f(fi)*c/V_inf,'%.2f') ', mode ' num2str(mi) ', ' ...
            '\lambda=' num2str(L(fi,mi),'%.2g')])
        xlim([min(x) max(x)]); ylim([min(y) max(y)])
        count=count+1;
    end
end

%Note (WAKE) Modes 1 and 2 are a wave pair. This is because POD is a real-valued function.
%So this means that, in order to produce a propagating wave, a pair of two phase-shifted modes are required,
%similar to a sine wave - cosine wave pair.

%% Animate the same modes.

%   Note how all wavepackets travel at approximately the same phase
%   speed c_ph. The reason is that their streamwise wavenumber k_x changes 
%   with frequency such that c_ph = omega/k_x is approximately constant.

figure
nt=100;
T=1/f(2); % period 
time=linspace(0,3*T,nt); % animate over 3 period
count=1;
n_f=2;
n_m=2;
for ti=1:nt
    for mi=[1 2]
        for fi=[2 21]
            subplot(n_m,n_f,count)
            pcolor(X,Y,real(squeeze(P(fi,:,:,mi)*exp(2i*pi*f(fi)*time(ti))))); 
            shading interp; axis equal tight; 
            caxis(max(abs(caxis))*[-1 1]); colorbar; colormap(gray(8192));
            set(gca,'FontSize',11);
            xlab=xlabel("$\frac{x}{c}$"); ylab=ylabel("$\frac{y}{c} \;\;\;\;$"); 
            set(xlab,'Interpreter', 'latex','FontSize',20); 
            set(ylab,'Interpreter','latex','Rotation',0,'FontSize',20); 
            title(['St=' num2str(f(fi)*c/V_inf,'%.2f') ', mode ' num2str(mi) ', ' ...
                '\lambda=' num2str(L(fi,mi),'%.2g')]);
            xlim([min(x) max(x)]); ylim([min(y) max(y)])
            count=count+1;
            hold on
        end
    end
    drawnow
    hold off
    count=1;
end

%% Plots the POD m-th mode for visualization comparison purposes
m=1; % mode index
figure
modeShape=real(squeeze(U_POD(m,:,:))); 
imagesc(x,flip(y),modeShape); axis equal tight; shading interp;
caxis([-max(abs(modeShape(:))) max(abs(modeShape(:)))]); 
colormap; colorbar
set(gca,'FontSize',11);
xlab=xlabel("$\frac{x}{c}$"); ylab=ylabel("$\frac{y}{c}$"); 
set(xlab,'Interpreter', 'latex','FontSize',20); 
set(ylab,'Interpreter','latex','Rotation',0,'FontSize',20); 
title(['Mode ' num2str(m,'%0.0f')]);

%% Plots the time coefficient matrix V for modes 1 and 2, for visualization purposes:

TimeCoefficients1=V_POD(:,1); %1st mode
TimeCoefficients_m=V_POD(:,m); %m-th mode

% Plot
figure
subplot(1,2,1);
plot(t,TimeCoefficients1,'k-'); hold on; axis equal; axis square;
plot(t,TimeCoefficients_m,'r-');
xlim([t(1) t(end)]);
title('Time Coefficients POD');
legend('Mode 1',['Mode ' num2str(m,'%0.0f')]);

% Zoom
subplot(1,2,2);
plot(t(1:200),TimeCoefficients1(1:200),'k-'); hold on; axis equal; axis square;
plot(t(1:200),TimeCoefficients_m(1:200),'r-');
title('Time Coefficients POD - Zoom');
legend('Mode 1',['Mode ' num2str(m,'%0.0f')]);

%% Frequency-time analysis.

% Computing the mode expansion coefficients of the first nModes modes
% using a windowing function consistent with the SPOD.

nFFT=1024; % To change coherently with the number of snapshots
window=hann(nFFT);
% Number of modes.
nModes=2;

% Expansion coefficient function tcoeffs.
a=tcoeffs(field_Snapshots,P,hann(nFFT),[],nModes);

%% Visualize the results. Frequency-time plot (dynamic spectrum).

% Plot
figure('Renderer','painters','Position',[10 10 900 600])
% First mode
pcolor(t*V_inf/c,f*c/V_inf,abs(squeeze(a(:,:,1)))); shading interp; colorbar;
title(strcat(['First mode frequency-time diagram,' ...
    ' '], num2str(sum(L(:,1))/sum(L(:))*100),'%3.1f'), '% of energy)');
xlab=xlabel("$t \frac{V_{\infty}{c}$"); 
ylab=ylabel("$St$"); caxis([-0.75 0.75].*caxis);
set(xlab,'Interpreter', 'latex','FontSize',20); 
set(ylab,'Interpreter','latex','Rotation',0,'FontSize',20); 
set(gca,'FontSize',11); 

% nModes modes (sum of the first nModes modes)
% subplot(2,1,2)
% pcolor(t,f,squeeze(abs(sum(a,nModes)))); shading interp; colorbar
% %daspect([100 1 1])
% title(['frequency-time diagram (sum of first ' num2str(nModes) ' modes, ' ...
%     '' num2str(sum(sum(L(:,1:nModes)))/sum(L(:))*100,'%3.1f') '% of energy)'])
% xlabel('time'), ylabel('frequency'), caxis([0 0.75].*caxis)

