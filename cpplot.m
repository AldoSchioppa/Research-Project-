function [Cp] = cpplot(filename,p_inf,c,rho_inf,u_inf,dim,surfaces)

%%%%%%%%%% Mean pressure coefficient %%%%%%%%%%
% The function allows the user to clean and plot the pressure coefficient
% data from Star-CCM+ XY plot tree. The code has been tailored for airfoils
% geometry being split into LE-Central Body-TE on both suction and pressure
% surface.
%%%%%%%%%%%%%%%%%
% INPUT:
% filename                 : file path.
% airfoil_coordinates      : Airfoil coordinates.
% Airfoil                  : Airfoil's name.
% c                        : Reference length.
% u_inf                    : Free-stream velocity.
% rho_inf                  : Free-stream density.
% dim                      : True or False. True if the signal is dimensional.
%%%%%%%%%%%%%%%%%%
% OUTPUT:
% Cp                       : A cell array containg the pressure coefficient both on
%                           pressure and suction surfaces, along with the
%                           corrisponding abscissa.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the data
data=readmatrix(filename);

% Due to the data format, in this step different vectors of data are
% extracted. 
p_LE_ss=data(:,6);  p_C_ss=data(:,8);  p_TE_ss=data(:,2); 
p_LE_ps=data(:,10); p_C_ps=data(:,12); p_TE_ps=data(:,4);

% Abscissa extraction. The data need to be sorted out. To obtain corrispondances with the signal, the function "sort" has
% been employed, which outputs the old indexes.
[dx_LE_ss,I_LE_ss]=sort(data(:,5)); [dx_C_ss,I_C_ss]=sort(data(:,7));  [dx_TE_ss,I_TE_ss]=sort(data(:,1));
[dx_LE_ps,I_LE_ps]=sort(data(:,9)); [dx_C_ps,I_C_ps]=sort(data(:,11)); [dx_TE_ps,I_TE_ps]=sort(data(:,3));

% Correspondances
p_LE_ss=p_LE_ss(I_LE_ss); p_C_ss=p_C_ss(I_C_ss); p_TE_ss=p_TE_ss(I_TE_ss);
p_LE_ps=p_LE_ps(I_LE_ps); p_C_ps=p_C_ps(I_C_ps); p_TE_ps=p_TE_ps(I_TE_ps);

% Clean the arrays from NaN values. It is necessary to keep track of the
% Nan being removed because there might be no NaN for the corresponding
% data when plotting.
[dx_LE_ss,I_dx_LE_ss]=rmmissing(dx_LE_ss); [dx_C_ss,I_dx_C_ss]=rmmissing(dx_C_ss); [dx_TE_ss,I_dx_TE_ss]=rmmissing(dx_TE_ss);
[dx_LE_ps,I_dx_LE_ps]=rmmissing(dx_LE_ps); [dx_C_ps,I_dx_C_ps]=rmmissing(dx_C_ps); [dx_TE_ps,I_dx_TE_ps]=rmmissing(dx_TE_ps);

[p_LE_ss,I_LE_ss]=rmmissing(p_LE_ss); [p_C_ss,I_C_ss]=rmmissing(p_C_ss); [p_TE_ss,I_TE_ss]=rmmissing(p_TE_ss);
[p_LE_ps,I_LE_ps]=rmmissing(p_LE_ps); [p_C_ps,I_C_ps]=rmmissing(p_C_ps); [p_TE_ps,I_TE_ps]=rmmissing(p_TE_ps);

% Coherence with eliminating NaN
dx_LE_ss(I_LE_ss(1:length(dx_LE_ss))==1)=[]; dx_C_ss(I_C_ss(1:length(dx_C_ss))==1)=[]; dx_TE_ss(I_TE_ss(1:length(dx_TE_ss))==1)=[]; 
dx_LE_ps(I_LE_ps(1:length(dx_LE_ps))==1)=[]; dx_C_ps(I_C_ps(1:length(dx_C_ps))==1)=[]; dx_TE_ps(I_TE_ps(1:length(dx_TE_ps))==1)=[]; 

p_LE_ss(I_dx_LE_ss(1:length(p_LE_ss))==1)=[]; p_C_ss(I_dx_C_ss(1:length(p_C_ss))==1)=[]; p_TE_ss(I_dx_TE_ss(1:length(p_TE_ss))==1)=[];
p_LE_ps(I_dx_LE_ps(1:length(p_LE_ps))==1)=[]; p_C_ps(I_dx_C_ps(1:length(p_C_ps))==1)=[]; p_TE_ps(I_dx_TE_ps(1:length(p_TE_ps))==1)=[];

% Assembling data
% Abscissa
dx_ss=[dx_LE_ss;dx_C_ss;dx_TE_ss];
dx_ps=[dx_LE_ps;dx_C_ps;dx_TE_ps];
% Presure 
p_ss = [p_LE_ss; p_C_ss; p_TE_ss];
p_ps = [p_LE_ps; p_C_ps; p_TE_ps];

% Dimensional coherence
min_len=min(length(dx_ss),length(p_ss));
dx_ss=dx_ss(1:min_len); p_ss=p_ss(1:min_len);
min_len=min(length(dx_ps),length(p_ps));
dx_ps=dx_ps(1:min_len); p_ps=p_ps(1:min_len);

% Nondimensionalization w.r.t. free-stream condition.
if dim==true
    p_ss=(p_ss-p_inf)./(0.5*rho_inf*u_inf^2);
    p_ps=(p_ps-p_inf)./(0.5*rho_inf*u_inf^2);
    dx_ss=dx_ss./c;
    dx_ps=dx_ps./c;
else
end 

% Output
Cp_plus=[dx_ps,p_ps];
Cp_minus=[dx_ss,p_ss];
Cp={Cp_minus;Cp_plus};

% Set the plot 
%Title=strcat('Time-Average Cp.');
if strcmp(surfaces,'both')
    plot(dx_ss,-p_ss,'*','color','k'); hold on; 
    plot(dx_ps,-p_ps,'o','color','r'); grid on; xlim([0 1]); axis tight; hold off
    xlab=xlabel("$\frac{x}{c} [-]$"); ylab = ylabel("$ -C_{p} [-] \;\;\;\;\;\;\;\;\;\;$"); legend_Cp=legend('Suction surface','Pressure surface','Location','southeast');
    set(xlab,'Interpreter','latex','FontSize',16); set(ylab,'Interpreter','latex','Rotation',90,'FontSize',12); set(legend_Cp,'Interpreter','latex','Fontsize',9);
    set(gca,'FontSize',10);
    ax=gca;
    ax.XLabel.FontSize=15;
    ax.YLabel.FontSize=15;
    title('Pressure Distribution');
elseif strcmp(surfaces,'upper') 
    plot(dx_ss,-p_ss,'o'); hold on;
    xlab=xlabel("$\frac{x}{c}$"); ylab = ylabel("$ -C_{p}^{-} \;\;\;\;$");  grid on; grid minor; axis tight;
    set(xlab,'Interpreter','latex','FontSize',16); set(ylab,'Interpreter','latex','Rotation',0,'FontSize',12); 
    set(gca,'FontSize',10);
    ax=gca;
    ax.XLabel.FontSize=15;
    ax.YLabel.FontSize=15;
else 
    plot(dx_ps,-p_ps,'.-','markersize', 12); hold on;
    xlab=xlabel("$\frac{x}{c}$"); ylab = ylabel("$ -C_{p}^{+} \;\;\;\;$"); 
    set(xlab,'Interpreter','latex','FontSize',16); set(ylab,'Interpreter','latex','Rotation',0,'FontSize',12);  
    set(gca,'FontSize',10);
    ax=gca;
    ax.XLabel.FontSize=15;
    ax.YLabel.FontSize=15;
end