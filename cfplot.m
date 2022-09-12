function [Cf]=cfplot(filename,c,rho_inf,u_inf,dim)

%%%%%%%%%% Mean Skin friction coefficient %%%%%%%%%%

% Data extraction 

% Read the data
data=readmatrix(filename);

% Due to the data format, in this step different vector of data are
% extracted.
tau_f_LE_ss=data(:,6);  tau_f_C_ss=data(:,8);  tau_f_TE_ss=data(:,2); 
tau_f_LE_ps=data(:,10); tau_f_C_ps=data(:,12); tau_f_TE_ps=data(:,4);

% Abscissa extraction. The data need to be sorted out. To obtain corrispondances with the signal, the function "sort" has
% been employed, which outputs the old indexes.
[dx_LE_ss,I_LE_ss]=sort(data(:, 5)); [dx_C_ss,I_C_ss]=sort(data(:, 7));  [dx_TE_ss,I_TE_ss]=sort(data(:, 1));
[dx_LE_ps,I_LE_ps]=sort(data(:, 9)); [dx_C_ps,I_C_ps]=sort(data(:, 11)); [dx_TE_ps,I_TE_ps]=sort(data(:, 3));

% Correspondances
tau_f_LE_ss=tau_f_LE_ss(I_LE_ss); tau_f_C_ss=tau_f_C_ss(I_C_ss); tau_f_TE_ss=tau_f_TE_ss(I_TE_ss);
tau_f_LE_ps=tau_f_LE_ps(I_LE_ps); tau_f_C_ps=tau_f_C_ps(I_C_ps); tau_f_TE_ps=tau_f_TE_ps(I_TE_ps);

% Clean the arrays from NaN values. It is necessary to keep track of the
% Nan being removed because there might be no NaN for the corresponding
% data when plotting. 
[dx_LE_ss,I_dx_LE_ss]=rmmissing(dx_LE_ss); [dx_C_ss,I_dx_C_ss]=rmmissing(dx_C_ss); [dx_TE_ss,I_dx_TE_ss]=rmmissing(dx_TE_ss);
[dx_LE_ps,I_dx_LE_ps]=rmmissing(dx_LE_ps); [dx_C_ps,I_dx_C_ps]=rmmissing(dx_C_ps); [dx_TE_ps,I_dx_TE_ps]=rmmissing(dx_TE_ps);

[tau_f_LE_ss, I_LE_ss] = rmmissing(tau_f_LE_ss); [tau_f_C_ss, I_C_ss] = rmmissing(tau_f_C_ss); [tau_f_TE_ss, I_TE_ss] = rmmissing(tau_f_TE_ss);
[tau_f_LE_ps, I_LE_ps] = rmmissing(tau_f_LE_ps); [tau_f_C_ps, I_C_ps] = rmmissing(tau_f_C_ps); [tau_f_TE_ps, I_TE_ps] = rmmissing(tau_f_TE_ps);

% Coherence with eliminating NaN
dx_LE_ss(I_LE_ss(1:length(dx_LE_ss))==1)=[]; dx_C_ss(I_C_ss(1:length(dx_C_ss))==1)=[]; dx_TE_ss(I_TE_ss(1:length(dx_TE_ss))==1)=[]; 
dx_LE_ps(I_LE_ps(1:length(dx_LE_ps))==1)=[]; dx_C_ps(I_C_ps(1:length(dx_C_ps))==1)=[]; dx_TE_ps(I_TE_ps(1:length(dx_TE_ps))==1)=[]; 

tau_f_LE_ss(I_dx_LE_ss(1:length(tau_f_LE_ss))==1)=[]; tau_f_C_ss(I_dx_C_ss(1:length(tau_f_C_ss))==1)=[]; tau_f_TE_ss(I_dx_TE_ss(1:length(tau_f_TE_ss))==1)=[];
tau_f_LE_ps(I_dx_LE_ps(1:length(tau_f_LE_ps))==1)=[]; tau_f_C_ps(I_dx_C_ps(1:length(tau_f_C_ps))==1)=[]; tau_f_TE_ps(I_dx_TE_ps(1:length(tau_f_TE_ps))==1)=[];

% Assembling data
% Abscissa
dx_ss=[dx_LE_ss;dx_C_ss;dx_TE_ss];
dx_ps=[dx_LE_ps;dx_C_ps;dx_TE_ps];
% Pressure 
tau_f_ss=[tau_f_LE_ss;tau_f_C_ss;tau_f_TE_ss];
tau_f_ps=[tau_f_LE_ps;tau_f_C_ps;tau_f_TE_ps];

% Dimensional coherence
min_len=min(length(dx_ss),length(tau_f_ss));
dx_ss=dx_ss(1:min_len); tau_f_ss=tau_f_ss(1:min_len);
min_len=min(length(dx_ps),length(tau_f_ps));
dx_ps=dx_ps(1:min_len); tau_f_ps=tau_f_ps(1:min_len);

% Nondimensionalization w.r.t. free-stream condition.
if dim==true
    tau_f_ss=tau_f_ss./(0.5*rho_inf*u_inf^2);
    tau_f_ps=tau_f_ps./(0.5*rho_inf*u_inf^2);
    dx_ss=dx_ss./c;
    dx_ps=dx_ps./c;
else
end 

% Output
Cf_plus=[dx_ps,tau_f_ps];
Cf_minus=[dx_ss,tau_f_ss];
Cf={Cf_minus;Cf_plus};

% Set the plot 
%Title=strcat('Time-average Cf.');
plot(dx_ss,tau_f_ss,'-','LineWidth',3); hold on;  % '#D95319'
plot(dx_ps,tau_f_ps,'-','LineWidth',3); hold on; grid on; grid minor; 
xlim([0 1]); %axis tight; 
text(0.487,0.1,sprintf('\\leftarrow @%.2f',0.487)); 
xline(0.487,'--','Color','k');
legend off;
xlab=xlabel("$\frac{x}{c} [-]$"); ylab=ylabel("$ C_{f} [-] \;\;\;\;\;\;\;\;$"); legend_Cp=legend('Suction surface','Pressure surface','Location', 'northeast');
set(xlab,'Interpreter','latex','FontSize',16); set(ylab,'Interpreter','latex','Rotation',90,'FontSize',14); set(legend_Cp,'Interpreter','latex','Fontsize',9);
set(gca,'FontSize',10);
ax=gca;
ax.XLabel.FontSize=15;
ax.YLabel.FontSize=15;
title('Skin Friction Coefficient Distribution');
end