function [] = statconvergenceplot(filename, n, transient,V_inf,rho_inf,c)
% Data extraction
data=readmatrix(filename);
signal=data(:,2); dt=data(:,1);

% Erase transient data
if transient~=0
    dt=dt(transient:end)*V_inf/c;
    signal=signal(transient:end)/(0.5*rho_inf*V_inf^2*c);
end

L=length(signal);
% Initialization
w_mean=zeros(1,L);
N = zeros(1,L);
jj=1;
while n<=L
    x=signal(1:n);
    w=dt(1:n);
    w_mean(jj)=wmean(x,w); % weighted average function
    N(jj)=n;
    n=1.1*n;
    jj=jj+1;
end
% The last iteration before breaking the while loop implies an extra index
% which has to be erased.
w_mean=w_mean(1:jj-1);
N=N(1:jj-1);
% Plot 
plot(N,w_mean,'.-','MarkerSize',16,'LineWidth',2,'Color','r'); grid on; axis tight; axis([min(N) max(N) min(abs(w_mean)) max(abs(w_mean))]); hold on;
plot(N,mean(signal)*ones(1,length(N)),'LineWidth',2,'color','k');
xlab=xlabel("$N [-]$"); ylab=ylabel("$\bar{C}_{l} [-] \;\;\;\;\;$"); legend1=legend('Weighted mean', 'Mean Signal','Location','southeast');
set(legend1,'Interpreter','latex','FontSize',9); 
set(xlab,'Interpreter','latex','FontSize',13); set(ylab,'Interpreter','latex','Rotation',90,'FontSize',16)
set(gca,'FontSize',10);
ax=gca;
ax.XLabel.FontSize=15;
ax.YLabel.FontSize=15;
title('Statistical Convergence');
%title('Statistical convergence.');
end
