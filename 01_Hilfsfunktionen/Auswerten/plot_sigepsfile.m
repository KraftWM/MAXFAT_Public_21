function [SIG,EPS,DLZ] = plot_sigepsfile(sigepsfile,savemode,varargin)

if nargin == 1
    savemode = 0;
end

if savemode
    jobname = varargin{1,1};
    outpath = varargin{1,2};
end
    
LW = 0.9;

% Daten Einlesen
[SIG,EPS,DLZ] = read_sigepsfile(sigepsfile);

% Plote Zeitverläufe
figure
subplot(2,1,1),grid on, hold on
latexinterpreter
plot(DLZ,SIG','LineWidth',LW);
ylabel('$\sigma$');
legend('xx','yy','xy')
subplot(2,1,2), grid on, hold on
latexinterpreter
plot(DLZ,EPS','LineWidth',LW);
legend('xx','yy','zz','xy')
xlabel('DLZ'); ylabel('$\varepsilon$');
latexinterpreter;
if savemode
    figname = [outpath,'/',jobname,'_ZeitVer.fig'];
    savefig(gcf,figname);
end

% Hyseresen
figure, grid on, hold on
latexinterpreter;
plot(EPS([1 2 4],:)',SIG','LineWidth',LW);
xlabel('$\varepsilon$'),ylabel('$\sigma$');
legend('xx','yy','xy')
latexinterpreter;
if savemode
    figname = [outpath,'/',jobname,'_Hyst.fig'];
    savefig(gcf,figname);
end

% Spannungsphase
figure, grid on, hold on
latexinterpreter;
plot(SIG(1,:),SIG(3,:),'LineWidth',LW);
xlabel('$\sigma_{xx}$'),ylabel('$\tau_{xy}$');
legend('Spannungsphase')
latexinterpreter;
if savemode
    figname = [outpath,'/',jobname,'_SpanPha.fig'];
    savefig(gcf,figname);
end


% Dehnungsphase
figure, grid on, hold on
latexinterpreter;
plot(EPS(1,:),EPS(4,:),'LineWidth',LW);
xlabel('$\varepsilon_{xx}$'),ylabel('$\gamma_{xy}$');
legend('Dehnungsphase')
latexinterpreter;
if savemode
    figname = [outpath,'/',jobname,'_DehnPha.fig'];
    savefig(gcf,figname);
end


% Amplituden
ASXX = hcm_origV2(SIG(1,:));
ASYY = hcm_origV2(SIG(2,:));
ASXY = hcm_origV2(SIG(3,:));
AEXX = hcm_origV2(EPS(1,:));
AEYY = hcm_origV2(EPS(2,:));
AGXY = hcm_origV2(EPS(4,:));
figure
subplot(2,1,1),grid on, hold on
latexinterpreter
plot(ASXX(1,:),ASXX(5,:),'d-','MarkerSize',7);
plot(ASYY(1,:),ASYY(5,:),'s-','MarkerSize',7);
plot(ASXY(1,:),ASXY(5,:),'<-','MarkerSize',7);
ylabel('$\sigma_a$');
legend('xx','yy','xy')
title('Amplituden')
subplot(2,1,2), grid on, hold on
latexinterpreter
plot(AEXX(1,:),AEXX(5,:),'d-','MarkerSize',7);
plot(AEYY(1,:),AEYY(5,:),'s-','MarkerSize',7);
plot(AGXY(1,:),AGXY(5,:),'<-','MarkerSize',7);
xlabel('N'); ylabel('$\varepsilon_a$');
legend('xx','yy','xy')
latexinterpreter;
if savemode
    figname = [outpath,'/',jobname,'_AmpVer.fig'];
    savefig(gcf,figname);
end

% Alle Figures Schließen
if savemode
    close all
end

end