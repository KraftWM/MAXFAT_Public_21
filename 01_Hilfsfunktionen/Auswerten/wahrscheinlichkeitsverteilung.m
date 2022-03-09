function [T9010,xm,s] = wahrscheinlichkeitsverteilung(Nexp,Ncalc,titlename,varargin)
% Berechnet auftretenswahrscheinlichkeit und plotet Sie auch
%
% ANNAHME: log10(Ncalc/Nexp) ist Normalverteilt
%
% INPUT:
% Nexp      - N Experiment
% Ncalc     - N Rechnung
% varargin  - plotoptions
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% filtern Dauerfeste Versuche
idx = (Ncalc < 1e9) & (Nexp < 1e9);

% -------------------------------------------------------------------------
% Werte
if length(Ncalc) == length(Nexp)
    x = Ncalc(idx)./Nexp(idx);                       % Verhältniss Lebensdauern
else
    T9010 = NaN;
    xm = NaN;
    s = NaN;
    ax = NaN;
    return;
end

% -------------------------------------------------------------------------
% Log Normalverteilung
% log werte
logx = log(x)';
% normalverteilung
pd = fitdist(logx,'Normal');
% Standartabweichung 
s = pd.sigma;
% Streuspanne (entlogaritmiert)
% T9010 = exp(pd.icdf(0.9))/exp(pd.icdf(0.1));
% T9010 = exp(2*(pd.icdf(0.9)-pd.mu));
% T9010 = exp(2*(pd.mu-pd.icdf(0.1)));
T9010 = exp(pd.icdf(0.9)-pd.icdf(0.1));
% Mittelwert (entlogaritmiert)
xm = exp(pd.mu);

% -------------------------------------------------------------------------
% Log Normalverteilung Plot
% Erstelle Figure
fig = figure; hold on
colors = lines(2);
% Einzeichnen Mittelwerte
% plot([pd.mu pd.mu],[-10000 10000],'Color',colors(2,:),'LineStyle',':','LineWidth',1.3,'handlevisibility','off')
% Einzeichnen 10% und 90% Überschreitungswahrscheinlichkeit
% plot([pd.icdf(0.9) pd.icdf(0.9)],[-10000 10000],'Color',colors(2,:),'LineWidth',1,'handlevisibility','off')
% plot([pd.icdf(0.1) pd.icdf(0.1)],[-10000 10000],'Color',colors(2,:),'LineWidth',1,'handlevisibility','off')
% Erstelle Normplot
h = normplot(logx);
% Ändere Beschriftungen
set(gca,'FontSize',16)  
xlabel('$N_{calc}/N_{exp}$','Interpreter','Latex')
ylabel('Wahrscheinlichkeit','Interpreter','Latex')
title(titlename,'Interpreter','Latex');
set(gcf,'Units','centimeter','Position',[2 2 15 20])
% Ändere Plot Ausehen
h(1).Color = [0 0 0];           % Alles Schwarz
h(1).Marker = 'o';              % Kreise als Marker für Daten
h(1).MarkerSize = 8;
h(1).MarkerFaceColor = [0 0 0];
h(2).Color = [0 0 0];     % Alles Schwarz
h(2).LineWidth = 1.3;
h(3).Color = [0 0 0];     % Alles Schwarz
h(3).LineWidth = 1.3;
% Legende
% legend('Test','Interpreter','Latex','Location','northwest')
% Annotation Mittelwert
dim = [7.5 4.0 .5 .5];
% dim = [0.5 0.3 .5 .5];
str = ['$x_m = ',num2str(xm),'$'];
annotation('textbox','Units','centimeters','Position',dim,'String',str,'FitBoxToText','on','Interpreter','Latex','FontSize',14,'EdgeColor','none');
% annotation('textbox','Position',dim,'String',str,'FitBoxToText','on','Interpreter','Latex','FontSize',14,'EdgeColor','none');
% Annotation Streuspanne
dim = [6.5 3.0 .5 .5];
% dim = [0.5 0.4 .5 .5];
str = ['$T_{90}/T_{10} = ',num2str(T9010),'$'];
annotation('textbox','Units','centimeters','Position',dim,'String',str,'FitBoxToText','on','Interpreter','Latex','FontSize',14,'EdgeColor','none');
% annotation('textbox','Position',dim,'String',str,'FitBoxToText','on','Interpreter','Latex','FontSize',14,'EdgeColor','none');
% da = annotation('doublearrow');
% da.Parent = fig.CurrentAxes;
% da.X = [pd.icdf(0.90) pd.icdf(0.10)];
% da.Units = 'centimeters';
% da.Y = [-42.0 -42.0];
% Umbenennen x Achse
base = 3;
n = max([ ceil(log(max(x))/log(base)), ceil(log(1/min(x))/log(base)) ]); % Maximale zweier Potenz
set(gca,'XLim',[log(base^-n),log(base^n)]);
N = -n:1:n;
XTicks = base.^N;
XLabels = cell(1,length(N));
for i = 1 : length(N)
    if N(i) < 0
        XLabels{i} = ['$\frac{1}{',num2str(base^(-N(i))),'}$'];
    else
        XLabels{i} = ['$',num2str(base^N(i)),'$'];
    end
end
fig.CurrentAxes.XTick = log(XTicks);
fig.CurrentAxes.TickLabelInterpreter = 'Latex';
fig.CurrentAxes.XTickLabel = XLabels;
end