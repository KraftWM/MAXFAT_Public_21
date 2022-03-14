function ax = plot_BWL(ax,N,Amp,messflag,titleName,legendentry,YAxisLabel,varargin)
% Plotet Bauteilwoehlerlinie
% INPUT:
% ax          - axis handle
% Nexp        - Experimentelle Ergebnisse
% Ncalc       - Rechenergebnisse 
% messflag    - 0 keine Messwerte 1 Messwerte
% titlename   - Name der figure
% legendentry - Eintrag fuer Legende
% YAxisLabel  - Label fuer y achse
% varargin    - plot optionen
%
% OUTPUT:
% ax          - axis handle
%--------------------------------------------------------------------------

% oeffne Figure
if ax == 0
    figure; grid on, hold on
    ax = gca;
    latexinterpreter;
    xmin = 1; xmax = 1e7;
    ymin = 1; ymax = 1e7;
    xmin = max( [xmin 10^floor(log10(min(N)))]);
    xmax = min( [xmax 10^ceil(log10(max(N)))]);
    ymin = max( [ymin 10^floor(log10(min(Amp)))]);
    ymax = min( [ymax 10^ceil(log10(max(Amp)))]);
else
    xmin = ax.XLim(1); xmax = ax.XLim(2);
    ymin = ax.YLim(1); ymax = ax.YLim(2);
    xmin = min( [xmin 10^floor(log10(min(N)))]);
    xmax = min( [1e7 max( [xmax 10^ceil(log10(max(N)))])]);
    ymin = min( [ymin 10^floor(log10(min(Amp)))]);
    ymax = max( [ymax 10^ceil(log10(max(Amp)))]);
end

% Sortiere Ergebnisse
[Amp,I] = sort(Amp);
N = N(I);


% Durchlaeufer
if messflag
    % Einfuegen dieses Codes am besten zwischen "figure(..." und dem Plotten der
    % Veruchspunkte, damit der Anfangspunkt des Pfeils hinter dem
    % Versuchspunkt verschwindet.

    % Die Variable N entspricht hier einem Vektor mit allen Anrisslebensdauern
    % aus den Versuchen. Die Variable eps_a entspricht einem Vektor mit allen
    % zugehoerigen, vorgegebenen Dehnungsamplituden.

    Dauerfest = 1999999;                                                       % Durchlaeufer bei 2.000.000 Schwingspiele (automatischer Stop der Pruefmaschine)

    N_D = N(N>Dauerfest);
    Amp_D = Amp(N>Dauerfest);

    pfeil = char(8599);                                                                  

    Anzahl_Durchlaeufer = length(N_D);

    for i = 1:Anzahl_Durchlaeufer
    text(N_D(i),Amp_D(i),pfeil,'HorizontalAlignment','left','VerticalAlignment','baseline','FontSize',20,'Interpreter','tex');
    end
else
    Dauerfest = 5999999;  
    N_D = N(N>Dauerfest);
    Amp_D = Amp(N>Dauerfest);
    Anzahl_Durchlaeufer = length(N_D);
    N(N>Dauerfest) = 2e40;
end

% Plote N-N
plot(ax,N,Amp,varargin{:});



% logScale
set(gca,'XScale','log','YScale','log')

% namen
xlabel('$N_{f}$'), ylabel(YAxisLabel)
title(titleName);

% Achsen beschneiden
axis([xmin xmax ymin ymax]);

% Legende Hinzufuegen
if strcmp(ax.Legend.String,'data1')
    ax.Legend.String = {legendentry};
else
    ax.Legend.String{length(ax.Legend.String)} = legendentry;
end

end