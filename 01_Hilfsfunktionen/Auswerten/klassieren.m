function [H,maxW] = klassieren(Data, ncl, plotopt,varargin)
% Funktion gibt Haeufigkeitsverteilung an
% Input :
% Data     - Datenreihe
% ncl      - Anzahl von Klassen
% plotopt  - 0 - keinen plot erstellen
%            1 - plot erstellen
% varargin - beliebige plot optionen
%
% Output:
% H    - Haeufigkeitssumme
% maxW - Maximalwert
%__________________________________________________________________________

% Speicher
H = zeros(1,ncl);     % Speicher fuer Haeufigkeitssumme
% Maximalwert
maxW = max(Data);
% Klassieren
dummy = int64(Data./maxW .* (ncl-1) + 1);
% Haefigkeitssumme
for i = 1:ncl
    idx = dummy >= i;
    H(i) = sum(idx);
end
% ... Barplot
if plotopt
    figure, grid on, hold on
    barh(1:ncl,H,'EdgeColor','k','LineWidth',1.3,'BarWidth',1,varargin{:})
    set(gca,'XScale','log')
    xlabel('Ueberschreitungshaeufigkeit')
    ylabel('Klasse')
end

end % Ende Funktion