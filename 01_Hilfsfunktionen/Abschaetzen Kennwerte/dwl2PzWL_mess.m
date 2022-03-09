function [PZ,Nlin,mJ,Q] = dwl2PzWL_mess(EA,SA,N,Ks,ns,E,Rm,ND)
% Funktion errechnet Pz Wöhlerlinie aus gemessenen Daten der
% Dehnungswöhlerlinie
%
% INPUT:
% EA - Dehnungsamplituden
% SA - Spannungsamplituden
% N  - Schwingspiele
%
% E     - EModul
% Ks,ns - Parameter RamOsg
% Rm    - Zugfestigkeit
%
% PZD0  - Dauerfestigkeit PZ
% 
%
% OUTPUT:
% PZ,Nlin   - Werte, N_out = N_in falls keine NaN in N
%  mJ       - Steigung
%  Q        - Lebensdauer Pz = 1
% PZD0,ND0  - Dauerfestigkeit
%
% _________________________________________________________________________

% Rausschmeisen aller Punkte mit NaN
notnans = ~isnan(N);
EA = EA(notnans);
SA = SA(notnans);
N = N(notnans);
% Rausschmeißen aller Durchläufer
notdl = N < ND;
EA = EA(notdl);
SA = SA(notdl);
N = N(notdl);

% Anzahl Datenpunkte
numdata = length(EA);
% Fließspannung
Rp02 = Ks * 0.002^ns;
sigF = (Rm+Rp02)/2;

% -------------------------------------------------------------------------
% Schleife über Datenpunkte
Weff = zeros(numdata,1);
for i = 1 : numdata
    % ... aktuelle Amplituden
    ea = EA(i);
    sa = SA(i);
    % ... Rissöffnungsspannung
    sop = sa*( 0.535*cos(pi/2*sa/sigF) - 0.344* sa/sigF);
    sop = max( [-sa min( [sa, sop] )] ); 
    % ... Rissschließdehnung
    ecl = -ea + (sop+sa)/E + 2 * ((sa+sop)/(2*Ks))^(1/ns);   
    % ... Rissschließspannung
    scl = closesig(sa,ea,ecl,E,Ks,ns);
    scl = max( [-sa min( [sa scl] )] );
    % ... effektive Verzerrungsenergiedichte
    Weff(i) =  (sa-scl)^2/(2*E) + (sa-scl)*((ea-ecl)-(sa-scl)/E)/(1+ns);
end

% -------------------------------------------------------------------------
% Steigung der Weff WL aus lin. Regression der log. Werte;
% ... log Werte
Nlin = N;
N = log10(N);
Wefflin = Weff;
Weff = log10(Weff);
% ... lineare Regression
p = polyfit(Weff,N,1);
% ... Auslesen Werte
mJ = -p(1);                      % Steigung
Q = 10^p(2);                     % Lebendsdauer bei Weff = 1;

% -------------------------------------------------------------------------
% a/c Verhältnis bestimmen sodass a/c = konst, d(a/c)/dn = 0

% ... Init a/c
azc = 1;
% ... Geometriefaktoren
[Y1A,Y1B] = Y_OFR_TYPA_ZD(azc,0);
% ... Iteration a/c
f = abs((Y1A/Y1B)^(2*mJ) - azc);
while f > 1e-4
    azc = 0.5 * ( azc + (Y1A/Y1B)^(2*mJ)); % Neuer Wert = mittelwert aus alt+a/c aus Rissfortschritt
    [Y1A,Y1B] = Y_OFR_TYPA_ZD(azc,0);
    f = abs((Y1A/Y1B)^(2*mJ) - azc);
end
% Hier ist Bisektion vielleicht stabiler

% -------------------------------------------------------------------------
% Q bestimmen mit Geometriefunktion (Lebensdauer bei PZ = 1)
Q = Q * (2*pi*Y1A^2)^mJ;

% -------------------------------------------------------------------------
% PZ werte bestimmen
PZ = 2* pi * Y1A^2 * Wefflin;

% -------------------------------------------------------------------------
% Kurzrissinitiierungsschwellwert
% % ... PZ Dauerfestigkeit
% PZD0 = 2*pi*Y1A^2*10^Weff(ndp);
% % ... Neuer Abknickpunkt
% ND0 = Q * PZD0^(-mJ);

% -------------------------------------------------------------------------
% Plot PZ WL
% plotPzWL(PZ,Nlin,mJ,Q,1,1e6);

end



function scl = closesig(s,e,ecl,E,Ks,ns)
% Berechnung Rissschließdehnung aus RO Parameter

% Sartwert
scl = s - E*(e-ecl);
% Fehler
f = e - (s-scl)/E - 2 *((s-scl)/(2*Ks))^(1/ns)- ecl;
% Newton
while abs(f) > 1e-7
    % Ableitung
    df = 1/E - 2/ns * ((s-scl)/(2*Ks))^(1/ns-1)*(-1/(2*Ks));
    % neuer Wert
    scl = scl - 1/df * f;
    % Fehler
    f = e - (s-scl)/E - 2 *((s-scl)/(2*Ks))^(1/ns)- ecl;
end
end

function [YIA,YIB] = Y_OFR_TYPA_ZD(azc,azr)
% Funktion gibt risslängenabhängige Geometriefunktionen in abhängigkeit des
% Risslängenverhältnisses a/c und Geometrieverhältniss a/rho
%
% Geometriefunktionen sind für Oberflächenriss unter Zugbeanspruchung 
% Risstiefe a=1, Zugspannung sx=1
% Punkt A - tiefster Punkt
% Punkt B - bei theta = 10°
% 
% Werte werden interpoliert aus Stützstellen berechnet von Hertel.
% Siehe Diss Hertl.
%
% INPUT:
% azc        - Verhältniss a zu c (Stützstellen 1,3/4,2/3,1/2,1/4)
% azr        - Verhältniss a zu rho (fiktiver Kerbradius) (Stützstellen 0,1/4,1/2,1,2,4)
%
% OUTPUT:
%  Y1A,Y1B   - Geometriefunktionen Mode I Stelle A und B
%
%__________________________________________________________________________


% -------------------------------------------------------------------------
% Stützstellen (aus Diss Hertel)
% % ... a zu c
ac = [ 0.25 0.5 0.667 0.75 1 ];
% ... a zu rho
ar = [ 0 0.25 0.5 1 2 4];

% ... Geometriefunktionen aus Abaqus
%  a/c= 0.25    0.5     0.667    0.75    1            a/rho
f1a = [1.0270, 0.8846, 0.7995, 0.7605, 0.6601;...   %  0
       0.6629, 0.5704, 0.5143, 0.4889, 0.4256;...   %  0.25
       0.5223, 0.4497, 0.4054, 0.3853, 0.3359;...   %  0.5
       0.3947, 0.3402, 0.3065, 0.2913, 0.2543;...   %  1
       0.2900, 0.2503, 0.2254, 0.2142, 0.1871;...   %  2
       0.2095, 0.1810, 0.1629, 0.1548, 0.1353];     %  4
f1b = [0.5789, 0.6883, 0.7165, 0.7227, 0.7244;...   %  0
       0.3867, 0.4596, 0.4782, 0.4823, 0.4857;...   %  0.25
       0.3053, 0.3635, 0.3782, 0.3816, 0.3848;...   %  0.5
       0.2311, 0.2753, 0.2866, 0.2891, 0.2920;...   %  1
       0.1704, 0.2024, 0.2109, 0.2128, 0.2152;...   %  2
       0.1236, 0.1462, 0.1524, 0.1538, 0.1557];     %  4

% -------------------------------------------------------------------------
% Input außerhalb des Interpolationsbereichs abfangen
% ... a zu c
azc(azc > 1) = 1;
azc(azc<0.25) = 0.25;
azr(azr>4) = 4;
azr(azr<0) = 0;
% if azc > 1
%     azc = 1;
% elseif azc < 0.25
%     azc = 0.25;
% end
% % ... a zu rho
% if azr > 4
%     azr = 4;
% elseif azr < 0
%     azr = 0;
% end

% -------------------------------------------------------------------------
% Input hat verschiedene Länge
lr = length(azr);
lc = length(azc);
if lr ~= lc
    if lr > lc
        azc = [azc,repelem(azc(end),lr-lc)];
    else
        azr = [azr,repelem(azr(end),lc-lr)];
        lr = lc;
    end
end

% -------------------------------------------------------------------------
% Interpolation von zwischenwerten

YIA = akima2d(6,5,ar,ac,f1a,lr,azr,azc);
YIB = akima2d(6,5,ar,ac,f1b,lr,azr,azc);

end % Ende Funktion