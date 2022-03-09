function [Weff_stuetz,d] = dwl2weffwl(E,ef,sf,bf,cf,ND,sigF,plotopt)
% -----------------------------------------------------------------
% Funktion bestimmt Weff Wöhlerlinie aus DWL(glatte Probe R=-1)
% für Pz
% INPUT:
% E             - E Modul
% ef,sf,b,c     - DWL
% OUTPUT:
% Weff_stuetz   - Stützstelle bei N = 1
% d             - Neigung
% -----------------------------------------------------------------
% ... Ramberg Osgood Parameter aus Konsostenzbedingung
ns = bf/cf;
Ks = sf/ef^ns;

%         ns = nstrich;
%         Ks = Kstrich;

% ... Speicher
ndp = 100;
N = zeros(1,ndp);   % Speicher Lebensdauer Wöhlerlinie
Weff = zeros(1,ndp);% Speicher effektive Verzerrungserngiedichte WL

% ... Datenpunkte erzeugen
for i = 1 : ndp
    % ... Lebensdauern von 100 bis ND
    N(i) = 100 * (ND/100)^((i-1)/(ndp-1));

    % ... Spannungs- und Dehnungsamplitude aus MBC Gleichung (ohne nst)
    sa = sf*(2*N(i))^bf;
    ea = sa/E + ef*(2*N(i))^cf;

    % ... Rissöffnungsspannung
    sop = sa*( 0.535*cos(pi/2*sa/sigF) - 0.344* sa/sigF);
    sop = max( [-sa min( [sa, sop] )] );

    % ... Rissschließdehnung
    ecl = -ea + (sop+sa)/E + 2 * ((sa+sop)/(2*Ks))^(1/ns);

    % ... Rissschließspannung
    scl = closesig(sa,ea,ecl,E,Ks,ns);
    scl = max( [-sa min( [sa scl] )] );

    % ... effektive Verzerrungsenergiedichte
    %             Weffe(i) = (sa-scl)^2/(2*E);
    %             Weffp(i) = (sa-scl)*((ea-ecl)-(sa-scl)/E)/(1+ns);
    Weff(i) =  (sa-scl)^2/(2*E) + (sa-scl)*((ea-ecl)-(sa-scl)/E)/(1+ns);

end % Ende Datenpunkte WL

% -----------------------------------------------------------------
% Steigung der Weff WL aus lin. Regression der log. Werte;
% ... log Werte
N = log10(N);
Weff = log10(Weff);
% ... lineare Regression
p = polyfit(Weff,N,1);
% ... Auslesen Werte
d = 1/p(1);                      
Weff_stuetz = 10^(-p(2)/p(1));   

% -----------------------------------------------------------------
% Falls geplotet werden soll
if plotopt
    figure,grid on, hold on
    plot(10.^N,10.^Weff,'x')
    plot([1 1e6],[Weff_stuetz,1e6^d*Weff_stuetz])
    set(gca,'XScale','log','YScale','log')
end
end % Ende Bestimmen Weff Wöhlerlinie


% ... bestimme schließspannung für einfaches ramberg osgood
function scl = closesig(s,e,ecl,E,Ks,ns)
% -----------------------------------------------------------------
% Berechnung Rissschließdehnung aus RO Parameter
% !!! mit verdopplungsgesetz auf halbästen
% INPUT:
%      s       - Spannung
%      e       - Dehnung
%      ecl     - Schließdehnung
%      E       - E-Modul
%    Ks,ns     - Parameter Ram. Osgood
% OUTPUT:
%     scl      - Schließspannung
% -----------------------------------------------------------------

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