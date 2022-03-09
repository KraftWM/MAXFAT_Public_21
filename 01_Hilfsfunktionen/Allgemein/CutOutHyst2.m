function [SUBDATA,In,Jn,Kn] = CutOutHyst2(DATA,I,J,K,Nf)
% Funktion schneidet schon geschlossene Hysteresen aus Zeitfolge aus
%
% INPUT:
% DATA     - Kompleter Datensatz
%  I       - Start Hysterese
%  J       - Umkehrpunkt der Hysterese
%  K       - Endpunkt
% NF       - Zeiger auf Nachfolger
%
% OUTPUT:
% SUBDATA  - Ausgeschnittene Hysteres
% In,Jn,Kn - Neue Indices der Umkehrpunkte
%__________________________________________________________________________
% Reserviere Speicher für Subdata (Annahme: Hysteresen kleiner als 10001
% Datenpunkte)
S = zeros(1,10000);
% ... Startwert Leseindex in Data
i = I;
% ... Speicherindex in Subdata
s = 1;
% ... Schleife
while i ~= K    
    % ... merke Umkehrpunkt
    if i == J
        Jn = s;
    end    
    % ... Speichern Data
    S(s) = i;    
    % ... Inkrementiere Speicherindex Subdata
    s = s + 1;   
    % ... Inkrementiere Leseindex Data
    i = Nf(i);    
end % Ende Schleife
% K dazufügen
S(s) = i;
% ... lösche Nullen
S(S==0) = [];
% Ausgeschnittener Datenteil
SUBDATA = DATA(:,S); 
% ... Neue Indices der Umkehrpunkte
In = 1;
Kn = s;

end % Ende Funktion