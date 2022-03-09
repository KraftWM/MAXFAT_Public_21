function [SUBDATA,In,Jn,Kn] = CutOutHyst(DATA,I,J,K,NF)
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

% Reserviere Speicher für Subdata
SUBDATA = zeros(size(DATA,1),K-I+1);
% ... Startwert Leseindex in Data
i = I;
% ... Speicherindex in Subdata
s = 1;
% ... Schleife
while i <= K
    
    % ... merke Umkehrpunkt
    if i == J
        Jn = s;
    end
    
    % ... Speichern Data
    SUBDATA(:,s) = DATA(:,i);
    
    % ... Inkrementiere Speicherindex Subdata
    s = s + 1;
    
    % ... Inkrementiere Leseindex Data
    i = NF(i);
    
end % Ende Schleife

% ... lösche Nullen
SUBDATA = SUBDATA(:,1:s-1); 

% ... Neue Indices der Umkehrpunkte
In = 1;
Kn = s-1;
% Jn = J-I+1;

end % Ende Funktion