function [L,DLZ] = zusammenfuegenLast(l,ndl)
% Funktion erzeugt mehre Durchläufe durch die Lastfolge
%
% INPUT:
% l       - Lastfolge (double array, size: ntens x ndata)
% ndl     - Anzahl Durchläufe durch die Lastfolge
%
% OUTPUT:
% L       - Zusammegesetzte Durchläufe
%           (double array, size: ntens x ndata*ndl)
% DLZ     - Durchlaufzähler
% 
%
% -------------------------------------------------------------------------
% Autor: Jan Kraft (TU Darmstadt)
% Stand: Februar 22
% -------------------------------------------------------------------------


%% Initialisieren der Lastfolge
ndata = size(l,2);                                      % Anzahl Datenpunkte
nkana = size(l,1);                                      % Anzahl Komponenten
[LAF,LZW] = init_lastfolge(l,ndata);
ndata_af = size(LAF,2);                                 % Anzahl Daten Anfahrt
ndata_zw = size(LZW,2);                                 % Anzahl Daten Zwischenwerte


%% Zusammenfügen Durchläufe und erzeugen Durchlaufzähler
% ... abfangen fehler bei ndl
if ndl < 1, ndl = 1; end                                % weniger als 1 Durchlauf
if ~floor(ndl)==ndl, ndl = floor(ndl); end              % ndl kein integer

% ... Datenmenge
ndata_gesamt = ndata_af + ndata + (ndl-1)*(ndata_zw+ndata);
L = zeros(nkana,ndata_gesamt);
DLZ = zeros(1,ndata_gesamt);
% ... erster Durchlauf
L(:,1:ndata_af+ndata) = [LAF,l];
DLZ(1:ndata_af+ndata) = 1/(ndata_af+ndata): 1/(ndata_af+ndata) : 1;
% ... alle anderen Durchläufe
first = ndata_af+ndata+1;
last = first + ndata_zw + ndata - 1;
if ndl > 1
    for dl = 2 : ndl
        L(:,first:last) = [LZW,l];
        DLZ(first:last) = (dl-1)+1/(ndata_zw + ndata):1/(ndata_zw + ndata):dl;
        first = last + 1;
        last = first + ndata_zw + ndata - 1;
    end
end

end % Ende superposition


% -------------------------------------------------------------------------
%                             Hilfsfunktionen:
%                              init_lastfolge
% -------------------------------------------------------------------------
function [LAF,LZW] = init_lastfolge(L,ndata)
% Initialisieren der Lastfolge. Bestimmen der Anfahrt und der Zwischenwerte
% zwischen Durchläufen.
% 
% Es werden automatisch 9 Zwischenwerte gesetzt
%
% INPUT:
%    L  -> Lastfolge e R^(numkana,numdata)
% ndata -> anzahl Datenpunkte
%
% OUTPUT
%   LAF  -> Anfahrtwert des Lastfolge
%   LZW  -> Zwischenwerte zum übergang auf neuen durchlauf
%
%__________________________________________________________________________


% ... Bestimme Anfahrt 
t = 0:0.1:0.9;                  % Pseudo Zeit
if sum(L(:,1).*L(:,1)) == 0     % erster Punkt startet in 0 ?
    LAF = [];%t .* L(:,1);
else
    LAF = t .* L(:,1);
end


% ... Bestimme Zwischenwerte
dL = L(:,1) - L(:,ndata);       % Zu überwindender Sprung
t = 0.1:0.1:0.9;                % Pseudo Zeit
if sum(dL.*dL) == 0
    LZW = [];%t .* dL + L(:,ndata);
else
    LZW = t .* dL + L(:,ndata);
end

end % Ende init_lastfolge



