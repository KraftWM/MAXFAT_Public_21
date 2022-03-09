function [ESIG,DLZ] = superposition(L,c,ndl)
% Funktion erzeugen pseudo elastischer spannungsverlauf im Ebenen
% Spannungszustand
%
% INPUT:
% L       - Lastfolge (double array, size: nkana x ndata)
% c       - �bertragungsfaktoren (double array, size: 3 x nkana)
% ndl     - Anzahl Durchl�ufe durch die Lastfolge, die in der
%           Kerbsimulation berechnet werden sollen
%
% OUTPUT:
% ESIG    - Verlauf pseudo elastische Spannungen 
%           (double array, size: 3 x ndata)
% DLZ     - Durchlaufz�hler
% 
%
% -------------------------------------------------------------------------
% Autor: Jan Kraft (TU Darmstadt)
% Stand: Januar 21
% -------------------------------------------------------------------------

%% Check Input
[ndata, nkana] = checkInput(L,c,ndl);

%% Initialisieren der Lastfolge
[LAF,LZW] = init_lastfolge(L,ndata);
ndata_af = size(LAF,2);                                 % Anzahl Daten Anfahrt
ndata_zw = size(LZW,2);                                 % Anzahl Daten Zwischenwerte

%% Superpostion 
ES = zeros(3,ndata);                                    % Speicher f�r Spannungen
ESAF = zeros(3,ndata_af);                               % Speicher f�r Anfahrt
ESZW = zeros(3,ndata_zw);                               % Speicher f�r Zwischenwerte
% ... Schleife �ber alle Kan�le
for i = 1:nkana
    ES = ES + c(:,i) .* L(i,:);
    if ndata_af > 0
        ESAF = ESAF + c(:,i) .* LAF(i,:);
    end
    if ndata_zw > 0
        ESZW = ESZW + c(:,i) .* LZW(i,:);
    end
end

%% Zusammenf�gen Durchl�ufe und erzeugen Durchlaufz�hler
% ... abfangen fehler bei ndl
if ndl < 1, ndl = 1; end                                % weniger als 1 Durchlauf
if ~floor(ndl)==ndl, ndl = floor(ndl); end              % ndl kein integer

% % ... erzeuge Spannungsverlauf
% ndata_dl = ndata_af + ndl*ndata + (ndl-1) * ndata_zw;   % Anzahl aller Datenpunkte
% ESIG = zeros(3,ndata_dl);                               % Speicher
% ESIG(:,1:ndata_af) = ESAF;                              % Startwerte
% % ... erzeuge Durchlaufz�hler
% DLZ = zeros(1,ndata_dl);                                % Speicher Durchlaufz�hler
% DLZ(1:ndata_af) = -1;                                   % Startwerte
% % ... schleife �ber Durchl�ufe bis auf den letzten
% first = ndata_af + 1;                                   % pointer erster nichtgesetzter Wert
% for dl = 1:ndl
%     last = first + ndata + ndata_zw - 1;                % pointer letzter zu setztender Wert
%     % ... Speichern Spannungen
%     ESIG(:,first:last) = [ES,ESZW];
%     % ... Speichern Durchlaufz�hler
%     DLZ(first:last) = (dl-1) + 1/(ndata_zw+ndata) : ...
%                                1/(ndata_zw+ndata) : ...
%                                dl;
%     first = last + 1;
% end
% ... Datenmenge
ndata_gesamt = ndata_af + ndata + (ndl-1)*(ndata_zw+ndata);
ESIG = zeros(3,ndata_gesamt);
DLZ = zeros(1,ndata_gesamt);
% ... erster Durchlauf
ESIG(:,1:ndata_af+ndata) = [ESAF,ES];
DLZ(1:ndata_af+ndata) = 1/(ndata_af+ndata): 1/(ndata_af+ndata) : 1;
% ... alle anderen Durchl�ufe
first = ndata_af+ndata+1;
last = first + ndata_zw + ndata - 1;
if ndl > 1
    for dl = 2 : ndl
        ESIG(:,first:last) = [ESZW,ES];
        DLZ(first:last) = (dl-1)+1/(ndata_zw + ndata):1/(ndata_zw + ndata):dl;
        first = last + 1;
        last = first + ndata_zw + ndata - 1;
    end
end

end % Ende superposition






% -------------------------------------------------------------------------
%                             Hilfsfunktionen:
%                                Inputcheck
% -------------------------------------------------------------------------
function [ndata, nkana] = checkInput(L,c,ndl)
% Bedingungen:
%  - Lastfolge und �bertragungsfaktoren sind numerischer Input
%  - gleiche Anzahl an Lastkan�len
%  - 3 �bertragungsfaktoren f�r jeden Lastkanal (ESZ)
%  - ndl ist positver integer
% -------------------------------------------------------------------------


% ... numerischer Input 
if ~isnumeric(L)
    eid = 'Input:notNumeric';
    msg = 'Die gegebene Lastfolge ist kein numerischer Input';
    throwAsCaller(MException(eid,msg))
elseif ~isnumeric(c)
    eid = 'Input:notNumeric';
    msg = 'Die gegebenen �bertragungsfaktoren sind kein numerischer Input';
    throwAsCaller(MException(eid,msg))
end

% ... keine NaN Werte
if sum(isnan(L),'all')
    eid = 'Input:isNaN';
    msg = 'Die gegebene Lastfolge enth�lt NaN Werte';
    throwAsCaller(MException(eid,msg))
elseif sum(isnan(c),'all')
    eid = 'Input:isNaN';
    msg = 'Die gegebenen �bertragungsfakoren enthalten NaN Werte';
    throwAsCaller(MException(eid,msg))
end

% ... 3 �bertragungsfaktoren
if size(c,1) ~= 3
    eid = 'Size:notEqual';
    msg = 'Falscher Spannungszustand in �bertragungsfaktoren';
    throwAsCaller(MException(eid,msg))
end

% ... Richtige Anzahl an Kan�len
if size(c,2) ~= size(L,1)
    eid = 'Size:notEqual';
    msg = 'Anzahl Lastkan�le f�r �bertragungsfaktoren und Last passen nicht zusammen';
    throwAsCaller(MException(eid,msg))
end

% ... ndl ist positiver interger
if ~isnumeric(ndl) || isnan(ndl) % kein numerischer Input
    eid = 'Input:notNumeric';
    msg = 'Die gegebene Anzahl an Durchl�ufen ist kein numerischer Input';
    throwAsCaller(MException(eid,msg))
end

% ... Alles richtig dann ermittle Arraygr��en
nkana = size(L,1);  % Anzahl Kan�le
ndata = size(L,2);  % Anzahl Datenpunkte

end % Ende checkInput



% -------------------------------------------------------------------------
%                             Hilfsfunktionen:
%                              init_lastfolge
% -------------------------------------------------------------------------
function [LAF,LZW] = init_lastfolge(L,ndata)
% Initialisieren der Lastfolge. Bestimmen der Anfahrt und der Zwischenwerte
% zwischen Durchl�ufen.
% 
% Es werden automatisch 9 Zwischenwerte gesetzt
%
% INPUT:
%    L  -> Lastfolge e R^(numkana,numdata)
% ndata -> anzahl Datenpunkte
%
% OUTPUT
%   LAF  -> Anfahrtwert des Lastfolge
%   LZW  -> Zwischenwerte zum �bergang auf neuen durchlauf
%
%__________________________________________________________________________


% ... Bestimme Anfahrt 
t = 0:0.1:0.9;                  % Pseudo Zeit
if sum(L(:,1).*L(:,1)) <= 1     % erster Punkt startet in 0 ?
    LAF = [];%t .* L(:,1);
else
    LAF = t .* L(:,1);
end


% ... Bestimme Zwischenwerte
dL = L(:,1) - L(:,ndata);       % Zu �berwindender Sprung
t = 0.1:0.1:0.9;                  % Pseudo Zeit
if sum(dL.*dL) <= 1
    LZW = [];%t .* dL + L(:,ndata);
else
    LZW = t .* dL + L(:,ndata);
end

end % Ende init_lastfolge



