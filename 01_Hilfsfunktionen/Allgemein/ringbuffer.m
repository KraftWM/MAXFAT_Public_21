classdef ringbuffer < handle
% Ringpuffer für double arrays der Größe  BufSize x VecSize 
%
%
% Erstellen:
% r = ringbuffer(Vektorgröße, 
%                Buffergröße (optional)
%                )
%
% 
% Vektorhinzufügen:
% bool = r.addVec(Vec);       1  erfolgreich, 0 nicht erfolgreich
%
% mehrere Vektoren hinzufügen:
% bool = r.addMultVec(Vec);   1  erfolgreich, 0 nicht erfolgreich
%
% Lesen Ersten (noch nicht gelesenen) Wert:
% vec = r.readVec();
%
% Lesen mehrere (noch nicht gelesene) Werte:
% vec = r.readMultVec();
% 
% Zurückgeben Wert im Leseindex (ohne Lesezugriff):
% vec = r.firstvalue()
% 
% Zurückgeben Wert im Schreibindex (ohne Lesezugriff):
% vec = r.lastvalue()
%
% Hilfsfunkionen:
% isEmpty()     - Testet ob Buffer leer ist
% isFull()      - Testet ob BUffer voll ist
% freeSpace()   - Gibt Anzahl noch freier Zellen zurück
% usedSpace()   - Gibt Anzahl der belegten Zellen zurück
%
% -------------------------------------------------------------------------


% -------------- Eigenschaften --------------------------------------------
properties (SetAccess = private, GetAccess = public)
    VecSize  (1,1) {mustBeInteger} = int64(1);                             % Vektorgröße
    BufSize (1,1) {mustBeInteger} = int64(5e5);                            % Buffergröße  
    fst (1,1) {mustBeInteger} = int64(0);                                  % Zeiger ersten Index (Leseindex)
    lst (1,1) {mustBeInteger} = int64(0);                                  % Zeiger letzten Index (Schreibindex-1) 
    data double {mustBeNumeric}                                            % Array für Daten    
end % Ende Eigenschaften


% -------------- Kontruktor -----------------------------------------------
methods (Access = public)
    function obj = ringbuffer(VSize,BSize)
        % ... Zeilen
        obj.VecSize = VSize;
        % ... Spalten
        if nargin == 2
            obj.BufSize = BSize;
        end
        % ... Speicher für Daten
        obj.data = NaN(obj.BufSize,obj.VecSize);        
    end % Ende Konstruktor
    
end % Ende Public Methoden
    
% -------------- Public Methoden ------------------------------------------
methods (Access = public)
    
    % Neuen Vektor hinzufügen
    function flag = addVec(obj,Vec)
        flag = 1;                                                          % Schreiben Erfolgreich
        % ... Teste ob Buffer voll
        if obj.isFull()
            flag = 0;                                                      % Schreiben nicht erfolgreich
            return;
        elseif obj.isEmpty()
            % inkrementiere ersten Index (index bei leerem Buffer sollte 0
            % sein)
            obj.fst = obj.fst + 1;
        end
        % ... Inkrementiere Zeiger auf letzten Wert
        obj.lst = mod(obj.lst,obj.BufSize)+1;
        l = obj.lst;
        % ... Speicher Vektor
        obj.data(l,:) = Vec;        
    end % Ende appendVec
    
    % Fügt mehrere Vektoren (array e R^(ndata,VecSize)) hinzu
    function flag = addMultVec(obj,Vec)
        flag = 1;                                                          % Schreiben Erfolgreich
        ndata = size(Vec,1);                                               % Größe neue Daten
        space = obj.freeSpace();                                           % Freier Speicher im Buffer
        
        % Berechne Indices um Vectoren dazuzufügen
        if space < ndata
            flag = 0;                                                      % Speichern nicht erfolgreich kein Speicher mehr da
        else
            % Buffer war Leer
            if obj.isEmpty()
                obj.fst = obj.fst + 1;                                     % init Leseposition
            end
            % Rechts dazufügen
            obj.data(mod(obj.lst:obj.lst+ndata-1,obj.BufSize)+1,:) = Vec;           
            % Neuer Schreibindex
            obj.lst = mod(obj.lst+ndata-1,obj.BufSize)+1;            
        end
        
    end
    
    % Ließt Werte
    function vec = readVec(obj)
        % ... Teste ob Buffer leer is
        if obj.isEmpty()
            vec = NaN;
        else
            % ... Lese Wert
            vec = obj.data(obj.fst,:);
            % ... Teste ob Buffer jetzt Leer is
            if obj.lst == obj.fst
                obj.lst = 0;
                obj.fst = 0;
            else
                % ... Inkrementiere Leseindex
                obj.fst = mod(obj.fst,obj.BufSize)+1;
            end
        end
            
    end % Ende readVec
    
    % Ließt mehrere Vektoren (nämlich die nächsten numrows Zeilen)
    function vec = readMultVec(obj,numrows)
        % ... Teste ob Buffer leer is
        if obj.isEmpty()
            vec = NaN;
        else
            % ... Wieviele Werte sind im Buffer
            used = obj.usedSpace();
            % ... Weniger Werte Im Buffer als gelesen werden soll
            if used < numrows
                % ... Lese einfach alle Werte die noch da sind
                vec = obj.data(mod(obj.fst-1:obj.fst+used-2,obj.BufSize)+1,:);
                % ... Inkrementiere Leseindex
                obj.fst = mod(obj.fst+used-2,obj.BufSize)+1;
            else
                % ... Lese Werte
                vec = obj.data(mod(obj.fst-1:obj.fst+numrows-2,obj.BufSize)+1,:);
                % ... Inkrementiere Leseindex
                obj.fst = mod(obj.fst+numrows-2,obj.BufSize)+1;
            end
            % ... Teste ob Buffer jetzt Leer is
            if obj.lst == obj.fst
                obj.lst = 0;
                obj.fst = 0;
            else
                % ... Inkrementiere Leseindex
                obj.fst = mod(obj.fst,obj.BufSize)+1;
            end
        end
        
    end
 
    % Teste ob Buffer voll ist
    function full = isFull(obj)
        full = mod(obj.lst,obj.BufSize) + 1 == obj.fst;
    end
    
    % Teste ob Buffer leer ist
    function empty = isEmpty(obj)
        empty = obj.fst == 0;
    end
    
    % Freier Speicher im Buffer
    function space = freeSpace(obj)
        if obj.isFull()
            space = 0;
        elseif obj.isEmpty()
            space = obj.BufSize;
        elseif obj.fst > obj.lst
            space = obj.fst - obj.lst - 1;
        elseif obj.fst < obj.lst
            space = obj.BufSize - obj.lst + obj.fst - 1;
        end        
    end 
    
    % Belegter Speicher im Buffer
    function used = usedSpace(obj)
        if obj.isFull()
            used = obj.BufSize;
        elseif obj.isEmpty()
            used = 0;
        elseif obj.fst > obj.lst
            used = obj.BufSize - obj.fst + 1 + obj.lst;
        elseif obj.fst < obj.lst
            used = obj.lst - obj.fst + 1;
        end 
    end
    
    % Gibt Wert im Schreibindex zurück (letzten Wert, ohne manipulation der
    % Indices, wird nicht als Lesezugriff gewertet)
    function vec = lastvalue(obj)
        % ... Testen ob Buffer leer ist
        if obj.isEmpty()
            vec = NaN;
        else
        % ... Lese Letzten Wert
            vec = obj.data(obj.lst,:);
        end
    end
    
    % Gibt Wert im Leseindex zurück (ersten Wert, ohne manipulation der
    % Indices, wird nicht als Lesezugriff gewertet)
    function vec = firstvalue(obj)
        % ... Testen ob Buffer leer ist
        if obj.isEmpty()
            vec = NaN;
        else
        % ... Lese Letzten Wert
            vec = obj.data(obj.fst,:);
        end
    end
    
end % Ende Public Methoden


end % Ende Klassendefinition Ringbuffer