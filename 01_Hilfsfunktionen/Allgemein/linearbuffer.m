classdef linearbuffer < handle
% linearer puffer für double arrays der Größe BufSize x VecSize 
%
% -------------------------------------------------------------------------

% -------------- Eigenschaften --------------------------------------------
properties (SetAccess = private, GetAccess = public)
    VecSize  (1,1) {mustBeInteger} = int64(1);                             % Vektorgröße
    BufSize (1,1) {mustBeInteger} = int64(5e5);                            % Buffergröße  
    data double {mustBeNumeric}                                            % Array für Daten    
    fst (1,1) {mustBeInteger} = int64(0);                                  % Zeiger ersten Index (Leseindex)
    lst (1,1) {mustBeInteger} = 0;                                         % Zeiger auf zuletzt gefüllten Wert
end % Ende Eigenschaften
% -------------- Kontruktor -----------------------------------------------
methods (Access = public)
    function obj = linearbuffer(VSize,BSize)
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
    
    % Hinzufügen von Daten
    function flag = addVec(obj,Vec)
        flag = 1;                                                          % Schreiben Erfolgreich
        ndata = size(Vec,1);                                               % Größe neue Daten
        space = obj.freeSpace();                                           % Freier Speicher
        % Berechne Indices um Vectoren dazuzufügen
        if space < ndata           
            flag = 0;                                                      % Speichern nicht erfolgreich kein Speicher mehr da
            if space > 0
                % hinzufügen
                obj.data(obj.lst+1:obj.lst+space,:) = Vec(1:space,:);           
                % Neuer Schreibindex
                obj.lst = obj.lst+space;
            end
        else
            if obj.isEmpty()
                obj.fst = 1;
            end
            % hinzufügen
            obj.data(obj.lst+1:obj.lst+ndata,:) = Vec;           
            % Neuer Schreibindex
            obj.lst = obj.lst+ndata;            
        end
    end
    
    % Lesen von Daten
    function vec = readVec(obj,numrows)
       % ... Teste ob Buffer leer is
        if obj.isEmpty()
            vec = NaN;
        else
            % ... Wieviele Werte sind im Buffer
            used = obj.usedSpace();
            % ... Weniger Werte Im Buffer als gelesen werden soll
            if used < numrows
                % ... Lese einfach alle Werte die noch da sind
                vec = obj.data(obj.fst:obj.fst+used-1,:);
                % ... Inkrementiere Leseindex
                obj.fst = obj.fst+used-1;
            else
                % ... Lese Werte
                vec = obj.data(obj.fst:obj.fst+numrows-1,:);
                % ... Inkrementiere Leseindex
                obj.fst = obj.fst+numrows-1;
            end
            % ... Teste ob Buffer jetzt Leer is
            if obj.lst == obj.fst
                obj.lst = 0;
                obj.fst = 0;
            else
                % ... Inkrementiere Leseindex
                obj.fst = obj.fst+1;
            end
        end 
    end
    
    % Ermitteln anzahl freie Zellen
    function space = freeSpace(obj)
        space = obj.BufSize - obj.lst;
    end
    
    % Ermittelt anzahl noch zu lesender Werte
    function used = usedSpace(obj)
        used = obj.lst - obj.fst + 1;
    end
    
    % Ermittelt ob Buffer voll ist
    function full = isFull(obj)
        full = obj.lst >= obj.BufSize;
    end
    
    % Ermittelt ob Buffer Leer ist
    function empty = isEmpty(obj)
        empty = obj.fst == 0;
    end
    
end % Ende Methoden

end % Ende Klassendefinition