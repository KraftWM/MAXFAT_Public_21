% Skript fügt alle Ordner und Unterordner dem Matlabpfad hinzu 


% alle elemente in Ordner
elem = dir;
% anzahl elemente
numelem = size(elem,1);
% schleife über elemente
for i = 1:numelem
    % aktuelles element
    e = elem(i);
    % element ist ordner
    if e.isdir
        % element ist nicht vor,zurück oder git
        if ~strcmp(e.name,'.') && ~strcmp(e.name,'..') && ~strcmp(e.name,'.git')
            % ordner und unterordner pfad
            p = genpath(e.name);
            % hinzufügen pfad
            addpath(p);
        end
    end
end

% Lösche Variablen
clear e elem i numelem p


