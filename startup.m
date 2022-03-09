% Skript f�gt alle Ordner und Unterordner dem Matlabpfad hinzu 


% alle elemente in Ordner
elem = dir;
% anzahl elemente
numelem = size(elem,1);
% schleife �ber elemente
for i = 1:numelem
    % aktuelles element
    e = elem(i);
    % element ist ordner
    if e.isdir
        % element ist nicht vor,zur�ck oder git
        if ~strcmp(e.name,'.') && ~strcmp(e.name,'..') && ~strcmp(e.name,'.git')
            % ordner und unterordner pfad
            p = genpath(e.name);
            % hinzuf�gen pfad
            addpath(p);
        end
    end
end

% L�sche Variablen
clear e elem i numelem p


