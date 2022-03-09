function [sigepsfile,numdata] = solve_matmodel_file(jobname,outpath,...
                                          matmodell,ntens,ndi,para,M,inkflag,load,ndl)
% Funktion Integriert ein Materialmodell und schreibt Ergebnisse in *.sig
% Datei
% INPUT:
%    jobname       - (str) name der Rechnung
%    outpath       - (str) Name des Speicherorts
%    matmodell     - (str) "Jiang","OhnoWang","Chaboche","KarimOhno",
%                          "Doring" -> bei Döring Modell gut aufpassen mit
%                          den Parametern, keine Entfestigung bei
%                          Spannungssteurung und nur 3D & ESZ
%    ntens         - (int) Anzahl Tensor komponenten
%    ndi           - (int) Anzahl Tensor komponenten auf diagonalen
%    para          - (double,array) Parameter des Modells
%    M             - (int) Anazhl Backstresstensoren
%   inkflag        - (int/bool) 0 spannungssteurung 1 dehnungssteuerung
%   load           - (double array) lastpfad e R ^(ntens,numdata)
%   ndl            - (int) Durchlaufe durch load, die simuliert werden
%                    sollen
% varargin         - optional, Startparameter
%
% OUTPUT:
% sigepsfile       - *.sig binary file mit Spannungen,Dehungen und
%                    Durchlaufzähler
% numdata          - Länge Lastfolge
%__________________________________________________________________________

% ... Init Zustandsvariablen
ZVAR0 = init_matmodel(matmodell,ntens,para,M);

% ... Materialfunktion
switch matmodell
    case 'Chaboche'
        if ntens == 6 || ntens == 3
            matfun = @chabocheRR2_mex;
%             matfun = @chabocheRR2;
        else
            matfun = @chaboche;
        end
    case 'OhnoWang'
        if ntens == 6 || ntens == 3
            matfun = @ohnowangRR2_mex;
%             matfun = @ohnowangRR2;
        else
            matfun = @ohnowang;
        end
    case 'Jiang'
        matfun = @jiang;
    case 'KarimOhno'
        if ntens == 6 || ntens == 3
            matfun = @karimohnoRR2_mex;
%             matfun = @karimohnoRR2;
        else
            matfun = @karimohno;
        end
    case 'Doring'
        matfun = @doring;
    case 'OWT'
        matfun = @OWT;
    otherwise
        msg = 'Falsches Materialmodell';
        error(msg)
end % Ende Unterscheidung Materialmodelle

% ... Zusammensetzten aller Durchläufe durch load
[load,DLZ] = zusammenfuegenLast(load,ndl);

% ... Daten der Lastfolge/Speicher Zustandsvariablen
numdata = size(load,2);
ZVAR = zeros(size(ZVAR0,1),numdata);
ZVAR(:,1) = ZVAR0;

% ... Schleife über Lastfolge
for i = 2:numdata
    % ... Inkrement
    ink = load(:,i) - load(:,i-1);
    % ... Integration
    ZVAR(:,i) = matfun(ntens,ndi,ink,ZVAR(:,i-1),inkflag,para);
end

% ... Öffne Output Datei
sigepsfile = [outpath,'/',jobname,'.sig']; % Dateiname
fout = fopen(sigepsfile,'w');              % Öffne Datei

% ... Unterscheide Spannungszustände 
if ntens == 6 && ndi == 3 % 3D
    if inkflag                                 % Schreibe Datei
        fwrite(fout,[DLZ;ZVAR(1:ntens,:);load],'double');
    else
        fwrite(fout,[DLZ;load;ZVAR(1:ntens,:)],'double');
    end
elseif ntens == 3 && ndi == 2 % ESZ
    if inkflag                                                             % Schreibe Datei
        [EPS, ~] = dehnungZZ(load,ZVAR(ntens+1:2*ntens,:),para(2));
        fwrite(fout,[DLZ;ZVAR(1:ntens,:);EPS],'double');
    else
        [EPS, ~] = dehnungZZ(ZVAR(1:ntens,:),ZVAR(ntens+1:2*ntens,:),para(2));
        fwrite(fout,[DLZ;load;EPS],'double');
    end
elseif ntens == 2 && ndi == 1 % Sigma - Tau
    if inkflag                                                             % Schreibe Datei
        [EPS, ~] = dehnungYYZZ(load,ZVAR(ntens+1:2*ntens,:),para(2));
        fwrite(fout,[DLZ;ZVAR(1:ntens,:);EPS],'double');        
    else
        [EPS, ~] = dehnungYYZZ(ZVAR(1:ntens,:),ZVAR(ntens+1:2*ntens,:),para(2));
        fwrite(fout,[DLZ;load;EPS],'double');
    end
elseif ntens == 1 && ndi == 1 % reiner Zug
    if inkflag                                                             % Schreibe Datei
        [EPS, ~] = dehnungYYZZ(load,ZVAR(ntens+1:2*ntens,:),para(2));
        fwrite(fout,[DLZ;ZVAR(1:ntens,:);EPS],'double');        
    else
        [EPS, ~] = dehnungYYZZ(ZVAR(1:ntens,:),ZVAR(ntens+1:2*ntens,:),para(2));
        fwrite(fout,[DLZ;load;EPS],'double');
    end
else % Spannungszustand nicht erkannt
    msg = 'Spannungszustand nicht erkannt';
    error(msg)
end

% ... schließe Datei
fclose(fout);                             

end % Ende Funktion