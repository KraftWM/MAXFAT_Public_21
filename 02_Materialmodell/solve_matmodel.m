function [OUT,EPSP,ZVAR1] = solve_matmodel(matmodell,ntens,ndi,para,M,inkflag,load,varargin)
% Funktion Integriert ein Materialmodell und gibt Spannungen und Dehnungen
% zurück
% INPUT:
%    matmodell     - (str) "Jiang","OhnoWang","Chaboche","KarimOhno","OWT"
%                          "Doring" -> bei Döring Modell gut aufpassen mit
%                          den Parametern, keine Entfestigung bei
%                          Spannungssteurung und nur 3D & ESZ
%    ntens         - (int) Anzahl Tensor komponenten
%    ndi           - (int) Anzahl Tensor komponenten auf diagonalen
%    para          - (double,array) Parameter des Modells
%    M             - (int) Anazhl Backstresstensoren
%   inkflag        - (int/bool) 0 spannungssteurung 1 dehnungssteuerung
%   load           - (double array) lastpfad e R ^(ntens,numdata)
% varargin         - optional, Startparameter
%
% OUTPUT:
% OUT              - Spannungen oder Dehnngen, je nach Lastvorgabe
% EPSP             - Plastische Dehungen
% ZVAR1            - Zustandsvariablen am Ende
%__________________________________________________________________________

% ... Init Zustandsvariablen
if nargin == 8
    ZVAR0 = varargin{1};
else
    ZVAR0 = init_matmodel(matmodell,ntens,para,M);
end

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

% ... auslesen Spannungen oder Dehnungen
OUT = ZVAR(1:ntens,:);
EPSP = ZVAR(ntens+1:2*ntens,:);
ZVAR1 = ZVAR(:,numdata);
end % Ende Funktion