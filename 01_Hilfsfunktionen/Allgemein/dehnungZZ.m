function [EPS, EPSP] = dehnungZZ(EPS,EPSP,nu)
% Funktion berechnet die nicht abgespeicherte Dehnungskomponte in ZZ 
% Richtung und fuegt sie der Dehnung hinzu
%
% INPUT:
%   EPS, EPSP -> Gesamt- & plastische Dehnung  aus ESZ,
%                wie folgt abgespeichert:
%             epsXX
%      eps =  epsYY
%            2epsXY
%
% OUTPUT:
%
%   EPS, EPSP -> Gesamt- & plastische Dehnung  aus ESZ,
%                wie folgt abgespeichert:
%             epsXX
%      eps =  epsYY
%             epsZZ
%            2epsXY
% _________________________________________________________________________

% elastische Dehnung
EPSE = EPS - EPSP;  

% plastische Komponente in ZZ Richtung
ezzp = - EPSP(1,:) - EPSP(2,:);                                            % plastischer Anteil der zz Dehnung aus 
                                                                           % inkompressiobilitaet plastischer Def.
% elastische Komponente in ZZ Richtung
ezze = -nu/(1-nu) .* ( EPSE(1,:) + EPSE(2,:) );                            % Elastischer Anteil aus Elastizitaetsgesetz

% gesamte Komponente in ZZ Richtung
ezz = ezze + ezzp;

% Zusammenschustern der Gesamtdehnungen
EPS = [EPS(1:2,:);ezz;EPS(3,:)];

% Zusammenschustern der plastischen Dehnungen
EPSP = [EPSP(1:2,:);ezzp;EPSP(3,:)];

end % Ende Funktion