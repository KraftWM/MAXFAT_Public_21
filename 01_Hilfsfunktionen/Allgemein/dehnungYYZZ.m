function [EPS, EPSP] = dehnungYYZZ(EPS,EPSP,nu)
% Funktion berechnet die nicht abgespeicherte Dehnungskomponten in YY und 
% ZZ Richtung und fügt sie der Dehnung hinzu
%
% INPUT:
%   EPS, EPSP -> Gesamt- & plastische Dehnung  aus sigma-tau,
%                wie folgt abgespeichert:
%             epsXX
%      eps =  2epsXY
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
% eyyp = - EPSP(1,:)./2 ;                                                    % plastischer Anteil der yy Dehnung aus 
ezzp = - EPSP(1,:)./2 ;                                                    % plastischer Anteil der zz Dehnung aus 
                                                                           % inkompressiobilität plastischer Def.
% elastische Komponente in ZZ Richtung
% eyye = -nu * EPSE(1,:);
ezze = -nu * EPSE(1,:);                                                    % Elastischer Anteil aus elastizitätsgesetz

% gesamte Komponente in ZZ Richtung
% eyy = eyye + eyyp;
ezz = ezze + ezzp;

% Zusammenschustern der gesamt- und plastischen dehnungen
if size(EPS,1) == 2
    EPS = [EPS(1,:);ezz;ezz;EPS(2,:)];
    EPSP = [EPSP(1,:);ezzp;ezzp;EPSP(2,:)];
else
    EPS = [EPS(1,:);ezz;ezz];
    EPSP = [EPSP(1,:);ezzp;ezzp];
end

end % Ende Funktion