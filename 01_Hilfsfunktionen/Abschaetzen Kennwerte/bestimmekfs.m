function [kds,k] = bestimmekfs(sf,ef,b,c,tf,gf,b0,c0,sigF,nu,E, ...
    flag)
% Berechnet k für Fatemi Socie Parameter
%
% INPUT:
% sf,ef,b,c     - Parameter DWL
% tf,gf,b0,c0   - Parameter GWL
% sigF          - Fließspannung
% nu            - Querdehnzahl (elastisch)
% flag          - Bestimmen des Parameters 
%                 0 - Integrales Mittel
%                 1 - logarithmisches Mittel
%
% OUTPUT:
% kds           - k/sigF
% k             - k Parameter Fatemi Socie Parameter
%
% ANMERKUNG:
% Vorgehen aus D. McClaflin, A. Fatemi 2003 - Torsional deformation
%              and fatigue of hardened steel including mean stress
%              and stress gradient effects
%
% Bestimme k für N = 10^2...10^7 und nehme Mittelwert
%
%______________________________________________________________

% ... automatische flagge
if nargin == 11
    flag = 0;
end


% ... berechne Schubmodul
G = E/(2*(1+nu));

% ... Bestimme Punkte und k über Mittelwert (log Verteilt)
if flag
    Nmin = 3;
    Nmax = 6;
    ndp = Nmax-Nmin+1;
    N = logspace(Nmin,Nmax,ndp);
    fak = ((tf/G*(2*N).^b0 + gf*(2*N).^c0) ./ ...
        ( (1+nu)*sf/E*(2*N).^b + 1.5*ef*(2*N).^c ) - 1);
    k = fak * ...
        2*sigF./(sf * (2*N).^b);
    k = 1/length(k) * sum(k);

else
    % ... bestimme Integralen Mittelwert über Simpson Integration
    N1 = 1e3; N2 = 1e7;
    nip = 1001;         % Anzahl Inkremente + 1
    h = (N2-N1)/(nip-1);
    N = N1:h:N2;
    ksimpson = ((tf/G*(2*N).^b0 + gf*(2*N).^c0) ./ ...
        ((1+nu)*sf/E*(2*N).^b + 1.5*ef*(2*N).^c ) - 1)* ...
        2*sigF./(sf * (2*N).^b);
    F = ksimpson(1) + ksimpson(nip);
    for i = 2:nip-1
        if mod(i,2) == 0
            F = F + 2*ksimpson(i);
        else
            F = F + 4*ksimpson(i);
        end
    end
    k = h/3 * F/(N2-N1);
end
% ... k/sigF
kds = k/sigF;

end % Ende bestimme kfs