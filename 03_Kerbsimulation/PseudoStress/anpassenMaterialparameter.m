function [para,q,ep_M] = anpassenMaterialparameter(maxsigv,...
                                          typ, eindkerb,modell,M,...
                                          E,nu,verfahren_flag,...
                                          varargin)
% Funktion bestimmt neue Parameter wenn Spannungsgrenzwert des Modells
% kleiner ist als die Maximale Spannung die in der Lastfolge aufgebracht
% werden soll. Nur für Werkstoff & Pseudo Stress Ansatz 
%
% INPUT:
% maxsigv        - maximale Vergleichsspannung der Lastfolge
% typ            - Werkstoff oder Strukturfließflächen (siehe ro2para $ 
%                  bfk2para)
% eindkerb       - 1D Kerbsimulation für BFK ('selbst' wenn Fließkurve 
%                  vorgegeben wird)
% modell         - Name der Plastimodells
% M              - Anzahl der Backstresstensoren
% E,nu           - Elastizität
% verfahren_flag - Flagge für Verfahren zum bestimmen der Parameter (siehe 
%                  ro2para $ bfk2para)
% varargin       - Variabler Input je nach Verfahren
%                  eindkerb = 'selbst' varargin(1) = bfk
%                                      varargin(2) = dummy
%                           = sonst    varargin(1) = Ks
%                                      varargin(2) = ns
%                 
%                  verfahren_flag = 2  varargin(3) = q 
%                                 = 3  varargin(3) = q 
%                                      varargin(4) = ep_M 
%                             
%                  ggf.                varargin(5) = Kp          
% OUTPUT:
% para    - Angepasste Parameter
% q,ep_M  - Werte zum anpassen der Parameter
% -------------------------------------------------------------------------
tol = maxsigv/10;                  % Toleranz um die die maximale Spannung hochgesetzt wird um numerische Fehler zu vermeiden
% Handle Variablen Input
switch eindkerb
    case 'selbst'
        bfk = varargin{1};
        zz = 2;
    otherwise
        Ks = varargin{1};
        ns = varargin{2};
        zz = 2;
end
switch verfahren_flag
    case {1,2,3}
        q = varargin{zz+1};
    otherwise
        msg = 'Nimm ein anderes Verfahren';
        error(msg);
end
if length(varargin) == 5
    Kp = varargin{5};
else
    Kp = NaN;
end
% Bestimme Stützstelle für maximale plastische Dehnung aus maximaler
% Spannung
switch eindkerb
    % ... aus Fließkurve
    case 'selbst'
        idx = find(bfk(:,1) >= maxsigv+tol,1);
        if isempty(idx)
            msg = 'Fließkurve reicht nicht weit genug';
            error(msg);
        else
            epmax = bfk(idx,2);
        end
    % ... aus Ramberg osgood
    otherwise
        % ... Werkstoff
        if strcmp(typ,'werkstoff')
            epmax = ((maxsigv+tol)/Ks)^(1/ns);
        % ... Strukturff (muss PseudoStress)
        else
            % ... Kerbnäherung
            sig = 0:1:2000;
            epsp = (sig./Ks).^(1/ns);
            switch eindkerb
                case 'Neuber'
                    esig = uniax_neuber(E,sig,epsp);
                case 'ESED'
                    esig = uniax_esed(E,ns,sig,epsp);
                case 'Seeger Beste'
                    esig = uniax_seegerbeste(E,sig,epsp,Ks,ns,Kp);
                case 'Seeger Heuler'
                    esig = uniax_neuberstern(E,sig,epsp,Ks,ns,Kp);
                case 'Neuber Stern'
                    esig = uniax_neuberstern(E,sig,epsp,Ks,ns,Kp);
                otherwise
                    msg = 'angegebenes Verfahren nicht implementiert';
                    error(msg)
            end
            epmax = interp1(esig,epsp,maxsigv+tol);
        end
end

% Bestimme angepasste Werte für Verfahren
switch verfahren_flag
    % ... Äquidistant in Spannungen
    case 1
        q = maxsigv + tol;
        ep_M = NaN;
        varinput{1} = q;
        varinput{2} = ep_M;
        varinput{3} = eindkerb;
        varinput{4} = Kp;
    % ... Geometrische Reihe
    case 2
        q = (0.0001/epmax)^(1/M);
        ep_M = NaN;
        varinput = cell(1,1);
        varinput{1} = q;
        varinput{2} = ep_M;
        varinput{3} = eindkerb;
        varinput{4} = Kp;
    % ... Verfahren von Simon
    case 3
        ep_M = epmax;
        varinput = cell(1,2);
        varinput{1} = q;
        varinput{2} = ep_M;
        varinput{3} = eindkerb;
        varinput{4} = Kp;
end
% Unterscheide aus Ramberg Osgood oder aus Fließkurve
switch eindkerb
    % ... Aus Fließkurve
    case 'selbst'
        para = bfk2paraV2(typ,bfk,M,modell,E,nu,verfahren_flag,varinput{:});
    % ... Aus Ramberg Osgood
    otherwise
        para = ro2paraV2(typ,...
                          E,nu,Ks,ns,M,...
                          modell,verfahren_flag,...
                          varinput{:});
end
end