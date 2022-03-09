function [ZVARneu] = chabochelang(desig,ZVAR,para,epara)
% Implementierung des Materialmodells nach chaboche, für pseudo
% stress aproach nach Lang 
%
% ! Nur Spannungssteuerung im ESZ
%
% QUELLE DEHNUNGSGESTEURTER IMPLIZITER EULER
% Kobayashi, Ohno 2002 - Implementation of cyclic plasticity models based
% on a general form of kinematic hardening
%
% QUELLE SPANNUNGSGESTEUERTER IMPLIZITER EULER:
% Eigener Algo.
%
%
% TODO: ersetzte zentrale Differenzen in Newtonverfahren durch analytische
% Ableitung
%
%   INPUT:
%         desig     -> Belastungsinkrement
%         ZVAR      -> Zustandsvariablen 
%         para      -> Parameter des Materialmodells
%        epara      -> Parameter des Pseudo Modells 
%
%
%    OUTPUT:
%        ZVARneu    -> neue zustandsvariablen nach Lastinkrement
%
%__________________________________________________________________________
%
% Zustandsvariablen:
% Zuerst Zustandsvariablen des Struckturmodells (wie bei Spannungssteuerung)
% dann zusätzliche Variablen des Materialmodells 
% ZVAR = [eeps    -> Pseudo Dehnungen (! nicht pseudo elastische Dehnungen)
%         epsp    -> "reale" Plastische Dehnungen 
%         ealphai -> pseudo Backstress
%         p       -> "reale" plastische Bogenlänge
%         alphai  -> "reale" Backstresstensoren
%__________________________________________________________________________    
%    NOTATIONEN:
%    Spannungen                 Verzerrungen
%
%           sig_11                     eps_11
%    sig =  sig_22               eps = eps_22
%           sig_12                     2eps_12
%
%__________________________________________________________________________
% Parameter:
% zuerst Material- dann Struckturmodell 
%     para = [E, nu,q,b,  c_i,  r_i,  r0,...
%     epara = [E, nu, eq,eb, ec_i, er_i, er0]
%       M = (length(para) - 8)/4
%     M = (length(para)-5)/2;
%     eM = (length(epara)-5)/2;
%
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
% |  Stand: November 2020                                                     |
%  ----------------------------------------------------------------------

% -------------------------------------------------------------------------
% !!!!! AKTUELL NUR ESZ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Für weitere Spannungszustände müssten "nur" die Modellgleichung
% implementiert werden
% -------------------------------------------------------------------------
ntens = 3;
ndi = 2;


%--------------------------------------------------------------------------
%                  Materialparamter ermitteln                                  
%--------------------------------------------------------------------------

% Anzahl Backstress
% M = (length(para)-8)/4;
M = (length(para)-5)/2;
eM = (length(epara)-5)/2;
% Elastizitätskonstanten
E = para(1);
nu = para(2);
er0 = epara(2*eM+5);
eq = epara(3);
eb = epara(4); if eb==0, eb=1e-40; end 
eQinf = eq/eb;
%--------------------------------------------------------------------------
%                        Hilfsmatrizen                                    
%--------------------------------------------------------------------------

% Steifigkeit
C = elast_steifigkeit(E,nu,ntens,ndi);
D = elast_nachgiebigkeit(E,nu,ntens,ndi);
% statische Matrizen
[P, P_line] = set_maps(ntens,ndi);

%--------------------------------------------------------------------------
%                  Zustandsvariablen auslesen                                 
%--------------------------------------------------------------------------

% pseudo gesamtdehnung 
eeps =  ZVAR(1:ntens);
% plastische dehnung
epsp = ZVAR(ntens+1:2*ntens);
% Pseudo (elastische9 Spannungen
esig = C * (eeps - epsp);
% pseudo backstress
ealphai = reshape(ZVAR(2*ntens+1:(eM+2)*ntens),ntens,eM);
% Backstress
alphai = reshape(ZVAR((eM+2)*ntens+2:(eM+M+2)*ntens+1),ntens,M);
% plastische bogenlänge
p = ZVAR((eM+2)*ntens+1);
% Größe der Fließfläche des Strukturmodells
eY = er0 + eQinf * ( 1 - exp(-eb*p));

%--------------------------------------------------------------------------
%                       Elastischer Trial Step
%--------------------------------------------------------------------------

es = P * esig;                                                             % Pseudo Spannungsdeviator
ea = sum(ealphai,2);                                                       % Pseudo Gesamtbackstress
des = P * desig;                                                           % Inkrement Pseudo Spannungsdeviator
es_tr = es + des;                                                          % Pseudo Versuchsspannungsdeviator
ebeta = es_tr - ea;                                                        % Pseudo effektive Spannung
F_tr = ebeta' * P_line * ebeta - 2/3 * eY^2;                          % Pseudo Überspannung
FTOL = 1e-7;


%--------------------------------------------------------------------------
%                   Trial Step wird angenommen
%--------------------------------------------------------------------------
if F_tr < FTOL   % Trial step völlig elastische
    
    % Updaten der pseudo gesamtdehnung
    eeps = eeps + D * desig;      
    % Updaten Zustandsvariablen
    ZVARneu = [eeps;ZVAR(ntens+1:end)];
    
%--------------------------------------------------------------------------
%                   Trial Step wird abgelehnt
%-------------------------------------------------------------------------- 
else
    
    % statische Matrizen
    [P, P_line, P_hat, A] = set_maps(ntens,ndi);
    
    % Integration des Inkrements
    
    % Matlab Version
    ZVARneu = chabocheESZlang(esig, eeps, epsp, ealphai,alphai,p,...
        esig+desig,...
        para,...
        epara,...
        D,P,P_line,P_hat,A);
    
    % Mex Version
    %         ZVARneu = CoderVersion_ESZ_KarimOhnoLang_mex(esig, eeps, epsp, ealphai,alphai,p,...
    %                            esig+desig,...
    %                            para,...
    %                            D,P,P_line,P_hat,A);
    


end % Ende Unterscheidung annehmen/ablehnen trial step
end % Ende Hauptfunktion











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen ebener Spannungszustand             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZVARneu = chabocheESZlang(esig_n, eeps_n, eepsp_n, ealpha_n,...
                                  alpha_n,p_n,...
                                  esig,...
                                  para,epara,...
                                  D,P,P_line,P_hat,A)
% Konkretes Materialmodell für ESZ Spannungszustände, gibt bei vorgabe 
% eines Lastinkrementes die geupdatetn Zustandsvariablen mit implizitem 
% Euler zurück .
%
%   INPUT:
%       esig_n       -> Pseudo Spannungen bei tn
%  eeps_n,eepsp_n    -> Pseudo Dehnungen bei tn
%       ealpha_n     -> Pseudo Backstresstensoren bei tn
%        alpha_n     -> Backstress bei tn
%       esig        -> Spannung bei tn+1
%         parameter -> Modellparameter
%         C,D       -> Elastische Steifigkeit und Nachgiebigkeit   
%   P, P_hat ...    -> Diverse Abbildungen
%
%
%   OUTPUT:
%         ZVARneu -> Zustandsvariablen am Ende des Inkrements
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
%  ----------------------------------------------------------------------

% -------------------------------------------------------------------------
%                        Materialparameter
% -------------------------------------------------------------------------
% Zuweisen der Materialparameter
ntens = 3;                                                                 % Tensorkomponenten
% M = (length(para)-8)/4;
M = (length(para)-5)/2;
eM = (length(epara)-5)/2;
% elastische 
% E = parameter(1);
% nu = parameter(2);
% kinematische Verfestigung Materialmodell
c_i = para(5:4+M);
r_i = para(5+M:4+2*M);
h_i = c_i .* r_i;
% Isotrope Verfestigung Materialmodell
r0 = para(5+2*M);
q = para(3);
b = para(4); if b==0, b=1e-40; end 
Qinf = q/b;
% kinematische Verfestigung Struckturmodell
ec_i = epara(5:4+eM);
er_i = epara(5+eM:4+2*eM);
eh_i = ec_i .* er_i;
% Isotrope Verfestigung Materialmodell
er0 = epara(5+2*eM);
eq = epara(3);
eb = epara(4); if eb==0, eb=1e-40; end 
eQinf = eq/eb;
% oft verwendete Konstanden
w3d2 = sqrt(3/2);
% Radien der Fließflächen
eY_n = er0 + eQinf*(1-exp(-eb*p_n));
Y_n = r0 + Qinf*(1-exp(-b*p_n));

% -------------------------------------------------------------------------
%                        Iterationsschleife
% ------------------------------------------------------------------------- 

% Speicher und Abbruchbedingungen
maxiter = 10000;                % Maximale Iterationsschleifen
dpiter = zeros(1,maxiter+1);  % Speicher plastische dehnungsinkremente
iter = 1;                     % Schleifenzähler 
tol = 1e-8;                   % toleranz Abbruchbedingung

% Startwerte
etheta_i = ones(1,eM);         % Faktoren radial return Strukturmodell
egamma = 1;                   % Faktor radial return Strukturmodell

% deviatorischer Anteilspannung
es = P * esig;

% Zielfunktion
phi = zielfunktion(es,ealpha_n,eY_n,...
                  dpiter(iter),etheta_i,egamma,...
                  eb,er0,eQinf,eh_i,...
                  P_line);

% Iterationsschleife
while abs(phi) > tol
    
    % inkremetiere iterationszähler
    iter = iter + 1;
    
    % verhindere endlosschleife
    if iter > maxiter + 1
        msg = ['Keine Konvergenz in Fixpunktiteration nach ', num2str(iter), ...
               ' Iterationen'];
        error(msg)
    end
        
    % Ableitung phi
    dphidp = centraldiff(es,ealpha_n,eY_n,...
                         dpiter(iter),etheta_i,egamma,...
                         eb,er0,eQinf,eh_i,...
                         P_line);
    
    % Neues Inkrement
    dpiter(iter) = dpiter(iter-1) - phi/dphidp;
    etheta_i = 1./(1+ec_i.*dpiter(iter));
    egamma = 1/(1+eb*dpiter(iter));
    
    % Neue Zielfunktion
    phi = zielfunktion(es,ealpha_n,eY_n,...
                  dpiter(iter),etheta_i,egamma,...
                  eb,er0,eQinf,eh_i,...
                  P_line);  
    
    
end % Ende Iterationsschleife

% disp(iter)
% -------------------------------------------------------------------------
%                        update Zustandsvariablen
% ------------------------------------------------------------------------- 
% neuer Radius 
eY_np1 = er0 + eQinf*(1-exp(-eb*(p_n+dpiter(iter))));
% normale
var1 = sum(etheta_i .* eh_i);
var2 = es - sum(etheta_i .* ealpha_n,2);
ebeta_np1 = eY_np1 * var2 / (eY_np1 + var1 * dpiter(iter));
n_np1 = w3d2 * P_hat * ebeta_np1./eY_np1;

% Inkrement plastische Dehnung
depsp = w3d2 * dpiter(iter) * n_np1;

% Backstresstensoren radial return Strukturmodell
ealpha_star = ealpha_n + 2/3 * eh_i .* (A * depsp);
ealpha_np1 = etheta_i .* ealpha_star;

% Backstress radial return Materialmodell
theta_i = 1./(1+c_i.*dpiter(iter));
alpha_star = alpha_n + 2/3 * h_i .* (A * depsp);
alpha_np1 = theta_i .* alpha_star;

% Plastische Dehnungen
eepsp_np1 = eepsp_n + depsp;

% plastische Bogenlänge
p_np1 = p_n + dpiter(iter);



% Inkrement der pseudo Dehnungen
eeps_np1 = eeps_n + D * (esig- esig_n) + depsp;

% Zusammenfassen der Zustandsvarablen
ZVARneu = [eeps_np1; ...
           eepsp_np1; ...
           reshape(ealpha_np1,ntens*eM,1); ...
           p_np1; ...
           reshape(alpha_np1,ntens*M,1)];   

end % Ende ESZ Materialmodell






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  zielfunktion                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi = zielfunktion(es,ealpha_n,eY_n,...
                            dpiter,etheta_i,egamma,...
                            eb,er0,eQinf,eh_i,...
                            P_line)

    % effektive Spannung
    ebeta_np1 = es - sum(etheta_i.*ealpha_n,2);

    % Fließfläche
    F_np1 = sqrt(3/2*ebeta_np1'*P_line*ebeta_np1)-egamma*eY_n;

    % Hilfsvariable
    hvar = (egamma*eb*(er0+eQinf) + sum(etheta_i.*eh_i));

    % Zielfunktion
    phi = dpiter - F_np1/hvar;

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Numerische Ableitung                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dphidp = centraldiff(es,ealpha_n,eY_n,...
                         dpiter,etheta_i,egamma,...
                         eb,er0,eQinf,eh_i,...
                         P_line)
    % Schrittweite
    epsilon = 1e-8;
    % vowärts
    phi1 = zielfunktion(es,ealpha_n,eY_n,...
                            dpiter+epsilon,etheta_i,egamma,...
                            eb,er0,eQinf,eh_i,...
                            P_line);
    % rückwärts
    phi2 = zielfunktion(es,ealpha_n,eY_n,...
                            dpiter-epsilon,etheta_i,egamma,...
                            eb,er0,eQinf,eh_i,...
                            P_line);
    % Ableitung
    dphidp = (phi1-phi2)/(2*epsilon);
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Verfahrensfunktion Spannungssteuerung                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dp = verfahren(s_tr,alpha_n,Y_np1,h_i,theta_i,P_line)
% Verfahrensfunktion für Spannungssteurung
% INPUT:
% s_tr     -> deviatorischer Anteil versuchsspannung
% alpha_n  -> aktueller Zustand Backstresstensoren
% Y_np1    -> neuer Zustand Fließfläche
% h_i      -> Materialparameter h_i = c_i * r_i
% theta_i  -> interne Variable
% P_line   -> Abbildung für skalarprodukt zwischen Spannungen (Variert je
%             nach spannungszustand
%
% OUTPUT:
% dp       -> neues Inkrement plastische Bogenlänge
%__________________________________________________________________________
    % trial effektive spannung
    a_n = sum(theta_i .* alpha_n,2);
    beta_tr = s_tr - a_n;
    % Versuchsfließspannunge
    F_tr = sqrt(3/2 * beta_tr' * P_line * beta_tr) - Y_np1;
    % Hilfsvariable 
    var1 = sum(h_i .* theta_i);
    % neues plastisches Dehnungsinkrement
    dp = F_tr/(var1);
end % Ende Verfahrensfunktion dehnungssteuerung








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Aitkens Delta Quadrat Verfahren                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dphat = deltaquadrat(dpiter,iter)
% Funktion zum beschleunigen der Konvergenz
    dphat = dpiter(iter) - ( dpiter(iter) - dpiter(iter-1) ).^2 / ...
        ( dpiter(iter) - 2*dpiter(iter-1) + dpiter(iter-2) );
    if dphat <= 0
        dphat = dpiter(iter);
    end
end






