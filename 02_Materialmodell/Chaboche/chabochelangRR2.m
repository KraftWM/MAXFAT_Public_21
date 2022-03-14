function [ZVARneu] = chabochelangRR2(desig,ZVAR,para,epara)
% Implementierung des Materialmodells nach chaboche, für pseudo
% stress aproach nach Lang 
%
% ! Nur Spannungssteuerung im ESZ
%
%
% QUELLE SPANNUNGSGESTEUERTER IMPLIZITER EULER:
% Eigener Algo.
%
%
%
%   INPUT:
%         desig     -> Belastungsinkrement
%         ZVAR      -> Zustandsvariablen 
%         para      -> Parameter des Pseudo Modells und des Materialmodells
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
%         eY      -> Radius FF Strukturmodell 
%         p       -> "reale" plastische Bogenlänge
%         alphai  -> "reale" Backstresstensoren
%         Y       -> "realer" Radius FF
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
%     para = [E, nu,b,Qinf,  c_i,  r_i,  r0,...
%     epara = [E, nu, eb,eQinf, ec_i, er_i, er0]
%     M = (length(para)-5)/2;
%     eM = (length(epara)-5)/2;
%
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
% |  Stand: Oktober 2021                                                     |
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
% er0 = epara(2*M+5);
% eb = epara(3);
% eQinf = epara(4);

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
alphai = reshape(ZVAR((eM+2)*ntens+3:(eM+M+2)*ntens+2),ntens,M);
% plastische bogenlänge
p = ZVAR((eM+2)*ntens+2);
% Größe der Fließfläche des Strukturmodells
eY = ZVAR((eM+2)*ntens+1);
% Größe der Fließfläche des Materialmodells
Y = ZVAR((eM+M+2)*ntens+3);

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
    ZVARneu = chabocheESZlang(esig, eeps, epsp,...
        ealphai,alphai,eY,Y,p,...
        esig+desig,...
        para,epara,...
        D,P,P_line,P_hat,A);   


end % Ende Unterscheidung annehmen/ablehnen trial step
end % Ende Hauptfunktion











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen ebener Spannungszustand             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZVARneu = chabocheESZlang(...
                                  esig_n,eeps_n,epsp_n,ealpha_n,alpha_n,...
                                    eY_n,Y_n,p_n,...                        % Zustandsvariablen
                                    esig_tr,...                                          % Versuchsspannung
                                    para,epara,...                                     % Parameter
                                    D,P,P_line,P_hat,A)
% -------------------------------------------------------------------------
%                        Materialparameter
% -------------------------------------------------------------------------
% Zuweisen der Materialparameter
ntens = 3;                                                                 % Tensorkomponenten
% M = (length(para)-8)/4;                                                  % Anzahl Backstresstensoren
M = (length(para)-5)/2;
eM = (length(epara)-5)/2;
% kinematische Verfestigung Materialmodell
c_i = para(5:4+M);
r_i = para(5+M:4+2*M);
h_i = c_i .* r_i;
% Isotrope Verfestigung Materialmodell
r0 = para(5+2*M);
b = para(3);
Qinf = para(4);
% kinematische Verfestigung Struckturmodell
ec_i = epara(5:4+eM);
er_i = epara(5+eM:4+2*eM);
eh_i = ec_i .* er_i;
% Isotrope Verfestigung Materialmodell
er0 = epara(5+2*eM);
eb = epara(3);
eQinf = epara(4); 

% oft verwendete Konstanden
w3d2 = sqrt(3/2);
w2d3 = sqrt(2/3);
% oft verwendete Konstanden
mynull = 1e-40;


% -------------------------------------------------------------------------
%                        Definitionen für Iteration
% ------------------------------------------------------------------------- 

% Speicher und Abbruchbedingungen
maxiter = 1000;                   % Maximale Iterationsschleifen
dpiter = zeros(1,maxiter+2);     % Speicher plastische dehnungsinkremente
normpiter = zeros(1,maxiter+2);  % Norm Abbruchbedingung (für Konvergnz Analyse)
iter = 1;                        % Schleifenzähler 
tol = 1e-4;                      % toleranz Abbruchbedingung
normp = 1;                       % Abbruchbedingung

% -------------------------------------------------------------------------
%                   Initialisiere startwerte Iteration
% ------------------------------------------------------------------------- 

% ... Korrekturfaktoren
% dp_npe = 0;
etheta_i = ones(1,eM);
eGamma = 1;
% ... Isotrope
eY_npe = eY_n;
% ... trial Spannungsdeviator
es_tr = P * esig_tr;

% -------------------------------------------------------------------------
%                   erster Iterationsschritt
% ------------------------------------------------------------------------- 
% ... Inkrement plastische Bogenlänge
ealpha_hat = sum(etheta_i.*ealpha_n,2);
ebeta_hat = es_tr - ealpha_hat;

% ---------------------
% Fixpunktiteration 
% ---------------------
dp_npe = (sqrt(1.5*(ebeta_hat'*P_line*ebeta_hat)) - eGamma*eY_n) ...
         /(eGamma*eb*(eQinf+er0)+sum(eh_i.*etheta_i));    
dpiter(iter) = dp_npe;
% ... Backstress
etheta_i = 1./(1+ec_i*dp_npe);
% ... Radius Fließfläche 
eGamma = 1/(1+eb*dp_npe);

% -------------------------------------------------------------------------
%                  Fixpunktiteration
% ------------------------------------------------------------------------- 
% ... Schleife
while normp > tol
    % ... Inkrementieren Schleifenzähler
    iter = iter + 1;
    % ... Inkrement plastische Bogenlänge
    ealpha_hat = sum(etheta_i.*ealpha_n,2);
    ebeta_hat = es_tr - ealpha_hat;
    
    % ---------------------
    % Fixpunktiteration
    % ---------------------
    dp_npe = (sqrt(1.5*(ebeta_hat'*P_line*ebeta_hat)) - eGamma*Y_n) ...
        /(eGamma*eb*(eQinf+er0)+sum(eh_i.*etheta_i));
    % Aitkens Delta^2 Verfahren
    if mod(iter,3) == 0
        dp_hat = dp_npe - ( dp_npe - dpiter(iter-1) )^2 / ...
        ( dp_npe - 2*dpiter(iter-1) + dpiter(iter-2) );
        if dp_hat > 0
            dp_npe = dp_hat;
        end
    end
    dpiter(iter) = dp_npe;
    % ... Backstress
    etheta_i = 1./(1+ec_i*dp_npe);
    % ... Radius Fließfläche
    eGamma = 1/(1+eb*dp_npe);
    eY_npe = eGamma * (eY_n + eb*(eQinf+er0)*dp_npe );
    % ... effektive Spannung & Fließflächen Normale
    beta_npe = eY_npe * ebeta_hat / ( eY_npe + dp_npe * sum(eh_i.*etheta_i));
%     n_npe = w3d2 * (P_hat * beta_npe)/ r0;
    n_npe = sqrt(beta_npe' * P_line * beta_npe); if n_npe == 0; n_npe = mynull; end
    n_npe = (P_hat * beta_npe)/n_npe;
%     fprintf('||n|| = %.15f\n',n_npe' * A * P_line * A * n_npe)
    % ... Prüfe Konvergenz
    normp = abs(1-dpiter(iter)/dpiter(iter-1));
    normpiter(iter) = normp;
    % ... Verhindere Endlosschleifen
    if iter > maxiter + 1
        msg = ['Keine Konvergenz in Fixpunktiteration im Chaboche Modell ',...
               'nach ', num2str(iter),' Iterationen. Aktuelle Fehlernorm: '...
               num2str(normp)];
        warning(msg)
        break;
    end
        
end % Ende Äußere Iterationsschleife
% fprintf('iter = %i normp = %.5f\n',iter,normp)
% figure(1)
% semilogy(normpiter(1:iter),'d-')
% figure(3)
% plot(dpiter(1:iter)), pause(0.5)
% -------------------------------------------------------------------------
%                        update Zustandsvariablen
% -------------------------------------------------------------------------
% ... effektive Spannung & Fließflächen Normale (Falls konvergenz ohne
% iteration)
beta_npe = eY_npe * ebeta_hat / ( eY_npe + dp_npe * sum(eh_i.*etheta_i));
%     n_npe = w3d2 * (P_hat * beta_npe)/ r0;
n_npe = sqrt(beta_npe' * P_line * beta_npe); if n_npe == 0; n_npe = mynull; end
n_npe = (P_hat * beta_npe)/n_npe;
% ... Backstress Material
theta_i = 1./(1+c_i*dp_npe);
alpha_npe = theta_i .* (alpha_n + w2d3 * dp_npe * (A*n_npe) * h_i);
% ... Radius FF Material
Gamma = 1/(1+b*dp_npe);
Y_npe = Gamma * (Y_n + b*(Qinf+r0)*dp_npe);
% ... Backstress Strukturmodell
ealpha_npe = etheta_i .* (ealpha_n + w2d3*(A*n_npe)*(ec_i.*er_i)*dp_npe);
% ... plastische Dehnungen
depsp = w3d2*dp_npe * n_npe;
epsp_npe = epsp_n + depsp;
% ... plastische Bogenlänge
p_npe = p_n + dp_npe;
% ... Dehnungen
eeps_npe = eeps_n + D * (esig_tr-esig_n) + depsp;
ZVARneu = [eeps_npe; epsp_npe ; reshape(ealpha_npe,ntens*eM,1);eY_npe; ...
           p_npe; reshape(alpha_npe,ntens*M,1);Y_npe];   

end % Ende Modellgleichung in ESZ