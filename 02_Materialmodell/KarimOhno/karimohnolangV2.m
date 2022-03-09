function [ZVARneu] = karimohnolangV2(desig,ZVAR,para)
% Implementierung des Materialmodells nach Karim und Ohno, für pseudo
% stress aproach nach Lang 
%
% ! Nur Spannungssteuerung im ESZ
% 
% Unterschied zu 'karimohnolang.m':
% hier wird auch isotrope Verfestigung mit berücksichtigt
%
% QUELLE MODELL:
% Aus Ohno,Karim et al. 2000 Uniaxial Ratchetting of 316FR steel at room
% tempreture part II
%
% Karim,Ohno 2000 Kinematic hardeningmodel suitable for ratchetting with
% steady state
%
% ABAQUS MANUEL (isotrope Verfestigung wie bei Chaboche Modell in Abaqus)
%
% QUELLE DEHNUNGSGESTEURTER IMPLIZITER EULER
% Kobayashi, Ohno 2002 - Implementation of cyclic plasticity models based
% on a general form of kinematic hardening
%
% QUELLE SPANNUNGSGESTEUERTER IMPLIZITER EULER:
% Eigener Algo.
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
%         p       -> "reale" plastische Bogenlänge
%         alphai  -> "reale" Backstresstensoren
%
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
%     para = [E, nu,  c_i,  r_i,  mu_i,  r0,  Qinf,  b...
%                    ec_i, er_i, emu_i, er0, eQinf, eb]
%       M = (length(para) - 8)/6
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
M = (length(para)-8)/6;
% Elastizitätskonstanten
E = para(1);
nu = para(2);
er0 = para(6*M+6);
eQinf = para(6*M+7);
eb = para(6*M+8);
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
ealphai = reshape(ZVAR(2*ntens+1:(M+2)*ntens),ntens,M);
% Backstress
alphai = reshape(ZVAR((M+2)*ntens+2:(2*M+2)*ntens+1),ntens,M);
% plastische bogenlänge
p = ZVAR((M+2)*ntens+1);
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
    ZVARneu = karimohnoESZlang(esig, eeps, epsp, ealphai,alphai,p,...
        esig+desig,...
        para,...
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
function ZVARneu = karimohnoESZlang(esig_n, eeps_n, eepsp_n, ealpha_n,...
                                  alpha_n,p_n,...
                                  esig,...
                                  parameter,...
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
M = (length(parameter)-8)/6;                                               % Anzahl TBST

% elastische 
E = parameter(1);
nu = parameter(2);
G = E/(2*(1+nu));
% kinematische Verfestigung Materialmodell
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
mu_i = parameter(3+2*M:2+3*M);
h_i = c_i .* r_i;
% Isotrope Verfestigung Materialmodell
r0 = parameter(3+3*M);
Qinf = parameter(4+3*M);
b = parameter(5+3*M);
% kinematische Verfestigung Struckturmodell
ec_i = parameter(6+3*M:5+4*M);
er_i = parameter(6+4*M:5+5*M);
emu_i = parameter(6+5*M:5+6*M);
eh_i = ec_i .* er_i;
% Isotrope Verfestigung Materialmodell
er0 = parameter(6+6*M);
eQinf = parameter(7+6*M);
eb = parameter(8+6*M);
% oft verwendete Konstanden
w3d2 = sqrt(3/2);
% Radien der Fließflächen
eY_n = er0 + eQinf*(1-exp(-eb*p_n));
Y_n = r0 + Qinf*(1-exp(-b*p_n));

% -------------------------------------------------------------------------
%                        Iterationsschleife
% ------------------------------------------------------------------------- 

% Speicher und Abbruchbedingungen
maxiter = 100;                % Maximale Iterationsschleifen
dpiter = zeros(1,maxiter+1);  % Speicher plastische dehnungsinkremente
iter = 1;                     % Schleifenzähler 
tol = 1e-4;                   % toleranz Abbruchbedingung

% Startwerte
normp = 1;                    % Abbruchbedingung
etheta_i = ones(1,M);         % Faktoren radial return Strukturmodell
etheta_i_tilde = ones(1,M);

% deviatorischer Anteilspannung
es = P * esig;

% neuer Radius der Fließfläche des Strukturmodells (für trial step)
eY_np1 = eY_n;

% erster schritt
dpiter(iter) = verfahren(es,ealpha_n,eY_np1,eh_i,etheta_i,P_line);

% bestimme neue effektive Spannung
var1 = sum(etheta_i .* eh_i);               % Hilfsvariable
var2 = es - sum(etheta_i .* ealpha_n,2);    % Hilfsvariable
ebeta_np1 = eY_np1 * var2 / (eY_np1 + var1 * dpiter(iter));

% Update Normale
n_np1 = w3d2 * P_hat * ebeta_np1./eY_np1;

% Inkrement plastische Dehnung
depsp = w3d2 * dpiter(iter) * n_np1;

% Prädiktor pseudo Backstress
ealpha_star = ealpha_n + 2/3 * eh_i .* (A * depsp);
ej_i = er_i./(er_i + emu_i.*eh_i*dpiter(iter));
ealpha_hash = ej_i .* ealpha_star;

% Iterationsschleife
while normp > tol
    
    % effektive Backstresstensoren
    ealpha_eff = sqrt( 3/2 * sum( (P_line*ealpha_hash) .* ealpha_hash) );
    idx = ealpha_eff <= er_i;
    etheta_i_tilde( idx ) = 1;
    etheta_i_tilde( ~idx ) = er_i(~idx)./ealpha_eff(~idx);
    etheta_i = ej_i .* etheta_i_tilde;
    
    % update radius der FF des Strukturmodells
    eY_np1 = er0 + eQinf * ( 1 - exp( -eb* (p_n + dpiter(iter)) ) );
    
    % Inkrementieren Schleifenzähler
    iter = iter + 1;
    % neues plastisches Dehnungsinkrement
    dpiter(iter) = verfahren(es,ealpha_n,eY_np1,eh_i,etheta_i,P_line);
    % Aitkens Delta^2 Verfahren
    if mod(iter,3) == 0
        dpiter(iter) = deltaquadrat(dpiter,iter);
    end
    % Update effektive spannung
    var1 = sum(etheta_i .* eh_i);
    var2 = es - sum(etheta_i .* ealpha_n,2);
    ebeta_np1 = eY_np1 * var2 / (eY_np1 + var1 * dpiter(iter));
    
    % Update Normale
    n_np1 = w3d2 * P_hat * ebeta_np1./eY_np1;
    
    % Inkrement plastische Dehnung
    depsp = w3d2 * dpiter(iter) * n_np1;
    
    % Prädiktor Backstress
    ealpha_star = ealpha_n + 2/3 * eh_i .* (A * depsp);
    ej_i = er_i./(er_i + emu_i.*eh_i*dpiter(iter));
    ealpha_hash = ej_i .* ealpha_star;
    
    % Check Konvergenz
    normp = abs(1-dpiter(iter)/dpiter(iter-1));
    
    % verhindere endlosschleife
    if iter > maxiter + 1
        msg = ['Keine Konvergenz in Fixpunktiteration nach ', num2str(iter), ...
               ' Iterationen'];
        error(msg)
    end
    
end % Ende Iterationsschleife

% Eventuelle Ausgabe um Konvergenz zu testen
% fid = fopen('Iter.dat','a');
% fprintf(fid,'%10.4f;%10.4f;%10.4f;%i\n',esig,iter);
% fclose(fid);

% -------------------------------------------------------------------------
%                        update Zustandsvariablen
% ------------------------------------------------------------------------- 
% Zustandsvariablen:
% Zuerst Zustandsvariablen des Struckturmodells (wie bei Spannungssteuerung)
% dann zusätzliche Variablen des Materialmodells 
% ZVAR = [eeps    -> Pseudo Dehnungen (! nicht pseudo elastische Dehnungen)
%         epsp    -> "reale" Plastische Dehnungen 
%         ealphai -> pseudo Backstress
%         p       -> "reale" plastische Bogenlänge
%         alphai  -> "reale" Backstresstensoren
%     para = [E, nu,  c_i,  r_i,  mu_i,  r0,...
%                    ec_i, er_i, emu_i, er0]

% Backstress radial return Materialmodell
alpha_star = alpha_n + 2/3 * h_i .* (A * depsp);
j_i = r_i./(r_i + mu_i.*h_i*dpiter(iter));
alpha_hash = j_i .* alpha_star;
alpha_eff = sqrt( 3/2 * sum( (P_line*alpha_hash) .* alpha_hash) ); 
idx1 = alpha_eff <= r_i;
idx2 = ~idx1;
theta_i_tilde( idx1 ) = 1;
theta_i_tilde( idx2 ) = r_i( idx2 ) ./ alpha_eff( idx2 );
theta_i = j_i .* theta_i_tilde;
alpha_np1 = theta_i .* alpha_star;

% Plastische Dehnungen
eepsp_np1 = eepsp_n + depsp;

% plastische Bogenlänge
p_np1 = p_n + dpiter(iter);

% Backstresstensoren radial return Strukturmodell
ealpha_np1 = etheta_i .* ealpha_star;

% Inkrement der pseudo Dehnungen
eeps_np1 = eeps_n + D * (esig- esig_n) + depsp;

% Zusammenfassen der Zustandsvarablen
ZVARneu = [eeps_np1; ...
           eepsp_np1; ...
           reshape(ealpha_np1,ntens*M,1); ...
           p_np1; ...
           reshape(alpha_np1,ntens*M,1)];   

end % Ende ESZ Materialmodell





















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






