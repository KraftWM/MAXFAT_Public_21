function [X_neu] = impohnowangI(ntens, ndi, ink, X, ink_flag, parameter)
% Materialmodell des Ohno Wang modells I:
%
% QUELLE MODELL:
% Aus Ohno et al. 1993 KINEMATIC HARDENING RULES WITH CRITICAL
% STATE OF DYNAMIC RECOVERY, PART I
%
% QUELLE DEHNUNGSGESTEURTER IMPLIZITER EULER
% Kobayashi, Ohno 2002 - Implementation of cyclic plasticity models based
% on a general form of kinematic hardening
%
% QUELLE SPANNUNGSGESTEUERTER IMPLIZITER EULER:
% Eigener Algo. Implementation zu Testzwecken
%
%   INPUT:
%         ntens     -> Anzahl Tensorkomponten (in Spannungen)
%         ndi       -> Anzahl Diagonalelemente, mehr kkommentare
%         ink       -> Belastungsinkrement
%         X         -> Zustandsvariablen [eps;epsp;alphai;p] bein spansteu
%                                        [sig;epsp;alphai;p] bein dehnsteu
%         ink_flag  -> (0) Spannungsgesteuerter Prozess 
%                      (1) Dehnungsgesteuerter Prozess
%         parameter -> Materialparameter des Modells
%                      [E,nu,c_i,r_i,r0]
%                      % chi_i = Inf im vergleich zu Ohno Wang Model II
%
%
%    OUTPUT:
%        X_neu -> neue zustandsvariablen nach Lastinkrement
%    
%    NOTATIONEN:
%    Spannungen                 Verzerrungen
%
%           sig_11                     eps_11
%    sig =  sig_22               eps = eps_22
%           sig_12                     2eps_12
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
% |  Stand: April 2019                                                  |
%  ----------------------------------------------------------------------

%--------------------------------------------------------------------------
%                  Materialparamter ermitteln                                  
%--------------------------------------------------------------------------

% Anzahl Backstress
M = (length(parameter)-3)/2;
% Elastizitätskonstanten
E = parameter(1);
nu = parameter(2);
r0 = parameter(end);
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

% Spannungen und Dehnungen
if ink_flag == 0 % spansteu
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = C * (eps - epsp);
    sig_tr = sig + ink;
    dsig_tr = ink;    % elastischer trial step
elseif ink_flag == 1 % dehnsteu
    sig = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    eps = D * sig + epsp;
    dsig_tr = C * ink; % elastischer trial step
    sig_tr = sig + dsig_tr;
else
    msg = '      keine oder flasche angabe der Laststeuerung';
    error(msg)
end
% Backstresstensoren
alpha = reshape( X(2*ntens+1:(M+2)*ntens) , ntens, M);
% plastische Bogenlänge
p = X(end);

%--------------------------------------------------------------------------
%                       Elastischer Trial Step
%--------------------------------------------------------------------------

% Trial Inkrement der Spannungen 
s = P * sig;                % Span dev
ds = P * dsig_tr;           % dev des Spannungsinkrements
s_tr = s + ds;              % Trial Spannungsdeviator
a = sum(alpha,2);           % Backstress;
beta = s_tr - a;                    % Relativ Spannung
if ntens == 1
    F_tr = abs(beta) - r0;
else
    F_tr = beta' * P_line * beta - 2/3 * r0^2; % Trial Fließfunktion
end
FTOL = 1e-10;                        % Toleranz für Abweichungen F ~= 0


%--------------------------------------------------------------------------
%                   Trial Step wird angenommen
%--------------------------------------------------------------------------
if F_tr < FTOL   % Trial step völlig elastische
    
    if ink_flag == 0 % spansteu;
        eps = eps + D * ink;
        X_neu = [eps;X(ntens+1:end)];
    elseif ink_flag == 1 % dehnsteu
        sig = sig + dsig_tr;
        X_neu = [sig;X(ntens+1:end)];
    end
    
%--------------------------------------------------------------------------
%                   Trial Step wird abgelehnt
%-------------------------------------------------------------------------- 
else
    
    % Integration je nach Spannungszustand
    if ntens == 6 % 3D

        % abbildungen
        [P, P_hat, P_tilde] = set_maps(ntens,ndi);
        % Funktion
        X_neu = impohno3D( sig, eps, epsp, alpha,p,...                     % Zustand
                           sig_tr,ink,ink_flag,...                        % Inkrement
                           parameter,...                                   % parameter
                           C,D,P,P_hat,P_tilde);                           % Abbildungen

    elseif ntens == 3 && ndi == 2 % ESZ
        
        % statische Matrizen
        [P, P_line, P_hat, A, P_check] = set_maps(ntens,ndi);
        % Funktion
        X_neu = impohnoESZ(sig, eps, epsp, alpha,p,...
                           sig_tr,ink,ink_flag,...
                           parameter,...
                           C,D,P,P_line,P_hat,A,P_check);
    elseif ntens == 1 % 1D

        msg = 'net implementiert';
        error(msg)
        
    end

end % Ende Unterscheidung annehmen/ablehnen trial step
end % Ende Hauptfunktion



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Verfahrensfunktion Dehnungssteuerung                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dp = verfahren_eps(s_tr,alpha_n,Y_np1,G,h_i,theta_i,P_line)
% Verfahrensfunktion für Dehnungssteurung
% INPUT:
% s_tr     -> deviatorischer Anteil versuchsspannung
% alpha_n  -> aktueller Zustand Backstresstensoren
% Y_np1    -> neuer Zustand Fließfläche
% G        -> Schubmodul
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
    dp = F_tr/(3*G+var1);
end % Ende Verfahrensfunktion dehnungssteuerung

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Verfahrensfunktion Spannungssteuerung                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dp = verfahren_sig(s_tr,alpha_n,Y_np1,G,h_i,theta_i,P_line)
% Verfahrensfunktion für Spannungssteurung
% INPUT:
% s_tr     -> deviatorischer Anteil versuchsspannung
% alpha_n  -> aktueller Zustand Backstresstensoren
% Y_np1    -> neuer Zustand Fließfläche
% G        -> Schubmodul
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
%                  Modellgleichungen 3D Spannungszustand                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X_neu = impohno3D(sig_n, eps_n, epsp_n, alpha_n,p_n,...
                           sig_tr,ink,ink_flag,...
                           parameter,...
                           C,D,P,P_hat,P_tilde)
% Konkretes Materialmodell für 3D Spannungszustände, gibt bei vorgabe eines
% Lastinkrementes die geupdatetn Zustandsvariablen mit implizitem Euler 
% zurück .
%
% Materialmodell:
% kinematische Verfestigung -> durch Summe von Teilbackstresstensoren
% jeder TBST verhält sich nach Ohno Wang Model I (Heaviside)
%   INPUT:
%       sig_n       -> Spannungen
%  eps_n,epsp_n     -> Dehnungen
%       alpha_n     -> Backstresstensoren
%       dsig_tr     -> Versuchsspannung ( richtige spannung bei
%                      spannungssteuerung )
%         ink       -> Lastinkrement
%         ink_flag  -> Belastungsart
%         parameter -> Modellparameter
%         C,D       -> Elastische Steifigkeit und Nachgiebigkeit   
%   P, P_hat ...    -> Diverse Abbildungen
%
%
%   OUTPUT:
%         dX -> Inkremente der Zustandsvariablen
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
%  ----------------------------------------------------------------------

% -------------------------------------------------------------------------
%                        Materialparameter
% -------------------------------------------------------------------------
% Zuweisen der Materialparameter
ntens = 6;                                                                 % Tensorkomponenten
M = (length(parameter)-3)/2;                                               % Anzahl TBST

% elastische 
E = parameter(1);
nu = parameter(2);
G = E/(2*(1+nu));
r0 = parameter(end);                                                       % startradius fliessfläche
% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
h_i = c_i .* r_i;
% oft verwendete Konstanden
w3d2 = sqrt(3/2);

% -------------------------------------------------------------------------
%                        Iterationsschleife
% ------------------------------------------------------------------------- 

% Speicher und Abbruchbedingungen
maxiter = 100;              % Maximale Iterationsschleifen
dpiter = zeros(1,maxiter);  % Speicher plastische dehnungsinkremente
iter = 1;                   % Schleifenzähler 
tol = 1e-4;                 % toleranz Abbruchbedingung
% Startwerte
normp = 1;
theta_i = ones(1,M);
% Setze verfahrensvorschrift je nach steuerung
if ink_flag == 0
    verfahren = @verfahren_sig;
elseif ink_flag == 1
    verfahren = @verfahren_eps;
end
% deviatorischer Anteil Trial spannung
s_tr = P * sig_tr;
% neuer radius FF = alter radius FF
Y_np1 = r0;
% erster schritt
dpiter(iter) = verfahren(s_tr,alpha_n,Y_np1,G,h_i,theta_i,P_hat);
% Update Hilfvariable
if ink_flag == 0
    var1 = sum(theta_i .* h_i);
elseif ink_flag == 1
    var1 = 3*G + sum(theta_i .* h_i);
end
% Update Hilfsvariable 2
var2 = s_tr - sum(theta_i .* alpha_n,2);
% Update effektive spannung
beta_np1 = Y_np1 * var2 / (Y_np1 + var1 * dpiter(iter));
% Update Normale
n_np1 = w3d2 * P_hat * beta_np1./Y_np1;
% Inkrement plastische Dehnung
depsp = w3d2 * dpiter(iter) * n_np1;
% Prädiktor Backstress
alpha_star = alpha_n + 2/3 * h_i .* (P_tilde * depsp);
% Iterationsschleife
while normp > tol
    % effektive Backstresstensoren
    alpha_eff = sqrt( 3/2 * sum( (P_hat*alpha_star) .* alpha_star) );
    idx1 = alpha_eff <= r_i;
    idx2 = alpha_eff > r_i;
    theta_i( idx1 ) = 1;
    theta_i( idx2 ) = r_i(idx2)./alpha_eff(idx2);
    % Inkrementieren Schleifenzähler
    iter = iter + 1;
    % neues plastisches Dehnungsinkrement
    dpiter(iter) = verfahren(s_tr,alpha_n,Y_np1,G,h_i,theta_i,P_hat);
    % Update Hilfvariable
    if ink_flag == 0
        var1 = sum(theta_i .* h_i);
    elseif ink_flag == 1
        var1 = 3*G + sum(theta_i .* h_i);
    end
    % Update Hilfsvariable 2
    var2 = s_tr - sum(theta_i .* alpha_n,2);
    % Update effektive spannung
    beta_np1 = Y_np1 * var2 / (Y_np1 + var1 * dpiter(iter));
    % Update Normale
    n_np1 = w3d2 * P_hat * beta_np1./Y_np1;
    % Inkrement plastische Dehnung
    depsp = w3d2 * dpiter(iter) * n_np1;
    % Prädiktor Backstress
    alpha_star = alpha_n + 2/3 * h_i .* (P_tilde * depsp);
    % Check Konvergenz
    normp = abs(1-dpiter(iter)/dpiter(iter-1));
    % verhindere endlosschleife
    if iter > maxiter + 1
        msg = ['Keine Konvergenz in Fixpunktiteration nach ', num2str(iter), ...
               ' Iterationen'];
        error(msg)
    end
end % Ende Iterationsschleife

% -------------------------------------------------------------------------
%                        update Zustandsvariablen
% ------------------------------------------------------------------------- 
%         X         -> Zustandsvariablen [eps;epsp;alphai;p] bein spansteu
%                                        [sig;epsp;alphai;p] bein dehnsteu
%         ink_flag  -> (0) Spannungsgesteuerter Prozess 
%                      (1) Dehnungsgesteuerter Prozess
%         parameter -> Materialparameter des Modells
%                      [E,nu,c_i,r_i,r0]
% Plastische Dehnungen
epsp_np1 = epsp_n + depsp;
% plastische Bogenlänge
p_np1 = p_n + dpiter(iter);
% Backstresstensoren
alpha_np1 = theta_i .* alpha_star;
% Inkrement in Spannungen oder Dehnungen
if ink_flag == 0
    eps_np1 = eps_n + D * ink + depsp;
    out = [eps_np1; epsp_np1];
elseif ink_flag == 1
    sig_np1 = sig_tr - C * depsp;
    out = [sig_np1; epsp_np1];
end
% beta_np1 = P * sig_np1 - sum(alpha_np1,2);
% F_np1 = sqrt(3/2 * beta_np1' * P_hat * beta_np1) - r0;
% fprintf('F = %.2d\n',F_np1);
X_neu = [out ; reshape(alpha_np1,ntens*M,1); p_np1];   

end % Ende 3D Materialmodell































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen 3D Spannungszustand                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X_neu = impohnoESZ(sig_n, eps_n, epsp_n, alpha_n,p_n,...
                           sig_tr,ink,ink_flag,...
                           parameter,...
                           C,D,P,P_line,P_hat,A,P_check)
% Konkretes Materialmodell für ESZ Spannungszustände, gibt bei vorgabe 
% eines Lastinkrementes die geupdatetn Zustandsvariablen mit implizitem 
% Euler zurück .
%
% Materialmodell:
% kinematische Verfestigung -> durch Summe von Teilbackstresstensoren
% jeder TBST verhält sich nach Ohno Wang Model I (Heaviside)
%   INPUT:
%       sig_n       -> Spannungen
%  eps_n,epsp_n     -> Dehnungen
%       alpha_n     -> Backstresstensoren
%       dsig_tr     -> Versuchsspannung ( richtige spannung bei
%                      spannungssteuerung )
%         ink       -> Lastinkrement
%         ink_flag  -> Belastungsart
%         parameter -> Modellparameter
%         C,D       -> Elastische Steifigkeit und Nachgiebigkeit   
%   P, P_hat ...    -> Diverse Abbildungen
%
%
%   OUTPUT:
%         dX -> Inkremente der Zustandsvariablen
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
%  ----------------------------------------------------------------------

% -------------------------------------------------------------------------
%                        Materialparameter
% -------------------------------------------------------------------------
% Zuweisen der Materialparameter
ntens = 3;                                                                 % Tensorkomponenten
M = (length(parameter)-3)/2;                                               % Anzahl TBST

% elastische 
E = parameter(1);
nu = parameter(2);
G = E/(2*(1+nu));
r0 = parameter(end);                                                       % startradius fliessfläche
% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
h_i = c_i .* r_i;
% oft verwendete Konstanden
w3d2 = sqrt(3/2);

% -------------------------------------------------------------------------
%                        Iterationsschleife
% ------------------------------------------------------------------------- 

% Speicher und Abbruchbedingungen
maxiter = 100;              % Maximale Iterationsschleifen
dpiter = zeros(1,maxiter);  % Speicher plastische dehnungsinkremente
iter = 1;                   % Schleifenzähler 
tol = 1e-4;                 % toleranz Abbruchbedingung
% Startwerte
normp = 1;
theta_i = ones(1,M);
% Setze verfahrensvorschrift je nach steuerung
if ink_flag == 0
    verfahren = @verfahren_sig;
elseif ink_flag == 1
    verfahren = @verfahren_eps;
end
% deviatorischer Anteil Trial spannung
s_tr = P * sig_tr;
% neuer radius FF = alter radius FF
Y_np1 = r0;
% erster schritt
dpiter(iter) = verfahren(s_tr,alpha_n,Y_np1,G,h_i,theta_i,P_line);
% Update Hilfvariable
if ink_flag == 0
    var1 = sum(theta_i .* h_i);
elseif ink_flag == 1
    var1 = 3*G + sum(theta_i .* h_i);
end
% Update Hilfsvariable 2
var2 = s_tr - sum(theta_i .* alpha_n,2);
% Update effektive spannung
beta_np1 = Y_np1 * var2 / (Y_np1 + var1 * dpiter(iter));
% Update Normale
n_np1 = w3d2 * P_hat * beta_np1./Y_np1;
% Inkrement plastische Dehnung
depsp = w3d2 * dpiter(iter) * n_np1;
% Prädiktor Backstress
alpha_star = alpha_n + 2/3 * h_i .* (A * depsp);
% Iterationsschleife
while normp > tol
    % effektive Backstresstensoren
    alpha_eff = sqrt( 3/2 * sum( (P_line*alpha_star) .* alpha_star) );
    idx1 = alpha_eff <= r_i;
    idx2 = alpha_eff > r_i;
    theta_i( idx1 ) = 1;
    theta_i( idx2 ) = r_i(idx2)./alpha_eff(idx2);
    % Inkrementieren Schleifenzähler
    iter = iter + 1;
    % neues plastisches Dehnungsinkrement
    dpiter(iter) = verfahren(s_tr,alpha_n,Y_np1,G,h_i,theta_i,P_line);
    % Update Hilfvariable
    if ink_flag == 0
        var1 = sum(theta_i .* h_i);
    elseif ink_flag == 1
        var1 = 3*G + sum(theta_i .* h_i);
    end
    % Update Hilfsvariable 2
    var2 = s_tr - sum(theta_i .* alpha_n,2);
    % Update effektive spannung
    beta_np1 = Y_np1 * var2 / (Y_np1 + var1 * dpiter(iter));
    % Update Normale
    n_np1 = w3d2 * P_hat * beta_np1./Y_np1;
    % Inkrement plastische Dehnung
    depsp = w3d2 * dpiter(iter) * n_np1;
    % Prädiktor Backstress
    alpha_star = alpha_n + 2/3 * h_i .* (A * depsp);
    % Check Konvergenz
    normp = abs(1-dpiter(iter)/dpiter(iter-1));
    % verhindere endlosschleife
    if iter > maxiter + 1
        msg = ['Keine Konvergenz in Fixpunktiteration nach ', num2str(iter), ...
               ' Iterationen'];
        error(msg)
    end
end % Ende Iterationsschleife

% -------------------------------------------------------------------------
%                        update Zustandsvariablen
% ------------------------------------------------------------------------- 
%         X         -> Zustandsvariablen [eps;epsp;alphai;p] bein spansteu
%                                        [sig;epsp;alphai;p] bein dehnsteu
%         ink_flag  -> (0) Spannungsgesteuerter Prozess 
%                      (1) Dehnungsgesteuerter Prozess
%         parameter -> Materialparameter des Modells
%                      [E,nu,c_i,r_i,r0]
% Plastische Dehnungen
epsp_np1 = epsp_n + depsp;
% plastische Bogenlänge
p_np1 = p_n + dpiter(iter);
% Backstresstensoren
alpha_np1 = theta_i .* alpha_star;
% Inkrement in Spannungen oder Dehnungen
if ink_flag == 0
    eps_np1 = eps_n + D * ink + depsp;
    out = [eps_np1; epsp_np1];
elseif ink_flag == 1
    sig_np1 = P\(beta_np1 + sum(alpha_np1,2));
    out = [sig_np1; epsp_np1];
end
% F_np1 = sqrt(3/2 * beta_np1' * P_line * beta_np1) - r0;
% fprintf('F = %.2d\n',F_np1);
X_neu = [out ; reshape(alpha_np1,ntens*M,1); p_np1];   

end % Ende 3D Materialmodell

