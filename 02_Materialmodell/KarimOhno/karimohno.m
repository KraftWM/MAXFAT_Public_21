function [X_neu, DEP , CEP] = karimohno(ntens, ndi, ink, X, ink_flag, parameter)
% Implementierung des Materialmodells nach Karim und Ohno,
% Kinematische Verfestigung ist lineare Interpolation zwischen Armstrong
% Fredrick und OhnoWang Modell I Kinematik
%
% QUELLE MODELL:
% Aus Ohno,Karim et al. 2000 Uniaxial Ratchetting of 316FR steel at room
% tempreture part II
%
% Karim,Ohno 2000 Kinematic hardeningmodel suitable for ratchetting with
% steady state
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
%                      [E,nu,c_i,r_i,mu_i,r0]
%                       mu_i(i) = 0 -> OhnoWang Modell 1 Kinematik
%                       mu_i(i) = 1 -> Armstrong Frederick Kinematik
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
% |  Stand: Mai 2020                                                     |
%  ----------------------------------------------------------------------

%--------------------------------------------------------------------------
%                  Materialparamter ermitteln                                  
%--------------------------------------------------------------------------

% Anzahl Backstress
M = (length(parameter)-3)/3;
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
FTOL = 1e-7;                        % Toleranz für Abweichungen F ~= 0


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
    
    if nargout == 2 % berechne Elastisch plastische steifigkeitsmatrix
        DEP = D;
    elseif nargout == 3
        DEP = D;
        CEP = C;
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
        % Matlab version
%         X_neu = karimohno3D(eps, epsp, alpha,p,...                       % Zustand
%                            sig_tr,ink,ink_flag,...                       % Inkrement
%                            parameter,...                                 % parameter
%                            C,D,P,P_hat,P_tilde);                         % Abbildungen
        % mex version
        X_neu = CoderVersion_3D_KarimOhno_mex(eps, epsp, alpha,p,...       % Zustand
                           sig_tr,ink,ink_flag,...                         % Inkrement
                           parameter,...                                   % parameter
                           C,D,P,P_hat,P_tilde);                           % Abbildungen


    elseif ntens == 3 && ndi == 2 % ESZ
        
        % statische Matrizen
        [P, P_line, P_hat, A, P_check] = set_maps(ntens,ndi);
        % Trenne Spannungs und Dehnungssteuerung
        if ink_flag == 0 % Spannungssteuerung
            % Matlab version
%             X_neu = karimohnoESZ(eps, epsp, alpha,p,...
%                                sig_tr,ink,ink_flag,...
%                                parameter,...
%                                D,P,P_line,P_hat,A);
            % Mex version
            X_neu = CoderVersion_ESZ_KarimOhno_mex(eps, epsp, alpha,p,...
                               sig_tr,ink,ink_flag,...
                               parameter,...
                               D,P,P_line,P_hat,A);
        else % Dehnungssteuerung, umschalten zu expliziter Integration
            
            % 1. ermittle elastischen Anteil
%             xel = elastink(s,a,r0,ds,ndi);
            xel = elastink2(s,a,P_line,r0,ds,FTOL);  
            % 2. ausführen elastischen Teil
            if ink_flag == 0 % spansteu;
                D = elast_nachgiebigkeit(E,nu,ntens,ndi);
                eps = eps + D * (xel*ink);
                X = [eps;X(ntens+1:end)];
            elseif ink_flag == 1 % dehnsteu
                sig = sig + xel*dsig_tr;
                X = [sig;X(ntens+1:end)];
            end

            % 3. plastischer Teil des inkrements
            ink = (1-xel)*ink;
            
            % Integration je nach Spannungszustand
            options = [];
            
            % Explizite Integration Matlab Version
            [~,X_neu] = rk87(@karimohnoESZ_exp,[0,1], X, options, ink,...
                            ink_flag, parameter, C, D, P, P_line, P_hat,...
                            A, P_check);
                        
            % nur letzten Zeitschritt ausgeben
             X_neu = X_neu(end,:)';
             
        end % Ende Verzweigung Spannungs/Dehnungssteuerung
        
    elseif ntens == 1 % 1D
        
        % Matlab Version
%         X_neu = karimohno1D(eps, epsp, alpha,p,...
%                            sig_tr,ink,ink_flag,...
%                            parameter...
%                            );
        
        % Mex Version
        X_neu = CoderVersion_1D_KarimOhno_mex(eps, epsp, alpha,p,...
                           sig_tr,ink,ink_flag,...
                           parameter...
                           );
        
    end
    
    % Falls tangentiale Steifigkeit gebraucht wird
    if nargout == 2 % berechne Elastisch plastische steifigkeitsmatrix
        DEP = tangnach_karim(X_neu,parameter,ink_flag,ntens,ndi,C,D);
    elseif nargout == 3
        DEP = tangnach_karim(X_neu,parameter,ink_flag,ntens,ndi,C,D);
        CEP = tangsteifigkeit_karim(X_neu,parameter,ink_flag,ntens,ndi);
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

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen 3D Spannungszustand                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X_neu = karimohno3D(eps_n, epsp_n, alpha_n,p_n,...
                           sig_tr,ink,ink_flag,...
                           parameter,...
                           C,D,P,P_hat,P_tilde)
% Konkretes Materialmodell für 3D Spannungszustände, gibt bei vorgabe eines
% Lastinkrementes die geupdatetn Zustandsvariablen mit implizitem Euler 
% zurück .
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
M = (length(parameter)-3)/3;                                               % Anzahl TBST

% elastische 
E = parameter(1);
nu = parameter(2);
G = E/(2*(1+nu));
r0 = parameter(end);                                                       % startradius fliessfläche
% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
mu_i = parameter(3+2*M:2+3*M);
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
theta_i_tilde = ones(1,M);
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
else %ink_flag == 1
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
j_i = r_i./(r_i + mu_i.*h_i*dpiter(iter));
alpha_hash = j_i .* alpha_star;
% Iterationsschleife
while normp > tol
    % effektive Backstresstensoren
    alpha_eff = sqrt( 3/2 * sum( (P_hat*alpha_hash) .* alpha_hash) );
    idx1 = alpha_eff <= r_i;
    idx2 = alpha_eff > r_i;
    theta_i_tilde( idx1 ) = 1;
    theta_i_tilde( idx2 ) = r_i(idx2)./alpha_eff(idx2);
    theta_i = j_i .* theta_i_tilde;
    % Inkrementieren Schleifenzähler
    iter = iter + 1;
    % neues plastisches Dehnungsinkrement
    dpiter(iter) = verfahren(s_tr,alpha_n,Y_np1,G,h_i,theta_i,P_hat);
    % Aitkens Delta^2 Verfahren
    if mod(iter,3) == 0
        dpiter(iter) = deltaquadrat(dpiter,iter);
    end
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
    j_i = r_i./(r_i + mu_i.*h_i*dpiter(iter));
    alpha_hash = j_i .* alpha_star;
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
%                      [E,nu,c_i,r_i,mu_i,r0]
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
    sig_np1 = sig_tr - 2*G * depsp;
    out = [sig_np1; epsp_np1];
end
% beta_np1 = P * sig_np1 - sum(alpha_np1,2);
% F_np1 = sqrt(3/2 * beta_np1' * P_hat * beta_np1) - r0;
% fprintf('F = %.2d\n',F_np1);
% fprintf('n = %i\n',iter-1);
X_neu = [out ; reshape(alpha_np1,ntens*M,1); p_np1];   

end % Ende 3D Materialmodell































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen ebener Spannungszustand             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X_neu = karimohnoESZ(eps_n, epsp_n, alpha_n,p_n,...
                           sig_tr,ink,ink_flag,...
                           parameter,...
                           D,P,P_line,P_hat,A)
% Konkretes Materialmodell für ESZ Spannungszustände, gibt bei vorgabe 
% eines Lastinkrementes die geupdatetn Zustandsvariablen mit implizitem 
% Euler zurück .
%
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
M = (length(parameter)-3)/3;                                               % Anzahl TBST

% elastische 
E = parameter(1);
nu = parameter(2);
G = E/(2*(1+nu));
r0 = parameter(end);                                                       % startradius fliessfläche
% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
mu_i = parameter(3+2*M:2+3*M);
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
theta_i_tilde = ones(1,M);
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
j_i = r_i./(r_i + mu_i.*h_i*dpiter(iter));
alpha_hash = j_i .* alpha_star;
% Iterationsschleife
while normp > tol
    % effektive Backstresstensoren
    alpha_eff = sqrt( 3/2 * sum( (P_line*alpha_hash) .* alpha_hash) );
    idx1 = alpha_eff <= r_i;
    idx2 = alpha_eff > r_i;
    theta_i_tilde( idx1 ) = 1;
    theta_i_tilde( idx2 ) = r_i(idx2)./alpha_eff(idx2);
    theta_i = j_i .* theta_i_tilde;
    % Inkrementieren Schleifenzähler
    iter = iter + 1;
    % neues plastisches Dehnungsinkrement
    dpiter(iter) = verfahren(s_tr,alpha_n,Y_np1,G,h_i,theta_i,P_line);
    % Aitkens Delta^2 Verfahren
    if mod(iter,3) == 0
        dpiter(iter) = deltaquadrat(dpiter,iter);
    end
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
    j_i = r_i./(r_i + mu_i.*h_i*dpiter(iter));
    alpha_hash = j_i .* alpha_star;
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
% fprintf('n = %i\n',iter-1);
X_neu = [out ; reshape(alpha_np1,ntens*M,1); p_np1];   

end % Ende ESZ Materialmodell

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Modellgleichungen Ebener Spannungszustand (explizite Integration)   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dX = ...
          karimohnoESZ_exp(~, X, ink, ink_flag, parameter, C, D, ...
                            P, P_line, P_hat, A, P_check)
% Konkretes Materialmodell für ESZ Spannungszustände, gibt bei vorgabe eines
% Lastinkrementes die Inkremente der inneren Variablen zurück. Dabei wird
% angenommen, das jedes übergebene Inkrement elastisch-plastische 
% Deformationen hervorruft.
%
%   INPUT:
%         t         -> Zeit (auch wenn nicht explizit gebraucht)
%         X         -> Zustand
%         ink       -> Lastinkrement
%         ink_flag  -> Belastungsart
%         parameter -> Modellparameter
%                      = 2 + num_alpha * 2 + 1
%         C,D       -> Elastische Steifigkeit und Nachgiebigkeit   
%   P, P_line ...   -> Diverse Abbildungen
%
%
%   OUTPUT:
%         dX -> Inkremente der Zustandsvariablen
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
%  ----------------------------------------------------------------------


%--------------------------------------------------------------------------
%                   Materialparameter
%--------------------------------------------------------------------------
% Zuweisen der Materialparameter
ntens = 3;                                                                 % Tensorkomponenten
M = (length(parameter)-3)/3;                                               % Anzahl TBST

% elastische 
r0 = parameter(end);                                                       % startradius fliessfläche
% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
mu_i = parameter(3+2*M:2+3*M);
h_i = c_i .* r_i;
% oft verwendete Konstanden
w3d2 = sqrt(3/2);
w2d3 = sqrt(2/3);
delta = 1e-40;
hohezahl = 1000000;
%--------------------------------------------------------------------------
%                   Zustandsvariablen
%--------------------------------------------------------------------------

% Spannungen und Dehnungen
if ink_flag == 0 % Spansteu.
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = C * (eps - epsp);
elseif ink_flag == 1 % Dehnsteu
    sig = X(1:ntens);
    % epsp = X(ntens+1:2*ntens);
end

% backstresstensoren
alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);

%--------------------------------------------------------------------------
%                   Normale an die Fließfläche
%--------------------------------------------------------------------------

% Spannungsdeviator
s = P * sig;
% Backstress
a = sum(alpha,2);
% normale an Fließfläche
n = (w3d2/r0).* P_hat * (s-a);
nTrans = n.';

%--------------------------------------------------------------------------
%                   Ableitung der Teilbackstresstensoren
%--------------------------------------------------------------------------

% Normen der Teilbackstresstensoren
norm_ai = sqrt( sum( (P_line*alpha) .* alpha ) );
norm_ai(norm_ai == 0) = delta; 

% Hilfvariablen
% Li = alpha./norm_ai;
% var1 = c_i .* r_i;
Li = w2d3*alpha./norm_ai;
var2 = A * n;
var3 = w3d2*norm_ai./r_i;
var3(var3 >= 1) = 1;
var3(var3 < 1) = 0;
% idx1 = var3 >= 1;
% idx0 = var3 < 1;
% var3(idx1) = 1;
% var3(idx0) = 0;
var4 = w3d2 * nTrans * P_check * Li - mu_i;
var4 = 0.5 * (var4 + abs(var4));

% Ableitung Teilbackstresstensoren
dalpha_dp = zeros(ntens,M);
for ii = 1 : M

    dalpha_dp(:,ii) = h_i(ii) * ( ...
                                w2d3 * var2 - ...
                                mu_i(ii)/r_i(ii) .* alpha(:,ii) - ...
                                var3(ii) * var4(ii) .* Li(:,ii) ...
                                );
    
end
% Ableitung gesamtbackstresstensor
da_dp = sum(dalpha_dp,2);


%--------------------------------------------------------------------------
%                   plastischer Tangentenmodul
%--------------------------------------------------------------------------

h = w3d2 * nTrans * P_check * da_dp;

%--------------------------------------------------------------------------
%                   inkremente plastische bogenlänge
%--------------------------------------------------------------------------

if ink_flag == 0 % spannungssteuerung
    dp = (w3d2/h) * (nTrans * ink); 
elseif ink_flag == 1 % dehnungssteuerung
    dp = (w3d2 .* nTrans * (C * ink)) / ...
         (h + 3/2.* nTrans * ( C * n));
else % Fehler
    msg = '      keine oder flasche angabe der Laststeuerung';
    error(msg)
end

%--------------------------------------------------------------------------
%                   inkremente der zustandsvariablen
%--------------------------------------------------------------------------

depsp=dp.*w3d2.*n;
dalpha=dalpha_dp.*dp;

%--------------------------------------------------------------------------
%                   Zusammenfassen der Inkremente
%--------------------------------------------------------------------------

if ink_flag == 0 % spansteu
    dout = D * ink + depsp;
elseif ink_flag == 1 % dehnsteu
    dout = C * (ink - depsp);
end
dX = [dout; depsp ; reshape(dalpha,ntens*M,1); dp]; 
    
end % Ende Modell ESZ explizite Integration

































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen 1D Spannungszustand                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X_neu = karimohno1D(eps_n, epsp_n, alpha_n,p_n,...
                           sig_tr,ink,ink_flag,...
                           parameter...
                           )
% Konkretes Materialmodell für 1D Spannungszustände, gibt bei vorgabe 
% eines Lastinkrementes die geupdatetn Zustandsvariablen mit implizitem 
% Euler zurück .
%
%   INPUT:
%       sig_n       -> Spannungen
%  eps_n,epsp_n     -> Dehnungen
%       alpha_n     -> Backstresstensoren
%       dsig_tr     -> Versuchsspannung ( richtige spannung bei
%                      spannungssteuerung )
%         ink       -> Lastinkrement
%         ink_flag  -> Belastungsart
%         parameter -> Modellparameter
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
ntens = 1;                                                                 % Tensorkomponenten
M = (length(parameter)-3)/3;                                               % Anzahl TBST

% elastische 
E = parameter(1);
nu = parameter(2);
G = E/(2*(1+nu));
r0 = parameter(end);                                                       % startradius fliessfläche
% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
mu_i = parameter(3+2*M:2+3*M);
h_i = c_i .* r_i;
% oft verwendete Konstanden
w3d2 = sqrt(3/2);
w2d3 = sqrt(2/3);

% -------------------------------------------------------------------------
%                        Iterationsschleife
% ------------------------------------------------------------------------- 

% Speicher und Abbruchbedingungen
maxiter = 100;                % Maximale Iterationsschleifen
dpiter = zeros(1,maxiter+3);  % Speicher plastische dehnungsinkremente
iter = 1;                     % Schleifenzähler 
tol = 1e-4;                   % toleranz Abbruchbedingung
% Startwerte
normp = 1;
theta_i = ones(1,M);
theta_i_tilde = ones(1,M);
% neuer radius FF = alter radius FF
Y_np1 = r0;
% erster schritt
if ink_flag == 0 %Spansteu
    dummy = sum(theta_i .* h_i);
    dpiter(iter) = w2d3 * (abs(sig_tr-sum(theta_i.*alpha_n)) - Y_np1)/dummy;
else
    dummy = sum(theta_i .* h_i) + E;
    dpiter(iter) = w2d3 * (abs(sig_tr-sum(theta_i.*alpha_n)) - Y_np1)/dummy;
end
% Update Hilfvariable
if ink_flag == 0
    var1 = sum(theta_i .* h_i)*w3d2;
else% ink_flag == 1
    var1 = (E + sum(theta_i .* h_i))*w3d2;
end
% Update Hilfsvariable 2
var2 = sig_tr - sum(theta_i .* alpha_n);
% Update effektive spannung
beta_np1 = Y_np1 * var2 / (Y_np1 + var1 * dpiter(iter));
% Update Normale
n_np1 = beta_np1/Y_np1;
% Inkrement plastische Dehnung
depsp = w3d2 * dpiter(iter) * n_np1;
% Prädiktor Backstress
alpha_star = alpha_n + h_i .* depsp;
j_i = r_i./(r_i + mu_i.*h_i*dpiter(iter));
alpha_hash = j_i .* alpha_star;
% Iterationsschleife
while normp > tol
    % effektiver Backstresstensor
    alpha_eff = abs(alpha_hash);
    idx1 = alpha_eff <= r_i;
    idx2 = alpha_eff > r_i;
    theta_i_tilde( idx1 ) = 1;
    theta_i_tilde( idx2 ) = r_i(idx2)./alpha_eff(idx2);
    theta_i = j_i .* theta_i_tilde;
    % Inkrementieren Schleifenzähler
    iter = iter + 1;
    % neues plastisches Dehnungsinkrement
    if ink_flag == 0 %Spansteu
        dummy = sum(theta_i .* h_i);
        dpiter(iter) = w2d3 * (abs(sig_tr-sum(theta_i.*alpha_n)) - Y_np1)/dummy;
    else
        dummy = sum(theta_i .* h_i) + E;
        dpiter(iter) = w2d3 * (abs(sig_tr-sum(theta_i.*alpha_n)) - Y_np1)/dummy;
    end
    % Aitkens Delta^2 Verfahren
    if mod(iter,3) == 0
        dpiter(iter) = deltaquadrat(dpiter,iter);
    end
    % Update Hilfvariable
    if ink_flag == 0
        var1 = sum(theta_i .* h_i)*w3d2;
    else% ink_flag == 1
        var1 = (E + sum(theta_i .* h_i))*w3d2;
    end
    % Update Hilfsvariable 2
    var2 = sig_tr - sum(theta_i .* alpha_n);
    % Update effektive spannung
    beta_np1 = Y_np1 * var2 / (Y_np1 + var1 * dpiter(iter));
    % Update Normale
    n_np1 = beta_np1/Y_np1;
    % Inkrement plastische Dehnung
    depsp = w3d2 * dpiter(iter) * n_np1;
    % Prädiktor Backstress
    alpha_star = alpha_n + h_i .* depsp;
    j_i = r_i./(r_i + mu_i.*h_i*dpiter(iter));
    alpha_hash = j_i .* alpha_star;
    % Check Konvergenz
    normp = abs(1-dpiter(iter)/dpiter(iter-1));
    % verhindere endlosschleife
    if iter > maxiter + 1
        msg = ['Keine Konvergenz in Fixpunktiteration nach 101 Iterationen'];
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
    eps_np1 = eps_n + ink/E + depsp;
    out = [eps_np1; epsp_np1];
else% ink_flag == 1
    sig_np1 = sig_tr - E * depsp;
    out = [sig_np1; epsp_np1];
end
% F_np1 = sqrt(3/2 * beta_np1' * P_line * beta_np1) - r0;
% fprintf('F = %.2d\n',F_np1);
% fprintf('n = %i\n',iter-1);
X_neu = [out ; reshape(alpha_np1,ntens*M,1); p_np1];   

end % Ende 1D Spannungszustand



% -------------------------------------------------------------------------
% Hilfsfunktionen
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Tangentialer elastisch-plastischer Steifigkeitstensor
function CEP = tangsteifigkeit_karim(X_neu,parameter,ink_flag,ntens,ndi)
        % Matrizen je nach Spannungszustand
        P = set_maps(ntens,ndi);
        % Parameter 
        M = (length(parameter)-3)/3;
        E = parameter(1);
        nu = parameter(2);
        r0 = parameter(end);
        % kinematische Verfestigung
        c_i = parameter(3:2+M);
        r_i = parameter(3+M:2+2*M);
        mu_i = parameter(3+2*M:2+3*M);
        h_i = c_i .* r_i;
        CE = elast_steifigkeit(E,nu,ntens,ndi);
        % Neue Zustände
        if ink_flag == 0 % spansteu
            eps = X_neu(1:ntens);
            epsp = X_neu(ntens+1:2*ntens);
            sig = CE * (eps - epsp);
        elseif ink_flag == 1 % dehnsteu
            sig = X_neu(1:ntens);
        end
        alpha = zeros(ntens,M);
        for i = 1:M
            alpha(:,i) = X_neu((i+1)*ntens+1:(i+2)*ntens);
        end
        % berechne normale
        s = P * sig;
        a = sum(alpha,2);
        beta = s - a;
        n = (sqrt(3/2)/r0).*diag([ones(1,ndi),2*ones(1,ntens-ndi)])*beta;
        % plastischer Modul 
        h = tangmod_karim(n,alpha,h_i,c_i,mu_i,r_i,ntens,ndi);
        % tangentiale Steifigkeit
        dummy = CE * n;
        CP = 3/2 .*(dummy * dummy')/(3/2 * n' * dummy + h) ;
        CEP = CE - CP;
end


% -------------------------------------------------------------------------
% Tangentialer elastisch-plastischer Nachgibigkeitstensor
function DEP = tangnach_karim(X,parameter,ink_flag,ntens,ndi,CEL,DEL)
% Funktion berechnet nachgiebigkeit im elastisch-plastischen

    % Abbildungen
    if ntens == 6
        [P,Paa] = set_maps(ntens,ndi);
        Pna = diag(ones(1,6));
    elseif ntens == 3
        [P,Paa,~,~,Pna] = set_maps(ntens,ndi);
    else
        [P,Paa,Pna] = set_maps(ntens,ndi);
    end
    % Parameter
    M = (length(parameter)-3)/3;                                               % Anzahl TBST
    r0 = parameter(end);                                                       % startradius fliessfläche
    c_i = parameter(3:2+M);
    r_i = parameter(3+M:2+2*M);
    mu_i = parameter(3+2*M:2+3*M);
    h_i = c_i .* r_i;

    % Zustandsvariablen
    if ink_flag == 0 % Spansteu.
        eps = X(1:ntens);
        epsp = X(ntens+1:2*ntens);
        sig = CEL * (eps - epsp);
    elseif ink_flag == 1 % Dehnsteu
        sig = X(1:ntens);
        % epsp = X(ntens+1:2*ntens);
    end
    alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);
    % Spannungsdeviator
    s = P * sig;
    % Backstress
    a = sum(alpha,2);
    % Normale an die Fließfläche
    n = sqrt(3/2) * diag([1,1,2])*(s - a) / r0;
    % tangentenmodul
    h = tangmod_karim(n,alpha,h_i,c_i,mu_i,r_i,ntens,ndi);
    % Tangetiale Nachgiebigkeit
    DEP = DEL + sqrt(3/2) * (n * n')/h;
end



% -------------------------------------------------------------------------
% Tangenten Modul Karim Ohno 
function h = tangmod_karim(n,alpha,h_i,zeta_i,mu_i,r_i,ntens,ndi)
% Anzahl Backstresstensoren
M = length(h_i);
% Spannungsdev je nach Spannungszustand
if ntens == 1
    P_line = 1;
elseif ntens == 3 && ndi == 2
    P_line = [2,1,0;...
             1,2,0;...
             0,0,1];
elseif ntens == 6
    P_line = diag([1,1,1,1,1,1]);
end   
h = sum(h_i);
for i = 1:M
    na = n' * P_line * alpha(:,i);
    h = h - sqrt(3/2) * zeta_i(i) * mu_i(i)* na ;
    norma = sqrt(3/2) * sqrt( sum( (P_line * alpha(:,i) .* alpha(:,i)) ));
    if norma >= r_i(i)
        dummy = sqrt(3/2)/r_i(i)*na-mu_i(i);
        dummy = 0.5 * (dummy + abs(dummy));
        h = h + sqrt(3/2) * zeta_i(i) * na * dummy;
    end
end
end
