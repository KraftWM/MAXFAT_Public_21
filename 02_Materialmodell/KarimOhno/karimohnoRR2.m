function [ZVARneu] = karimohnoRR2(ntens, ndi, ink, ZVAR, ink_flag, parameter)
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
% Integration mit Radial Return Methode nach eigenem Algorithmuss
%
% Implementiert sind:
%        3D Spannungszustand
%        ESZ
%
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
% |  Stand: Oktober 2021                                                     |
%  ----------------------------------------------------------------------


%--------------------------------------------------------------------------
%                  Materialparamter ermitteln                                  
%--------------------------------------------------------------------------
% mynull = 1e-40;
% Anzahl Backstress
M = (length(parameter)-3)/3;                                               % anzahl Backstresstensoren
% Elastizitätskonstanten
E = parameter(1);
nu = parameter(2);
G = E/(2*(1+nu));
r0 = parameter(end);
% kinematische Verfestigung
% c_i = parameter(3:2+M);
% r_i = parameter(3+M:2+2*M);
% h_i = c_i.*r_i;

%--------------------------------------------------------------------------
%                        Hilfsmatrizen                                    
%--------------------------------------------------------------------------

% Steifigkeit
C = elast_steifigkeit(E,nu,ntens,ndi);
D = elast_nachgiebigkeit(E,nu,ntens,ndi);
% statische Matrizen
[P, P_line, P_tilde, A, P_check] = set_maps(ntens,ndi);

%--------------------------------------------------------------------------
%                  Zustandsvariablen auslesen                                 
%--------------------------------------------------------------------------

% Spannungen und Dehnungen
if ink_flag == 0 % spansteu
    eps = ZVAR(1:ntens);
    epsp = ZVAR(ntens+1:2*ntens);
    sig = C * (eps - epsp);
    dsig_tr = ink;    % elastischer trial step
    sig_tr = sig + ink;
elseif ink_flag == 1 % dehnsteu
    sig = ZVAR(1:ntens);
    dsig_tr = C * ink; % elastischer trial step
    epsp = ZVAR(ntens+1:2*ntens);
    eps = D*sig + epsp;
    sig_tr = sig + dsig_tr;
else
    msg = '      keine oder flasche angabe der Laststeuerung';
    error(msg)
end
% Backstresstensoren
alpha = reshape( ZVAR(2*ntens+1:(M+2)*ntens) , ntens, M);
p = ZVAR((M+2)*ntens+1);

%--------------------------------------------------------------------------
%                       Elastischer Trial Step (für ganzes Inkrement)
%--------------------------------------------------------------------------

% Trial Inkrement der Spannungen 
s = P * sig;                                               % Span dev
ds = P * dsig_tr;                                          % dev des Spannungsinkrements
s_tr = s + ds;                                             % Trial Spannungsdeviator
a = sum(alpha,2);                                          % Backstress;
beta = s_tr - a;                                           % Relativ Spannung
F_tr = 3/2 * (beta' * (P_line * beta)) - r0^2;             % Trial Fließfunktion

FTOL = 1e-7;                                               % Toleranz für Abweichungen F ~= 0

%--------------------------------------------------------------------------
%                   Trial Step wird angenommen
%--------------------------------------------------------------------------

if F_tr < FTOL   % Trial step völlig elastische
    
    ZVARneu = NaN;    % Dummy Wert für Compiler
    if ink_flag == 0 % spansteu;
        eps = eps + D * ink;
        ZVARneu = [eps;ZVAR(ntens+1:end)];
    elseif ink_flag == 1 % dehnsteu
        sig = sig + dsig_tr;
        ZVARneu = [sig;ZVAR(ntens+1:end)];
    end
    
    
%--------------------------------------------------------------------------
%                   Trial Step wird abgelehnt
%--------------------------------------------------------------------------
else
    
    
    % Integration je nach Spannungszustand
    if ntens == 6 % 3D
        
        % abbildungen
        %             [P, P_hat, P_tilde] = set_maps(ntens,ndi);
        % Integration radial return
        % Schleife über alle subinkremente
        ZVARneu = karimohno3D(...
            eps,epsp,alpha,p,...                                   % Zustandsvariablen
            sig_tr,...                                                 % Versuchsspannung
            ink,ink_flag,...                                           % skaliertes Inkrement und Steuerung
            parameter,G,...                                            % Parameter
            D,P,P_line,P_tilde);                                        % Abbildungen
        
    else
        
        % abbildungen
        %             [P, P_line, P_tilde, A, P_check] = set_maps(ntens,ndi);
        % Integration radial return
        ZVARneu = karimohnoESZ(...
            sig,eps,epsp,alpha,p,...                                       % Zustandsvariablen
            sig_tr,...                                                     % Versuchsspannung
            ink,ink_flag,...                                               % skaliertes Inkrement und Steuerung
            parameter,G,...                                                % Parameter
            D,P,P_line,P_tilde,A,P_check);                                 % Abbildungen
        
    end
    
end % Ende Unterscheidung Trial Step
    
   
end % Ende Hauptfunktion




% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
% Materialmodell 3D
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
function ZVARneu = karimohno3D(...
                       eps_n,epsp_n,alpha_n,p_n,...                  % Zustandsvariablen
                       sig_tr,...                                          % Versuchsspannung
                       ink,ink_flag,...                                    % Inkrement und Steuerung
                       parameter,G,...                                     % Parameter
                       D,P,P_hat,P_tilde)
                   
                   
% -------------------------------------------------------------------------
%                        Materialparameter
% -------------------------------------------------------------------------
% Zuweisen der Materialparameter
ntens = 6;                                                                 % Tensorkomponenten
M = (length(parameter)-3)/3;                                               % anzahl Backstresstensoren

% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
mu_i = parameter(3+2*M:2+3*M);
r0 = parameter(end);                                                       % startradius fliessfläche
h_i = c_i.*r_i;

% oft verwendete Konstanden
w3d2 = sqrt(3/2);
w2d3 = sqrt(2/3);
% oft verwendete Konstanden
mynull = 1e-40;

% -------------------------------------------------------------------------
%                        Korrektor
% -------------------------------------------------------------------------
if ink_flag % Dehnungssteuerung
    delta = 1;
else        % Spannungssteuerung
    delta = 0;
end


% -------------------------------------------------------------------------
%                        Definitionen für Iteration
% ------------------------------------------------------------------------- 

% Speicher und Abbruchbedingungen
maxiter = 1000;                   % Maximale Iterationsschleifen
dpiter = zeros(1,maxiter+2);     % Speicher plastische dehnungsinkremente
normpiter = zeros(1,maxiter+2);  % Norm Abbruchbedingung (für Konvergnz Analyse)
iter = 1;                        % Schleifenzähler 
tol = 1e-4;                      % toleranz Abbruchbedingung
tolbs = 1e-6;                    % toleranz iteration Backstresstensor
normp = 1;                       % Abbruchbedingung


% -------------------------------------------------------------------------
%                   Initialisiere startwerte Iteration
% ------------------------------------------------------------------------- 

% ... Korrekturfaktoren
theta_i = ones(1,M);

% ... Backstress
alpha_npe = alpha_n;

% ... trial Spannungsdeviator
s_tr = P * sig_tr;


% -------------------------------------------------------------------------
%                   erster Iterationsschritt
% ------------------------------------------------------------------------- 
% ... Inkrement plastische Bogenlänge
alpha_hat = sum(theta_i.*alpha_n,2);
beta_hat = s_tr - alpha_hat;

% ---------------------
% Fixpunktiteration 
% ---------------------
% dp_npe = (sqrt(1.5*(beta_hat'*P_hat*beta_hat)) - r0) ...
%          /(r0+(3*G*delta + sum(h_i.*theta_i)));
dp_npe = (sqrt(1.5*(beta_hat'*P_hat*beta_hat)) - r0) ...
         /((3*G*delta + sum(h_i.*theta_i)));    
dpiter(iter) = dp_npe;

% ... effektive Spannung & Fließflächen Normale
beta_npe = r0 * beta_hat / ( r0 + dp_npe * (3*G*delta + sum(h_i.*theta_i)));
n_npe = sqrt(beta_npe' * P_hat * beta_npe); if n_npe == 0; n_npe = mynull; end
n_npe = (P_hat * beta_npe)/n_npe;
% ... Update Backstress
theta_i_hash = 1./(1+c_i.*mu_i*dp_npe);
alpha_star = alpha_n + w2d3 * dp_npe * (P_tilde*n_npe) * h_i; 
alpha_hash = theta_i_hash .* alpha_star;
norm_ai_hash = w3d2 * sqrt( sum( (P_hat*alpha_hash) .* alpha_hash) );
fi_hash = norm_ai_hash - r_i;
idx = fi_hash > 0;
fi_hash(idx) = 1;
fi_hash(~idx) = 0;
theta_i = theta_i_hash .* ( 1 +  fi_hash .* (r_i./norm_ai_hash - 1));
% -------------------------------------------------------------------------
%                        Iteration
% ------------------------------------------------------------------------- 
% ... Schleife
while normp > tol
    
    % ... Inkrementiere Iterationszähler
    iter = iter + 1;
    % ... Inkrement plastische Bogenlänge
    alpha_hat = sum(theta_i.*alpha_n,2);
    beta_hat = s_tr - alpha_hat;
    % --------------------
    % Fixpunktiteration
    % --------------------
    % Aitkens Delta^2 Verfahren
%     dp_npe = (sqrt(1.5*(beta_hat'*P_hat*beta_hat)) - r0) ...
%          /(r0+(3*G*delta + sum(h_i.*theta_i)));  
    dp_npe = (sqrt(1.5*(beta_hat'*P_hat*beta_hat)) - r0) ...
         /((3*G*delta + sum(h_i.*theta_i)));  
    if mod(iter,3) == 0
        dp_hat = dp_npe - ( dp_npe - dpiter(iter-1) )^2 / ...
        ( dp_npe - 2*dpiter(iter-1) + dpiter(iter-2) );
        if dp_hat > 0
            dp_npe = dp_hat;
        end
    end
    dpiter(iter) = dp_npe;
    % ... effektive Spannung & Fließflächen Normale
    beta_npe = r0 * beta_hat / ( r0 + dp_npe * (3*G*delta + sum(h_i.*theta_i)));
    n_npe = sqrt(beta_npe' * P_hat * beta_npe); if n_npe == 0; n_npe = mynull; end
    n_npe = (P_hat * beta_npe)/n_npe;
    % ... Update Backstress
    theta_i_hash = 1./(1+c_i.*mu_i*dp_npe);
    alpha_star = alpha_n + w2d3 * dp_npe * (P_tilde*n_npe) * h_i;
    alpha_hash = theta_i_hash .* alpha_star;
    norm_ai_hash = w3d2 * sqrt( sum( (P_hat*alpha_hash) .* alpha_hash) );
    fi_hash = norm_ai_hash - r_i;
    idx = fi_hash > 0;
    fi_hash(idx) = 1;
    fi_hash(~idx) = 0;
    theta_i = theta_i_hash .* ( 1 +  fi_hash .* (r_i./norm_ai_hash - 1));
    % ... Prüfe Konvergenz
    normp = abs(1-dpiter(iter)/dpiter(iter-1));
    normpiter(iter) = normp;
    % ... Verhindere Endlosschleifen
    if iter > maxiter + 1
        msg = ['Keine Konvergenz in Fixpunktiteration im Ohno Wang Modell ',...
               'nach ', num2str(iter),' Iterationen. Aktuelle Fehlernorm: '...
               num2str(dpiter(iter))];
        warning(msg)
        break;
    end
        
end % Ende Äußere Iterationsschleife
% fprintf('iter = %i normp = %.5f\n',iter,normp)
% figure(1)
% semilogy(normpiter(1:iter),'d-') 
% axis([1e0 iter+10 1e-6 max(normpiter)])
% figure(3)
% plot(dpiter(1:iter)), pause(0.5)


% -------------------------------------------------------------------------
%                        update Zustandsvariablen
% ------------------------------------------------------------------------- 
% ... Backstress Tensoren
alpha_npe = theta_i .* alpha_star;
% ... plastische Dehnung
depsp = w3d2 * dp_npe * n_npe;
epsp_npe = epsp_n + depsp;
% ... plastische Bogenlänge
p_npe = p_n + dp_npe;
% .. Rest
out = [];              % Dummy init für Compiler
if ink_flag == 0
    eps_npe = eps_n + D * ink + depsp;
    out = [eps_npe; epsp_npe];
elseif ink_flag == 1
    sig_npe = sig_tr - 2*G * P_tilde * depsp;
    out = [sig_npe; epsp_npe];
end
ZVARneu = [out ; reshape(alpha_npe,ntens*M,1); p_npe];   

end % Ende Modellgleichung in 3D




% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
% Materialmodell ESZ
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
function ZVARneu = karimohnoESZ(...
                       sig_n,eps_n,epsp_n,alpha_n,p_n,...                  % Zustandsvariablen
                       sig_tr,...                                          % Versuchsspannung
                       ink,ink_flag,...                                    % Inkrement und Steuerung
                       parameter,G,...                                     % Parameter
                       D,P,P_line,P_hat,A,~)                       % Abblidungen
                   
% -------------------------------------------------------------------------
%                        Materialparameter
% -------------------------------------------------------------------------
% Zuweisen der Materialparameter
ntens = 3;                                                                 % Tensorkomponenten
M = (length(parameter)-3)/3;                                               % anzahl Backstresstensoren
nu = parameter(2);
% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
mu_i = parameter(3+2*M:2+3*M);
r0 = parameter(end);                                                       % startradius fliessfläche
h_i = c_i.*r_i;

% oft verwendete Konstanden
w3d2 = sqrt(3/2);
w2d3 = sqrt(2/3);
% oft verwendete Konstanden
mynull = 1e-40;
DELTA3 = [-1/3;-1/3;0];  % Hilfstensor
fak = (1-2*nu)/(1-nu);

% -------------------------------------------------------------------------
%                        Korrektor
% -------------------------------------------------------------------------
if ink_flag % Dehnungssteuerung
    delta = 1;
else        % Spannungssteuerung
    delta = 0;
end


% -------------------------------------------------------------------------
%                        Definitionen für Iteration
% ------------------------------------------------------------------------- 

% Speicher und Abbruchbedingungen
maxiter = 1000;                   % Maximale Iterationsschleifen
dpiter = zeros(1,maxiter+2);     % Speicher plastische dehnungsinkremente
normpiter = zeros(1,maxiter+2);  % Norm Abbruchbedingung (für Konvergnz Analyse)
iter = 1;                        % Schleifenzähler 
tol = 1e-4;                      % toleranz Abbruchbedingung
tolbs = 1e-6;                    % toleranz iteration Backstresstensor
normp = 1;                       % Abbruchbedingung


% -------------------------------------------------------------------------
%                   Initialisiere startwerte Iteration
% ------------------------------------------------------------------------- 

% ... Korrekturfaktoren
dp_npe = 0;
theta_i = ones(1,M);
% ... Backstress
alpha_npe = alpha_n;
% ... trial Spannungsdeviator
s_tr = P * sig_tr;
% ... trial Fließflächennormale
n_npe = P * sig_n - sum(alpha_npe,2);
var1 = sqrt(n_npe' * P_line * n_npe); if var1 == 0, var1 = mynull; end
n_npe = (P_hat * n_npe)./var1; 

% -------------------------------------------------------------------------
%                   erster Iterationsschritt
% ------------------------------------------------------------------------- 
% ... Inkrement plastische Bogenlänge
alpha_hat = sum(theta_i.*alpha_n,2) ...
            - w3d2*delta *2*G*fak*dp_npe*(-n_npe(1) - n_npe(2))*DELTA3;
beta_hat = s_tr - alpha_hat;

% ---------------------
% Fixpunktiteration 
% ---------------------
% dp_npe = (sqrt(1.5*(beta_hat'*P_line*beta_hat)) - r0) ...
%          /(r0+(3*G*delta + sum(h_i.*theta_i)));    
dp_npe = (sqrt(1.5*(beta_hat'*P_line*beta_hat)) - r0) ...
         /((3*G*delta + sum(h_i.*theta_i)));
dpiter(iter) = dp_npe;

% ... effektive Spannung & Fließflächen Normale
beta_npe = r0 * beta_hat / ( r0 + dp_npe * (3*G*delta + sum(h_i.*theta_i)));
n_npe = sqrt(beta_npe' * P_line * beta_npe); if n_npe == 0; n_npe = mynull; end
n_npe = (P_hat * beta_npe)/n_npe;
% ... Update Backstress
theta_i_hash = 1./(1+c_i.*mu_i*dp_npe);
alpha_star = alpha_n + w2d3 * dp_npe * (A*n_npe) * h_i; 
alpha_hash = theta_i_hash .* alpha_star;
norm_ai_hash = w3d2 * sqrt( sum( (P_line*alpha_hash) .* alpha_hash) );
fi_hash = norm_ai_hash - r_i;
idx = fi_hash > 0;
fi_hash(idx) = 1;
fi_hash(~idx) = 0;
theta_i = theta_i_hash .* ( 1 +  fi_hash .* (r_i./norm_ai_hash - 1));
% -------------------------------------------------------------------------
%                        Iteration
% ------------------------------------------------------------------------- 
% ... Schleife
while normp > tol
    % ... Inkrementieren Schleifenzähler
    iter = iter + 1;
    % ... Inkrement plastische Bogenlänge
    alpha_hat = sum(theta_i.*alpha_n,2) ...
        - w3d2*delta *2*G*fak*dp_npe*(-n_npe(1) - n_npe(2))*DELTA3;
    beta_hat = s_tr - alpha_hat;
    
    % ---------------------
    % Fixpunktiteration
    % ---------------------
%     dp_npe = (sqrt(1.5*(beta_hat'*P_line*beta_hat)) - r0) ...
%          /(r0+(3*G*delta + sum(h_i.*theta_i))); 
    dp_npe = (sqrt(1.5*(beta_hat'*P_line*beta_hat)) - r0) ...
         /((3*G*delta + sum(h_i.*theta_i)));
    % Aitkens Delta^2 Verfahren
    if mod(iter,3) == 0
        dp_hat = dp_npe - ( dp_npe - dpiter(iter-1) )^2 / ...
        ( dp_npe - 2*dpiter(iter-1) + dpiter(iter-2) );
        if dp_hat > 0
            dp_npe = dp_hat;
        end
    end
    dpiter(iter) = dp_npe;
    % ... effektive Spannung & Fließflächen Normale
    beta_npe = r0 * beta_hat / ( r0 + dp_npe * (3*G*delta + sum(h_i.*theta_i)));
    n_npe = sqrt(beta_npe' * P_line * beta_npe); if n_npe == 0; n_npe = mynull; end
    n_npe = (P_hat * beta_npe)/n_npe;
    % ... Update Backstress
    theta_i_hash = 1./(1+c_i.*mu_i*dp_npe);
    alpha_star = alpha_n + w2d3 * dp_npe * (A*n_npe) * h_i;
    alpha_hash = theta_i_hash .* alpha_star;
    norm_ai_hash = w3d2 * sqrt( sum( (P_line*alpha_hash) .* alpha_hash) );
    fi_hash = norm_ai_hash - r_i;
    idx = fi_hash > 0;
    fi_hash(idx) = 1;
    fi_hash(~idx) = 0;
    theta_i = theta_i_hash .* ( 1 +  fi_hash .* (r_i./norm_ai_hash - 1));
%     fprintf('||n|| = %.15f\n',n_npe' * A * P_line * A * n_npe)
    % ... Prüfe Konvergenz
    normp = abs(1-dpiter(iter)/dpiter(iter-1));
    normpiter(iter) = normp;
    % ... Verhindere Endlosschleifen
    if iter > maxiter + 1
        msg = ['Keine Konvergenz in Fixpunktiteration im Ohno Wang Modell ',...
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
alpha_npe = theta_i .* alpha_star;
depsp = w3d2*dp_npe * n_npe;
epsp_npe = epsp_n + depsp;
p_npe = p_n + dp_npe;
out = [];              % Dummy init für Compiler
if ink_flag == 0
    eps_npe = eps_n + D * ink + depsp;
    out = [eps_npe; epsp_npe];
elseif ink_flag == 1
%     sig_npe = sig_tr - 2*G * A * depsp ...
%     + (parameter(1)/(3*(1-2*nu))*[1;1;0] +2*G*DELTA3)*fak*(-depsp(1)-depsp(2));
    s_npe = beta_npe + sum(alpha_npe,2);
    sig_npe = s_npe + (s_npe(1)+s_npe(2)) * [1 1 0]';
    out = [sig_npe; epsp_npe];
end
ZVARneu = [out ; reshape(alpha_npe,ntens*M,1); p_npe];   

end % Ende Modellgleichung in ESZ
