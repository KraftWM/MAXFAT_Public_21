%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Matlab Coder Version von:                                            %
%    Modellgleichung 3D Karim Ohno Modell                                 %
%                                                                         %
%    Aufgerufen in:                                                       %
%    karimohno.m                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen 3D Spannungszustand                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X_neu = CoderVersion_3D_KarimOhno(eps_n, epsp_n, alpha_n,p_n,...
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
dpiter = zeros(1,maxiter+3);  % Speicher plastische dehnungsinkremente
iter = 1;                   % Schleifenzähler 
tol = 1e-4;                 % toleranz Abbruchbedingung
% Startwerte
normp = 1;
theta_i = ones(1,M);
theta_i_tilde = ones(1,M);
% Setze verfahrensvorschrift je nach steuerung
% if ink_flag == 0
%     verfahren = @verfahren_sig;
% elseif ink_flag == 1
%     verfahren = @verfahren_eps;
% end
% deviatorischer Anteil Trial spannung
s_tr = P * sig_tr;
% neuer radius FF = alter radius FF
Y_np1 = r0;
% erster schritt
if ink_flag == 0
    dpiter(iter) = verfahren_sig(s_tr,alpha_n,Y_np1,G,h_i,theta_i,P_hat);
else
    dpiter(iter) = verfahren_eps(s_tr,alpha_n,Y_np1,G,h_i,theta_i,P_hat);
end
% dpiter(iter) = verfahren(s_tr,alpha_n,Y_np1,G,h_i,theta_i,P_hat);
% Update Hilfvariable
if ink_flag == 0
    var1 = sum(theta_i .* h_i);
else % ink_flag == 1
    var1 = 3*G + sum(theta_i .* h_i);
end
% Update Hilfsvariable 2
var2 = s_tr;
for ii = 1:M
    var2 = var2 - theta_i(ii) * alpha_n(:,ii);
end
% var2 = s_tr - sum(theta_i .* alpha_n,2);
% Update effektive spannung
beta_np1 = Y_np1 * var2 / (Y_np1 + var1 * dpiter(iter));
% Update Normale
n_np1 = w3d2 * P_hat * beta_np1./Y_np1;
% Inkrement plastische Dehnung
depsp = w3d2 * dpiter(iter) * n_np1;
% Prädiktor Backstress
alpha_hash = alpha_n;
j_i = r_i./(r_i + mu_i.*h_i*dpiter(iter));
dummy = 2/3 * P_tilde * depsp;
for ii = 1 : M
    alpha_hash(:,ii) = j_i(ii) * ( alpha_hash(:,ii) + h_i(ii) * dummy );
end
% alpha_star = alpha_n + 2/3 * h_i .* (P_tilde * depsp);
% j_i = r_i./(r_i + mu_i.*h_i*dpiter(iter));
% alpha_hash = j_i .* alpha_star;
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
    if ink_flag == 0
        dpiter(iter) = verfahren_sig(s_tr,alpha_n,Y_np1,G,h_i,theta_i,P_hat);
    else
        dpiter(iter) = verfahren_eps(s_tr,alpha_n,Y_np1,G,h_i,theta_i,P_hat);
    end
%     dpiter(iter) = verfahren(s_tr,alpha_n,Y_np1,G,h_i,theta_i,P_hat);
    % Aitkens Delta^2 Verfahren
    if mod(iter,3) == 0
        dpiter(iter) = deltaquadrat(dpiter,iter);
    end
    % Update Hilfvariable
    if ink_flag == 0
        var1 = sum(theta_i .* h_i);
    else % ink_flag == 1
        var1 = 3*G + sum(theta_i .* h_i);
    end
    % Update Hilfsvariable 2
    var2 = s_tr;
    for ii = 1:M
        var2 = var2 - theta_i(ii) * alpha_n(:,ii);
    end
%     var2 = s_tr - sum(theta_i .* alpha_n,2);
    % Update effektive spannung
    beta_np1 = Y_np1 * var2 / (Y_np1 + var1 * dpiter(iter));
    % Update Normale
    n_np1 = w3d2 * P_hat * beta_np1./Y_np1;
    % Inkrement plastische Dehnung
    depsp = w3d2 * dpiter(iter) * n_np1;
    % Prädiktor Backstress
    alpha_hash = alpha_n;
    j_i = r_i./(r_i + mu_i.*h_i*dpiter(iter));
    dummy = 2/3 * P_tilde * depsp;
    for ii = 1 : M
        alpha_hash(:,ii) = j_i(ii) * ( alpha_hash(:,ii) + h_i(ii) * dummy );
    end
%     alpha_star = alpha_n + 2/3 * h_i .* (P_tilde * depsp);
%     j_i = r_i./(r_i + mu_i.*h_i*dpiter(iter));
%     alpha_hash = j_i .* alpha_star;
    % Check Konvergenz
    normp = abs(1-dpiter(iter)/dpiter(iter-1));
    % verhindere endlosschleife
    if iter > maxiter + 1
        msg = 'Keine Konvergenz in Fixpunktiteration nach 101 Iterationen';
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
alpha_np1 = zeros(ntens,M);
for ii = 1:M
    alpha_np1(:,ii) = theta_i_tilde(ii) * alpha_hash(:,ii);
end
% alpha_np1 = theta_i .* alpha_star;
% Inkrement in Spannungen oder Dehnungen
if ink_flag == 0
    eps_np1 = eps_n + D * ink + depsp;
    out = [eps_np1; epsp_np1];
else % ink_flag == 1
    sig_np1 = sig_tr - 2*G * depsp;
    out = [sig_np1; epsp_np1];
end
% beta_np1 = P * sig_np1 - sum(alpha_np1,2);
% F_np1 = sqrt(3/2 * beta_np1' * P_hat * beta_np1) - r0;
% fprintf('F = %.2d\n',F_np1);
% fprintf('n = %i\n',iter-1);
X_neu = [out ; reshape(alpha_np1,ntens*M,1); p_np1];   

end % Ende 3D Materialmodell








