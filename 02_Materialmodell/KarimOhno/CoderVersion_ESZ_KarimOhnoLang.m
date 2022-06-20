%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Matlab Coder Version von:                                            %
%    Modellgleichung ESZ Karim Ohno Modell für Kerbnäherungsverfahren     %
%    nach Lang                                                            %
%                                                                         %
%    Aufgerufen in:                                                       %
%    karimohnolang.m                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZVARneu = CoderVersion_ESZ_KarimOhnoLang(esig_n, eeps_n, eepsp_n, ealpha_n,...
                                  alpha_n,p_n,...
                                  esig,...
                                  parameter,...
                                  D,P,P_line,P_hat,A)
% Konkretes Materialmodell für ESZ Spannungszustände, gibt bei vorgabe 
% eines Lastinkrementes die geupdatetn Zustandsvariablen mit implizitem 
% Euler zurück .
%
%   INPUT:
%       esig_n       -> Pseudo Spannungen
%  eeps_n,eepsp_n    -> Pseudo Dehnungen
%       ealpha_n     -> Pseudo Backstresstensoren
%        alpha_n     -> Backstress
%       esig        -> Spannung bei np1
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
M = (length(parameter)-4)/6;                                               % Anzahl TBST

% elastische 
E = parameter(1);
nu = parameter(2);
G = E/(2*(1+nu));
r0 = parameter(3+3*M);                                                     % startradius fliessfläche
er0 = parameter(end);                                                     % startradius fliessfläche
% kinematische Verfestigung Materialmodell
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
mu_i = parameter(3+2*M:2+3*M);
h_i = c_i .* r_i;
% kinematische Verfestigung Struckturmodell
ec_i = parameter(4+3*M:3+4*M);
er_i = parameter(4+4*M:3+5*M);
emu_i = parameter(4+5*M:3+6*M);
eh_i = ec_i .* er_i;
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
etheta_i = ones(1,M);
theta_i = ones(1,M);
etheta_i_tilde = ones(1,M);
% deviatorischer Anteilspannung
es = P * esig;
% neuer radius FF = alter radius FF
Y_np1 = er0;
% erster schritt
dpiter(iter) = verfahren_sig(es,ealpha_n,Y_np1,G,eh_i,etheta_i,P_line);
% Update Hilfvariable
var1 = sum(etheta_i .* eh_i);
% Update Hilfsvariable 2
var2 = es;
for ii = 1 : M
    var2 = var2 - etheta_i(ii) * ealpha_n(:,ii);
end
% var2 = es - sum(etheta_i .* ealpha_n,2);
% Update effektive spannung
beta_np1 = Y_np1 * var2 / (Y_np1 + var1 * dpiter(iter));
% Update Normale
n_np1 = w3d2 * P_hat * beta_np1./Y_np1;
% Inkrement plastische Dehnung
depsp = w3d2 * dpiter(iter) * n_np1;
% Prädiktor pseudo Backstress
ealpha_hash = ealpha_n;
ej_i = er_i./(er_i + emu_i.*eh_i*dpiter(iter));
dummy = 2/3 * A * depsp;
for ii = 1 : M
    ealpha_hash(:,ii) = ej_i(ii) * ( ealpha_hash(:,ii) + eh_i(ii) * dummy );
end
% ealpha_star = ealpha_n + 2/3 * eh_i .* (A * depsp);
% ej_i = er_i./(er_i + emu_i.*eh_i*dpiter(iter));
% ealpha_hash = ej_i .* ealpha_star;
% Iterationsschleife
while normp > tol
    % effektive Backstresstensoren
    ealpha_eff = sqrt( 3/2 * sum( (P_line*ealpha_hash) .* ealpha_hash) );
    idx1 = ealpha_eff <= er_i;
    idx2 = ealpha_eff > er_i;
    etheta_i_tilde( idx1 ) = 1;
    etheta_i_tilde( idx2 ) = er_i(idx2)./ealpha_eff(idx2);
    etheta_i = ej_i .* etheta_i_tilde;
    % Inkrementieren Schleifenzähler
    iter = iter + 1;
    % neues plastisches Dehnungsinkrement
    dpiter(iter) = verfahren_sig(es,ealpha_n,Y_np1,G,eh_i,etheta_i,P_line);
    % Aitkens Delta^2 Verfahren
    if mod(iter,3) == 0
        dpiter(iter) = deltaquadrat(dpiter,iter);
    end
    % Update Hilfvariable
    var1 = sum(etheta_i .* eh_i);
    % Update Hilfsvariable 2
    var2 = es;
    for ii = 1 : M
        var2 = var2 - etheta_i(ii) * ealpha_n(:,ii);
    end
%     var2 = es - sum(etheta_i .* ealpha_n,2);
    % Update effektive spannung
    beta_np1 = Y_np1 * var2 / (Y_np1 + var1 * dpiter(iter));
    % Update Normale
    n_np1 = w3d2 * P_hat * beta_np1./Y_np1;
    % Inkrement plastische Dehnung
    depsp = w3d2 * dpiter(iter) * n_np1;
    % Prädiktor Backstress
    ealpha_hash = ealpha_n;
    ej_i = er_i./(er_i + emu_i.*eh_i*dpiter(iter));
    dummy = 2/3 * A * depsp;
    for ii = 1 : M
        ealpha_hash(:,ii) = ej_i(ii) * ( ealpha_hash(:,ii) + eh_i(ii) * dummy );
    end
%     ealpha_star = ealpha_n + 2/3 * eh_i .* (A * depsp);
%     ej_i = er_i./(er_i + emu_i.*eh_i*dpiter(iter));
%     ealpha_hash = ej_i .* ealpha_star;
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

% Backstress radial return
alpha_np1 = zeros(ntens,M);
alpha_star = alpha_n;
j_i = r_i./(r_i + mu_i.*h_i*dpiter(iter));
dummy = 2/3 * A * depsp;
for ii = 1 : M 
    alpha_star(:,ii) = alpha_star(:,ii) + h_i(ii) * dummy;
end
alpha_eff = sqrt( 3/2 * sum( (P_line*alpha_star) .* alpha_star) ); 
idx1 = alpha_eff <= r_i;
idx2 = alpha_eff > r_i;
theta_i( idx1 ) = j_i( idx1 );
theta_i( idx2 ) = r_i( idx2 ) ./ alpha_eff( idx2 );
for ii = 1 : M
    alpha_np1(:,ii) = theta_i(ii) * alpha_star(:,ii);
end
% 
% alpha_star = alpha_n + 2/3 * h_i .* (A * depsp);
% j_i = r_i./(r_i + mu_i.*h_i*dpiter(iter));
% alpha_eff = sqrt( 3/2 * sum( (P_line*alpha_star) .* alpha_star) ); 
% idx1 = alpha_eff <= r_i;
% idx2 = alpha_eff > r_i;
% theta_i( idx1 ) = j_i( idx1 );
% theta_i( idx2 ) = r_i( idx2 ) ./ alpha_eff( idx2 );
% alpha_np1 = theta_i .* alpha_star;
% Plastische Dehnungen
eepsp_np1 = eepsp_n + depsp;
% plastische Bogenlänge
p_np1 = p_n + dpiter(iter);
% Backstresstensoren
ealpha_np1 = zeros(ntens,M);
for ii = 1:M
    ealpha_np1(:,ii) = etheta_i_tilde(ii) * ealpha_hash(:,ii);
end 
% ealpha_np1 = etheta_i .* ealpha_star;
% Inkrement in Spannungen oder Dehnungen
eeps_np1 = eeps_n + D * (esig- esig_n) + depsp;
out = [eeps_np1; eepsp_np1];
% Zusammenfassen der Zustandsvarablen
% F_np1 = sqrt(3/2 * beta_np1' * P_line * beta_np1) - r0;
% fprintf('F = %.2d\n',F_np1);
% fprintf('n = %i\n',iter-1);
ZVARneu = [out; reshape(ealpha_np1,ntens*M,1); p_np1;...
                reshape(alpha_np1,ntens*M,1)];   

end % Ende ESZ Materialmodell

