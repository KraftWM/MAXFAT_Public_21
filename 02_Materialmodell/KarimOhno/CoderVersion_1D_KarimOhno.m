%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Matlab Coder Version von:                                            %
%    Modellgleichung 1D Karim Ohno Modell                                 %
%                                                                         %
%    Aufgerufen in:                                                       %
%    karimohno.m                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X_neu = CoderVersion_1D_KarimOhno(eps_n, epsp_n, alpha_n,p_n,...
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