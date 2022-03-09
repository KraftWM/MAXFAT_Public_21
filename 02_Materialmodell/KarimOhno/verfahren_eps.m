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
    M = length(theta_i);
    ntens = size(s_tr,1);
    a_n = zeros(ntens,1);
    for ii = 1:M
    a_n = a_n + theta_i(ii) * alpha_n(:,ii);
    end
%     a_n = sum(theta_i .* alpha_n,2);
    beta_tr = s_tr - a_n;
    % Versuchsfließspannunge
    F_tr = sqrt(3/2 * beta_tr' * P_line * beta_tr) - Y_np1;
    % Hilfsvariable 
    var1 = sum(h_i .* theta_i);
    % neues plastisches Dehnungsinkrement
    dp = F_tr/(3*G+var1);
end % Ende Verfahrensfunktion dehnungssteuerung













