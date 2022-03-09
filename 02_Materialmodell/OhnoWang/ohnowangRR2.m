function [ZVARneu] = ohnowangRR2(ntens, ndi, ink, ZVAR, ink_flag, parameter)
% Materialmodell:
% Aus Ohno et al. 1993 KINEMATIC HARDENING RULES WITH CRITICAL
% STATE OF DYNAMIC RECOVERY, PART I
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
%         ndi       -> Anzahl Diagonalelemente
%         ink       -> Belastungsinkrement
%         ZVAR      -> Zustandsvariablen [eps;epsp;alphai;p] bein spansteu
%                                        [sig;epsp;alphai;p] bein dehnsteu
%         ink_flag  -> (0) Spannungsgesteuerter Prozess 
%                      (1) Dehnungsgesteuerter Prozess
%         parameter -> Materialparameter des Modells
%                      [E,nu,c_i,r_i,chi_i,r0]
%
%
%    OUTPUT:
%        ZVARneu    -> neue zustandsvariablen nach Lastinkrement
%    
%    NOTATIONEN:
%    Spannungen                 Verzerrungen
%
%           sig_11                     eps_11
%    sig =  sig_22               eps = eps_22
%           sig_12                     2eps_12
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
% |  Stand: August 2021                                                 |
%  ----------------------------------------------------------------------


%--------------------------------------------------------------------------
%                  Materialparamter ermitteln                                  
%--------------------------------------------------------------------------
% mynull = 1e-40;
% Anzahl Backstress
M = (length(parameter)-3)/3;
% Elastizitätskonstanten
E = parameter(1);
nu = parameter(2);
G = E/(2*(1+nu));
r0 = parameter(end);
% kinematische Verfestigung
% c_i = parameter(3:2+M);
% r_i = parameter(3+M:2+2*M);
% h_i = c_i.*r_i;
% chi_i = parameter(3+2*M:2+3*M);

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
F_tr = 3/2 * (beta' * (P_line * beta)) - r0^2;                 % Trial Fließfunktion

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
        ZVARneu = ow3D(...
            eps,epsp,alpha,p,...                                   % Zustandsvariablen
            sig_tr,...                                                 % Versuchsspannung
            ink,ink_flag,...                                           % skaliertes Inkrement und Steuerung
            parameter,G,...                                            % Parameter
            D,P,P_line,P_tilde);                                        % Abbildungen
        
    else
        
        % abbildungen
        %             [P, P_line, P_tilde, A, P_check] = set_maps(ntens,ndi);
        % Integration radial return
        ZVARneu = owESZ(...
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
function ZVARneu = ow3D(...
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
M = (length(parameter)-3)/3;                                               % Anzahl TBST

% elastische 

r0 = parameter(end);                                                       % startradius fliessfläche
% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
h_i = c_i.*r_i;
chi_i = parameter(3+2*M:2+3*M);
% oft verwendete Konstanden
w3d2 = sqrt(3/2);
mynull = 1e-40;
EINS4 = diag([1 1 1 1 1 1]);  % Einheitstensor

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
dp_npe = (sqrt(1.5*(beta_hat'*P_hat*beta_hat)) - r0) ...
         /(w3d2*(2*G*delta + sum(h_i.*theta_i)));    
dpiter(iter) = dp_npe;
% ... effektive Spannung & Fließflächen Normale
beta_npe = r0 * beta_hat / ( r0 + w3d2 * dp_npe * (2*G*delta + sum(h_i.*theta_i)));
%     n_npe = w3d2 * (P_hat * beta_npe)/ r0;
n_npe = sqrt(beta_npe' * P_hat * beta_npe); if n_npe == 0; n_npe = mynull; end
n_npe = (P_hat * beta_npe)/n_npe;
% ... Backstress Tensoren
% Normen Backstress
norm_ai = sqrt( sum( (P_hat*alpha_npe) .* alpha_npe) );
norm_ai(norm_ai == 0) = mynull;
Li = alpha_npe./norm_ai;
% Wichtefunktion
var1 = norm_ai./r_i; var1(var1>1) = 1;
var1 = var1.^chi_i;
var2 = n_npe' * Li;
var2 = 0.5 * (var2 + abs(var2));
w_i_npe = var1 .* var2;
% Newtonverfahren für alle Backstresstensoren
alpha_tr = alpha_n + dp_npe * h_i .* repmat(P_tilde * n_npe,1,M);
theta_i = 1./(1+w_i_npe.*c_i*dp_npe);
transn = n_npe';
for i = 1 : M
    
    % Iteration
%     [ai,wi,norma] = NewtonBackstress(i,...                                 % Index Backstress
%         theta_i(i) * alpha_tr(:,i),...                                     % Startwert Newtoniteration
%         alpha_tr(:,i),...                                                  % trial Backstress
%         dp_npe,n_npe',...                                                  % Inkrement plastische Bogenlänge und FFNormale
%         c_i(i),r_i(i),chi_i(i),...                                         % Parameter
%         P_hat,EINS4,EINS4,...                                              % Abblidungen
%         tolbs,maxiter);                                                    % Iterationsoptionen
    % Startwerte (aus Explizitem Euler)
    ai_tr = alpha_tr(:,i);
    wi = w_i_npe(i);
    ci = c_i(i);
    ri = r_i(i);
    chii = chi_i(i);
    ai = 1/( 1 + wi * ci * dp_npe) * ai_tr;
    % Berechne neue Werte
    norma = sqrt(ai' * P_hat * ai); if norma == 0, norma = 1e-40; end
    Li = ai./norma;
    var1 = norma/ri; if var1 > 1, var1 = 1; end
    var1 = var1^chii;
    var2 = transn * Li;
    var2 = 0.5 * (var2 + abs(var2));
    wi= var1 .* var2;
    % Aktueller Fehler
    f = ai * ( 1 + wi * ci * dp_npe) - ai_tr;
    normf = sqrt(f' * P_hat * f);
    iterbs = 1;
    % Newtoniteration
    while normf > tolbs
        % Inkrementiere iter
        iterbs = iterbs + 1;
        % berechne Ableitung
        k1 = 1 + wi * ci * dp_npe;
        if var2 == 0
            k2 = 0;
            k3 = 0;
            var6 = 0;
            var7 = 0;
        else
            k2 = ci*dp_npe * var1 *  var2 * (chii -1);  ...
                k3 = ci*dp_npe * var1;
            var6 = Li' * P_hat * Li;
            var7 = transn * Li;
        end
        invdRdA = 1/k1 * EINS4 - k2/(k1*(k1+k2*var6+k3*var7)) * (Li*(P_hat*Li)') - k3/(k1*(k1+k2*var6+k3*var7)) * (Li * (transn));
%         invdRdA = 1/k1 * EINS4 - k2/(k1*(k1+k2*var6+k3*var7)) * (Li*(P_hat*Li)') - k3/(k1*(k1+k2*var6+k3*var7)) * (n_npe * Li');
        % neuer Backstress
        ai = ai - invdRdA * f;
        % Update alles mögliche
        norma = sqrt(ai' * P_hat * ai); if norma == 0, norma = mynull; end
        Li = ai./norma;
        % Wichtefunktion
        var1 = norma/ri; if var1 > 1, var1 = 1; end
        var1 = var1^chii;
        var2 = transn * Li;
        var2 = 0.5 * (var2 + abs(var2));
        wi = var1 .* var2;
        % berechne Fehler
        f = ai * ( 1 + wi * ci * dp_npe) - ai_tr;
        normf = sqrt(f' * P_hat * f);
        %             normfiter(iterbs) = normf;
        % keine Endlossschleifen
        if iterbs > maxiter
            msg = ['Keine Konvergenz in Backstress ',num2str(bs),' im Ohno Wang Modell ',...
                'nach ', num2str(iterbs),' Iterationen. Aktuelle Fehlernorm: '...
                num2str(normf)];
            warning(msg);
            %                 fprintf('dp = %.10f n = [%.4f %.4f %.4f %.4f %.4f %.4f] ',dp_npe,n_npe)
            %                 figure, grid on, hold on
            %                 plot(normfiter,'.-');
            %                 set(gca,'YScale','log');
            break;
        end
    end % Ende Iterationsschleife
    % Speichern konvergierte Werte
    alpha_npe(:,i) = ai;
    w_i_npe(i) = wi;
    norm_ai(i) = norma;
    theta_i(i) = 1/(1+wi*c_i(i)*dp_npe);
    
end % Ende Schleife Backstress Tensoren
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
%     if mod(iter-1,3) == 0
%         dp_hat = dpiter(iter-1) - ( dpiter(iter-1) - dpiter(iter-2) )^2 / ...
%         ( (dpiter(iter-1) - dpiter(iter-2)) - (dpiter(iter-2)- dpiter(iter-3)));
%         if dp_hat > 0
%             dp_npe = dp_hat;
%         else
%             dp_npe = (sqrt(1.5*(beta_hat'*P_hat*beta_hat)) - r0) ...
%                  /(w3d2*(2*G*delta + sum(h_i.*theta_i)));
%         end
%     else
%         dp_npe = (sqrt(1.5*(beta_hat'*P_hat*beta_hat)) - r0) ...
%                  /(w3d2*(2*G*delta + sum(h_i.*theta_i)));
%     end
%     % Aitkens Delta^2 Verfahren
    dp_npe = (sqrt(1.5*(beta_hat'*P_hat*beta_hat)) - r0) ...
                 /(w3d2*(2*G*delta + sum(h_i.*theta_i)));
    if mod(iter,3) == 0
        dp_hat = dp_npe - ( dp_npe - dpiter(iter-1) )^2 / ...
        ( dp_npe - 2*dpiter(iter-1) + dpiter(iter-2) );
        if dp_hat > 0
            dp_npe = dp_hat;
        end
    end

    dpiter(iter) = dp_npe;
    % ... effektive Spannung & Fließflächen Normale
    beta_npe = r0 * beta_hat / ( r0 + w3d2 * dp_npe * (2*G*delta + sum(h_i.*theta_i)));
%     n_npe = w3d2 * (P_hat * beta_npe)/ r0;
    n_npe = sqrt(beta_npe' * P_hat * beta_npe); if n_npe == 0; n_npe = mynull; end
    n_npe = (P_hat * beta_npe)/n_npe;
%     fprintf('||n|| = %.15f\n',n_npe' * P_tilde * n_npe)
    % ... Backstress Tensoren
    % Normen Backstress
    norm_ai = sqrt( sum( (P_hat*alpha_npe) .* alpha_npe) );
    norm_ai(norm_ai == 0) = mynull;
    Li = alpha_npe./norm_ai;
    % Wichtefunktion
    var1 = norm_ai./r_i; var1(var1>1) = 1;
    var1 = var1.^chi_i;
    var2 = n_npe' * Li;
    var2 = 0.5 * (var2 + abs(var2));
    w_i_npe = var1 .* var2;
    % Newtonverfahren für alle Backstresstensoren
    alpha_tr = alpha_n + dp_npe * h_i .* repmat(P_tilde * n_npe,1,M);
    theta_i = 1./(1+w_i_npe.*c_i*dp_npe);
    transn = n_npe';
    for i = 1 : M
        
        % Iteration
        %     [ai,wi,norma] = NewtonBackstress(i,...                                 % Index Backstress
        %         theta_i(i) * alpha_tr(:,i),...                                     % Startwert Newtoniteration
        %         alpha_tr(:,i),...                                                  % trial Backstress
        %         dp_npe,n_npe',...                                                  % Inkrement plastische Bogenlänge und FFNormale
        %         c_i(i),r_i(i),chi_i(i),...                                         % Parameter
        %         P_hat,EINS4,EINS4,...                                              % Abblidungen
        %         tolbs,maxiter);                                                    % Iterationsoptionen
        % Startwerte (aus Explizitem Euler)
        ai_tr = alpha_tr(:,i);
        wi = w_i_npe(i);
        ci = c_i(i);
        ri = r_i(i);
        chii = chi_i(i);
        ai = 1/( 1 + wi * ci * dp_npe) * ai_tr;
        % Berechne neue Werte
        norma = sqrt(ai' * P_hat * ai); if norma == 0, norma = 1e-40; end
        Li = ai./norma;
        var1 = norma/ri; if var1 > 1, var1 = 1; end
        var1 = var1^chii;
        var2 = transn * Li;
        var2 = 0.5 * (var2 + abs(var2));
        wi= var1 .* var2;
        % Aktueller Fehler
        f = ai * ( 1 + wi * ci * dp_npe) - ai_tr;
        normf = sqrt(f' * P_hat * f);
        iterbs = 1;
        % Newtoniteration
        while normf > tolbs
            % Inkrementiere iter
            iterbs = iterbs + 1;
            % berechne Ableitung
            k1 = 1 + wi * ci * dp_npe;
            if var2 == 0
                k2 = 0;
                k3 = 0;
                var6 = 0;
                var7 = 0;
            else
                k2 = ci*dp_npe * var1 *  var2 * (chii -1);  ...
                    k3 = ci*dp_npe * var1;
                var6 = Li' * P_hat * Li;
                var7 = transn * Li;
            end
            invdRdA = 1/k1 * EINS4 - k2/(k1*(k1+k2*var6+k3*var7)) * (Li*(P_hat*Li)') - k3/(k1*(k1+k2*var6+k3*var7)) * (Li * (transn));
            % neuer Backstress
            ai = ai - invdRdA * f;
            % Update alles mögliche
            norma = sqrt(ai' * P_hat * ai); if norma == 0, norma = mynull; end
            Li = ai./norma;
            % Wichtefunktion
            var1 = norma/ri; if var1 > 1, var1 = 1; end
            var1 = var1^chii;
            var2 = transn * Li;
            var2 = 0.5 * (var2 + abs(var2));
            wi = var1 .* var2;
            % berechne Fehler
            f = ai * ( 1 + wi * ci * dp_npe) - ai_tr;
            normf = sqrt(f' * P_hat * f);
            %             normfiter(iterbs) = normf;
            % keine Endlossschleifen
            if iterbs > maxiter
                msg = ['Keine Konvergenz in Backstress ',num2str(i),' im Ohno Wang Modell ',...
                    'nach ', num2str(iterbs),' Iterationen. Aktuelle Fehlernorm: '...
                    num2str(normf)];
                warning(msg);
                %                 fprintf('dp = %.10f n = [%.4f %.4f %.4f %.4f %.4f %.4f] ',dp_npe,n_npe)
                %                 figure, grid on, hold on
                %                 plot(normfiter,'.-');
                %                 set(gca,'YScale','log');
                break;
            end
        end % Ende Iterationsschleife
        % Speichern konvergierte Werte
        alpha_npe(:,i) = ai;
        w_i_npe(i) = wi;
        norm_ai(i) = norma;
        theta_i(i) = 1/(1+wi*c_i(i)*dp_npe);
        
    end % Ende Schleife Backstress Tensoren
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
% figure(200), hold on
% semilogy(1:iter-1,normpiter(2:iter),'k--s','LineWidth',1) 
% axis([1e0 100 1e-6 max(normpiter)])
% figure(3)
% plot(dpiter(1:iter)), pause(0.5)


% -------------------------------------------------------------------------
%                        update Zustandsvariablen
% ------------------------------------------------------------------------- 
depsp = dp_npe * n_npe;
epsp_npe = epsp_n + depsp;
p_npe = p_n + dp_npe;
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
function ZVARneu = owESZ(...
                       sig_n,eps_n,epsp_n,alpha_n,p_n,...                  % Zustandsvariablen
                       sig_tr,...                                          % Versuchsspannung
                       ink,ink_flag,...                                    % Inkrement und Steuerung
                       parameter,G,...                                     % Parameter
                       D,P,P_line,P_hat,A,P_check)                       % Abblidungen
                   
% -------------------------------------------------------------------------
%                        Materialparameter
% -------------------------------------------------------------------------
% Zuweisen der Materialparameter
ntens = 3;                                                                 % Tensorkomponenten
M = (length(parameter)-3)/3;                                               % Anzahl TBST

% elastische 
% E = parameter(1);
nu = parameter(2);
r0 = parameter(end);                                                       % startradius fliessfläche
% K = E/(3*(1-2*nu));
% kinematische Verfestigung
c_i = parameter(3:2+M);
r_i = parameter(3+M:2+2*M);
h_i = c_i.*r_i;
chi_i = parameter(3+2*M:2+3*M);
% oft verwendete Konstanden
w3d2 = sqrt(3/2);
mynull = 1e-40;
EINS = diag([1 1 1]);    % Einheitstensor
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
% ... trial Inkrement
dp_npe = 0;
% ... Korrekturfaktoren
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
          - delta * 2*G*fak*dp_npe*(-n_npe(1) - n_npe(2))*DELTA3;
beta_hat = s_tr - alpha_hat;
dp_npe = (sqrt(1.5*(beta_hat'*P_line*beta_hat)) - r0) ...
         /(w3d2*(2*G*delta + sum(h_i.*theta_i)));
dpiter(iter) = dp_npe;
% ... effektive Spannung & Fließflächen Normale
beta_npe = r0 * beta_hat / ( r0 + w3d2 * dp_npe * (2*G*delta + sum(h_i.*theta_i)));
%     n_npe = w3d2 * (P_hat * beta_npe)/ r0;
n_npe = sqrt(beta_npe' * P_line * beta_npe); if n_npe == 0; n_npe = mynull; end
n_npe = (P_hat * beta_npe)/n_npe;
% ... Backstress Tensoren
% Normen Backstress
norm_ai = sqrt( sum( (P_line*alpha_npe) .* alpha_npe) );
norm_ai(norm_ai == 0) = mynull;
Li = alpha_npe./norm_ai;
% Wichtefunktion
var1 = norm_ai./r_i; var1(var1>1) = 1;
var1 = var1.^chi_i;
var2 = n_npe' * P_check * Li;
var2 = 0.5 * (var2 + abs(var2));
w_i_npe = var1 .* var2;
% oft verwendeter Kram
alpha_tr = alpha_n + dp_npe * h_i .* repmat(A * n_npe,1,M);            % trial Backstress
theta_i = 1./(1+w_i_npe.*c_i*dp_npe);
transn = n_npe';                                                       % Transponierte Normale
% Schleife über Backstress
for i = 1 : M
    
    % aktueller backstress
    %         ai = alpha_npe(:,i);
    ai_tr = alpha_tr(:,i);
    wi = w_i_npe(i);
    ci = c_i(i);
    ri = r_i(i);
    chii = chi_i(i);
    ai = 1/( 1 + wi * ci * dp_npe) * ai_tr;
    norma = sqrt(ai' * P_line * ai); if norma == 0, norma = mynull; end
    Li = ai./norma;
    var1 = norma/ri; if var1 > 1, var1 = 1; end
    var1 = var1^chii;
    var2 = transn * P_check * Li;
    var2 = 0.5 * (var2 + abs(var2));
    wi= var1 .* var2;
    f = ai * ( 1 + wi * ci * dp_npe) - ai_tr;
    normf = sqrt(f' * P_line * f);
    %         normfiter = zeros(1,maxiter+1);
    %         normf(1) = normf;
    %         fprintf('i= %i ||f0|| = %.10f ',i,normf);
    % Newtonverfahren
    iterbs = 1;
    tolbs = 1e-6;
    % Iterationsschleife newtonverfahren
    while normf > tolbs
        % Inkrementiere iter
        iterbs = iterbs + 1;
        % berechne Ableitung
        k1 = 1 + wi * ci * dp_npe;
        if var2 == 0
            k2 = 0;
            k3 = 0;
            var6 = 0;
            var7 = 0;
        else
            k2 = ci*dp_npe * var1 *  var2 * (chii -1);  ...
                k3 = ci*dp_npe * var1;
            var6 = Li'*P_line*Li;
            var7 = n_npe'*P_check*Li;
        end
        %             dRdA = k1 * EINS4 + k2 * (Li*Li') + k3 * (Li * (A*n_npe)');
        invdRdA = 1/k1 * EINS - k2/(k1*(k1+k2*var6+k3*var7)) * (Li*(P_line*Li)') - k3/(k1*(k1+k2*var6+k3*var7)) * (Li * (P_check*n_npe)');
        % neuer Backstress
        %             ai = ai - dRdA\f;
        ai = ai - invdRdA * f;
        % Update alles mögliche
        norma = sqrt(ai' * P_line * ai); if norma == 0, norma = mynull; end
        Li = ai./norma;
        % Wichtefunktion
        var1 = norma/ri; if var1 > 1, var1 = 1; end
        var1 = var1^chii;
        var2 = transn * P_check * Li;
        var2 = 0.5 * (var2 + abs(var2));
        wi= var1 .* var2;
        % berechne Fehler
        f = ai * ( 1 + wi * ci * dp_npe) - ai_tr;
        normf = sqrt(f' * P_line * f);
        %             normfiter(iterbs) = normf;
        % keine Endlossschleifen
        if iterbs > 10*maxiter
            msg = ['Keine Konvergenz in Backstress',num2str(i),' im Ohno Wang Modell ',...
                'nach ', num2str(iterbs),' Iterationen. Aktuelle Fehlernorm: '...
                num2str(normf)];
            warning(msg);
            %                 fprintf('dp = %.10f n = [%.4f %.4f %.4f %.4f %.4f %.4f] ',dp_npe,n_npe)
            %                 figure, grid on, hold on
            %                 plot(normfiter,'.-');
            %                 set(gca,'YScale','log');
            break;
        end
    end % Ende Newtonverfahren
    %         fprintf('iterbs: %i  ||f|| = %.10f\n',iterbs,normf);
    % Speichern konvergierte Ergebnisse
    norm_ai(i) = norma;
    w_i_npe(i) = wi;
    alpha_npe(:,i) = ai;
    theta_i(i) = 1/(1+wi*c_i(i)*dp_npe);
    
end % Ende Schleife Über Backstress

% -------------------------------------------------------------------------
%                  Fixpunktiteration
% ------------------------------------------------------------------------- 
% ... Schleife
while normp > tol
    % ... Inkrementieren Schleifenzähler
    iter = iter + 1;
    % ... Inkrement plastische Bogenlänge
    alpha_hat = sum(theta_i.*alpha_n,2) ...
        - delta *2*G*fak*dp_npe*(-n_npe(1) - n_npe(2))*DELTA3;
    beta_hat = s_tr - alpha_hat;
    % Aitkens Delta^2 Verfahren
%     if mod(iter-1,3) == 0
%         dp_hat = dpiter(iter-1) - ( dpiter(iter-1) - dpiter(iter-2) )^2 / ...
%         ( (dpiter(iter-1) - dpiter(iter-2)) - (dpiter(iter-2)- dpiter(iter-3)));
%         if dp_hat > 0
%             dp_npe = dp_hat;
%         else
%             dp_npe = (sqrt(1.5*(beta_hat'*P_line*beta_hat)) - r0) ...
%                      /(w3d2*(2*G*delta + sum(h_i.*theta_i)));
%         end
%     else
%         dp_npe = (sqrt(1.5*(beta_hat'*P_line*beta_hat)) - r0) ...
%         /(w3d2*(2*G*delta + sum(h_i.*theta_i)));
%     end
    % Aitkens Delta^2 Verfahren
    dp_npe = (sqrt(1.5*(beta_hat'*P_line*beta_hat)) - r0) ...
        /(w3d2*(2*G*delta + sum(h_i.*theta_i)));
    if mod(iter,3) == 0
        dp_hat = dp_npe - ( dp_npe - dpiter(iter-1) )^2 / ...
        ( dp_npe - 2*dpiter(iter-1) + dpiter(iter-2) );
        if dp_hat > 0
            dp_npe = dp_hat;
        end
    end
    dpiter(iter) = dp_npe;
    % ... effektive Spannung & Fließflächen Normale
    beta_npe = r0 * beta_hat / ( r0 + w3d2 * dp_npe * (2*G*delta + sum(h_i.*theta_i)));
%     n_npe = w3d2 * (P_hat * beta_npe)/ r0;
    n_npe = sqrt(beta_npe' * P_line * beta_npe); if n_npe == 0; n_npe = mynull; end
    n_npe = (P_hat * beta_npe)/n_npe;
%     fprintf('||n|| = %.15f\n',n_npe' * A * P_line * A * n_npe)
    % ... Backstress Tensoren
    % oft verwendeter Kram
    alpha_tr = alpha_n + dp_npe * h_i .* repmat(A * n_npe,1,M);            % trial Backstress
%     theta_i = 1./(1+w_i_npe.*c_i*dp_npe);
    transn = n_npe';                                                       % Transponierte Normale
    % Schleife über Backstress
    for i = 1 : M
        
        % aktueller backstress
        %         ai = alpha_npe(:,i);
        ai_tr = alpha_tr(:,i);
        wi = w_i_npe(i);
        ci = c_i(i);
        ri = r_i(i);
        chii = chi_i(i);
        ai = 1/( 1 + wi * ci * dp_npe) * ai_tr;
        norma = sqrt(ai' * P_line * ai); if norma == 0, norma = mynull; end
        Li = ai./norma;
        var1 = norma/ri; if var1 > 1, var1 = 1; end
        var1 = var1^chii;
        var2 = transn * P_check * Li;
        var2 = 0.5 * (var2 + abs(var2));
        wi= var1 .* var2;
        f = ai * ( 1 + wi * ci * dp_npe) - ai_tr;
        normf = sqrt(f' * P_line * f);
        %         normfiter = zeros(1,maxiter+1);
        %         normf(1) = normf;
        %         fprintf('i= %i ||f0|| = %.10f ',i,normf);
        % Newtonverfahren
        iterbs = 1;
        tolbs = 1e-6;
        % Iterationsschleife newtonverfahren
        while normf > tolbs
            % Inkrementiere iter
            iterbs = iterbs + 1;
            % berechne Ableitung
            k1 = 1 + wi * ci * dp_npe;
            if var2 == 0
                k2 = 0;
                k3 = 0;
                var6 = 0;
                var7 = 0;
            else
                k2 = ci*dp_npe * var1 *  var2 * (chii -1);  ...
                    k3 = ci*dp_npe * var1;
                var6 = Li'*P_line*Li;
                var7 = transn*P_check*Li;
            end
            %             dRdA = k1 * EINS4 + k2 * (Li*Li') + k3 * (Li * (A*n_npe)');
            invdRdA = 1/k1 * EINS - k2/(k1*(k1+k2*var6+k3*var7)) * (Li*(P_line*Li)') - k3/(k1*(k1+k2*var6+k3*var7)) * (Li * (transn*P_check));
            % neuer Backstress
            %             ai = ai - dRdA\f;
            ai = ai - invdRdA * f;
            % Update alles mögliche
            norma = sqrt(ai' * P_line * ai); if norma == 0, norma = mynull; end
            Li = ai./norma;
            % Wichtefunktion
            var1 = norma/ri; if var1 > 1, var1 = 1; end
            var1 = var1^chii;
            var2 = transn * P_check * Li;
            var2 = 0.5 * (var2 + abs(var2));
            wi= var1 .* var2;
            % berechne Fehler
            f = ai * ( 1 + wi * ci * dp_npe) - ai_tr;
            normf = sqrt(f' * P_line * f);
            %             normfiter(iterbs) = normf;
            % keine Endlossschleifen
            if iterbs > 10*maxiter
                msg = ['Keine Konvergenz in Backstress',num2str(i),' im Ohno Wang Modell ',...
                    'nach ', num2str(iterbs),' Iterationen. Aktuelle Fehlernorm: '...
                    num2str(normf)];
                warning(msg);
                %                 fprintf('dp = %.10f n = [%.4f %.4f %.4f %.4f %.4f %.4f] ',dp_npe,n_npe)
                %                 figure, grid on, hold on
                %                 plot(normfiter,'.-');
                %                 set(gca,'YScale','log');
                break;
            end
        end % Ende Newtonverfahren
%         fprintf('                    iterbs: %i  ||f|| = %.10f\n',iterbs,normf);
        % Speichern konvergierte Ergebnisse
        norm_ai(i) = norma;
        w_i_npe(i) = wi;
        alpha_npe(:,i) = ai;
        theta_i(i) = 1/(1+wi*c_i(i)*dp_npe);
    end % Ende Schleife Über Backstress
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
depsp = dp_npe * n_npe;
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




% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
% Newtonverfahren Backstress
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
function [ai,wi,norma] = NewtonBackstress(bs,...                           % Index Backstress
    ai,...                                      % Startwert Newtoniteration
    ai_tr,...                                   % trial Backstress
    dp,transn,...                       % Inkrement plastische Bogenlänge und FFNormale
    ci,ri,chii,...                              % Parameter
    P_line,P_check,EINS,...                     % Abblidungen
    tol,maxiter)                          % Iterationsoptionen

% 1. Upadate Versuch
%         ai = 1/( 1 + wi * ci * dp) * ai_tr;
% Berechne neue Werte
norma = sqrt(ai' * P_line * ai); if norma == 0, norma = 1e-40; end
Li = ai./norma;
var1 = norma/ri; if var1 > 1, var1 = 1; end
var1 = var1^chii;
var2 = transn * P_check * Li;
var2 = 0.5 * (var2 + abs(var2));
wi= var1 .* var2;
% Aktueller Fehler
f = ai * ( 1 + wi * ci * dp) - ai_tr;
normf = sqrt(f' * P_line * f);
iterbs = 1;
% Newtoniteration
while normf > tol
    % Inkrementiere iter
    iterbs = iterbs + 1;
    % berechne Ableitung
    k1 = 1 + wi * ci * dp;
    if var2 == 0
        k2 = 0;
        k3 = 0;
        var6 = 0;
        var7 = 0;
    else
        k2 = ci*dp * var1 *  var2 * (chii -1);  ...
            k3 = ci*dp * var1;
        var6 = Li' * P_line * Li;
        var7 = transn * P_check * Li;
    end
    invdRdA = 1/k1 * EINS - k2/(k1*(k1+k2*var6+k3*var7)) * (Li*(P_line*Li)') - k3/(k1*(k1+k2*var6+k3*var7)) * (Li * (transn*P_check));
    % neuer Backstress
    ai = ai - invdRdA * f;
    % Update alles mögliche
    norma = sqrt(ai' * P_line * ai); if norma == 0, norma = mynull; end
    Li = ai./norma;
    % Wichtefunktion
    var1 = norma/ri; if var1 > 1, var1 = 1; end
    var1 = var1^chii;
    var2 = transn * P_check * Li;
    var2 = 0.5 * (var2 + abs(var2));
    wi = var1 .* var2;
    % berechne Fehler
    f = ai * ( 1 + wi * ci * dp) - ai_tr;
    normf = sqrt(f' * P_line * f);
    %             normfiter(iterbs) = normf;
    % keine Endlossschleifen
    if iterbs > maxiter
        msg = ['Keine Konvergenz in Backstress ',num2str(bs),' im Ohno Wang Modell ',...
            'nach ', num2str(iterbs),' Iterationen. Aktuelle Fehlernorm: '...
            num2str(normf)];
        warning(msg);
        %                 fprintf('dp = %.10f n = [%.4f %.4f %.4f %.4f %.4f %.4f] ',dp_npe,n_npe)
        %                 figure, grid on, hold on
        %                 plot(normfiter,'.-');
        %                 set(gca,'YScale','log');
        break;
    end
end % Ende Newtoniteration
end % Ende Newton Verfahren Backstress