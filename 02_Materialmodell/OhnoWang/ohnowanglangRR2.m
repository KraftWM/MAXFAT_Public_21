function [ZVARneu] = ohnowanglangRR2(desig, ZVAR, para, epara)
% Ohno Wang Plastizitätsmodell für pseudo Stress approach nach Asatz
% von Lang, implementierung nach Döring
%
% ! Nur Spannungssteuerung im ESZ
%
% QUELLE:
%       Ohno et al. 1993 -  KINEMATIC HARDENING RULES WITH CRITICAL
%                           STATE OF DYNAMIC RECOVERY, PART I
%       Lang et al       - A multiaxial stress-strain correction scheme
%
% Integration mit Radial Return Methode nach eigenem Algorithmus
%
% INPUT:
%  dESIG -> Inkrement in pseudo spannungen
%  ZVAR  -> Zustandsvariablen
%  para  -> Parameter des Pseudo Modells und des Materialmodells
%
% OUTPUT:
%  ZVARneu -> neuer zustand nach inkrement
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
% Darstellung von Tensoren
%         sig11              eps11
%  sig =  sig22       eps =  eps22
%         sig12             2eps12
%__________________________________________________________________________
%
% Parameter:
% zuerst Material- dann Struckturmodell
% Parameter:
% zuerst Material- dann Struckturmodell 
%     para = [E, nu,  c_i,  r_i,  chi_i,  r0]
%       M = (length(para)-3)/3;
%     epara = [E, nu,  ec_i, er_i, echi_i, er0]
%       eM = (length(epara)-3)/3; 
%
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik            |
% |  Stand: August 2021                                                 |
%  ----------------------------------------------------------------------
ntens = 3;
ndi = 2;


%--------------------------------------------------------------------------
%                  Materialparamter ermitteln
%--------------------------------------------------------------------------
M = (length(para)-3)/3;                                               % Anzahl Backstresstensoren
eM = (length(epara)-3)/3;                                               % Anzahl Backstresstensoren

E = para(1);                                                          % Elastizitätsmodul
nu = para(2);                                                         % Querdehnzahl
er0 = epara(3*eM+3);                                                    % Radius pseudo Fließfläche
% G = E/(2*(1+nu));                                                        % Schubmodul


%--------------------------------------------------------------------------
%                        Hilfsmatrizen
%--------------------------------------------------------------------------

% Steifigkeit
CEL = elast_steifigkeit(E,nu,ntens,ndi);
DEL = elast_nachgiebigkeit(E,nu,ntens,ndi);
% statische Matrizen
[P, P_line] = set_maps(ntens,ndi);

%--------------------------------------------------------------------------
%            Identifikation der Zustandsvariablen
%--------------------------------------------------------------------------

% pseudo gesamtdehnung
eeps =  ZVAR(1:ntens);
% plastische dehnung
epsp = ZVAR(ntens+1:2*ntens);
% Pseudo (elastische9 Spannungen
esig = CEL * (eeps - epsp);
% pseudo backstress
ealphai = reshape(ZVAR(2*ntens+1:(eM+2)*ntens),ntens,eM);
% Backstress
alphai = reshape(ZVAR((eM+2)*ntens+2:(eM+M+2)*ntens+1),ntens,M);
% plastische bogenlänge
p = ZVAR((eM+2)*ntens+1);

%--------------------------------------------------------------------------
%               elastischer Trial Step
%--------------------------------------------------------------------------

es = P * esig;                                                             % Pseudo Spannungsdeviator
ea = sum(ealphai,2);                                                       % Pseudo Gesamtbackstress
des = P * desig;                                                           % Inkrement Pseudo Spannungsdeviator
es_tr = es + des;                                                          % Pseudo Versuchsspannungsdeviator
ebeta = es_tr - ea;                                                        % Pseudo effektive Spannung
F_tr = ebeta' * P_line * ebeta - 2/3 * er0^2;                              % Pseudo Überspannung
FTOL = 1e-7;
%--------------------------------------------------------------------------
%                   Trial Step wird angenommen
%--------------------------------------------------------------------------

if F_tr < FTOL   % Trial step völlig elastische
    
    % Updaten der pseudo gesamtdehnung
    eeps = eeps + DEL * desig;
    % Updaten Zustandsvariablen
    ZVARneu = [eeps;ZVAR(ntens+1:end)];
    
    
    %--------------------------------------------------------------------------
    %                   Trial Step wird abgelehnt
    %--------------------------------------------------------------------------
    
else
    
    % abbildungen
    [P, P_line, P_hat, A, P_check] = set_maps(ntens,ndi);
    % Integration radial return
    ZVARneu = owlangESZ(...
        esig,eeps,epsp,ealphai,alphai,p,...                                % Zustandsvariablen
        esig+desig,...                                                  % Versuchsspannung
        para,epara,...                                                      % Parameter
        DEL,P,P_line,P_hat,A,P_check);                                     % Abbildungen
    
end

end % Ende Hauptfunktion










% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
% Materialmodell ESZ
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
function ZVARneu = owlangESZ(...
    esig_n,eeps_n,epsp_n,ealpha_n,alpha_n,p_n,...       % Zustandsvariablen
    esig_tr,...                                          % Versuchsspannung
    para,epara,...                                       % Parameter
    D,P,P_line,P_hat,A,P_check)                         % Abblidungen

% -------------------------------------------------------------------------
%    Parameter
% -------------------------------------------------------------------------
% Spannungszustand
ntens = 3;
% ndi = 2;
% Anzahl Backstresstensoren
M = (length(para)-3)/3;  
eM = (length(epara)-3)/3;                                                                                                               % Querdehnzahl

% Material
c_i = para(3:2+M);
r_i = para(3+M:2+2*M);
h_i = c_i .* r_i;
chi_i = para(3+2*M:2+3*M);
% Strucktur
ec_i = epara(3:2+eM);
er_i = epara(3+eM:2+2*eM);
eh_i = ec_i .* er_i;
echi_i = epara(3+2*eM:2+3*eM);
er0 = epara(3+3*eM);
% oft verwendete Konstanden
w3d2 = sqrt(3/2);
mynull = 1e-40;
EINS = diag([1 1 1]);    % Einheitstensor
% DELTA3 = [-1/3;-1/3;0];  % Hilfstensor
% fak = (1-2*nu)/(1-nu);

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
% dp_npe = 0;
% ... Korrekturfaktoren
etheta_i = ones(1,eM);
% ... Backstress
ealpha_npe = ealpha_n;
% ... trial Spannungsdeviator
es_tr = P * esig_tr;
% ... trial Fließflächennormale
% n_npe = P * esig_n - sum(ealpha_npe,2);
% var1 = sqrt(n_npe' * P_line * n_npe); if var1 == 0, var1 = mynull; end
% n_npe = (P_hat * n_npe)./var1; 

% -------------------------------------------------------------------------
%                   erster Iterationsschritt
% ------------------------------------------------------------------------- 
% ... Inkrement plastische Bogenlänge
ealpha_hat = sum(etheta_i.*ealpha_n,2);
ebeta_hat = es_tr - ealpha_hat;
dp_npe = (sqrt(1.5*(ebeta_hat'*P_line*ebeta_hat)) - er0) ...
         /(w3d2*sum(eh_i.*etheta_i));
dpiter(iter) = dp_npe;
% ... effektive Spannung & Fließflächen Normale
ebeta_npe = er0 * ebeta_hat / ( er0 + w3d2 * dp_npe * sum(eh_i.*etheta_i));
%     n_npe = w3d2 * (P_hat * ebeta_npe)/ er0;
n_npe = sqrt(ebeta_npe' * P_line * ebeta_npe); if n_npe == 0; n_npe = mynull; end
n_npe = (P_hat * ebeta_npe)/n_npe;
% ... Backstress Tensoren
% Normen Backstress
norm_aei = sqrt( sum( (P_line*ealpha_npe) .* ealpha_npe) );
norm_aei(norm_aei == 0) = mynull;
eLi = ealpha_npe./norm_aei;
% Wichtefunktion
var1 = norm_aei./er_i; var1(var1>1) = 1;
var1 = var1.^echi_i;
var2 = n_npe' * P_check * eLi;
var2 = 0.5 * (var2 + abs(var2));
ew_i_npe = var1 .* var2;
% oft verwendeter Kram
ealpha_tr = ealpha_n + dp_npe * eh_i .* repmat(A * n_npe,1,eM);            % trial Backstress
etheta_i = 1./(1+ew_i_npe.*ec_i*dp_npe);
transn = n_npe';                                                       % Transponierte Normale
% Schleife über Backstress
for i = 1 : eM
    
    % Startwert
    ai = 1/( 1 + ew_i_npe(i) * ec_i(i) * dp_npe) * ealpha_tr(:,i);
    
    % Newtonverfahren
%     [ai,wi,norma] = NewtonBackstress(i,...                             % Index Backstress
%         ai,...                                                         % Startwert Newton
%         ealpha_tr(:,i),...                                             % trial Backstress
%         dp_npe,transn,...                                              % Inkrement plastische Bogenlänge und FFNormale
%         ec_i(i),er_i(i),echi_i(i),...                                  % Parameter
%         P_line,P_check,EINS,...                                        % Abblidungen
%         tolbs,10*maxiter,' ');                                         % Iterationsoptionen
    ci = ec_i(i);
    ri = er_i(i);
    chii = echi_i(i);
    ai_tr = ealpha_tr(:,i);
    norma = sqrt(ai' * P_line * ai); if norma == 0, norma = 1e-40; end
    Li = ai./norma;
    var1 = norma/ri; if var1 > 1, var1 = 1; end
    var1 = var1^chii;
    var2 = transn * P_check * Li;
    var2 = 0.5 * (var2 + abs(var2));
    wi= var1 .* var2;
    % Aktueller Fehler
    f = ai * ( 1 + wi * ci * dp_npe) - ai_tr;
    normf = sqrt(f' * P_line * f);
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
            var6 = Li' * P_line * Li;
            var7 = transn * P_check * Li;
        end
        invdRdA = 1/k1 * EINS - k2/(k1*(k1+k2*var6+k3*var7)) * (Li*(P_line*Li)') - k3/(k1*(k1+k2*var6+k3*var7)) * (Li * (transn*P_check));
        % neuer Backstress
        ai = ai - invdRdA * f;
        % Update alles mögliche
        norma = sqrt(ai' * P_line * ai); if norma == 0, norma = 1e-40; end
        Li = ai./norma;
        % Wichtefunktion
        var1 = norma/ri; if var1 > 1, var1 = 1; end
        var1 = var1^chii;
        var2 = transn * P_check * Li;
        var2 = 0.5 * (var2 + abs(var2));
        wi = var1 .* var2;
        % berechne Fehler
        f = ai * ( 1 + wi * ci * dp_npe) - ai_tr;
        normf = sqrt(f' * P_line * f);
        %             normfiter(iterbs) = normf;
        % keine Endlossschleifen
        if iterbs > maxiter
            msg = ['Keine Konvergenz in ',mystr,' Backstress',num2str(bs),' im Ohno Wang Modell ',...
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
    % Speichern konvergierte Ergebnisse
    norm_aei(i) = norma;
    ew_i_npe(i) = wi;
    ealpha_npe(:,i) = ai;
    etheta_i(i) = 1/(1+wi*ec_i(i)*dp_npe);
    
end % Ende Schleife Über Backstress

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
    % Aitkens Delta^2 Verfahren
%     if mod(iter-1,3) == 0
%         dp_hat = dpiter(iter-1) - ( dpiter(iter-1) - dpiter(iter-2) )^2 / ...
%         ( (dpiter(iter-1) - dpiter(iter-2)) - (dpiter(iter-2)- dpiter(iter-3)));
%         if dp_hat > 0
%             dp_npe = dp_hat;
%         else
%             dp_npe = (sqrt(1.5*(ebeta_hat'*P_line*ebeta_hat)) - er0) ...
%                       /(w3d2*sum(eh_i.*etheta_i));
%         end
%     else
%         dp_npe = (sqrt(1.5*(ebeta_hat'*P_line*ebeta_hat)) - er0) ...
%         /(w3d2*sum(eh_i.*etheta_i));
%     end
    % Aitkens Delta^2 Verfahren
    dp_npe = (sqrt(1.5*(ebeta_hat'*P_line*ebeta_hat)) - er0) ...
        ./(w3d2*sum(eh_i.*etheta_i));
    if mod(iter,3) == 0
        dp_hat = dp_npe - ( dp_npe - dpiter(iter-1) )^2 / ...
        ( dp_npe - 2*dpiter(iter-1) + dpiter(iter-2) );
        if dp_hat > 0
            dp_npe = dp_hat;
        end
    end
    dpiter(iter) = dp_npe;
    % ... effektive Spannung & Fließflächen Normale
    ebeta_npe = er0 * ebeta_hat / ( er0 + w3d2 * dp_npe * sum(eh_i.*etheta_i));
%     n_npe = w3d2 * (P_hat * beta_npe)/ r0;
    n_npe = sqrt(ebeta_npe' * P_line * ebeta_npe); if n_npe == 0; n_npe = mynull; end
    n_npe = (P_hat * ebeta_npe)/n_npe;
%     fprintf('||n|| = %.15f\n',n_npe' * A * P_line * A * n_npe)
    % ... Backstress Tensoren
    % oft verwendeter Kram
    ealpha_tr = ealpha_n + dp_npe * eh_i .* repmat(A * n_npe,1,eM);            % trial Backstress
%     theta_i = 1./(1+w_i_npe.*c_i*dp_npe);
    transn = n_npe';                                                       % Transponierte Normale
    % Schleife über Backstress
    for i = 1 : eM
        
        % Startwert
        ai = 1/( 1 + ew_i_npe(i) * ec_i(i) * dp_npe) * ealpha_tr(:,i);
        
        % Newtonverfahren
%         [ai,wi,norma] = NewtonBackstress(i,...                             % Index Backstress
%             ai,...                                                         % Startwert Newton
%             ealpha_tr(:,i),...                                             % trial Backstress
%             dp_npe,transn,...                                              % Inkrement plastische Bogenlänge und FFNormale
%             ec_i(i),er_i(i),echi_i(i),...                                  % Parameter
%             P_line,P_check,EINS,...                                        % Abblidungen
%             tolbs,10*maxiter,' ');                                         % Iterationsoptionen
        ci = ec_i(i);
        ri = er_i(i);
        chii = echi_i(i);
        ai_tr = ealpha_tr(:,i);
        norma = sqrt(ai' * P_line * ai); if norma == 0, norma = 1e-40; end
        Li = ai./norma;
        var1 = norma/ri; if var1 > 1, var1 = 1; end
        var1 = var1^chii;
        var2 = transn * P_check * Li;
        var2 = 0.5 * (var2 + abs(var2));
        wi= var1 .* var2;
        % Aktueller Fehler
        f = ai * ( 1 + wi * ci * dp_npe) - ai_tr;
        normf = sqrt(f' * P_line * f);
%         normfspeicher = zeros(1,maxiter);
%         normfspeicher(1) = normf;
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
                var6 = Li' * P_line * Li;
                var7 = transn * P_check * Li;
            end
%             dRdA = k1 * EINS + k2 * (Li*(P_line*Li)') + k3 * (Li * (P_check*n_npe)');
            invdRdA = 1/k1 * EINS - k2/(k1*(k1+k2*var6+k3*var7)) * (Li*(P_line*Li)') - k3/(k1*(k1+k2*var6+k3*var7)) * (Li * (transn*P_check));
            % neuer Backstress
            ai = ai - invdRdA * f;
%             ai = ai - dRdA\f;
            % Update alles mögliche
            norma = sqrt(ai' * P_line * ai); if norma == 0, norma = 1e-40; end
            Li = ai./norma;
            % Wichtefunktion
            var1 = norma/ri; if var1 > 1, var1 = 1; end
            var1 = var1^chii;
            var2 = transn * P_check * Li;
            var2 = 0.5 * (var2 + abs(var2));
            wi = var1 .* var2;
            % berechne Fehler
            f = ai * ( 1 + wi * ci * dp_npe) - ai_tr;
            normf = sqrt(f' * P_line * f);
%             normfspeicher(iterbs) = normf;
            % keine Endlossschleifen
            if iterbs > maxiter
                msg = ['Keine Konvergenz in ',mystr,' Backstress',num2str(bs),' im Ohno Wang Modell ',...
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
        % Speichern konvergierte Ergebnisse
        norm_aei(i) = norma;
        ew_i_npe(i) = wi;
        ealpha_npe(:,i) = ai;
        etheta_i(i) = 1/(1+wi*ec_i(i)*dp_npe);
%         fprintf('%i\n',iterbs)
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
%                        Newton Backstress Materialmodell
% -------------------------------------------------------------------------
% ... Berechne Backstress
alpha_npe = alpha_n;
% ... Initialisiere Backstresstensoren
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
alpha_tr = alpha_n + dp_npe * h_i .* repmat(A * n_npe,1,M);                % trial Backstress
% Schleife über Backstress
for i = 1 : M
    
    
    % Startwert
    ai = 1/( 1 + w_i_npe(i) * c_i(i) * dp_npe) * alpha_tr(:,i);
    
    % Newtonverfahren
%     [ai,wi,norma] = NewtonBackstress(i,...                             % Index Backstress
%         ai,...                                                         % Startwert Newton
%         alpha_tr(:,i),...                                              % trial Backstress
%         dp_npe,transn,...                                              % Inkrement plastische Bogenlänge und FFNormale
%         c_i(i),r_i(i),chi_i(i),...                                     % Parameter
%         P_line,P_check,EINS,...                                        % Abblidungen
%         tolbs,10*maxiter,' ');                                          % Iterationsoptionen
    ci = c_i(i);
    ri = r_i(i);
    chii = chi_i(i);
    ai_tr = alpha_tr(:,i);
    norma = sqrt(ai' * P_line * ai); if norma == 0, norma = 1e-40; end
    Li = ai./norma;
    var1 = norma/ri; if var1 > 1, var1 = 1; end
    var1 = var1^chii;
    var2 = transn * P_check * Li;
    var2 = 0.5 * (var2 + abs(var2));
    wi= var1 .* var2;
    % Aktueller Fehler
    f = ai * ( 1 + wi * ci * dp_npe) - ai_tr;
    normf = sqrt(f' * P_line * f);
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
            var6 = Li' * P_line * Li;
            var7 = transn * P_check * Li;
        end
        invdRdA = 1/k1 * EINS - k2/(k1*(k1+k2*var6+k3*var7)) * (Li*(P_line*Li)') - k3/(k1*(k1+k2*var6+k3*var7)) * (Li * (transn*P_check));
        % neuer Backstress
        ai = ai - invdRdA * f;
        % Update alles mögliche
        norma = sqrt(ai' * P_line * ai); if norma == 0, norma = 1e-40; end
        Li = ai./norma;
        % Wichtefunktion
        var1 = norma/ri; if var1 > 1, var1 = 1; end
        var1 = var1^chii;
        var2 = transn * P_check * Li;
        var2 = 0.5 * (var2 + abs(var2));
        wi = var1 .* var2;
        % berechne Fehler
        f = ai * ( 1 + wi * ci * dp_npe) - ai_tr;
        normf = sqrt(f' * P_line * f);
        %             normfiter(iterbs) = normf;
        % keine Endlossschleifen
        if iterbs > maxiter
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
    end % Ende Newtoniteration
    % Speichern konvergierte Ergebnisse
    norm_ai(i) = norma;
    w_i_npe(i) = wi;
    alpha_npe(:,i) = ai;
%     fprintf('%i\n',iterbs)
    
end % Ende Schleife Über Backstress


% -------------------------------------------------------------------------
%                        update Zustandsvariablen
% -------------------------------------------------------------------------
depsp = dp_npe * n_npe;
epsp_npe = epsp_n + depsp;
p_npe = p_n + dp_npe;
eeps_npe = eeps_n + D * (esig_tr-esig_n) + depsp;
ZVARneu = [eeps_npe; epsp_npe ; reshape(ealpha_npe,ntens*eM,1); p_npe; reshape(alpha_npe,ntens*M,1)];
end % Ende Modellgleichung in ESZ





% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
% Newtonverfahren Backstress
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
% function [ai,wi,norma] = NewtonBackstress(bs,...                           % Index Backstress
%     ai,...                                      % Startwert Newtoniteration
%     ai_tr,...                                   % trial Backstress
%     dp,transn,...                       % Inkrement plastische Bogenlänge und FFNormale
%     ci,ri,chii,...                              % Parameter
%     P_line,P_check,EINS,...                     % Abblidungen
%     tol,maxiter,mystr)                          % Iterationsoptionen
% 
% % 1. Upadate Versuch
% %         ai = 1/( 1 + wi * ci * dp) * ai_tr;
% % Berechne neue Werte
% norma = sqrt(ai' * P_line * ai); if norma == 0, norma = 1e-40; end
% Li = ai./norma;
% var1 = norma/ri; if var1 > 1, var1 = 1; end
% var1 = var1^chii;
% var2 = transn * P_check * Li;
% var2 = 0.5 * (var2 + abs(var2));
% wi= var1 .* var2;
% % Aktueller Fehler
% f = ai * ( 1 + wi * ci * dp) - ai_tr;
% normf = sqrt(f' * P_line * f);
% iterbs = 1;
% % Newtoniteration
% while normf > tol
%     % Inkrementiere iter
%     iterbs = iterbs + 1;
%     % berechne Ableitung
%     k1 = 1 + wi * ci * dp;
%     if var2 == 0
%         k2 = 0;
%         k3 = 0;
%         var6 = 0;
%         var7 = 0;
%     else
%         k2 = ci*dp * var1 *  var2 * (chii -1);  ...
%         k3 = ci*dp * var1;
%         var6 = Li' * P_line * Li;
%         var7 = transn * P_check * Li;
%     end
%     invdRdA = 1/k1 * EINS - k2/(k1*(k1+k2*var6+k3*var7)) * (Li*(P_line*Li)') - k3/(k1*(k1+k2*var6+k3*var7)) * (Li * (transn*P_check));
%     % neuer Backstress
%     ai = ai - invdRdA * f;
%     % Update alles mögliche
%     norma = sqrt(ai' * P_line * ai); if norma == 0, norma = 1e-40; end
%     Li = ai./norma;
%     % Wichtefunktion
%     var1 = norma/ri; if var1 > 1, var1 = 1; end
%     var1 = var1^chii;
%     var2 = transn * P_check * Li;
%     var2 = 0.5 * (var2 + abs(var2));
%     wi = var1 .* var2;
%     % berechne Fehler
%     f = ai * ( 1 + wi * ci * dp) - ai_tr;
%     normf = sqrt(f' * P_line * f);
%     %             normfiter(iterbs) = normf;
%     % keine Endlossschleifen
%     if iterbs > maxiter
% %         msg = ['Keine Konvergenz in ',mystr,' Backstress',num2str(bs),' im Ohno Wang Modell ',...
% %             'nach ', num2str(iterbs),' Iterationen. Aktuelle Fehlernorm: '...
% %             num2str(normf)];
% %         warning(msg);
%         %                 fprintf('dp = %.10f n = [%.4f %.4f %.4f %.4f %.4f %.4f] ',dp_npe,n_npe)
%         %                 figure, grid on, hold on
%         %                 plot(normfiter,'.-');
%         %                 set(gca,'YScale','log');
%         break;
%     end
% end % Ende Newtoniteration
% end % Ende Newton Verfahren Backstress
