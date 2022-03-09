function [ZVARneu] = karimohnolangRR2(desig,ZVAR,para,epara)
% Implementierung des Materialmodells nach Karim und Ohno, für pseudo
% stress aproach nach Lang 
%
% ! Nur Spannungssteuerung im ESZ
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
% Parameter
% zuerst Material- dann Struckturmodell 
%     para = [E, nu,  c_i,  r_i,  mu_i,  r0]
%       M = (length(para)-3)/3;
%     epara = [E, nu,  ec_i, er_i, emu_i, er0]
%       eM = (length(epara)-3)/3; 
%
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
% |  Stand: Mai 2020                                                     |
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
% M = (length(para)-4)/6;
M = (length(para)-3)/3;
eM = (length(epara)-3)/3;

% Elastizitätskonstanten
E = para(1);
nu = para(2);
er0 = epara(3*eM+3);

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
ealphai = reshape(ZVAR(2*ntens+1:(eM+2)*ntens),ntens,eM);
% Backstress
alphai = reshape(ZVAR((eM+2)*ntens+2:(eM+M+2)*ntens+1),ntens,M);
% plastische bogenlänge
p = ZVAR((eM+2)*ntens+1);

%--------------------------------------------------------------------------
%                       Elastischer Trial Step
%--------------------------------------------------------------------------
esig_tr = esig + desig;
es = P * esig;                                                             % Pseudo Spannungsdeviator
ea = sum(ealphai,2);                                                       % Pseudo Gesamtbackstress
des = P * desig;                                                           % Inkrement Pseudo Spannungsdeviator
es_tr = es + des;                                                          % Pseudo Versuchsspannungsdeviator
ebeta = es_tr - ea;                                                        % Pseudo effektive Spannung
F_tr = ebeta' * P_line * ebeta - 2/3 * er0^2;                          % Pseudo Überspannung
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
    
    % Integration je nach Spannungszustand
    if ntens == 6 % 3D

        msg = 'nicht implementiert';
        error(msg)

    elseif ntens == 3 && ndi == 2 % ESZ
        
        % statische Matrizen
        [P, P_line, P_hat, A] = set_maps(ntens,ndi);
        % Funktion
        % Matlab Version
        ZVARneu = karimohnoESZlang(...
            esig,eeps,epsp,ealphai,alphai,p,...                              % Zustandsvariablen
            esig_tr,...                                                    % Versuchsspannung                                                           % skaliertes Inkrement und Steuerung
            para,epara,...                                                % Parameter
            D,P,P_line,P_hat,A);
        

    elseif ntens == 1 % 1D

        msg = 'net implementiert';
        error(msg)
        
    end

end % Ende Unterscheidung annehmen/ablehnen trial step
end % Ende Hauptfunktion


% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
% Materialmodell ESZ
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
function ZVARneu = karimohnoESZlang(...
                       esig_n,eeps_n,epsp_n,ealpha_n,alpha_n,p_n,...                  % Zustandsvariablen
                       esig_tr,...                                          % Versuchsspannung
                       para,epara,...                                     % Parameter
                       D,P,P_line,P_hat,A)                       % Abblidungen
                   
% -------------------------------------------------------------------------
%                        Materialparameter
% -------------------------------------------------------------------------
% Zuweisen der Materialparameter
ntens = 3;                                                                 % Tensorkomponenten
M = (length(para)-3)/3;
eM = (length(epara)-3)/3;
% Material
c_i = para(3:2+M);
r_i = para(3+M:2+2*M);
h_i = c_i .* r_i;
mu_i = para(3+2*M:2+3*M);
% Strucktur
ec_i = epara(3:2+eM);
er_i = epara(3+eM:2+2*eM);
eh_i = ec_i .* er_i;
emu_i = epara(3+2*eM:2+3*eM);
er0 = epara(3+3*eM);

% oft verwendete Konstanden
w3d2 = sqrt(3/2);
w2d3 = sqrt(2/3);
% oft verwendete Konstanden
mynull = 1e-40;



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
% dp_npe = 0;
etheta_i = ones(1,eM);
% ... Backstress
% ealpha_npe = ealpha_n;
% ... trial Spannungsdeviator
es_tr = P * esig_tr;

% -------------------------------------------------------------------------
%                   erster Iterationsschritt
% ------------------------------------------------------------------------- 
% ... Inkrement plastische Bogenlänge
ealpha_hat = sum(etheta_i.*ealpha_n,2); 
ebeta_hat = es_tr - ealpha_hat;

% ---------------------
% Fixpunktiteration 
% ---------------------
% dp_npe = (sqrt(1.5*(ebeta_hat'*P_line*ebeta_hat)) - er0) ...
%          /(er0+sum(eh_i.*etheta_i));    
dp_npe = (sqrt(1.5*(ebeta_hat'*P_line*ebeta_hat)) - er0) ...
         /(sum(eh_i.*etheta_i));  
dpiter(iter) = dp_npe;

% ... effektive Spannung & Fließflächen Normale
beta_npe = er0 * ebeta_hat / ( er0 + dp_npe * sum(eh_i.*etheta_i));
n_npe = sqrt(beta_npe' * P_line * beta_npe); if n_npe == 0; n_npe = mynull; end
n_npe = (P_hat * beta_npe)/n_npe;
% ... Update Backstress
etheta_i_hash = 1./(1+ec_i.*emu_i*dp_npe);
ealpha_star = ealpha_n + w2d3 * dp_npe * (A*n_npe) * eh_i; 
ealpha_hash = etheta_i_hash .* ealpha_star;
norm_aei_hash = w3d2 * sqrt( sum( (P_line*ealpha_hash) .* ealpha_hash) );
efi_hash = norm_aei_hash - er_i;
idx = efi_hash > 0;
efi_hash(idx) = 1;
efi_hash(~idx) = 0;
etheta_i = etheta_i_hash .* ( 1 +  efi_hash .* (er_i./norm_aei_hash - 1));
% -------------------------------------------------------------------------
%                        Iteration
% ------------------------------------------------------------------------- 
% ... Schleife
while normp > tol
    % ... Inkrementieren Schleifenzähler
    iter = iter + 1;
    % ... Inkrement plastische Bogenlänge
    ealpha_hat = sum(etheta_i.*ealpha_n,2);
    ebeta_hat = es_tr - ealpha_hat;
    
    % ---------------------
    % Fixpunktiteration
    % ---------------------
%     dp_npe = (sqrt(1.5*(ebeta_hat'*P_line*ebeta_hat)) - er0) ...
%          /(er0+sum(eh_i.*etheta_i));  
    dp_npe = (sqrt(1.5*(ebeta_hat'*P_line*ebeta_hat)) - er0) ...
         /(er0+sum(eh_i.*etheta_i));  
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
    beta_npe = er0 * ebeta_hat / ( er0 + dp_npe * sum(eh_i.*etheta_i));
    n_npe = sqrt(beta_npe' * P_line * beta_npe); if n_npe == 0; n_npe = mynull; end
    n_npe = (P_hat * beta_npe)/n_npe;
    % ... Update Backstress
    etheta_i_hash = 1./(1+ec_i.*emu_i*dp_npe);
    ealpha_star = ealpha_n + w2d3 * dp_npe * (A*n_npe) * eh_i;
    ealpha_hash = etheta_i_hash .* ealpha_star;
    norm_aei_hash = w3d2 * sqrt( sum( (P_line*ealpha_hash) .* ealpha_hash) );
    efi_hash = norm_aei_hash - er_i;
    idx = efi_hash > 0;
    efi_hash(idx) = 1;
    efi_hash(~idx) = 0;
    etheta_i = etheta_i_hash .* ( 1 +  efi_hash .* (er_i./norm_aei_hash - 1));
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
% ... Backstress Material
theta_i_hash = 1./(1+c_i.*mu_i*dp_npe);
alpha_star = alpha_n + w2d3 * dp_npe * (A*n_npe) * h_i;
alpha_hash = theta_i_hash .* alpha_star;
norm_ai_hash = w3d2 * sqrt( sum( (P_line*alpha_hash) .* alpha_hash) );
fi_hash = norm_ai_hash - r_i;
idx = fi_hash > 0;
fi_hash(idx) = 1;
fi_hash(~idx) = 0;
theta_i = theta_i_hash .* ( 1 +  fi_hash .* (r_i./norm_ai_hash - 1));
alpha_npe = theta_i .* alpha_star;
% ... pseudo backsress
ealpha_npe = etheta_i .* ealpha_star;
% ... plastische Dehnungen
depsp = w3d2*dp_npe * n_npe;
epsp_npe = epsp_n + depsp;
% ... plastische Bogenlänge
p_npe = p_n + dp_npe;
% ... Dehnungen
eeps_npe = eeps_n + D * (esig_tr-esig_n) + depsp;
ZVARneu = [eeps_npe; epsp_npe ; reshape(ealpha_npe,ntens*eM,1); p_npe; reshape(alpha_npe,ntens*M,1)];   

end % Ende Modellgleichung in ESZ
