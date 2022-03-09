function [SIG, EPS, EPSP, ALPHA, R, P] = neuberpropsig( ...
                                   ZVAR0, para, ESIG, EEPS, ...
                                   ntens, ndi, material, numink)
% Implentation der Neuber Methode für mehrachsige Kerbnäherung
%
%  !!!!!    Implementierung sodass Spannungen proportional zu 
%           Pseudo Spannungen sind !!!!!!!
%
% INPUT:
% ZVAR0     -> Startwerte der Zustandsvariablen, Zustandsvariblen für die
%              Materialmodelle wie bei dehnungsgesteuerter integration
% para      -> Parameter Materialmodel
% ESIG      -> Pseudo elastische Spannungen
% EEPS      -> Pseudo elastische Dehnungen
% ntens     -> Anzahl Tensorkomponenten
% ndi       -> Anzahl Hauptdiagonalelemente
% material  -> Definiert welches Materialmodell verwendet wird
% numink    -> Anzahl inkremente der Lastfolge
%
% OUTPUT:
% SIG    -> elastisch-plastischer Spannungsverlauf
% EPS    -> elastisch-plastischer Dehnungsverlauf
% EPSP   -> plastischer Dehnungsverlauf
% ALPHA  -> Verlauf der Backstresstensoren
% R      -> Radius FF
% P      -> plastische Bogenlänge
%
% ------------------------------------------------------------------------
% Autor: Jan Kraft                                                        |
% Stand: Januar 2020                                                      |
% ------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Prüfe ob ebener Spannungszustand

if ntens ~= 3 && ndi ~= 2
    msg = 'ESED akt. nur für ESZ gedacht';
    error(msg)
end

% -------------------------------------------------------------------------
% materialfunktion

switch material
    case 'Chaboche'
        
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
        
    case 'OhnoWang'
        
%         fprintf('noch ohnowang modell implementieren\n');
        matfun = @ohnowang;
        
    case 'Jiang'
        
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
        
    otherwise
        
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
end

% -------------------------------------------------------------------------
% Elastizitätskonstanten
E = para(1);                                                               % E-Modul
nu = para(2);                                                              % Querdehnzahl
CEL = elast_steifigkeit(E,nu,ntens,ndi);                                   % Steifigkeit
% DEL = elast_nachgiebigkeit(E,nu,ntens,ndi);                                % Nachgiebigkeit

% -------------------------------------------------------------------------
% Init Referenzzustände
% Annhame Lastfolge startet im Nullzustand -> Alle Refs mit 0
% inititalisieren
% REF(3,5) = [REF_ESIG, REF_EEPS, REF_SIG, REF_EPS]
REF = zeros(3,4);

% -------------------------------------------------------------------------
% Finde ersten plastischen Zustand
weiter = 1;
M = [2,-1,0;-1,2,0;0,0,3]./3;
ML = [2,1,0;1,2,0;0,0,2];
jj = 1;
r = para(end);
while weiter
    s = M * ESIG(:,jj);
    F = s'*ML*s - 2/3*r^2;
    if F < 0
        jj = jj + 1;
    else
        jj = jj - 1;
        weiter = 0;
    end
end

% -------------------------------------------------------------------------
% Speicher für Zustandsvariablen
ZVAR = zeros(size(ZVAR0,1),numink);
% ZVAR(:,1) = ZVAR0;
ZVAR(1:ntens,1:jj) = EEPS(:,1:jj);

% -------------------------------------------------------------------------
% Init altes Inkrement der pseudo Größen
dESIG0 = zeros(3,1);

% -------------------------------------------------------------------------
% Speichern der Energie (zum testen)
% EESED = zeros(3,numink);

% -------------------------------------------------------------------------
% Speicher Differenz aus soll Energieinkrement und ist Inkrement
% f = zeros(3,numink);
% dESED = zeros(3,1);

% -------------------------------------------------------------------------
% speicher für proportionalitätsfaktoren
rho = ones(1,numink);

% -------------------------------------------------------------------------
% Variablen für die Interation
maxiter = 100;
tol = 1e-16;

% -------------------------------------------------------------------------
% Hauptschleife über alle Inkremente der Lastfolge
for ii = jj + 1 : numink
    
    % =====================================================================
    % Inkremente in Pseudo elastischen Größen
    dESIG = ESIG(:,ii) - ESIG(:,ii-1);
    dEEPS = EEPS(:,ii) - EEPS(:,ii-1);
    
    % =====================================================================
    % Prüfe auf Umkehrpunkt
    ukp = dESIG .* dESIG0 < 0;
    
    % =====================================================================
    % Update die Referenzzustände
    SIG = CEL * (ZVAR(1:ntens,ii-1) - ZVAR(1+ntens:2*ntens,ii-1));
    REF(ukp,:) = [ESIG(ukp,ii-1), EEPS(ukp,ii-1), SIG(ukp), ZVAR(ukp,ii-1)];
    
    
    % =====================================================================
    % Prüfe auf erste Ecke 
    eck = [false,false,false];
    for ee = 1 : 3
        if dESIG(ee) ~= 0 && dESIG0(ee) == 0 && ii > jj + 1
            eck(ee) = true;
        end
    end
    
    % =====================================================================
    % Speichern pseudo Spannungsinkrement falls nicht gleich 0
    idx = dESIG ~= 0;
    dESIG0(idx) = dESIG(idx);
    
    % =====================================================================
    % Inkrement der totalen pseudo Verzerrungsenergiedichte
    dETSED = (ESIG(:,ii-1) - REF(:,1)).* dEEPS + ...
             (EEPS(:,ii-1) - REF(:,2)).* dESIG + ...
              dESIG .* dEEPS;
    detsed = sum(dETSED);

     % Speichern der Energie (zu testzwecken)
%     EESED(:,ii) = EESED(:,ii-1) + dEESED;
    
    % =====================================================================
    % Newtoniteration
    rat = rho(ii-1); % strain ratio
    ratiter = zeros(1,maxiter);
    ratiter(1) = rat;
    iter = 0;
    fnorm2 = 1;
    fiter=zeros(1,maxiter);
    verfahren = true; % true = newton, false = Bisektion
    
    while fnorm2 > tol
        
        % keine endlosschleifen
        iter = iter + 1;
        if iter > maxiter + 1
            msg = ['Keine Konvergenz in Newton Verfahren, Lastschritt ',...
				       num2str(ii), '. Fehlerquadrate:', num2str(fnorm2)];
			warning(msg)
            break
        end
        
        % speichern verhältnis
        ratiter(iter) = rat;
        
        % Dehnungsinkrement
        dSIG = rat .* dESIG;
        
        % integration Dehnungsgesteuert
        [ZVAR(:,ii),DEP] = matfun(3, 2, dSIG, ZVAR(:,ii-1), 0, para);
        
        % Spannungsinkrement
		dEPS = ZVAR(1:ntens, ii) - ZVAR(1:ntens, ii-1);
        
        % Energieinkrement
        G = ZVAR(1:ntens,ii-1) - REF(:,4);
        P = SIG - REF(:,3);
        dTSED = P .* dEPS + ...
                G .* dSIG + ...
                dSIG .* dEPS;
        dtsed = sum(dTSED);
        
        % Zielfunktion
        f = detsed - dtsed;
        fiter(iter) = f;
        
        % Abbruchkriterium
        fnorm2 = f^2;
        
        % Verzweigung Newton/Bisektion
        if fnorm2 > tol && verfahren % newton
        
            % MAtrizen an Stelle ii
            G = G + dEPS;
            P = P + dSIG;
            % Ableitung Zielfunktion
            dummy = DEP * dESIG;
            dfdrat = sum(G .* dESIG + P .* dummy);

            % inkrement verhältnis
            drat = f/dfdrat;
    %         rat = rat + drat;
            % Bisektion falls newton osziliert
            if iter > 2
                if fiter(iter) * fiter(iter-1) < 0
                    if fiter(iter) > 0
                        fplus = fiter(iter);
                        ratplus = rat;
                        fminus = fiter(iter-1);
                        ratminus = ratiter(iter-1);
                    else
                        fplus = fiter(iter-1);
                        ratplus = ratiter(iter-1);
                        fminus = fiter(iter);
                        ratminus = rat;
                    end
                    dis = -fiter(iter-1) / (fiter(iter)-fiter(iter-1));
                    rat = dis * rat + (1-dis)*ratiter(iter-1);
                    verfahren = false;
                else
                    % neues Verhältnis
                    rat = rat + drat;
                end
            else
                % neues Verhältnis
                rat = rat + drat;
            end
        
        elseif fnorm2 > tol % Bisektion
            
            % Neue Intervallgrenzen
            if fiter(iter) > 0
                fplus = fiter(iter);
                ratplus = rat;
            else
                fminus = fiter(iter);
                ratminus = rat;
            end
            
            % linearer Abstandsschätzer
            dis = -fminus/(fplus-fminus);
            
            % neues Verhältniss 
            rat = dis * ratplus + (1-dis) * ratminus;
            
        end % Ende Verzweigung Bisektion/Newton 
        
        
    end % Ende Iterationsschleife
    
    fprintf('%i %i\n',iter,ii)
    
    % Speichern verhältnis
    rho(ii) = rat;
    
end % Ende Schleife über Inkremente 


% -------------------------------------------------------------------------
% Herrauslesen der lokalen Größen je nach Material
switch material
    case 'Chaboche'
        
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)        
%         EPS = ZVAR(1:ntens,:);
%         EPSP = ZVAR(ntens+1:2*ntens,:);
%         ALPHA = zeros(ntens,M,numink);
%         for i = 1:M
%             ALPHA(:,i,:) = ZVAR((i+1)*ntens+1:(i+2)*ntens,:);
%         end
%         R = ZVAR(end-1,:);
%         P = ZVAR(end,:);
        
    case 'OhnoWang'
        
        M = (length(para)-3)/3;
        EPS = ZVAR(1:ntens,:);
        EPSP = ZVAR(ntens+1:2*ntens,:);
        SIG = CEL * ( EPS - EPSP);
        ALPHA = zeros(ntens,M,numink);
        for i = 1:M
            ALPHA(:,i,:) = ZVAR((i+1)*ntens+1:(i+2)*ntens,:);
        end
        P = ZVAR(end,:);
        R = para(end) * ones(1,numink);
        
    case 'Jiang'
        
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
%         EPS = ZVAR(1:ntens,:);
%         EPSP = ZVAR(ntens+1:2*ntens,:);
%         ALPHA = zeros(ntens,M,numink);
%         for i = 1:M
%             ALPHA(:,i,:) = ZVAR((i+1)*ntens+1:(i+2)*ntens,:);
%         end
%         P = ZVAR(end-1,:);
%         R = para(end) * ones(1,numink);
        
    otherwise
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
end




end