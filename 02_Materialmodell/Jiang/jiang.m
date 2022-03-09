function [X_neu,CEP] = jiang(ntens, ndi, ink, X, ink_flag, parameter)
% Jiang Plastizitätsmodell 
%
% QUELLE:
%        Jiang, Sehitoglu - Modeling of Cyclic Rachtetting Plasticity
%
% INPUT:
%         ntens     -> Anzahl Tensorkomponten (in Spannungen)
%         ndi       -> Anzahl Diagonalelemente, mehr kkommentare
%        ink:       -> Belastungsinkrement
%        X:         -> Zustandsvariablen
%        ink_flag:  -> (0) Spansteu oder (1) Dehnsteu
%        parameter: -> Materialparameter für das Modell
%
% OUTPUT:
%        X_neu:     -> Zustandsvariable nach Belastungsinkrement
%
% Spannungsgesteuerte Prozesse:
% X = [eps; epsp; alphai; p; rm]
% ink_flag = 0
%
% Dehnungsgesteuerte Prozesse:
% X = [sig; epsp; alphai ;p; rm]
% ink_flag = 1
%
% Materialparameter:
% paramter = [E, nu, 
%             a_chi, b_chi, 
%             ak, ck, k1 
%             cm,
%             c_i0, a_i1, b_i1, a_i2, b_i2,
%             r_i, 
%             Q_i]
%
% Elastische Parameter:
% E     ->      Elastizitätsmodul
% nu    ->      Querdehnzahl
% k1    ->      Schubfließspannung
%
% Pamaneter für nichtproportionalität und Spannungsleveleffekt:
% a_chi, b_chi, Q_i      
%
% Non-Masing-Verhalten:
% ak, ck
%
% Gedächtnisfläche:
% cm
%
% Zyklische Ver- und Entfestigung:
% c_i0, a_i1, b_i1, a_i2, b_i2
%
% Verfestigungsregel:
% r_i
%
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft                                                   |
% |  Stand: Januar 2020                                                 |
%  ----------------------------------------------------------------------

%--------------------------------------------------------------------------
%                  Materialparamter ermitteln                                  
%--------------------------------------------------------------------------
M = (length(parameter)-8)/7;                                               % Anzahl Backstresstensoren

E = parameter(1);                                                          % Elastizitätsmodul
nu = parameter(2);                                                         % Querdehnzahl   
C = elast_steifigkeit(E,nu,ntens,ndi);

ak = parameter(5);                                                         % parameter zum bestimmen der Größe der FF
ck = parameter(6);
k1 = parameter(7);

%--------------------------------------------------------------------------
%                        Hilfsmatrizen                                    
%--------------------------------------------------------------------------

% statische Matrizen
[P, P_line] = set_maps(ntens,ndi);

%--------------------------------------------------------------------------
%            Identifikation der Zustandsvariablen                         
%--------------------------------------------------------------------------

if ink_flag == 0                                                           % Spannungssteuerung
    eps = X(1:ntens);                                                      % Bestimmen der Dehnung
    epsp = X(ntens+1:2*ntens);                                             % Bestimmen der platischen Dehnung
    sig = C * (eps - epsp);
    dsig_tr = ink;                                                         % Versuchsspannung
    
elseif ink_flag == 1                                                       % Dehnungssteuerung
    sig = X(1:ntens);                                                      % Bestimmen der aktuellen Spannung
    dsig_tr = C * ink;                                                     % Berechnung des Versuchsspannungsinkrements
end

alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);                         % akt. Backstresstensoren ermitteln


%--------------------------------------------------------------------------
%               Radius der Fließfläche berechnen                          
%--------------------------------------------------------------------------
rm = X((M+2)*ntens+2);                                                     % Radius der Gedächtnisfläche bestimmen
k = k1*(1 + ak * exp(ck * rm));                                            % Radius der Fließfläche berechnen

%--------------------------------------------------------------------------
%                       Elastischer Trial Step
%--------------------------------------------------------------------------


if ntens == 1
    s = P * sig;                                                           % Berechnung des Spannungsdeviators
    a = sum(alpha,2);                                                      % Berechnung des Backstresstensors;
    ds = P * dsig_tr;                                                      % Deviator des Spannungsinkrements
    s_tr = s + ds;                                                         % Deviator der Versuchsspannung
    beta = s_tr - a;                                                       % Relativ Spannung
    F_tr = abs(beta) - sqrt(3) * k;
elseif ntens == 2 % sigma-tau
    s = P .* sig;                                                          % Berechnung des Spannungsdeviators
    a = sum(alpha,2);                                                      % Berechnung des Backstresstensors;
    ds = P .* dsig_tr;                                                     % Deviator des Spannungsinkrements
    s_tr = s + ds;                                                         % Deviator der Versuchsspannung
    beta = s_tr - a;                                                       % Relativ Spannung
    F_tr = beta' * (P_line .* beta) - 2 * k^2;                             % Trial Fließfunktion
else
    s = P * sig;                                                           % Berechnung des Spannungsdeviators
    a = sum(alpha,2);                                                      % Berechnung des Backstresstensors;
    ds = P * dsig_tr;                                                      % Deviator des Spannungsinkrements
    s_tr = s + ds;                                                         % Deviator der Versuchsspannung
    beta = s_tr - a;                                                       % Relativ Spannung
    F_tr = beta' * P_line * beta - 2 * k^2;                                % Trial Fließfunktion
end   
FTOL = 1e-7;                                                              % Toleranz für Abweichungen F ~= 0



%--------------------------------------------------------------------------
%                   Trial Step wird angenommen
%--------------------------------------------------------------------------
if F_tr < FTOL                                                             % Schritt ist vollkommen elastisch
        
       if ink_flag == 0    
            D = elast_nachgiebigkeit(E,nu,ntens,ndi);
            eps = eps + D * ink;                                           % Dehnung neu berechnen
            X_neu = [eps; X(ntens+1:end)];                                 % Zustandsvariablen neu bestimmen
            
       elseif ink_flag == 1                                                
            sig = sig + dsig_tr;                                           % Spannung neu berechnen
            X_neu = [sig; X(ntens+1:end)];                                 % Zustandsvariablen neu bestimmen
            
       end
       
       % Falls tangentiale Steifigkeit gebraucht wird
       if nargout == 2 % Berechne elastisch plastische nachgiebigkeit
            CEP = elast_steifigkeit(E,nu,ntens,ndi);
       end

       
%--------------------------------------------------------------------------
%                   Trial Step wird abgelehnt
%--------------------------------------------------------------------------  
else                                                                       % Schritt enhält einen platischen Anteil
        %------------------------------------------------------------------
        %     1. elastischen Anteil ermitteln
        %------------------------------------------------------------------
%             xel = elastink(s,a,k,ds,ndi,'jiang');
            xel = elastink2(s,a,P_line,k,ds,FTOL,'jiang');           
        %------------------------------------------------------------------
        %     2. elastischen Anteil aufbringen
        %------------------------------------------------------------------        
            if ink_flag == 0  
                D = elast_nachgiebigkeit(E,nu,ntens,ndi);
                eps = eps + D * (xel*ink);                                  % Dehnung nach Aufbringen des elastischen Anteils berechnen
                X = [eps; X(ntens+1:end)];                                  % Zustandsvariablen neu bestimmen
            
            elseif ink_flag == 1                                            
                sig = sig + xel*dsig_tr;                                    % Spannung nach aufbringen des elatischen Anteils berechnen
                X = [sig; X(ntens+1:end)];                                  % Zustandsvariablen neu bestimmen
            
            end
        %------------------------------------------------------------------
        %    3. platischen Anteil bestimmen
        %------------------------------------------------------------------  
            ink = (1-xel)*ink;                                              % Plastischer Anteil vom Inrement ermitteln
        
        %------------------------------------------------------------------
        %    4. Integration des platischen Anteils je nach
        %    Spannungszustands
        %------------------------------------------------------------------
            % Integration je nach Spannungszustand
            options = [];
            if ntens == 6 % 3D
                
                % elast nachgiebigkeit
                D = elast_nachgiebigkeit(E,nu,ntens,ndi);
                % Abblidungen
                [P, P_line, P_hat] = set_maps(ntens,ndi);
                
                % Integration 
                
                % Matlab Version
%                 [~,X_neu] = rk87(@materialmodell3d,[0,1], X, options, ink,...
%                         ink_flag, parameter, C, D, P, P_line, P_hat);

                % Mex Version
                [~,X_neu] = rk87(@CoderVersion_3D_Jiang_mex,[0,1], X, options,...
                            M,ink, ink_flag, parameter, C, D, P, P_line, P_hat);
                    
            elseif ntens == 3 && ndi == 2 % ESZ
                
                % elast nachgiebigkeit
                D = elast_nachgiebigkeit(E,nu,ntens,ndi);
                % abbildungen
                [P, P_line, P_hat, A, P_check] = set_maps(ntens,ndi);
                % Integration
                
                % Matlab Version
%                 [~,X_neu] = rk87(@materialmodellESZ,[0,1], X, options, ink,...
%                         ink_flag, parameter, C, D, P, P_line, P_hat, A, ...
%                         P_check);

                % Mex version
                [~,X_neu] = rk87(@CoderVersion_ESZ_Jiang_mex,[0,1], X, options, ...
                            M, ink, ink_flag, parameter, C, D, P, P_line, P_hat, A, ...
                            P_check);
            
            elseif ntens == 2 && ndi == 1 % sigma-tau
                
                % elast nachgiebigkeit
                D = elast_nachgiebigkeit(E,nu,ntens,ndi);
                % abbildungen
                [P, P_line, P_hat, A, P_check] = set_maps(ntens,ndi);
                % Integration
                
                % Matlab Version
%                 [~,X_neu] = rk87(@materialmodellST,[0,1], X, options, ink,...
%                         ink_flag, parameter, C, D, P, P_line, P_hat, A, ...
%                         P_check);
                
                % Mex version
                [~,X_neu] = rk87(@CoderVersion_Jiang_ST_mex,[0,1], X, options, ...
                            M, ink, ink_flag, parameter, C, D, P, P_line, P_hat, A, ...
                            P_check);    
                    
            elseif ntens == 1 % 1D
                
                D = 1/C;
                
                [~,X_neu] = rk87(@materialmodell1d,[0,1], X, options, ink,...
                        ink_flag, parameter);
                    
            end
    
            % nur letzten Schritt ausgeben
            X_neu = X_neu(end,:)'; 
            
            % Testen der Konsistenzbedingung
            X_neu = konsistenzbedingung(X_neu,M,P,P_line,ntens,C,D,ink_flag,k1,ak,ck);
            
            % Falls tangentiale Steifigkeit gebraucht wird
            if nargout == 2 % berechne Elastisch plastische steifigkeitsmatrix
                CEP = tangsteifigkeit_jiang(X_neu,parameter,ink_flag,ntens,ndi);
            end                                          % Letzten Zeitschritt ausgeben
    
end

end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Prüfe die konsistenzbedingung                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = konsistenzbedingung(X,M,P,P_line,ntens,C,D,ink_flag,k1,ak,ck)
% Window Ausgabe
% fprintf('Zustandsvariablen Original:\n')
% fprintf('%.32d \n',X)


% Identifikation der Zustandsvariablen
if ink_flag == 0                                                         % Spannungssteuerung
    eps = X(1:ntens);                                                   % Dehnungstensor
    epsp = X(ntens+1:2*ntens);                                          % Plastischer Anteil der Dehnung
    sig = C * (eps - epsp);                                             % Spannung
    
elseif ink_flag == 1                                                     % Dehnungssteuerung
    sig = X(1:ntens);                                                   % Spannung
end
alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);                             % akt. Backstresstensoren ermitteln
a = sum(alpha,2);
%  Radius der Fließfläche berechnen
delta = 1e-40;                                                          % abfangen der null
rm = X((M+2)*ntens+2);
if rm == 0                                                              % teilen durch Null Vermeiden
    rm = delta;
end
r = k1 * ( 1 + ak * exp( ck * rm ) );                                  % radius FF


% Relativspannung
% Fließfläche
if ntens == 1
    s = P * sig;
    beta = s - a;
    SV = 1/3 * beta^2;
elseif ntens == 2
    s = P .* sig;
    beta = s - a;
    SV = 0.5 * beta' .* (P_line .* beta);
else
    s = P * sig;
    beta = s - a;
    SV = 0.5 * beta' * P_line * beta;
end
% F1 = SV - k^2;
% fprintf('Überspannung Orig: %.4d \n',F1)

% Prüfe Konsistenzbedingung
if SV - r^2 > 1e-11
    
    % Fehler/Warnung wenn zu stark verletzt
    if SV - r^2 > 1e-1
        msg = 'Konsistenzbedinung verletzt';
        error(msg)
    elseif SV - r^2 > 1e-4
        msg = 'Konsistenzbedinung geringfügig verletzt';
        warning(msg)
    end

    % Korrigiere Spannungsdeviator (Rel. Span. radial auf FF zurückprojezieren)
    fak = r/sqrt(SV);
    if ntens == 6
        % ... 3D nur Deviator Korrigieren
        s = a + fak*beta;
        sh = (sig(1)+sig(2)+sig(3))/3;
%         eins = [1; 1; 1; 0; 0; 0];
        % Korrigierte Spannung 
%         sig = s + [sh; sh; sh; 0; 0; 0];
        sig = s;
        sig(1) = sig(1) + sh;
        sig(2) = sig(2) + sh;
        sig(3) = sig(3) + sh;
    elseif ntens == 3
        % ... ESZ nur Deviator korriegieren 
        s = a + fak*beta;
        % korrigiere Spannungen
        sig = [2 1 0; 1 2 0; 0 0 1] * s;
%         sh = (sig(1)+sig(2))/3;
%         eins = [1; 1; 0];
    elseif ntens == 2
        % ... sig-tau nur Deviator korrigieren
        s = a + fak*beta;
        % korrigiere Spannungen
        sig = [1.5*s(1);s(2)];
    elseif ntens == 1
        % ... 1D Spannung korrigieren
        s = a + fak*beta;
        % Korrigierte Spannung 
        sig =  s;
    end
    


    % Korrigierte Zustandsgrößen
    if ink_flag == 0 % Spansteu.
        epse = D*sig;
        eps = epse+epsp;
        X(1:ntens) = eps;

    elseif ink_flag == 1 % Dehnsteu
        X(1:ntens) = sig;
    end


end
end % Ende Prüfen Konsitenzbedingung

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Materialgleichungen 3D                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dX = materialmodell3d(~, X, ink, ink_flag, parameter, C, D, P,...
                               P_line, P_hat)
% Konkretes Materialmodell nach Jiang für 3d
%
% INPUT:
%         t         -> Zeit (auch wenn nicht explizit gebraucht)
%         X         -> Zustand
%         ink       -> Lastinkrement
%         ink_flag  -> Belastungsart
%         parameter -> Modellparameter
%         C,D       -> Elastische Steifigkeit und Nachgiebigkeit   
%   P, P_hat ...    -> Diverse Abbildungen
%
%   OUTPUT:
%         dX -> Inkremente der Zustandsvariablen
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
%  ----------------------------------------------------------------------   


    %-----------------------------------------------------------------------
    %            Identifikation Materialparameter                       
    %-----------------------------------------------------------------------
    ntens = 6;                                                             % Tensorkomponenten
    M = (length(parameter)-8)/7;                                          % Anzahl Backstresstensoren 
    
    a_chi = parameter(3);
    b_chi = parameter(4);
    ak = parameter(5);
    ck = parameter(6);
    k1 = parameter(7);
    cm = parameter(8);
    c_i0 = parameter(9:9+M-1);
    a_i1 = parameter(9+M:9+2*M-1);
    b_i1 = parameter(9+2*M:9+3*M-1);
    a_i2 = parameter(9+3*M:9+4*M-1);
    b_i2 = parameter(9+4*M:9+5*M-1);
    r_i = parameter(9+5*M:9+6*M-1);
    Q_i = parameter(9+6*M:9+7*M-1);
   

    %-----------------------------------------------------------------------
    %            Identifikation der Zustandsvariablen                        
    %-----------------------------------------------------------------------

   if ink_flag == 0                                                         % Spannungssteuerung
        eps = X(1:ntens);                                                   % Dehnungstensor                                    
        epsp = X(ntens+1:2*ntens);                                          % Plastischer Anteil der Dehnung
        sig = C * (eps - epsp);                                             % Spannung

   elseif ink_flag == 1                                                     % Dehnungssteuerung
        sig = X(1:ntens);                                                   % Spannung  
   end
   alpha = reshape(X(13:(M+2)*ntens),ntens,M);                             % akt. Backstresstensoren ermitteln

   %-----------------------------------------------------------------------
   %               Radius der Fließfläche berechnen                         
   %-----------------------------------------------------------------------
   delta = 1e-40;                                                          % abfangen der null 
   p = X((M+2)*ntens+1); 
   rm = X((M+2)*ntens+2);
   if rm == 0                                                              % teilen durch Null Vermeiden
       rm = delta;
   end
   k = k1 .* ( 1 + ak * exp( ck * rm ) );                                  % radius FF

   %-----------------------------------------------------------------------
   %               Normale an die Fließfläche berechnen                      
   %-----------------------------------------------------------------------    
   s = P * sig;                                                            % Spannnungs deviator
   a = sum(alpha,2);                                                       % Backstress
    
   beta = s - a;                                                           % effektive Spannung
 
   norm_a = sqrt( sum( (P_line*a) .* a) );                                 % norm backstress
   if norm_a == 0                                                          % Abfangen teilen durch null
        norm_a = delta;
   end

   n = P_line * beta./(sqrt(2)*k);                                         % Normale
   nTrans = n';                                                            % Tanspositon
   
   %-----------------------------------------------------------------------
   %              Ableitungen der Teilbackstresstensoren                    
   %-----------------------------------------------------------------------

   c_i = c_i0 .* (ones(1,M) + a_i1.*exp(-b_i1.*p) ...                      % Variable für zyklische Ver- bzw. Entfestigung
                                     + a_i2.*exp(-b_i2.*p));

   % normen der Teilbackstresstensoren 
   norm_ai = sqrt( sum( (P_line*alpha) .* alpha) );
   norm_ai(norm_ai == 0) = delta; 
   
   % Hilfsvariablen
   Li = alpha ./ norm_ai;                                                  % Normale des i-ten Backstresstensors
   var0 = nTrans * Li;                                                     % Hilfsvariable
   chi = Q_i .* (2 - var0) .* (1 + a_chi * exp(b_chi * rm));               % Exponent für Nichtproportionalität und Spannungsleveleffekt auf Rachetting
   var2 = P_hat * n * r_i;                                                 % Hilfsvariable
   var3 = norm_ai ./ r_i;
   var3(var3 > 1) = 1;
   var3 = var3.^chi.*alpha;
   
   % Ableitung Teilbackstresstensoren
   % !!!! IM ESZ UND IM OHNO ALS SCHLEIFE IMPLENTIERT; SCHLEIFE SCHEIN EIN
   % BISSCHEN ( ABER NICHT WESENTLICH SCHNELLER )
   dalpha_dp = c_i .* (var2 - var3);
   
   % Ableitung Gesamtbackstresstensor
   da_dp = sum(dalpha_dp,2); 


   %-----------------------------------------------------------------------
   %                  Ableitungen der Gedächtnissfläche                     
   %-----------------------------------------------------------------------
    
    % Hilfsvariablen
    L = a./norm_a;
    dummy = sign(norm_a - rm);
    hg = 0.5 * (dummy + abs(dummy));
    dummy = sum( (P_line * da_dp) .* L);
    mla = 0.5 * ( dummy + abs(dummy) );
    dummy = 1 - norm_a/rm;
    dummy = 0.5 * ( dummy + abs(dummy) );
    
    % Ableitung mem.surf
    drm_dp = hg * mla - cm * dummy;
    % Ableitung radius fließfläche
    dk_dp = k1 * ak * ck * exp(ck*rm) * drm_dp;
    % Abfangen von 0*Inf = NaN
    dk_dp(isnan(dk_dp)) = 0; 
    
   %-----------------------------------------------------------------------
   %                       Plastischer Modul                        
   %-----------------------------------------------------------------------

    h = nTrans * da_dp + sqrt(2) * dk_dp;                                   % Plastischer Modul

   %-----------------------------------------------------------------------
   %                Inkrement der plastischen Bogenlänge                     
   %-----------------------------------------------------------------------
        if ink_flag == 0                                                        
            
            dp = (nTrans * ink)/h;                                          % Plast. Inkrement bei Spannungssteuerung 
            
        elseif ink_flag == 1  
            
            dp = ( nTrans * (C * ink)) / (h + nTrans * (C * n));            % Plast. Inkrement bei Dehnungssteuerung

        end
    
   %-----------------------------------------------------------------------
   %                  Inkremente der Zustandsvariablen                      
   %-----------------------------------------------------------------------
    depsp = n * dp;
    dalpha = dalpha_dp * dp;
    drm = drm_dp * dp;

   %-----------------------------------------------------------------------
   %                   Zusammenfassen der Inkremente                        
   %-----------------------------------------------------------------------

    if ink_flag == 0                                                        % Spannungssteuerung
        deps = D*ink + depsp;
        out = vertcat(deps, depsp);
        
    elseif ink_flag == 1                                                    % Dehnungssteuerung
        dsig = C * (ink - depsp);
        out = vertcat(dsig, depsp);

    end

    tmp = reshape(dalpha,ntens*M,1);                                       % Umsortieren von dalpha zu einem Vektor der Größe ntens*num_alpha
    dX = vertcat(out, tmp, dp, drm);                                       % Inkremente der Zustandsvariablen    

end % Ende 3D

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Materialgleichungen ESZ                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dX = materialmodellESZ(~, X, ink, ink_flag, parameter, C, D, ...
                                 P, P_line, P_hat, A, P_check)
% Konkretes Materialmodell nach Jiang für ESZ
%
% INPUT:
%         t         -> Zeit (auch wenn nicht explizit gebraucht)
%         X         -> Zustand
%         ink       -> Lastinkrement
%         ink_flag  -> Belastungsart
%         parameter -> Modellparameter
%         C,D       -> Elastische Steifigkeit und Nachgiebigkeit   
%   P, P_line ...   -> Diverse Abbildungen
%
%
%   OUTPUT:
%         dX -> Inkremente der Zustandsvariablen
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
%  ----------------------------------------------------------------------

    %-----------------------------------------------------------------------
    %            Identifikation Materialparameter                       
    %-----------------------------------------------------------------------
    ntens = 3;                                                             % Tensorkomponenten
    M = (length(parameter)-8)/7;                                           % Anzahl Backstresstensoren 
    
    a_chi = parameter(3);
    b_chi = parameter(4);
    ak = parameter(5);
    ck = parameter(6);
    k1 = parameter(7);
    cm = parameter(8);
    c_i0 = parameter(9:9+M-1);
    a_i1 = parameter(9+M:9+2*M-1);
    b_i1 = parameter(9+2*M:9+3*M-1);
    a_i2 = parameter(9+3*M:9+4*M-1);
    b_i2 = parameter(9+4*M:9+5*M-1);
    r_i = parameter(9+5*M:9+6*M-1);
    Q_i = parameter(9+6*M:9+7*M-1);

   %-----------------------------------------------------------------------
   %            Identifikation der Zustandsvariablen                        
   %-----------------------------------------------------------------------

       if ink_flag == 0                                                     
           eps = X(1:ntens);                                               % Dehnung aus Zustandsvariable lesen                                
           epsp = X(ntens+1:2*ntens);                                      % Plastischer Anteil der Dehnung aus Zustandsvariable lesen
           sig = C * (eps - epsp);                                         % Spannung berechnen

       elseif ink_flag == 1                                                 
           sig = X(1:ntens);                                               % Spannung aus Zustandsvariable lesen
            
       end
   alpha = reshape(X(7:(M+2)*ntens),ntens,M);                              % akt. Backstresstensoren ermitteln

   %-----------------------------------------------------------------------
   %               Radius der Fließfläche berechnen                         
   %-----------------------------------------------------------------------
   delta = 1e-40;
   p = X((M+2)*ntens+1);                                                   % plastische Bogenlänge
   rm = X((M+2)*ntens+2);                                                  % Radius der Gedächtnisfläche
   if rm == 0                                                           
        rm = delta;                                                        % Durch 0 teilen abfangen
   end

   %-----------------------------------------------------------------------
   %               Normale an die Fließfläche berechnen                      
   %-----------------------------------------------------------------------    
   s = P * sig;                                                            % Deviator des Spannungstensors
   a = sum(alpha,2);                                                       % Backstress
   k = k1 .* ( 1 + ak * exp( ck * rm ) );                                  % radius FF
   beta = s - a;                                                           % effektive Spannung
   
   n = P_hat * (beta./(sqrt(2)*k));
   nTrans = n';                                                            % Transposition                             
   
   norm_a = sqrt( sum( (P_line*a) .* a) );                                 % Norm Backstress
   norm_a( norm_a == 0) = delta;

   %-----------------------------------------------------------------------
   %              Ableitungen der Teilbackstresstensoren                    
   %-----------------------------------------------------------------------

   c_i = c_i0 .* (ones(1,M) + a_i1.*exp(-b_i1.*p) ...                      % Variable für zyklische Ver- bzw. Entfestigung
                                    + a_i2.*exp(-b_i2.*p));
                                 
   % Norm von alpha_i ermitteln
   norm_ai = sqrt( sum( (P_line*alpha) .* alpha ) );
   norm_ai(norm_ai == 0) = delta; 

%    % Hilfsvariablen
%    Li = alpha ./ norm_ai;                                                  % Normale des i-ten Backstresstensors
%    var0 = nTrans * P_check * Li;                                           % Hilfsvariable
%    chi = Q_i .* (2 - var0) .* (1 + a_chi * exp(b_chi * rm));               % Exponent für Nichtproportionalität und Spannungsleveleffekt auf Rachetting
%    var0 = A * n * r_i;
%    var1 = norm_ai./r_i;
%    var1(var1 > 1) = 1;
%    var1 = var1.^(chi).*alpha;
%    % Ableitung Teilbackstress
%    dalpha_dp = c_i.* ( var0 - var1);
   
   % Hilfsvariablen ( als schleife is es irgendwie schneller )
   Li = alpha ./ norm_ai; 
   var0 = nTrans * P_check * Li;                                           % Hilfsvariable
   chi = Q_i .* (2 - var0) .* (1 + a_chi * exp(b_chi * rm)); 
   var0 = c_i .* r_i;
   var1 = norm_ai./r_i;
   var1( var1 > 1) = 1;
   var2 = A * n;
   
   % schleife über alle Backstresstensoren
   dalpha_dp = zeros(ntens,M);
   for i = 1 : M
       
       dalpha_dp(:,i) = var0(i) * ( var2 - var1(i).^(chi(i)+1).* Li(:,i) );
       
   end
   
   % Ableitung gesamtbackstress
   da_dp = sum(dalpha_dp,2);                                               % Ableitung: Backstresstensor nach plast. Bogenlänge 


   %-----------------------------------------------------------------------
   %                  Ableitungen der Gedächtnissfläche                     
   %-----------------------------------------------------------------------
    
    L = a./norm_a;                                                         % Normierter Backstresstensor
    dummy = sign(norm_a - rm);                                             % Hilfsvariable
    hg = 0.5 * ( dummy + abs(dummy));                                      % Heaviside
    dummy = sum( (P_line * da_dp) .* L);                                   % Hilfsvariable
    mla = 0.5 * ( dummy + abs(dummy) );                                    % maccaulay
    dummy = 1- norm_a/rm;                                                  % Hilfsvariable
    dummy = 0.5 * ( dummy + abs(dummy) );                                  % Hilfsvariable
    
    % Ableitung Mem. surf 
    drm_dp = hg * mla - cm * dummy;
    
    % Ableitung FF
    dk_dp = k1 * ak * ck * exp(ck*rm) * drm_dp;                            % Ableitung: Fließfläche nach plast. Bogenlänge
    % Abfangen von 0*Inf = NaN
    dk_dp(isnan(dk_dp)) = 0; 
    
   %-----------------------------------------------------------------------
   %                       Plastischer Modul                        
   %-----------------------------------------------------------------------

    h = nTrans * P_check * da_dp + sqrt(2) * dk_dp;                        % Plastischer Modul

   %-----------------------------------------------------------------------
   %                Inkrement der plastischen Bogenlänge                     
   %-----------------------------------------------------------------------
  
   if ink_flag == 0                                                        
            
       dp = (nTrans * ink)/h;                                              % Plast. Inkrement bei Spannungssteuerung 
            
   elseif ink_flag == 1  
            
       dp = ( nTrans * (C * ink)) / (h + (nTrans * (C * n)));              % Plast. Inkrement bei Dehnungssteuerung

   end
    
   %-----------------------------------------------------------------------
   %                  Inkremente der Zustandsvariablen                      
   %-----------------------------------------------------------------------
   
    depsp = n * dp;                                                        % Inkrement der plast. Dehnung
    dalpha = dalpha_dp * dp;                                               % Inkrement der Backstresstensoren
    drm = drm_dp * dp;                                                     % Inkrement Radius der Gedächtnisfläche

   %-----------------------------------------------------------------------
   %                   Zusammenfassen der Inkremente                        
   %-----------------------------------------------------------------------

    if ink_flag == 0                                                        	
        deps = D * ink + depsp;                                              % Inkrement der Dehnung berechnen
        out = vertcat(deps, depsp);                                        % Hilfsvariable
        
    elseif ink_flag == 1                                                    
        dsig = C * (ink - depsp);                                          % Inkrement der Spannung berechnen
        out = vertcat(dsig, depsp);                                        % Hilfsvariable

    end
    dX = vertcat(out, reshape(dalpha,ntens*M,1), dp, drm);                                       % Inkremente der Zustandsvariablen    

end % Ende ESZ





















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Materialgleichungen ESZ                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dX = materialmodellST(~, X, ink, ink_flag, parameter, C, D, ...
                                 P, P_line, P_hat, A, P_check)
% Konkretes Materialmodell nach Jiang für sigma-tau Spannungszustand
%
% INPUT:
%         t         -> Zeit (auch wenn nicht explizit gebraucht)
%         X         -> Zustand
%         ink       -> Lastinkrement
%         ink_flag  -> Belastungsart
%         parameter -> Modellparameter
%         C,D       -> Elastische Steifigkeit und Nachgiebigkeit   
%   P, P_line ...   -> Diverse Abbildungen
%
%
%   OUTPUT:
%         dX -> Inkremente der Zustandsvariablen
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2021 Tu Darmstadt FG Werkstoffmechanik             |
%  ----------------------------------------------------------------------

%-----------------------------------------------------------------------
%            Identifikation Materialparameter
%-----------------------------------------------------------------------
ntens = 2;                                                             % Tensorkomponenten
M = (length(parameter)-8)/7;                                           % Anzahl Backstresstensoren

a_chi = parameter(3);
b_chi = parameter(4);
ak = parameter(5);
ck = parameter(6);
k1 = parameter(7);
cm = parameter(8);
c_i0 = parameter(9:9+M-1);
a_i1 = parameter(9+M:9+2*M-1);
b_i1 = parameter(9+2*M:9+3*M-1);
a_i2 = parameter(9+3*M:9+4*M-1);
b_i2 = parameter(9+4*M:9+5*M-1);
r_i = parameter(9+5*M:9+6*M-1);
Q_i = parameter(9+6*M:9+7*M-1);

%-----------------------------------------------------------------------
%            Identifikation der Zustandsvariablen
%-----------------------------------------------------------------------

if ink_flag == 0
    eps = X(1:ntens);                                               % Dehnung aus Zustandsvariable lesen
    epsp = X(ntens+1:2*ntens);                                      % Plastischer Anteil der Dehnung aus Zustandsvariable lesen
    sig = C * (eps - epsp);                                         % Spannung berechnen
    
elseif ink_flag == 1
    sig = X(1:ntens);                                               % Spannung aus Zustandsvariable lesen
    
end
alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);                              % akt. Backstresstensoren ermitteln

%-----------------------------------------------------------------------
%               Radius der Fließfläche berechnen
%-----------------------------------------------------------------------
delta = 1e-40;
p = X((M+2)*ntens+1);                                                   % plastische Bogenlänge
rm = X((M+2)*ntens+2);                                                  % Radius der Gedächtnisfläche
if rm == 0
    rm = delta;                                                        % Durch 0 teilen abfangen
end

%-----------------------------------------------------------------------
%               Normale an die Fließfläche berechnen
%-----------------------------------------------------------------------
s = P .* sig;                                                            % Deviator des Spannungstensors
a = sum(alpha,2);                                                       % Backstress
k = k1 .* ( 1 + ak * exp( ck * rm ) );                                  % radius FF
beta = s - a;                                                           % effektive Spannung

n = P_hat .* (beta./(sqrt(2)*k));
nTrans = n';                                                            % Transposition

norm_a = sqrt( sum( (P_line.*a) .* a) );                                 % Norm Backstress
norm_a( norm_a == 0) = delta;

%-----------------------------------------------------------------------
%              Ableitungen der Teilbackstresstensoren
%-----------------------------------------------------------------------

c_i = c_i0 .* (ones(1,M) + a_i1.*exp(-b_i1.*p) ...                      % Variable für zyklische Ver- bzw. Entfestigung
    + a_i2.*exp(-b_i2.*p));

% Norm von alpha_i ermitteln
norm_ai = sqrt( sum( (P_line.*alpha) .* alpha ) );
norm_ai(norm_ai == 0) = delta;

%    % Hilfsvariablen

% Hilfsvariablen ( als schleife is es irgendwie schneller )
Li = alpha ./ norm_ai;
var0 = nTrans * (P_check .* Li);                                           % Hilfsvariable
chi = Q_i .* (2 - var0) .* (1 + a_chi * exp(b_chi * rm));
var0 = c_i .* r_i;
var1 = norm_ai./r_i;
var1( var1 > 1) = 1;
var2 = A .* n;

% schleife über alle Backstresstensoren
dalpha_dp = zeros(ntens,M);
for i = 1 : M
    
    dalpha_dp(:,i) = var0(i) * ( var2 - var1(i).^(chi(i)+1).* Li(:,i) );
    
end

% Ableitung gesamtbackstress
da_dp = sum(dalpha_dp,2);                                               % Ableitung: Backstresstensor nach plast. Bogenlänge


%-----------------------------------------------------------------------
%                  Ableitungen der Gedächtnissfläche
%-----------------------------------------------------------------------

L = a./norm_a;                                                         % Normierter Backstresstensor
dummy = sign(norm_a - rm);                                             % Hilfsvariable
hg = 0.5 * ( dummy + abs(dummy));                                      % Heaviside
dummy = sum( (P_line .* da_dp) .* L);                                   % Hilfsvariable
mla = 0.5 * ( dummy + abs(dummy) );                                    % maccaulay
dummy = 1- norm_a/rm;                                                  % Hilfsvariable
dummy = 0.5 * ( dummy + abs(dummy) );                                  % Hilfsvariable

% Ableitung Mem. surf
drm_dp = hg * mla - cm * dummy;

% Ableitung FF
dk_dp = k1 * ak * ck * exp(ck*rm) * drm_dp;                            % Ableitung: Fließfläche nach plast. Bogenlänge
% Abfangen von 0*Inf = NaN
dk_dp(isnan(dk_dp)) = 0;

%-----------------------------------------------------------------------
%                       Plastischer Modul
%-----------------------------------------------------------------------

h = nTrans * (P_check .* da_dp) + sqrt(2) * dk_dp;                        % Plastischer Modul

%-----------------------------------------------------------------------
%                Inkrement der plastischen Bogenlänge
%-----------------------------------------------------------------------

if ink_flag == 0
    
    dp = (nTrans * ink)/h;                                              % Plast. Inkrement bei Spannungssteuerung
    
elseif ink_flag == 1
    
    dp = ( nTrans * (C * ink)) / (h + (nTrans * (C * n)));              % Plast. Inkrement bei Dehnungssteuerung
    
end

%-----------------------------------------------------------------------
%                  Inkremente der Zustandsvariablen
%-----------------------------------------------------------------------

depsp = n * dp;                                                        % Inkrement der plast. Dehnung
dalpha = dalpha_dp * dp;                                               % Inkrement der Backstresstensoren
drm = drm_dp * dp;                                                     % Inkrement Radius der Gedächtnisfläche

%-----------------------------------------------------------------------
%                   Zusammenfassen der Inkremente
%-----------------------------------------------------------------------

if ink_flag == 0
    deps = D * ink + depsp;                                              % Inkrement der Dehnung berechnen
    out = vertcat(deps, depsp);                                        % Hilfsvariable
    
elseif ink_flag == 1
    dsig = C * (ink - depsp);                                          % Inkrement der Spannung berechnen
    out = vertcat(dsig, depsp);                                        % Hilfsvariable
    
end
dX = vertcat(out, reshape(dalpha,ntens*M,1), dp, drm);                                       % Inkremente der Zustandsvariablen

end % Ende sigma-tau

























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Materialgleichungen 1D                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dX = materialmodell1d(~, X, ink, ink_flag, parameter)
% Konkretes Materialmodell nach Jiang für ESZ
%
% INPUT:
%         t         -> Zeit (auch wenn nicht explizit gebraucht)
%         X         -> Zustand
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

    %-----------------------------------------------------------------------
    %            Identifikation Materialparameter                       
    %-----------------------------------------------------------------------
    ntens = 1;                                                             % Tensorkomponenten
    M = (length(parameter)-8)/7;                                           % Anzahl Backstresstensoren

    E = parameter(1);                                                      % Elastizitätsmodul
    
    a_chi = parameter(3);
    b_chi = parameter(4);
    ak = parameter(5);
    ck = parameter(6);
    k1 = parameter(7);
    cm = parameter(8);
    c_i0 = parameter(9:9+M-1);
    a_i1 = parameter(9+M:9+2*M-1);
    b_i1 = parameter(9+2*M:9+3*M-1);
    a_i2 = parameter(9+3*M:9+4*M-1);
    b_i2 = parameter(9+4*M:9+5*M-1);
    r_i = parameter(9+5*M:9+6*M-1);
    Q_i = parameter(9+6*M:9+7*M-1);
   
   
   
    %-----------------------------------------------------------------------
    %            statische Matrizen                       
    %-----------------------------------------------------------------------
%     C = elast_steifigkeit(E,nu,ntens,ndi);
    
    
    %----------------------------------------------------------------------
    %            Identifikation der Zustandsvariablen                        
    %----------------------------------------------------------------------

    if ink_flag == 0       
        eps = X(1);                                                        % Dehnung aus Zustandsvariable lesen                                
        epsp = X(2);                                                    % Plastischer Anteil der Dehnung aus Zustandsvariable lesen
        sig = E * (eps - epsp);                                         % Spannung berechnen

    elseif ink_flag == 1                                                 
        sig = X(1);                                                     % Spannung aus Zustandsvariable lesen
            
    end
       
    alpha =  reshape(X(3:(M+2)*ntens),ntens, M);          % akt. Backstresstensoren ermitteln

    %----------------------------------------------------------------------
    %               Radius der Fließfläche berechnen                         
    %----------------------------------------------------------------------
    delta = 1e-40;
    p = X((M+2)*ntens+1);                                           % plastische Bogenlänge
    rm = X((M+2)*ntens+2);                                          % Radius der Gedächtnisfläche
    if rm == 0                                                           
        rm = delta;                                                     % Rm = 0 abfangen
    end

    %----------------------------------------------------------------------
    %               Normale an die Fließfläche berechnen                      
    %----------------------------------------------------------------------  
    
    a =  sum(alpha,2);                                                      % Backstresstensor berechnen
    n =  sign(sig - a);                                                     % Normale zur Fließfläche

    norm_a = sqrt(2/3) * abs(a);                                            % Norm: Backstresstensor
    if norm_a == 0                                                   
        norm_a = delta;                                                    % Norm = 0 abfangen
    end

    %----------------------------------------------------------------------
    %              Ableitungen der Teilbackstresstensoren                    
    %----------------------------------------------------------------------

    c_i = c_i0 .* (ones(1,M) + a_i1.*exp(-b_i1.*p) ...              % Variable für zyklische Ver- bzw. Entfestigung
                                    + a_i2.*exp(-b_i2.*p));
                                 
    norm_ai = sqrt(2/3) * abs(alpha);                                       % Norm von alpha_i ermitteln
    for i = 1 : M   
       if norm_ai(i) == 0                                                     
           norm_ai(i) = delta;                                              % Abfangen teilen durch null   
       end
        
    end
    
    Li =  sign(alpha);
    var0 = ones(1,M) * sign( sig - a) .* sign(alpha);                      % Hilfsvariable
    var3 = norm_ai ./ r_i;
    var3(var3 > 1) = 1;
    chi = Q_i .* (2 - var0) .* (1 + a_chi * exp(b_chi * rm));               % Exponent für Nichtproportionalität und Spannungsleveleffekt auf Rachetting
    dalpha_dp = sqrt(6)/2 *c_i .* r_i .*...                                 
               (n - (var3).^ (chi+1) .* Li);
   
    
    da_dp = sum(dalpha_dp,2);                                               % Ableitung: Backstresstensor nach plast. Bogenlänge 

    %----------------------------------------------------------------------
    %                  Ableitungen der Gedächtnissfläche                     
    %----------------------------------------------------------------------
    
    L =  sign(a);                                                           % Normierter Backstresstensor
    hg = heaviside( (norm_a - rm) );                                        % Hilfsvariable
    mla = macaulay(sqrt(6)/3 * L * da_dp );                                 % Hilsfvariable
    drm_dp = hg * mla - cm * macaulay(1- norm_a/rm);                        % Ableitung: Gedächtnisfläche nach plast. Bogenlänge

    dk_dp = k1 * ak * ck * exp(ck*rm) * drm_dp;                             % Ableitung: Fließfläche nach plast. Bogenlänge
    % Abfangen von 0*Inf = NaN
    dk_dp(isnan(dk_dp)) = 0; 
    
    %----------------------------------------------------------------------
    %                       Plastischer Modul                        
    %----------------------------------------------------------------------

    h = sqrt(6)/3 * n * da_dp + sqrt(2) * dk_dp;                            % Plastischer Modul

    %----------------------------------------------------------------------
    %                Inkrement der plastischen Bogenlänge                     
    %----------------------------------------------------------------------
    
    if ink_flag == 0 
                    
        dp = (sqrt(6)/3 * n * ink)/ h;                                  % Plast. Inkrement bei Spannungssteuerung   
            
    elseif ink_flag == 1
            
        dp = (sqrt(6)/3 * n * E * ink) / (h + 2/3  * E) ;            % Plast. Inkrement bei Dehnungssteuerung
        
    end
    
    %----------------------------------------------------------------------
    %                  Inkremente der Zustandsvariablen                      
    %----------------------------------------------------------------------
    
    depsp = sqrt(6)/3 * n * dp;                                             % Inkrement der plast. Dehnung
    dalpha = dalpha_dp * dp;                                                % Inkrement der Backstresstensoren
    drm = drm_dp * dp;                                                      % Inkrement Radius der Gedächtnisfläche

    %----------------------------------------------------------------------
    %                   Zusammenfassen der Inkremente                        
    %----------------------------------------------------------------------

    if ink_flag == 0                                                        	
        out = ink/E + depsp;                                                % Inkrement der Dehnung berechnen
        
    elseif ink_flag == 1                                                    
        out = E * (ink - depsp);                                            % Inkrement der Spannung berechnen

    end
    tmp = reshape(dalpha,ntens*M,1);                                % Umsortieren von dalpha zu einem Vektor der Größe ntens*num_alpha
    dX = vertcat(out, depsp, tmp, dp, drm);                                 % Inkremente der Zustandsvariablen    

end % Ende 1D


