%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Matlab Coder Version von:                                            %
%    Modellgleichung 3D Jiang Modell                                      %
%                                                                         %
%    Aufgerufen in:                                                       %
%    jiang.m                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dX = CoderVersion_3D_Jiang(~, X, M, ink, ink_flag, parameter, C, D, P,...
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
%     M = (length(parameter)-8)/7;                                          % Anzahl Backstresstensoren 
    
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
   else
        error('Faslche Steuerung angegeben. (0) Spansteu. oder (1) Dehnsteu.')
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
   Li = zeros(ntens,M);
   var1 = zeros(1,M);
   for kk = 1 : M
       Li(:,kk) = alpha(:,kk)/norm_ai(kk);
       var1(kk) = norm_ai(kk)/r_i(kk);
   end
   var0 = nTrans * Li;                                                     % Hilfsvariable
   chi = Q_i .* (2 - var0) .* (1 + a_chi * exp(b_chi * rm));               % Exponent für Nichtproportionalität und Spannungsleveleffekt auf Rachetting
   var0 = c_i .* r_i;
   var2 = P_hat * n;                                                       % Hilfsvariable
   var1(var1 > 1) = 1;
   
   % schleife über alle Backstresstensoren
   dalpha_dp = zeros(ntens,M);
   for i = 1 : M
       
       dalpha_dp(:,i) = var0(i) * ( var2 - var1(i).^(chi(i)+1).* Li(:,i) );
       
   end
   
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
        
        else
            
            error('Faslche Steuerung angegeben. (0) Spansteu. oder (1) Dehnsteu.')

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
        
    else
        
        error('Faslche Steuerung angegeben. (0) Spansteu. oder (1) Dehnsteu.')

    end

    tmp = reshape(dalpha,ntens*M,1);                                       % Umsortieren von dalpha zu einem Vektor der Größe ntens*num_alpha
    dX = vertcat(out, tmp, dp, drm);                                       % Inkremente der Zustandsvariablen    

end % Ende 3D