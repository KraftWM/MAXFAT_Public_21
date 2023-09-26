function [esig,eepsp] = uniax_seegerbeste(E,sig,epsp,Kprime,nprime,Kp)
% einachsiger Seeger Beste
%  !!!!!  zur ermittlung von BAUTEILFLIEßKURVEN !!!!!
%
% INPUT:
% E    -> Emodul
% sig  -> ELATISCH PLASTISCHE Spannungen
% epsp -> PLASTISCHE Dehnungen
% Kp   -> Traglastformzahl
%
% OUTPUT:
%  esig   -> pseudo elastische Spannung für esig-approach
%  eepsp  -> pseudo plastische Dehnung für eeps-approach
%
% -------------------------------------------------------------------------

% Korrigiere Kp = 1
if Kp == 1
    Kp = 1.0001;
elseif Kp < 1
    Kp = 1.0001;
end 

% Neuber Regel zum erzeigen von Schätzwerten
esig = sqrt(sig .* (sig + E .* epsp));
eeps = esig/E;

% parameter für newton itreation
maxiter = 1000;
tol = 1e-7;
dL = 1e-5; % Schrittweite für zentrale Differenzen
fiter = zeros(1,maxiter+1); 

% Schleife über alle Stützstellen
for ii = 1 : length(sig)
    
    % aktuelle elastische Spannung (Startwerte aus Neuber)
    L = esig(ii);
    
    % Fange elastischen Bereich ab
    if epsp(ii) == 0      
        % Spannung ist Gleich 
        esig(ii) = sig(ii);
        % Abspeichern zugehörige Dehnung
        eeps(ii) = esig(ii)/E;
        
    % Elastisch - plastische Iteration 
    else
        % aktuelle totale Verzerrungsenergiedichte
        omega = sig(ii)*(sig(ii)/E + epsp(ii));
        %     omega = sig(ii)*(sig(ii)/E + (sig(ii)/Kprime)^(1/nprime));
        
        % Zielfunktion
        [f,d] = zielfun(L,sig(ii),omega,E,Kprime,nprime,Kp);
        
        
        % Starte Newton iteration
        iter = 0;
        while abs(f) > tol
            
            % Debugg Zielfunktion
            fiter(iter+1) = f;
            
            % keine endlosschleife
            if iter > maxiter
                msg = 'keine Konvergenz in Seeger/Beste';
                error(msg)
            end
            
            % Ableitung Zentrale Differenzen
            fp = zielfun(L+dL,sig(ii),omega,E,Kprime,nprime,Kp);
            fm = zielfun(L-dL,sig(ii),omega,E,Kprime,nprime,Kp);
            dfdL = (fp - fm)/(2*dL);
            
            % neue elastische Spannung
            L = L - f/dfdL;
            
            % Zielfunktion
            [f,d] = zielfun(L,sig(ii),omega,E,Kprime,nprime,Kp);
            
            % Inkrement iterationszähler
            iter = iter + 1;
            
        end % Ende Iteration
        
        % Abspeichern elastische Spannung
        esig(ii) = L;
        % Abspeichern zugehörige elastische Dehnung
        eeps(ii) = esig(ii)/E;%d * Kp * (esig(ii)/(Kp*E) + (esig(ii)/(Kp*Kprime))^(1/nprime));
        
    end % Ende Verzweigung elastische/elastisch plastisch
    
end % Ende Schleife über Inkremente

% pseudo plastische dehnung
eepsp = eeps - sig/E;
eepsp(eepsp<0) = 1e-40;

end % Ende Hauptfunktion

% Hilfsfunktion zielfunktion für seeger beste
function [f,d,u] = zielfun(L,sig,omega,E,Kprime,nprime,Kp)

    estar = L/(Kp*E) + (L/(Kp*Kprime))^(1/nprime);                         % Hilfsvariable
    fak = L/sig;                                                           % Hilfsvariable
    u = pi/2 * ( (fak-1)/(Kp-1) );                                         % Hilfsvariable
    u = min(u,pi/2);
    u = max(u,0);
    if u <= 1e-3
         d = 1 + u^2/6 + 2*u^3/45 + (1/fak)^2 - 1/fak;
    else
        d = 2/u^2 * log(1/cos(u)) + (1/fak)^2 - 1/fak;                     % Hilfsvariable
    end
%     d =  min(2/u^2,1e40) * abs(log(1/cos(u)))  + fak^2 - fak;
    f = d*L*Kp*estar - omega;
    
end % Ende Ableitung
