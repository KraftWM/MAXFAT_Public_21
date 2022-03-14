function [esig,eepsp] = uniax_neuberstern(E,sig,epsp,Kprime,nprime,Kp)
% einachsiger Neuber mit Sternterm aus Richtlinie Nichtlinear
% !!!!!  zur Ermittlung von BAUTEILFLIEßKURVEN !!!!!
%
% INPUT:
% E             -> Emodul
% sig           -> ELATISCH PLASTISCHE Spannungen
% epsp          -> PLASTISCHE Dehnungen
% Kprime,nprime -> Parameter Ramberg Osgood
% Kp            -> Traglastformzahl
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
maxiter = 100;
tol = 1e-8;

% Schleife über alle Stützstellen
for ii = 1 : length(sig)
    
    % aktuelle elastische Spannung (Startwerte aus Neuber)
    L = esig(ii);
    
    % aktuelle totale Verzerrungsenergiedichte
    omega = sig(ii)*(sig(ii)/E + epsp(ii));
    
    % Zielfunktion f = 0
    f = L^2/E + L * Kp * (L/(Kp*Kprime))^(1/nprime) - omega;
    
    % Starte Newton iteration
    iter = 0;
    while abs(f) > tol
        
        % keine endlosschleife
        if iter > maxiter
            msg = 'keine Konvergenz in Seeger/Heuler';
            error(msg)
        end
        
        % Ableitung
        dfdL = 2*L/E + Kp * (L/(Kp*Kprime))^(1/nprime) * ( 1 + 1/nprime);
        
        % neue elastische Spannung
        L = L - f/dfdL;
        
        % Zielfunktion f = 0
        f = L^2/E + L * Kp * (L/(Kp*Kprime))^(1/nprime) - omega;
        
        % Inkrement iterationszähler
        iter = iter + 1;
        
    end % Ende Schleife Newtoniteration
    
    % Abspeichern elastische Spannung
    esig(ii) = L;
    % Abspeichern zugehörige elastische Dehnung
    eeps(ii) = Kp * (L/(Kp*E) + (L/(Kp*Kprime))^(1/nprime));
    
end % Ende Schleife über Stützstellen

% pseudo plastische dehnung
eepsp = eeps - sig/E;

end % Ende Funktion
