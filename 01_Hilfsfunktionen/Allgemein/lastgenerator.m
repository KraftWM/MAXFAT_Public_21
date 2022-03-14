function load = lastgenerator(typ,amp,pha,fre,mit,ndp)
% Funktion generiert Lastfolge fuer die Zeit t = [0,1]
% 
% INPUT:
% typ    - 1 = Sinus, 2 = Dreieck, 3 = Trapez
% amp    - Amplitude
% pha    - Phasenverschiebung (in Grad)
% fre    - Frequenz
% mit    - Mittelwert
% ndp    - Anzahl Datenpunkte (insgesamt)
% 
% OUTPUT:
% load   - Lastfolge
% -------------------------------------------------------------------------

if typ == 2 % Dreieck
    % https://en.wikipedia.org/wiki/Triangle_wave
    % ... dummy Zeit
    t = 0:1/ndp:1-1/ndp;
    % ... Wellenlaenge
    T = 1/fre;
    % ... Phasenverschiebung normieren
    pha = pha/360*T;
    % ... Last
    load = 4*amp/T *( abs( mod((t-T/4+pha),T) - T/2) - T/4)+mit;
elseif typ == 3 % Trapez
    % ... dummy Zeit
    t = 0:1/ndp:1-1/ndp;
    % ... Wellenlaenge
    T = 1/fre;
    % ... Phasenverschiebung normieren
    pha = pha/360*T;
    % ... Last (Dreieck mit doppelter Amplitude)
    load = 8*amp/T *( abs( mod((t-T/4+pha),T) - T/2) - T/4)+ mit;
    % ... Abschneiden
    load(load > amp+mit) = amp+mit;
    load(load < mit-amp) = mit-amp;    
else % Sinus
    % ... Umrechnen der Phasenverschiebung 
    pha = pha * pi/180;
    % ... dummy Zeit
    T = 0:1/ndp:1-1/ndp;
    % ... Last
    load = amp * sin(2*pi*fre*T+pha) + mit;
end

end % Ende Funktion