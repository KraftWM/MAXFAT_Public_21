function [ZVARneu] = SpanIntDoring(ntens,ndi,ZVAR0,para,dSIG,dEPS0)
% Spannungsgesteuerte Integration des Döring Modells
% Wegen der Möglichkeit des Entfestigenden Materialverhaltens werden
% Spannungsinkremente mit den Döringmodell iterirt, das Materialmodell wird
% dabei dehnungsgesteuert integriert. Die Iteration des
% Spannungsinkrements erfolgt mit den Newton verfahren
%
%
% INPUT:
% ntens       - Anzahl tensorkomponenten
% ndi         - Anzahl Hauptdiagonalelemente
% ZVAR0       - Zustandsvariablen des (dehnungsgesteuerten) Döring-Modells
% para        - Materialparameter
% dSIG        - Aufzubringenden Spannungsinkrement
% dEPS0       - Startwert Dehnungsinkrement für Newton Iteration
% 
% OUTPUT
% ZVARneu     - Geupdatete Zustandsvariablen Döring Modell
%__________________________________________________________________________

% ... Anfangszustand
SIG0 = ZVAR0(1:ntens);
dEPS = dEPS0;
% ... Endzustand
SIG1 = SIG0 + dSIG;
% ... Nullter Iterationsschritt
[ZVARneu, ~,DEP] = doring(ntens,ndi,dEPS,ZVAR0,1,para);
SIG = ZVARneu(1:ntens);
% ... Fehler f = SIG1 - SIG(dEPS)
f = SIG1 - SIG;
fnorm2 = sum(f.*f);
% ... Newton Iteration
iter = 0;         % Iterationszähler
maxiter = 100;    % Abbruchbedingung keine Konvergenz
tol = 1e-6;       % Toleranz bei 3ter Nachkommastelle in Spannungen
while fnorm2 > tol
    % ... neues Dehnungsinkrement dEPS = dEPS + CEP^-1*f
    dEPS = dEPS + DEP*f;
    % ... Integration mit Neuem Dehnungsinkrement
    [ZVARneu,~,DEP] = doring(ntens,ndi,dEPS,ZVAR0,1,para);
    SIG = ZVARneu(1:ntens);
    % ... Fehler f = SIG1 - SIG(dEPS)
    f = SIG1 - SIG;
    fnorm2 = sum(f.*f);
    % ... keine Endlosschleifen
    iter = iter + 1;
    if iter > maxiter
        msg = ['Keine Konvergenz im Newtonverfahren für',...
               ' spannungsgesteurte Integration de Döring Modells'];
        error(msg);
    end
end % Ende Iterationsschleife
end % Ende Funktion