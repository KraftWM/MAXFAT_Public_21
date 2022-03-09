function signal = sinus_zwischenwerte(UKP,nzw)
% Funktion füllt eine Umkehrpunktfolge mit zwischenwerten mit Sinus auf
%
% INPUT:
% UKP     - Umkehrpunktfolge (nkana,nukp)
% nzw     - Anzahl zwischenwerte die zwischen jeden UKP eingefügt werden
% 
%
% OUTPUT:
% signal  - UKP mit zwischenwerten
%
% -------------------------------------------------------------------------

% Größe des Signals
nkana = size(UKP,1);
nukp = size(UKP,2);

% Speicher fürs Signal
signal = NaN(nkana,nukp+(nukp-1)*nzw);
signal(:,1) = UKP(:,1);
zsignal = 1;

% Zwischenwerte Hochsetzten
nzw = nzw + 1;

% Schleife über alle UKP
for i = 2:nukp
    
    for j = 1:nkana
       
        if UKP(j,i) < UKP(j,i-1)
            t = 0.25+0.5/nzw:0.5/nzw:0.75;
            X0 = (UKP(j,i) + UKP(j,i-1))/2;
            XA = X0 - UKP(j,i);
            X = X0 + XA * sin(2*pi*t);
        else
            t = 0.5+0.5/nzw:0.5/nzw:1;
            X0 = (UKP(j,i) + UKP(j,i-1))/2;
            XA = UKP(j,i) - X0;
            X = X0 + XA * cos(2*pi*t);
        end
        
        signal(j,zsignal + 1 : zsignal + nzw) = X;

    end
    
    zsignal = zsignal + nzw;
end
end