function [Lnew,verl] = smoothSignal(L,maxGrad,varargin)
% Glätte gefiltertes Signal.
% Vorgehen aus Guide to load analysis for durability in vehicle 
% engineering (2014, Wiley) Kap 4.3.1.3:
%    zwischen Zwei monotonen Punkten:
%       - Einfügen Geraden
%    zwischen monotonen & Umkehrpunkten
%       - Einfügen viertel Sinus
%    zwischen zwei Umkehrpunkten
%       - Einfügen halber Sinus
%
% QUELLE:
%
% 
% INPUT:
% L        - Lastsignal (nkana x ndata)
% maxGrad  - maximaler Gradient den Signal haben Soll
% varargin - variabler Input, Zeitpunkte
%
% OUTPUT:
% L      - Verlängertes Lastsignal
% verl   - Verlängerung des Signals 
%
% -------------------------------------------------------------------------


% Init
nkana = size(L,1);     % Anzahl Lastkanäle
ndata = size(L,2);     % Anzahl Datenpunkte
Ndata = ndata;         % Anzahl Datenpunkte nachdem Werte eigefügt sind
if nargin == 3
    t = varargin{1,1};
    dt = diff(t);
else
    dt = ones(1,Ndata-1);                % Sampling Rate (auf Standartwert, vllt später noch 
                           % als Input
end
% LZW = cell(1,Ndata-2); % Speicher für zwischenwerte
Lnew = L(:,1);             % neues L mit ZW !!!! Schlechte Lösung, vllt noch verbessern
% Erkenne Umkehrpunkte
dL = diff(L,[],2);
UKP = dL(:,1:end-1).*dL(:,2:end) < 0;
UKP = [true(nkana,1),UKP,true(nkana,1)]; % Ersten und letzten Punkt immer behalten

% Schleife über alle Punkte
for i = 2:ndata-1
    % Übergang von L(:,i-1) -> L(:,i)
    % Anzahl einzufügender Punkte
    if any(UKP(:,i-1)) || any(UKP(:,i)) % Mindestens einer der Punkte ist UKP
        n = max(ceil( abs(dL(:,i-1)*pi)/(2*maxGrad*dt(i-1)) - 1 ));
    else % Monotoner Part
        n = max(ceil( abs(dL(:,i-1))/(maxGrad*dt(i-1)) - 1 ));
    end
    % Inkrementiere Ndata
    Ndata = Ndata + n;
    % Schleife über Lastkanäle
    Lzw = zeros(nkana,n);                % Speicher Zwischenwerte
    for j = 1:nkana
        if UKP(j,i-1) && UKP(j,i)        % Beides UKP
            Lzw(j,:) = einfuegenHalberSinus(L(j,i-1),L(j,i),n);
        elseif UKP(j,i-1)                % Nur i-1 UKP
            Lzw(j,:)  = einfuegenViertelSinus(L(j,i-1),L(j,i),n,pi,L(j,i));
        elseif UKP(j,i)                  % Nur i UKP
            Lzw(j,:)  = einfuegenViertelSinus(L(j,i-1),L(j,i),n,3*pi/2,L(j,i-1));
        else                             % Keine UKP
            Lzw(j,:)  = einfuegenLinie(L(j,i-1),L(j,i),n);
        end
    end
    % Speicher Zwischenwerte
%     LZW{i-1} = Lzw;
    Lnew = [Lnew,Lzw,L(:,i)];
end % Ende Schleife über alle Datenpunkte

% Zusammenbauen

verl = Ndata/ndata;
end % Ende Hauptfunktion


% -------------------------------------------------------------------------
% Hilfsfunktionen Zwischenwerte einfügen
function L = einfuegenLinie(L0,L1,n)
    n = n + 1;
    t = 1/n:1/n:1-1/n;
    L = (L1-L0) .* t + L0;
end

function L = einfuegenHalberSinus(L0,L1,n)
    amp = (L1 - L0)/2;
    mean = (L1 + L0)/2;
    n = n + 1;
    t = 1/n:1/n:1-1/n;
    L = amp * cos(pi.*t+pi) + mean;
end

function L = einfuegenViertelSinus(L0,L1,n,pha,mean)
    amp = (L1 - L0);
%     mean = L1;
    n = n + 1;
    t = 1/n:1/n:1-1/n;
    L = amp * cos(pi/2.*t+pha) + mean;
end