function W = akima2d(LX,LY,X,Y,Z,N,U,V)
% Interpolation in x-y Ebene
%
% Quelle:
% Akima 1974 - A Method of Bivariate Interpolation and Smooth Surface
% Fitting Based on Local Procedures
%
% 
% !!! Interpolierte Werte müssen im Netz liegen. Es wird nicht extrapoliert
% 
%
% INPUT:
%  LX,LY   - (int) anzahl Datenpunkte in X oder Y
%   X,Y    - (double array 1 x nX) Datenpunkte des Netzes in x-y Ebene
%    Z     - (double array nX x nY) Werte der Funktion Z(i,j)=z(X(i),Y(j))
%  N       - (int) anzahl Datenpunkte für ausgabe
%  U,V     - (double array 1 x no) Koordinaten der Punkte für die die Werte
%            interpoliert werden sollen.
%
% OUPUT:
%  W       - (double array 1 x no) interpolierte Werte an den Koordinaten 
%            W(i) = z(x(i),y(i))
%__________________________________________________________________________

% ... Initialisieren Speicher
ZA = zeros(5,2);     % Dividierte Differenzen von Z in Richtung x
ZB = zeros(2,5);     % Dividierte Differenzen von Z in Richtung y
ZAB = zeros(3,3);    % Dividierte Differenzen 2.Ordnung von Z in Richtung x und y
ZX = zeros(4,4);     % Partielle Ableitung von Z in Richtung x
ZY = zeros(4,4);     % Partielle Ableitung von Z in Richtung y
ZXY = zeros(4,4);    % Partielle Ableitung 2. Ordnung von Z in Richtung x und y
W = zeros(1,N);      % Output Werte

% ... Initialisieren innerer Variablen
LX0 = LX;
LXM1 = LX0 - 1;
LXM2 = LXM1 - 1;
LXP1 = LX0 + 1;
LY0 = LY;
LYM1 = LY0 - 1;
LYM2 = LYM1 - 1;
LYP1 = LY0 + 1;
N0 = N;

% ... Error Checks
if LXM2 < 0
    msg = 'Zu wenige Stützpunkte in X. Mindestens 2 benötigt.';
    error(msg);
end
if LYM2 < 0
    msg = 'Zu wenige Stützpunkte in Y. Mindestens 2 benötigt.';
    error(msg);
end

% ... Hauptschleif über alle Paare (x(k),y(k))
for K = 1 : N0
    UK = U(K);
    VK = V(K);
    % ... Finde Punkt in X für den gilt
    % UK >= X(IX-1) $$ UK < X(IX)
    if LXM2 == 0
        IX = 2;
    elseif UK >= X(LX0)
        IX = LXP1;
    elseif UK < X(1)
        IX = 1;
    else
        IX = 2;
        while UK < X(IX-1) || UK >= X(IX)
            IX = IX + 1;
        end
    end
    % ... Finde Punkt in Y für den gilt
    % VK >= Y(IY-1) $$ VK< Y(IY)
    % UK >= X(IX-1) $$ UK < X(IX)
    if LYM2 == 0
        IY = 2;
    elseif VK >= Y(LY0)
        IY = LYP1;
    elseif VK < Y(1)
        IY = 1;
    else
        IY = 2;
        while VK < Y(IY-1) || VK >= Y(IY)
            IY = IY + 1;
        end
    end
    % ... Rauslesen Werte für dividierte differenzen
    JX = IX;
    if JX == 1
        JX = 2;
    elseif JX == LXP1
        JX = LX0;
    end
    JY = IY;
    if JY == 1
        JY = 2;
    elseif JY == LYP1
        JY = LY0;
    end
    JXM2 = JX - 2;
    JXML = JX - LX0;
    JYM2 = JY - 2;
    JYML = JY - LY0;
    % ... ... im Reckteck, dass den Punkt enthält
    X3 = X(JX-1);
    X4 = X(JX);
    A3 = 1.0/(X4-X3);
    Y3 = Y(JY-1);
    Y4 = Y(JY);
    B3 = 1.0/(Y4-Y3);
    Z33 = Z(JX-1,JY-1);
    Z43 = Z(JX,JY-1);
    Z34 = Z(JX-1,JY);
    Z44 = Z(JX,JY);
    Z3A3 = (Z43-Z33)*A3;
    Z4A3 = (Z44-Z34)*A3;
    Z3B3 = (Z34-Z33)*B3;
    Z4B3 = (Z44-Z43)*B3;
    ZA3B3 = (Z4B3-Z3B3)*A3;
    % ... ... in x Richtung
    if LXM2 == 0
        Z3A2 = Z3A3;
        Z4A2 = Z4A3;
        Z3A4 = Z3A3 + Z3A3 - Z3A2;
        Z4A4 = Z4A3 + Z4A3 - Z4A2;
        ZA2B3 = (Z4A2-Z3A2)*B3;
        ZA4B3 = (Z4A4-Z3A4)*B3;
        Z3A1 = Z3A2 + Z3A2 - Z3A3;
        Z4A1 = Z4A2 + Z4A2 - Z4A3;
        Z3A5 = Z3A4 + Z3A4 - Z3A3;
        Z4A5 = Z4A4 + Z4A4 - Z4A3;
    elseif JXM2 == 0
        X5 = X(JX+1);
        A4 = 1.0/(X5-X4);
        Z53 = Z(JX+1,JY-1);
        Z54 = Z(JX+1,JY);
        Z3A4 = (Z53-Z43)*A4;
        Z4A4 = (Z54-Z44)*A4;
        Z3A2 = Z3A3 + Z3A3 - Z3A4;
        Z4A2 = Z4A3 + Z4A3 - Z4A4;
        ZA2B3 = (Z4A2-Z3A2)*B3;
        ZA4B3 = (Z4A4-Z3A4)*B3;
        Z3A1 = Z3A2 + Z3A2 - Z3A3;
        Z4A1 = Z4A2 + Z4A2 - Z4A3;
        if JX >= LXM1
            Z3A5 = Z3A4 + Z3A4 - Z3A3;
            Z4A5 = Z4A4 + Z4A4 - Z4A3;
        else
            A5 = 1.0/(X(JX+2)-X5);
            Z3A5 = (Z(JX+2,JY-1)-Z53)*A5;
            Z4A5 = (Z(JX+2,JY)-Z54)*A5;
        end
    else
        X2 = X(JX-2);
        A2 = 1.0/(X3-X2);
        Z23 = Z(JX-2,JY-1);
        Z24 = Z(JX-2,JY);
        Z3A2 = (Z33-Z23)*A2;
        Z4A2 = (Z34-Z24)*A2;
        if JXML == 0
            Z3A4 = Z3A3 + Z3A3 - Z3A2;
            Z4A4 = Z4A3 + Z4A3 - Z4A2;
        else
            X5 = X(JX+1);
            A4 = 1.0/(X5-X4);
            Z53 = Z(JX+1,JY-1);
            Z54 = Z(JX+1,JY);
            Z3A4 = (Z53-Z43)*A4;
            Z4A4 = (Z54-Z44)*A4;
        end
        ZA2B3 = (Z4A2-Z3A2)*B3;
        ZA4B3 = (Z4A4-Z3A4)*B3;
        if JX <= 3
            Z3A1 = Z3A2 + Z3A2 - Z3A3;
            Z4A1 = Z4A2 + Z4A2 - Z4A3;
        else
            A1 = 1.0/(X2-X(JX-3));
            Z3A1 = (Z23-Z(JX-3,JY-1))*A1;
            Z4A1 = (Z24-Z(JX-3,JY))*A1;
        end
        if JX >= LXM1
            Z3A5 = Z3A4 + Z3A4 - Z3A3;
            Z4A5 = Z4A4 + Z4A4 - Z4A3;
        else
            A5 = 1.0/(X(JX+2)-X5);
            Z3A5 = (Z(JX+2,JY-1)-Z53)*A5;
            Z4A5 = (Z(JX+2,JY)-Z54)*A5;
        end
    end % Ende Werte X
    % ... ... in x Richtung
    if LYM2 == 0
        Z3B2 = Z3B3;
        Z4B2 = Z4B3;
        Z3B4 = Z3B3 + Z3B3 - Z3B2;
        Z4B4 = Z4B3 + Z4B3 - Z4B2;
        ZA3B2 = (Z4B2-Z3B2)*A3;
        ZA3B4 = (Z4B4-Z3B4)*A3;
        Z3B1 = Z3B2 + Z3B2 - Z3B3;
        Z4B1 = Z4B2 + Z4B2 - Z4B3;
        Z3B5 = Z3B4 + Z3B4 - Z3B3;
        Z4B5 = Z4B4 + Z4B4 - Z4B3;
    elseif JYM2 == 0
        Y5 = Y(JY+1);
        B4 = 1.0/(Y5-Y4);
        Z35 = Z(JX-1,JY+1);
        Z45 = Z(JX,JY+1);
        Z3B4 = (Z35-Z34)*B4;
        Z4B4 = (Z45-Z44)*B4;
        Z3B2 = Z3B3 + Z3B3 - Z3B4;
        Z4B2 = Z4B3 + Z4B3 - Z4B4;
        ZA3B2 = (Z4B2-Z3B2)*A3;
        ZA3B4 = (Z4B4-Z3B4)*A3;
        Z3B1 = Z3B2 + Z3B2 - Z3B3;
        Z4B1 = Z4B2 + Z4B2 - Z4B3;
        if JY >= LYM1
            Z3B5 = Z3B4 + Z3B4 - Z3B3;
            Z4B5 = Z4B4 + Z4B4 - Z4B3;
        else
            B5 = 1.0/(Y(JY+2)-Y5);
            Z3B5 = (Z(JX-1,JY+2)-Z35)*B5;
            Z4B5 = (Z(JX,JY+2)-Z45)*B5;
        end
    else
        Y2 = Y(JY-2);
        B2 = 1.0/(Y3-Y2);
        Z32 = Z(JX-1,JY-2);
        Z42 = Z(JX,JY-2);
        Z3B2 = (Z33-Z32)*B2;
        Z4B2 = (Z43-Z42)*B2;
        if JYML == 0
            Z3B4 = Z3B3 + Z3B3 - Z3B2;
            Z4B4 = Z4B3 + Z4B3 - Z4B2;
        else
            Y5 = Y(JY+1);
            B4 = 1.0/(Y5-Y4);
            Z35 = Z(JX-1,JY+1);
            Z45 = Z(JX,JY+1);
            Z3B4 = (Z35-Z34)*B4;
            Z4B4 = (Z45-Z44)*B4;
        end
        ZA3B2 = (Z4B2-Z3B2)*A3;
        ZA3B4 = (Z4B4-Z3B4)*A3;
        if JY <= 3
            Z3B1 = Z3B2 + Z3B2 - Z3B3;
            Z4B1 = Z4B2 + Z4B2 - Z4B3;
        else
            B1 = 1.0/(Y2-Y(JY-3));
            Z3B1 = (Z32-Z(JX-1,JY-3))*B1;
            Z4B1 = (Z42-Z(JX,JY-3))*B1;
        end
        if JY >= LYM1
            Z3B5 = Z3B4 + Z3B4 - Z3B3;
            Z4B5 = Z4B4 + Z4B4 - Z4B3;
        else
            B5 = 1.0/(Y(JY+2)-Y5);
            Z3B5 = (Z(JX-1,JY+2)-Z35)*B5;
            Z4B5 = (Z(JX,JY+2)-Z45)*B5;
        end
    end % Ende Werte in y
    % ... ... in diagonale Richtung
    if LXM2 == 0
        ZA2B2 = ZA3B2;
        ZA4B2 = ZA3B2;
        ZA2B4 = ZA3B4;
        ZA4B4 = ZA3B4;
    elseif LYM2 == 0
        ZA2B2 = ZA2B3;
        ZA2B4 = ZA2B3;
        ZA4B2 = ZA4B3;
        ZA4B4 = ZA4B3;
    elseif JXML == 0
        if JYM2 == 0
            ZA2B4 = (Z3B4-(Z(JX-2,JY+1)-Z24)*B4)*A2;
            ZA2B2 = ZA2B3 + ZA2B3 - ZA2B4;
        else
            ZA2B2 = (Z3B2-(Z23-Z(JX-2,JY-2))*B2)*A2;
            if JYML == 0
                ZA2B4 = ZA2B3 + ZA2B3 - ZA2B2;
            else
                ZA2B4 = (Z3B4-(Z(JX-2,JY+1)-Z24)*B4)*A2;
            end
        end
        ZA4B2 = ZA3B2 + ZA3B2 - ZA2B2;
        ZA4B4 = ZA3B4 + ZA3B4 - ZA2B4;
    elseif JYM2 == 0
        ZA4B4 = ((Z(JX+1,JY+1)-Z54)*B4-Z4B4)*A4;
        ZA4B2 = ZA4B3 + ZA4B3 - ZA4B4;
        if JXM2 == 0
            ZA2B2 = ZA3B2 + ZA3B2 - ZA4B2;
            ZA2B4 = ZA3B4 + ZA3B4 - ZA4B4;
        else
            ZA2B4 = (Z3B4-(Z(JX-2,JY+1)-Z24)*B4)*A2;
            ZA2B2 = ZA2B3 + ZA2B3 - ZA2B4;
        end
    else
        ZA4B2 = ((Z53-Z(JX+1,JY-2))*B2-Z4B2)*A4;
        if JYML == 0
            ZA4B4 = ZA4B3 + ZA4B3 - ZA4B2;
        else
            ZA4B4 = ((Z(JX+1,JY+1)-Z54)*B4-Z4B4)*A4;
        end
        if JXM2 == 0
            ZA2B2 = ZA3B2 + ZA3B2 - ZA4B2;
            ZA2B4 = ZA3B4 + ZA3B4 - ZA4B4;
        else
            ZA2B2 = (Z3B2-(Z23-Z(JX-2,JY-2))*B2)*A2;
            if JYML == 0
                ZA2B4 = ZA2B3 + ZA2B3 - ZA2B2;
            else
                ZA2B4 = (Z3B4-(Z(JX-2,JY+1)-Z24)*B4)*A2;
            end
        end
    end % Ende Diagonale Richtung
    % ... Array (ZA,ZB,ZAB) Elemente Füllen
    ZA(1) = Z3A1;
    ZA(2) = Z3A2;
    ZA(3) = Z3A3;
    ZA(4) = Z3A4;
    ZA(5) = Z3A5;
    ZA(6) = Z4A1;
    ZA(7) = Z4A2;
    ZA(8) = Z4A3;
    ZA(9) = Z4A4;
    ZA(10) = Z4A5;
    
    ZB(1) = Z3B1;
    ZB(2) = Z4B1;
    ZB(3) = Z3B2;
    ZB(4) = Z4B2;
    ZB(5) = Z3B3;
    ZB(6) = Z4B3;
    ZB(7) = Z3B4;
    ZB(8) = Z4B4;
    ZB(9) = Z3B5;
    ZB(10) = Z4B5;
    
    ZAB(1) = ZA2B2;
    ZAB(2) = ZA3B2;
    ZAB(3) = ZA4B2;
    ZAB(4) = ZA2B3;
    ZAB(5) = ZA3B3;
    ZAB(6) = ZA4B3;
    ZAB(7) = ZA2B4;
    ZAB(8) = ZA3B4;
    ZAB(9) = ZA4B4;
    
    % ... partielle Ableitung als gewichtete Mittelwerte der dividierten
    % differenzen
    for JY = 2 : 3
        for JX = 2 : 3
            W2 = abs(ZA(JX+2,JY-1)-ZA(JX+1,JY-1));
            W3 = abs(ZA(JX,JY-1)-ZA(JX-1,JY-1));
            SW = W2 + W3;
            if SW == 0
                WX2 = 0.5;
                WX3 = 0.5;
            else
                WX2 = W2/SW;
                WX3 = W3/SW;
            end
            ZX(JX,JY) = WX2*ZA(JX,JY-1) + WX3*ZA(JX+1,JY-1);
            W2 = abs(ZB(JX-1,JY+2)-ZB(JX-1,JY+1));
            W3 = abs(ZB(JX-1,JY)-ZB(JX-1,JY-1));
            SW = W2 + W3;
            if SW == 0
                WY2 = 0.5;
                WY3 = 0.5;
            else
                WY2 = W2/SW;
                WY3 = W3/SW;
            end
            ZY(JX,JY) = WY2*ZB(JX-1,JY) + WY3*ZB(JX-1,JY+1);
            ZXY(JX,JY) = ...
                         WY2*(WX2*ZAB(JX-1,JY-1)+WX3*ZAB(JX,JY-1)) + ...
                         WY3*(WX2*ZAB(JX-1,JY)+WX3*ZAB(JX,JY));
        end
    end

    % ... Werte Rausschreiben
    ZX33 = ZX(6);
    ZX43 = ZX(7);
    ZX34 = ZX(10);
    ZX44 = ZX(11);

    ZY33 = ZY(6);
    ZY43 = ZY(7);
    ZY34 = ZY(10);
    ZY44 = ZY(11);

    ZXY33 = ZXY(6);
    ZXY43 = ZXY(7);
    ZXY34 = ZXY(10);
    ZXY44 = ZXY(11);

    P00 = Z33;
    P01 = ZY33;
    P10 = ZX33;
    P11 = ZXY33;

    % ... Koeffizienten fürs Polynom
    ZX3B3 = (ZX34-ZX33)*B3;
    ZX4B3 = (ZX44-ZX43)*B3;
    ZY3A3 = (ZY43-ZY33)*A3;
    ZY4A3 = (ZY44-ZY34)*A3;
    A = ZA3B3 - ZX3B3 - ZY3A3 + ZXY33;
    B = ZX4B3 - ZX3B3 - ZXY43 + ZXY33;
    C = ZY4A3 - ZY3A3 - ZXY34 + ZXY33;
    D = ZXY44 - ZXY43 - ZXY34 + ZXY33;
    E = A + A - B - C;
    A3SQ = A3*A3;
    B3SQ = B3*B3;
    P02 = (2.0*(Z3B3-ZY33)+Z3B3-ZY34)*B3;
    P03 = (-2.0*Z3B3+ZY34+ZY33)*B3SQ;
    P12 = (2.0*(ZX3B3-ZXY33)+ZX3B3-ZXY34)*B3;
    P13 = (-2.0*ZX3B3+ZXY34+ZXY33)*B3SQ;
    P20 = (2.0*(Z3A3-ZX33)+Z3A3-ZX43)*A3;
    P21 = (2.0*(ZY3A3-ZXY33)+ZY3A3-ZXY43)*A3;
    P22 = (3.0*(A+E)+D)*A3*B3;
    P23 = (-3.0*E-B-D)*A3*B3SQ;
    P30 = (-2.0*Z3A3+ZX43+ZX33)*A3SQ;
    P31 = (-2.0*ZY3A3+ZXY43+ZXY33)*A3SQ;
    P32 = (-3.0*E-C-D)*B3*A3SQ;
    P33 = (D+E+E)*A3SQ*B3SQ;

    % ... Berechnen des Polynoms
    DY = VK - Y3;
    Q0 = P00 + DY*(P01+DY*(P02+DY*P03));
    Q1 = P10 + DY*(P11+DY*(P12+DY*P13));
    Q2 = P20 + DY*(P21+DY*(P22+DY*P23));
    Q3 = P30 + DY*(P31+DY*(P32+DY*P33));
    DX = UK - X3;
    W(K) = Q0 + DX*(Q1+DX*(Q2+DX*Q3));
end % Ende Schleife über alle Rückgabewerte

end % Ende Funktion
