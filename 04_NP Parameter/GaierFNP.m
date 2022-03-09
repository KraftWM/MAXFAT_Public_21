function [FNP] = GaierFNP(sig)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%--------------------------------------------------------------------------
%                      Eingabedaten auswerten
%--------------------------------------------------------------------------

if size(sig,2) == 3
    sig_xx = sig(:,1)';
    sig_yy = sig(:,2)';
    sig_zz = 0*sig_yy;
    tau_xy = sig(:,3)';
    tau_xz = 0*tau_xy;
    tau_yz = 0*tau_xy;
    
elseif size(sig,2) == 6
    sig_xx = sig(:,1)';
    sig_yy = sig(:,2)';
    sig_zz = sig(:,3)';
    tau_xy = sig(:,4)';
    tau_xz = sig(:,5)';
    tau_yz = sig(:,6)';
    
else
    msg = 'kein ebener oder 3D Spannungszustand';
    error(msg);
    
end

%--------------------------------------------------------------------------
%                  Berechnung Nichtproportionalitätskennzahl
%--------------------------------------------------------------------------

S = [sig_xx; sqrt(2)*tau_xy; sig_yy; sqrt(2)*tau_xz; sig_zz; sqrt(2)*tau_yz]';

I11 = sum( S(:,2).^2 + S(:,3).^2 + S(:,4).^2 + S(:,5).^2 + S(:,6).^2 );
I22 = sum( S(:,1).^2 + S(:,3).^2 + S(:,4).^2 + S(:,5).^2 + S(:,6).^2 );
I33 = sum( S(:,1).^2 + S(:,2).^2 + S(:,4).^2 + S(:,5).^2 + S(:,6).^2 );
I44 = sum( S(:,1).^2 + S(:,2).^2 + S(:,3).^2 + S(:,5).^2 + S(:,6).^2 );
I55 = sum( S(:,1).^2 + S(:,2).^2 + S(:,3).^2 + S(:,4).^2 + S(:,6).^2 );
I66 = sum( S(:,1).^2 + S(:,2).^2 + S(:,3).^2 + S(:,4).^2 + S(:,5).^2 );

I12 = -sum( S(:,1) .* S(:,2) );
I13 = -sum( S(:,1) .* S(:,3) );
I14 = -sum( S(:,1) .* S(:,4) );
I15 = -sum( S(:,1) .* S(:,5) );
I16 = -sum( S(:,1) .* S(:,6) );

I23 = -sum( S(:,2) .* S(:,3) );
I24 = -sum( S(:,2) .* S(:,4) );
I25 = -sum( S(:,2) .* S(:,5) );
I26 = -sum( S(:,2) .* S(:,6) );% I26 = sum( S(:,2) .* S(:,5) );

I34 = -sum( S(:,3) .* S(:,4) );
I35 = -sum( S(:,3) .* S(:,5) );
I36 = -sum( S(:,3) .* S(:,6) );

I45 = -sum( S(:,4) .* S(:,5) );
I46 = -sum( S(:,4) .* S(:,6) );

I56 = -sum( S(:,5) .* S(:,6) );

I = [I11 I12 I13 I14 I15 I16;...
     I12 I22 I23 I24 I25 I26;...
     I13 I23 I33 I34 I35 I36;...
     I14 I24 I34 I44 I45 I46;...
     I15 I25 I35 I45 I55 I56;...
     I16 I26 I36 I46 I56 I66;];
 
l = sort(eig(I),'descend');
 
FNP = sqrt(l(6)/l(5));

if abs(imag(FNP))>0
    msg = ['FNP Gaier komplex. Realteil verwenden: FNP = ',num2str(FNP)];
    warning(msg)
    FNP = real(FNP);
end
end

