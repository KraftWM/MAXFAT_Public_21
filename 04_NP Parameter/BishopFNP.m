function [FNP] = BishopFNP(sig)
%UNTITLED3 Summary of this function goes here
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
x = [sig_xx; sig_yy; sig_zz; sqrt(2)*tau_xy; sqrt(2)*tau_xz; sqrt(2)*tau_yz];

xtk = x(:,1:end-1);
xtk1 = x(:,2:end);

Lk = zeros(size(x,2)-1,1);
for k = 1 : size(x,2)-1
    Lk(k) = norm( xtk1(:,k) - xtk(:,k) );
end

xbar = 1/(2*sum(Lk)) * sum( Lk .* (xtk1' + xtk') ); 

y = x - xbar';
ytk = y(:,1:end-1)';
ytk1 = y(:,2:end)';


I11 = 1/6 * sum (Lk .* ( ytk(:,1).* ytk(:,1) + ytk1(:,1).* ytk1(:,1) + (ytk(:,1) + ytk1(:,1)).*(ytk(:,1) + ytk1(:,1)) ));
I22 = 1/6 * sum (Lk .* ( ytk(:,2).* ytk(:,2) + ytk1(:,2).* ytk1(:,2) + (ytk(:,2) + ytk1(:,2)).*(ytk(:,2) + ytk1(:,2)) ));
I33 = 1/6 * sum (Lk .* ( ytk(:,3).* ytk(:,3) + ytk1(:,3).* ytk1(:,3) + (ytk(:,3) + ytk1(:,3)).*(ytk(:,3) + ytk1(:,3)) ));
I44 = 1/6 * sum (Lk .* ( ytk(:,4).* ytk(:,4) + ytk1(:,4).* ytk1(:,4) + (ytk(:,4) + ytk1(:,4)).*(ytk(:,4) + ytk1(:,4)) ));
I55 = 1/6 * sum (Lk .* ( ytk(:,5).* ytk(:,5) + ytk1(:,5).* ytk1(:,5) + (ytk(:,5) + ytk1(:,5)).*(ytk(:,5) + ytk1(:,5)) ));
I66 = 1/6 * sum (Lk .* ( ytk(:,6).* ytk(:,6) + ytk1(:,6).* ytk1(:,6) + (ytk(:,6) + ytk1(:,6)).*(ytk(:,6) + ytk1(:,6)) ));

I12 = 1/6 * sum (Lk .* ( ytk(:,1).* ytk(:,2) + ytk1(:,1).* ytk1(:,2) + (ytk(:,1) + ytk1(:,1)).*(ytk(:,2) + ytk1(:,2)) ));
I13 = 1/6 * sum (Lk .* ( ytk(:,1).* ytk(:,3) + ytk1(:,1).* ytk1(:,3) + (ytk(:,1) + ytk1(:,1)).*(ytk(:,3) + ytk1(:,3)) ));
I14 = 1/6 * sum (Lk .* ( ytk(:,1).* ytk(:,4) + ytk1(:,1).* ytk1(:,4) + (ytk(:,1) + ytk1(:,1)).*(ytk(:,4) + ytk1(:,4)) ));
I15 = 1/6 * sum (Lk .* ( ytk(:,1).* ytk(:,5) + ytk1(:,1).* ytk1(:,5) + (ytk(:,1) + ytk1(:,1)).*(ytk(:,5) + ytk1(:,5)) ));
I16 = 1/6 * sum (Lk .* ( ytk(:,1).* ytk(:,6) + ytk1(:,1).* ytk1(:,6) + (ytk(:,1) + ytk1(:,1)).*(ytk(:,6) + ytk1(:,6)) ));

I23 = 1/6 * sum (Lk .* ( ytk(:,2).* ytk(:,3) + ytk1(:,2).* ytk1(:,3) + (ytk(:,2) + ytk1(:,2)).*(ytk(:,3) + ytk1(:,3)) ));
I24 = 1/6 * sum (Lk .* ( ytk(:,2).* ytk(:,4) + ytk1(:,2).* ytk1(:,4) + (ytk(:,2) + ytk1(:,2)).*(ytk(:,4) + ytk1(:,4)) ));
I25 = 1/6 * sum (Lk .* ( ytk(:,2).* ytk(:,5) + ytk1(:,2).* ytk1(:,5) + (ytk(:,2) + ytk1(:,2)).*(ytk(:,5) + ytk1(:,5)) ));
I26 = 1/6 * sum (Lk .* ( ytk(:,2).* ytk(:,6) + ytk1(:,2).* ytk1(:,6) + (ytk(:,2) + ytk1(:,2)).*(ytk(:,6) + ytk1(:,6)) ));

I34 = 1/6 * sum (Lk .* ( ytk(:,3).* ytk(:,4) + ytk1(:,3).* ytk1(:,4) + (ytk(:,3) + ytk1(:,3)).*(ytk(:,4) + ytk1(:,4)) ));
I35 = 1/6 * sum (Lk .* ( ytk(:,3).* ytk(:,5) + ytk1(:,3).* ytk1(:,5) + (ytk(:,3) + ytk1(:,3)).*(ytk(:,5) + ytk1(:,5)) ));
I36 = 1/6 * sum (Lk .* ( ytk(:,3).* ytk(:,6) + ytk1(:,3).* ytk1(:,6) + (ytk(:,3) + ytk1(:,3)).*(ytk(:,6) + ytk1(:,6)) ));

I45 =  1/6 * sum (Lk .* ( ytk(:,4).* ytk(:,5) + ytk1(:,4).* ytk1(:,5) + (ytk(:,4) + ytk1(:,4)).*(ytk(:,5) + ytk1(:,5)) ));
I46 =  1/6 * sum (Lk .* ( ytk(:,4).* ytk(:,6) + ytk1(:,4).* ytk1(:,6) + (ytk(:,4) + ytk1(:,4)).*(ytk(:,6) + ytk1(:,6)) ));

I56 = 1/6 * sum (Lk .* ( ytk(:,5).* ytk(:,6) + ytk1(:,5).* ytk1(:,6) + (ytk(:,5) + ytk1(:,5)).*(ytk(:,6) + ytk1(:,6)) ));

I = [I11 I12 I13 I14 I15 I16;...
     I12 I22 I23 I24 I25 I26;...
     I13 I23 I33 I34 I35 I36;...
     I14 I24 I34 I44 I45 I46;...
     I15 I25 I35 I45 I55 I56;...
     I16 I26 I36 I46 I56 I66;];
 

l = sort(eig(I),'descend');
 
FNP = sqrt(l(2)/l(1)); 
end

