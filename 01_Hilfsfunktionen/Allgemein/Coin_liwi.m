function [Sv,I1t,I2t,yfft_I1,yfft_I2]=Coin_liwi(t,Stress,ka,method)
% -------------------------------------------------------------------------
% Complexe Invarianten Hypothese mit Hilfe der lineare Wellen Interferenz
% 
% Aus dem mehrachsigen Spannungszustand werden 2 invariante Sinuswellen (I1,I2)
% ermittelt, wobei I1 repräsentativ für die Normalspannung und I2
% repräsentativ für die Schubspannung ist.
% CR Alex Schmidt
% -------------------------------------------------------------------------
%% Input
%
%   t                  und dehnungen die mit HCM gezÃ¤hlt wird
%                    default = 7 (normaldehnung)
%   t        - Zeitvektor
%   Stress   - (:,6) für den 3D Tensor [xx yy zz xy xz yz]
%              (:,3) für den 2D Tensor [xx yy xy]
%   ka       - Verhältnis von Zugwechsel- und Schubwechselfestigkeit
%              (default) sqrt(3) für von-Mises
%   method   - Umgang mit NP Beanspruchungen 'Preu','Bonte' oder Skalar [0.5..1]
%              (default) 'Preu'
%              0.5 (Bonte): Die Phase reduziert die Amplitude bei NP Beanspruchungen völlig
%              1 (Preumont): Die Phase hat keinen Einfluss auf NP Beanspruchungen
%% Output
%  Sv     - Vergleichsspannungssignal passend zu t
%  I1t       - Invariante Schwingung der Normalspannung
%  I2t       - Invariante Schwingung der Schubspannung
%
%% ------------------------------------------------------------------------
% check Input

    if nargin<4
       
        method='Preu';             
        
    end

    if nargin<3

        ka=sqrt(3); % mises

    else

        if isempty(ka)

            ka=sqrt(3); % mises

        end

    end
    
    Stress=Spaltenvektor(Stress);
    
% -------------------------------------------------------------------------    

    Fs=1/(t(2)-t(1)); % Abtastfrequenz

    ka2=ka.^2;

    if size(Stress,2)==3
        
        [~,yfft]=fourierFFT([Stress(:,1)+eps.*Stress(:,3) Stress(:,2:3)],Fs);
        
        Ax=abs(yfft(:,1));
        Ay=abs(yfft(:,2));
        Axy=abs(yfft(:,3));

        phiX=angle(yfft(:,1));
        phiY=angle(yfft(:,2));
        phiXY=angle(yfft(:,3));  
  
        
        if strcmpi(method,'Bonte')        
        
            AmpM=sqrt((1./2).*(Ax.^2 + ka2.*Axy.^2 + Ay.^2 + (2-ka2).* Ax.*Ay.*cos(phiX - phiY) + abs(Ax.^2.*exp(2i.*phiX) + ka2.*Axy.^2.*exp(2i.*phiXY) + Ay.^2.*exp(2i.*phiY) + (2-ka2).* Ax.*Ay.*exp(1i.*(phiX + phiY)))));
        
        elseif strcmpi(method,'Preu')   
            
            AmpM=sqrt(Ax.^2 + Ay.^2 + (2-ka2).*Ax.*Ay.*cos(phiX - phiY)  + ka2.*Axy.^2);
            
        else
            
            kb=method;            
            AmpM=sqrt(kb.*(Ax.^2 + ka2.*Axy.^2 + Ay.^2 + (2-ka2).* Ax.*Ay.*cos(phiX - phiY)) + (1-kb).*abs(Ax.^2.*exp(2i.*phiX) + ka2.*Axy.^2.*exp(2i.*phiXY) + Ay.^2.*exp(2i.*phiY) + (2-ka2).* Ax.*Ay.*exp(1i.*(phiX + phiY))));
        
        end        
        
        BpM=0.5*angle(Ax.^2.*exp(2i.*(phiX)) + Ay.^2.*exp(2i.*(phiY)) +(2-ka2).*Ax.*Ay.*exp(1i.*(phiY+phiX)) + ka2.*Axy.^2.*exp(2i.*(phiXY)));   
        
        
        I1vec=Ax.*exp(1i.*(phiX))+Ay.*exp(1i.*(phiY));   
        
        if nargout > 1
            I1_amp=sqrt((1./4).*((1 - kb).* abs(Ax.*exp(1i.*phiX) + Ay.*exp(1i.*phiY)).^2 + kb.*(Ax.^2 + Ay.^2 + 2.*Ax.*Ay.*cos(phiX - phiY))));

            I2_amp=sqrt((1./4).*(3.*(1 - kb).* abs(Ax.*exp(1i.*phiX) + Ay.*exp(1i.*phiY)).^2 + 2.*ka2.*(-1 + kb).*abs(2.*Axy.^2.*exp(2.*1i.*phiXY) - 2.*Ax.*Ay.*exp(1i.*(phiX + phiY))) + kb.*(3.*(Ax.^2 + Ay.^2) + 4.*Axy.^2.*ka2 + 2.*Ax.*Ay.*(3 - 2.*ka2).*cos(phiX - phiY))));
            I2s_vec=(1./4).*(3.*Ax.^2.*exp(2.*1i.*phiX) + 3.*Ay.^2.*exp(2.*1i.*phiY) + 4.*Axy.^2.*exp(2.*1i.*phiXY).*ka2 - 2.*Ax.*Ay.*exp(1i.*(phiX + phiY)).*(-3 + 2.*ka2));
 
        end
        
    elseif size(Stress,2)==6
        
%         [~,yfft]=fourierFFT([Stress(:,1:3)+eps.*Stress(:,4:6) Stress(:,4:6)],Fs);
        
        [~,yfft]=fourierFFT(Stress,Fs);
        
        Ax=abs(yfft(:,1));
        Ay=abs(yfft(:,2));
        Az=abs(yfft(:,3));
        Axy=abs(yfft(:,4));
        Axz=abs(yfft(:,5));
        Ayz=abs(yfft(:,6));

        phiX=angle(yfft(:,1));
        phiY=angle(yfft(:,2));
        phiZ=angle(yfft(:,3));
        phiXY=angle(yfft(:,4));  
        phiXZ=angle(yfft(:,5));    
        phiYZ=angle(yfft(:,6));    

        
        if strcmpi(method,'Bonte')
            
            AmpM=sqrt((1./2).*(Ax.^2 + Ay.^2 + Az.^2 + (2-ka2).* Ax.*Ay.*cos(phiX - phiY) + (2-ka2).* Ax.*Az.*cos(phiX - phiZ) + (2-ka2).* Ay.*Az.*cos(phiY - phiZ) + ka2.*Axy.^2 + ka2.*Axz.^2 + ka2.*Ayz.^2 + abs(Ax.^2.*exp(2i.*phiX) + Ay.^2.*exp(2i.*phiY) + Az.^2.*exp(2i.*phiZ) +(2-ka2).*  Ax.*Ay.*exp(1i.*(phiX + phiY)) +(2-ka2).* Ax.*Az.*exp(1i.*(phiX + phiZ)) +(2-ka2).* Ay.*Az.*exp(1i.*(phiY + phiZ)) + ka2.*Axy.^2.*exp(2i.*phiXY) + ka2.*Axz.^2.*exp(2i.*phiXZ)  + ka2.*Ayz.^2.*exp(2i.*phiYZ) )));
        
        elseif strcmpi(method,'Preu')
            
            AmpM=sqrt(Ax.^2 + Ay.^2 + Az.^2 + (2-ka2).* Ax.*Ay.*cos(phiX - phiY) + (2-ka2).* Ax.*Az.*cos(phiX - phiZ) + (2-ka2).* Ay.*Az.*cos(phiY - phiZ) + ka2.*Axy.^2 + ka2.*Axz.^2 + ka2.*Ayz.^2);
            
        else
            
            kb=method;            
            AmpM=sqrt(kb.*(Ax.^2 + Ay.^2 + Az.^2 + (2-ka2).* Ax.*Ay.*cos(phiX - phiY) + (2-ka2).* Ax.*Az.*cos(phiX - phiZ) + (2-ka2).* Ay.*Az.*cos(phiY - phiZ) + ka2.*Axy.^2 + ka2.*Axz.^2 + ka2.*Ayz.^2) + (1-kb).*abs(Ax.^2.*exp(2i.*phiX) + Ay.^2.*exp(2i.*phiY) + Az.^2.*exp(2i.*phiZ) +(2-ka2).*  Ax.*Ay.*exp(1i.*(phiX + phiY)) +(2-ka2).* Ax.*Az.*exp(1i.*(phiX + phiZ)) +(2-ka2).* Ay.*Az.*exp(1i.*(phiY + phiZ)) + ka2.*Axy.^2.*exp(2i.*phiXY) + ka2.*Axz.^2.*exp(2i.*phiXZ)  + ka2.*Ayz.^2.*exp(2i.*phiYZ) ));
        
        end
            
        BpM=angle(sqrt(Ax.^2.*exp(2i.*(phiX)) + Ay.^2.*exp(2i.*(phiY)) + Az.^2.*exp(2i.*(phiZ)) +(2-ka2).* Ax.*Ay.*exp(1i.*(phiY+phiX)) +(2-ka2).* Ax.*Az.*exp(1i.*(phiZ+phiX)) +(2-ka2).* Ay.*Az.*exp(1i.*(phiY+phiZ)) + ka2.*Axy.^2.*exp(2i.*(phiXY)) + ka2.*Axz.^2.*exp(2i.*(phiXZ)) + ka2.*Ayz.^2.*exp(2i.*(phiYZ)) ));

        
        I1vec=Ax.*exp(1i.*phiX) + Ay.*exp(1i.*phiY) + Az.*exp(1i.*phiZ);
        
        if nargout > 1

            I1_amp=sqrt((1./4).*((1 - kb).* abs(Ax.*exp(1i.*phiX) + Ay.*exp(1i.*phiY) + Az.*exp(1i.*phiZ)).^2 + kb.*(Ax.^2 + Ay.^2 + Az.^2 + 2.*Ax.*Ay.* cos(phiX - phiY) + 2.*Ax.*Az.*cos(phiX - phiZ) + 2.*Ay.*Az.*cos(phiY - phiZ))));

%             I1_amp=sqrt((1./4).*kb.*(Ax.^2 + Ay.^2 + Az.^2 + 2.*Ax.*Ay.*cos(phiX - phiY) + 2.*Ax.*Az.*cos(phiX - phiZ) + 2.*Ay.*Az.*cos(phiY - phiZ)));
            I2_amp=sqrt(AmpM.^2-I1_amp.^2);
            
%             I2_amp=sqrt((1./4).*(3.*(1 - kb).* abs(Ax.*exp(1i.*phiX) + Ay.*exp(1i.*phiY) + Az.*exp(1i.*phiZ)).^2 + 4.*ka2.*(-1 + kb).* abs(Axy.^2.*exp(2.*1i.*phiXY) + Axz.^2.*exp(2.*1i.*phiXZ) - Ax.*Ay.*exp(1i.*(phiX + phiY)) + Ayz.^2.*exp(2.*1i.*phiYZ) - Ax.*Az.*exp(1i.*(phiX + phiZ)) - Ay.*Az.*exp(1i.*(phiY + phiZ))) + kb.*(3.*(Ax.^2 + Ay.^2 + Az.^2) + 4.*(Axy.^2 + Axz.^2 + Ayz.^2).*ka2 - 2.*(-3 + 2.*ka2).*(Ax.*Ay.*cos(phiX - phiY) + Ax.*Az.*cos(phiX - phiZ) + Ay.*Az.*cos(phiY - phiZ)))));
%             I2_amp=sqrt((1./4).*kb.*(3.*(Ax.^2 + Ay.^2 + Az.^2) + 4.*(Axy.^2 + Axz.^2 + Ayz.^2).*ka2 - 2.*(-3 + 2.*ka2).*(Ax.*Ay.*cos(phiX - phiY) + Ax.*Az.*cos(phiX - phiZ) + Ay.*Az.*cos(phiY - phiZ))));
     
            I2s_vec=(3./4).*(Ax.*exp(1i.*phiX) + Ay.*exp(1i.*phiY) + Az.*exp(1i.*phiZ)).^2 + (Axy.^2.*exp(2.*1i.*phiXY) + Axz.^2.*exp(2.*1i.*phiXZ) - Ax.*Ay.*exp(1i.*(phiX + phiY)) + Ayz.^2.*exp(2.*1i.*phiYZ) - Ax.*Az.*exp(1i.*(phiX + phiZ)) - Ay.*Az.*exp(1i.*(phiY + phiZ))).*ka2;
            
        end
        

    else        
        
        warning('wrong number of tensor components')
        return
        
    end   
    
    

    I1_ph=angle(I1vec);
    
  
    

    Z1=exp(1i.*BpM);
    Z2=exp(1i.*(BpM+pi));

    Z=[Z1 Z2];

    Zi=exp(1i.*I1_ph);
    [~,idx]=min(abs(angle(Z./Zi)),[],2);
    idm = sub2ind(size(Z), 1:size(Z, 1), idx.').';

    PhiVM=angle(Z(idm));
    
    yfft_Sv=AmpM.*exp(1i.*(PhiVM));

    [Sv]=ifourierFFT(yfft_Sv,t);
    
    if nargout>1
       
        I2_ph=angle(sqrt(I2s_vec));    
        
        Z1=exp(1i.*I2_ph);
        Z2=exp(1i.*(I2_ph+pi));

        Z=[Z1 Z2];
        Zi=exp(1i.*PhiVM);
        [~,idx]=min(abs(angle(Z./Zi)),[],2);
        idm = sub2ind(size(Z), 1:size(Z, 1), idx.').';
   
        I2_ph=angle(Z(idm));   
        
        yfft_I1=I1_amp.*exp(1i.*(I1_ph));
        yfft_I2=I2_amp.*exp(1i.*(I2_ph));

        [I1t]=ifourierFFT(yfft_I1,t);  
        [I2t]=ifourierFFT(yfft_I2,t);  
        
    end
    

end



function [f,y_fft]=fourierFFT(y,Fs)


   NFFT=length(y);
   y_f=fft(y,NFFT)/NFFT;

   halbe=int64(NFFT/2+1);
    
   f = Fs/2*linspace(0,1,halbe);

   y_fft=2*y_f(1:halbe,:);
   y_fft(1,:)=y_fft(1,:)/2;
    
end

function [y]=ifourierFFT(y_fft,t)

    m=length(t);
    y=real(ifft(y_fft,m)*m);
    y(isnan(y))=0;
    
end

function [y]=Spaltenvektor(x)
% prüft, ob ein Zeilen oder Spaltenvektor vorliegt
% falls nicht, wird zu einem transformiert
% Annahme: Signallänge > Anzahl Signalwerte

[d1,d2]=size(x);


    if d1<d2
        % es soll IMMER über einzelne Signale transformiert werden
        % Signallänge > Anzahl
        x=transpose(x);

    end

y=x;

end

