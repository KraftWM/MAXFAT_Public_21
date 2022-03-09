function [ Ext ] = find_extrema( S,K )
%Funktion zum Erkennen von Extremstellen in einem Zeitverlauf
%   Funktion zum Erkennen von Extremstellen in einem Zeitverlauf, durch
%   überprüfen der nächsten Nachbarn.
%   S(:,1)->Zeit
%   S(:,2)->Zeitabhängiges Signal
%   K-> Schranke zum Filtern (für Versuchauswertung von Alex, kann auf
%   nullgestetzt werden wenn verwendet um vergleichsspannung auszurechnen)
%
% OUTPUT:
%    Ext = [Zeiten der Extreme,Extrema]
% 
%
%

% Abfangen flasche eingabe
if size(S,2) > size(S,1)
    S = S';
end

% Abfangen wenn für S keine Zeit übergeben wurde
if size(S,2) == 1
    ink = 1:length(S);
    S = [ink',S];
end

% figure(1)
% plot(S(:,1),S(:,2))
% hold on, grid on
n=2;
ende=length(S);
last_ext=2;    %was war das letzte Extremum ? 0 ein Maximum 1 ein Minimum 
Ext=zeros(1,2);
for i=6:ende-5
   %% Finde Maxima  
    if S(i,2)>=S(i-1,2)
     if S(i,2)>=S(i-2,2)
      if S(i,2)>=S(i-3,2)
       if S(i,2)>=S(i-5,2)
        if S(i,2)>=S(i+1,2)
         if S(i,2)>=S(i+2,2)
          if S(i,2)>=S(i+3,2)
            if S(i,2)>=S(i+5,2)                              %Hier ist jetzt ein Maximum gefunden 
%               plot(S(i,1),S(i,2),'x')
              if last_ext==0                                               % Das Gefundene Maximum ist das 2. Maximum in Folge -> gespeichert wird das größere
               if S(i,2)>Ext(n-1,2)
                 Ext(n-1,1)=S(i,1);
                 Ext(n-1,2)=S(i,2);
               end  
              else                                                         %Das gefundene Maximum folgt auf ein Minimum
               if n>4
               %K=schranke(Ext(n-1,2),Ext(n-2,2),Ext(n-3,2),Ext(n-4,2),typ );
               end
               if S(i,2)<Ext(n-1,2)                                        % das gefundene Maxima ist kleiner als das letzte Minimum wird also nicht gespeichert
                   
               elseif abs(S(i,2)-Ext(n-1,2))<K                                 %Abstand der Extrema ist kleiner als Schranke..
               
               else
                 Ext(n,1)=S(i,1);
                 Ext(n,2)=S(i,2);
                 n=n+1;
                 last_ext=0;

               end
              end
            end
          end
          end
         end
        end
      end
     end
   else
   %% Finde Minima 
    if S(i,2)<=S(i-1,2)
     if S(i,2)<=S(i-2,2)
      if S(i,2)<=S(i-3,2)
       if S(i,2)<=S(i-5,2)
        if S(i,2)<=S(i+1,2)
         if S(i,2)<=S(i+2,2)
          if S(i,2)<=S(i+3,2)
           if S(i,2)<=S(i+5,2)                                                                %Hier ist ein Minimum gefunden
%            plot(S(i,1),S(i,2),'o')
              if last_ext==1                                               %Wenn zwei Minima aufeinander folgen Speichere nur das kleinste
                 if S(i,2)<Ext(n-1,2)
                 Ext(n-1,1)=S(i,1);
                 Ext(n-1,2)=S(i,2); 
                 end     
              else                                                         %Hier ein Minimum auf ein Maximum 
                 if n>4
                  %K=schranke(Ext(n-1,2),Ext(n-2,2),Ext(n-3,2),Ext(n-4,2),typ );
                 end
                 if abs(S(i,2)-Ext(n-1,2))<K
                     
%                  elseif abs(S(i,2)-Ext(n-1,2))<K                               %Abstand der Extrema ist kleiner als Schranke..
                 else
                 Ext(n,1)=S(i,1);
                 Ext(n,2)=S(i,2);
                 n=n+1;
                 last_ext=1;
                 end

              end
           end
          end
          end
         end
        end
      end
     end
    end
   end
  
end
%disp('fertig')
% if size(Ext,1) >= 2
%     Ext=Ext(2:length(Ext),:);
% else 
%     Ext = [];
% end
Ext(1,:) = [];
% plot(S(:,1),S(:,2))
% hold on
% plot(Ext(:,1),Ext(:,2),'xr')


end

