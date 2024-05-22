hold on;
grid on;
box on;

global C0 C1 h Mgiven d b
C1 = 15103.41;
h = 0.010;
H = h/2;

C0 = 530.495;
b = 1; 
I = b*h^3/12;


% M0 = [];
% C_0 = [];
% He0 = [];
% 
% M1 = [];
% C_1 = [];
% He1 = [];

% M2 = [];
% C_2 = [];
% He2 = [];

% M3 = [];
% C_3 = [];
% He3 = [];
% 
% M4 = [];
% C_4 = [];
% He4 = [];
% 
% M5 = [];
% C_5 = [];
% He5 = [];


xr(1) = h/2; %he 
xr(2) = 0; %c
Say = @(y) C0 + C1*(y - xr(2)); %ноль на нейтральной плоскости
k = Say(xr(1))./xr(1); % плюс xr(2) не нужно, т к сигма игрек уже в новой
S_e = 0;
n = 0;
for d = 1.3:0.000005:1.3
    Say_min = min(Say(h/2),Say(-h/2)*d); 
    Say_max = max(Say(h/2),Say(-h/2)*d); 
    
  for Mgiven = 0:0.0005:0.05
      %Стадия 1 
      n = n+1;
         S_e_1 = xr(1)*Mgiven/I;
         S_e_2 = -(h/2)*Mgiven/I;
         if (S_e_1 <= Say_min) && (xr(1) == h/2 )
%              plot(Mgiven,  xr(1), '+')
%              plot(Mgiven,  xr(2), '+')
            plot(Mgiven, S_e_1, '+');
            plot(Mgiven, S_e_2, '+'); 
%               C_4(n) = xr(2);
%               He4(n) = xr(1);
%               M4(n) = Mgiven;
         else
             if Say(h/2 + xr(2))> Say(-h/2 + xr(2))*d && S_e < Say(h/2 - xr(2)) 
              display("Stage 2b");                       
              [xr, fr, ex] = fsolve(@function2B_hm, [h/2, 0], optimset('TolX', 1.0e-5));
               k = d*(C0 - C1*(xr(1)+xr(2)))./xr(1); %Поменять знак перед С1 на плюс. В уравнениях тоже, тк С1 глобальная. 
               S_e = k*(h/2+xr(2));
               plot(Mgiven,S_e, 'o')
               plot(Mgiven, -Say_min, 'o')
%                  plot(Mgiven,  xr(1), '*')
%                  plot(Mgiven,  xr(2), '*')
%               C_4(n) = xr(2);
%               He4(n) = xr(1);
%               M4(n) = Mgiven;
              end
       if Say(h/2 + xr(2)) < Say(-h/2 + xr(2))*d && S_e < d*Say(-h/2 + xr(2)) %здесь S_e тоже для сжатой                                                        
             
               display("Stage 2A");% предел текучести достигается на растянутой стороне
            % Может на самом деле работать и при положительном направлении
            % градиента при больших значениях d, например, при d=1.5;

               %S_e_2 = -(h/2)*Mgiven/I;
               [xr, fr, ex] = fsolve(@function2A, [h/2, 0], optimset('TolX', 1.0e-5));
              k = (C0 + C1*(xr(1)+xr(2)))./xr(1);
              S_e = k*(h/2 - xr(2)); %для сжатой, стоит обозначить переменные ????(напряжения на противотоложной стороне)
              plot(Mgiven, Say_min, '*');
              plot(Mgiven,-S_e,'*');
%              plot(Mgiven, xr(1), 'o')
%               plot(Mgiven, xr(2), 'o')
%               C_4(n) = xr(2);
%               He4(n) = xr(1);
%               M4(n) = Mgiven;
               
          end
          
          %начало третьей стадии
          %для третьей стадии должно быть два разных вывода для случая а и
          %b 
          if (S_e >= Say_max) || (S_e == 0 )
          display('Stage 3');
          [xr, fr, ex] = fsolve(@function3N, [h/2, 0], optimset('TolX', 1.0e-5));
           plot(Mgiven ,-Say_min,'^');
           plot(Mgiven,Say_max, '^')
             hem = d*xr(1)*(C0-C1*xr(2))/(C0-C1*xr(2) + C1*xr(1)*(d+1)); 
%            plot(Mgiven,  xr(1), '^')
%            plot(Mgiven,  xr(2), '^')
%               C_4(n) = xr(2);
%               He4(n) = xr(1);
%               M4(n) = Mgiven;
          end
         
       end
          if xr(1)<=10^(-6)
              break;
          end
      
  end
end



Mz = h^2*(3*C0*(d+1)+C1*h*(d-1))/24
% 
% M0;
% C_0;
% He0;
% % 
% M1;
% C_1;
% He1;
% % % % % % 
% M2;
% C_2;
% He2;
% % % % % % 
% M3;
% C_3;
% He3;
% % % % % % 
% M4;
% C_4;
% He4;
% % % 
% % M5;
% % C_5;
% % He5;
% plot(M0./0.00755, C_0./h,'r','LineWidth',3);
% plot(M0./0.00755, He0./h,'r','LineWidth',3);
% % 
% plot(M1./0.00835, C_1./h,'y','LineWidth',3);
% plot(M1./0.00835, He1./h,'y','LineWidth',3);
% % % % % 
% % % % 
% plot(M2./0.0091, C_2./h,'k','LineWidth',3);
% plot(M2./0.0091, He2./h,'k','LineWidth',3);
% % % % % 
% % % % % 
% plot(M3./0.00985, C_3./h,'g','LineWidth',3);
% plot(M3./0.00985, He3./h,'g','LineWidth',3);
% % % % % % 
% plot(M4./0.0101, C_4./h,'b','LineWidth',3);
% plot(M4./0.0101, He4./h,'b','LineWidth',3);
% % 
% % plot(M5./0.02075, C_5./h,'c','LineWidth',3);
% % plot(M5./0.02075, He5./h,'c','LineWidth',3);
% % % % % 
% % % % % 
%   legend({['d = 1.0'],[''],['d = 1.1'],[''],['d = 1.2'],[''],['d = 1.3'],[''],['d = 1.4'],['']},'Location', 'southwest');
% % % % 
%  xlabel('M/Mкр');
%  ylabel('c/h, he/h');





    