%function bend1ASimpleCorr2
%elastic (Stage 1) and elastic-plastic  stages
%numerical solution
box on;
grid on;
hold on;
global CC0 CC1 C0 C1 Mgiven d h
%CC0 CC1 are at t=[0: t1y] in Say = CC0 + CC1*(y);
%Say = C0 + C1*(y-c);  %yield stress at t > t1y
CC0 = 530.495;
CC1 = -15103.41
h0 = 0.014;
Mgiven = 0.012;
ap = 0.0005;  %corrosion rate at sigma > 0 in [m per year]
an = 0.0001  %corrosion rate at sigma < 0

Sap = 6*Mgiven./(h0.^2);  %max Sa at sigma > 0  at t = 0
Sayp = CC0 + CC1*(h0/2);  %max Say at sigma > 0  at t = 0
dif = Sap - Sayp;
dt = 0.01;
hp = 0; %da h positive
hn = 0; %da h negative
h = h0;
t = 0
xr(2) = 0;
xr(1) = h0/2;
C =[];
D =[];
H = [];
T =[];
p =0;

%%%%%%%%%%%Stage 1A
for d = 1.0:0.05:1.0
    p = p+1;

while dif < 0
    
    xr(1) = h/2;
    t = t + dt;
    hp = hp + ap*dt;
    hn = hn + an*dt;
    h = h0 - (hp + hn);
    Sap = 6*Mgiven./(h.^2);
    San =  -6*Mgiven./(h.^2)
    Sayp = CC0 + CC1*(h0/2-hp);
    Sayn = d*(CC0 + CC1*(-h0/2+hn));
    dif = Sap - Sayp;
    %plot(h, xr(1), 'd')
    %plot(t, xr(2), '*r')
    plot(t, Sap,'*')
    plot(t, San,'*')
     
end
   t1y = t;
   Sap1 = Sap;
   Sayp1= Sayp;
   hp1 = hp;
   hn1 = hn;
   h1y = h;
   
    
%%%%%%%%%Stage 2A

%k = (C0 + C1*(he-c))/he;
%plotting for finding initial approximation
% fimplicit(@(he,c) (he.*(C0 - C1*(c - he)))/2 + ((2*c + h - 2*he).*(4*C0 - 2*C1*c +C1*h + 2*C1*he))/8 - ((c - h/2).^2.*(C0 - C1*(c - he)))./(2*he), [-0.007 0.007]); %Sigma =0
% hold on
% fimplicit(@(he,c) (C0/2 - (C1*c)/2)*(c + h/2).^2 - he.^2.*(C0/2 - (C1*c)/2) - (C1*he.^3)/3 + (C1*(c + h/2).^3)/3 + (he.^2.*(C0 - C1*(c - he)))/3 - ((c - h/2).^3.*(C0 - C1*(c - he)))./(3*he) - Mgiven', [-0.007 0.007]); %M-Mgiven=0
 
%solution for stage 2A
%Say = C0 + C1*(y-c);
%xr = [he, c]
%fr = [Sigma, M - Mgiven]
Sayn = d*(CC0 + CC1*(-h0/2+hn)); %max Say at sigma < 0  at t = t1y %Проверить домножение на d
dif2 = Sap - Sayn; % San = Sap at t1y
while dif2 < 0
   % p = p+1
    t = t + dt;
    hp = hp + ap*dt;
    hn = hn + an*dt;
    h = h0 - (hp + hn);
    cCor = (hp-hn)/2;   %shifp of central axis 
    C0 = CC0-CC1*cCor;
    C1 = CC1;
    [xr, fr, ex] = fsolve(@function2A, [h/2, 0], optimset('TolX', 1.0e-15)); 
    %hold on
    %plot(xr(1), xr(2), '*r')
    k = (C0 + C1*(xr(1)-xr(2)))/xr(1);  %k = (C0 + C1*(he-c))/he;
    San1 = -6*Mgiven./(h.^2)
    San = (h/2-xr(2)).*k;   %San = (h/2-c).*k -- absolute value of Sa at -h/2
    Sayn = d*(CC0 + CC1*(-h0/2+hn)); %max Say at sigma < 0  at t = t1y %проверить про домножение на d 
    Sayp = CC0 + CC1*(h0/2-hp);
    dif2 = San - Sayn;
  % plot(h, xr(1), '*')
   %plot(t, xr(2), '*g')
    C(p) = xr(2);
    T(p) = t;
     H(p) = xr(1);
   plot(t,Sayp,'*r');
   plot(t, -San,'*r')
end
hp2 = hp;
hn2 = hn;

   if CC1>0 
%%%%%%%%%%%%%%%%%%%%Stage 2B
%!kc = (C0 - C1*(he+c))/he;
%Say = C0 + C1*(y-c);
%xr = [he, c]
%fr = [Sigma, M - Mgiven]
Sayp = CC0 + CC1*(h0/2-hp); %max Say at sigma > 0  at t = t1y 
dif2 = Sap - Sayp; % San = Sap at t1y
while dif2 < 0
    t = t + dt;
    hp = hp + ap*dt;
    hn = hn + an*dt;
    h = h0 - (hp + hn);
    cCor = (hp-hn)/2;   %shifp of central axis 
    C0 = CC0-CC1*cCor;
    C1 = CC1;
    [xr, fr, ex] = fsolve(@function2B_f, [h/2, 0], optimset('TolX', 1.0e-5)); 
    %hold on
    %plot(xr(1), xr(2), '*r')
    k = d*(C0 - C1*(xr(1)+xr(2)))/xr(1);  %k = (C0 + C1*(he-c))/he;
    San = (-h/2+xr(2)).*k;   %San = (h/2-c).*k -- absolute value of Sa at -h/2
    Sayp = CC0 + CC1*(h0/2-hp); %max Say at sigma < 0  at t = t1y 
    dif2 = Sap - Sayp;
    plot(t, xr(2), 'd')
end

   end;
   
%      t1y
%    Sap1
%    Sayp1
%    hp1
%    hn1
%    h1y
%    
%    t2y = t
%    %San2 = San
%    Sayn2= Sayn
   % hp
%    hn
%    h

 %%%%%%%%%%%%%Stage 3A
 
  while xr(1) >= 10^-16
      %p = p+1
    t = t + dt;
    hp = hp + ap*dt;
    hn = hn + an*dt;
    h = h0 - (hp + hn);
    cCor = (hp - hn)/2;   
    C0 = CC0-CC1*cCor;
    C1 = CC1;
    [xr, fr, ex] = fsolve(@function3N, [h/2, 0], optimset('TolX', 1.0e-15)); 
%     k = (C0 + C1*(xr(1)-xr(2)))/xr(1);  %k = (C0 + C1*(he-c))/he;
%     San = (h/2-xr(2)).*k;   %San = (h/2-c).*k -- absolute value of Sa at -h/2
%     Sayn = CC0 + CC1*(-h0/2+hn); %max Say at sigma < 0  at t = t1y
%     %hold on
     %plot(t, xr(2), '*b')
     %plot(h, xr(1), 'o');
%      C(p) = xr(2)
%      T(p) = t;
%      H(p) = xr(1);
    Sayp = CC0 + CC1*(h0/2-hp);
    Sayn = d*(CC0 + CC1*(-h0/2+hn));
    plot(t, Sayp,'o');
    plot(t, -Sayn, 'o')
  end
 T(p) = t;
 D(p) = d;
end

% legend({'an = ap = 0.0001 м/год', 'ap = 0.0001, an = 0.0002 м/год','ap = 0.0003, an = 0.0001 м/год' },'Location','northwest' )
%  xlabel('d')
%   ylabel('t*(Долговечность лет)');
% Cn;
% Tn; 
% C2 ; 
% T2 ;
% C3 = C;
% T3 = T;

% plot(Cn, Tn);
% plot(C2, T2);
%plot(C3, T3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stage2A: tensed plastic zone 
%!Sigma = int(k*y, y, -(h/2-c), he) + int(Say, y, he, h/2+c);
%cc = solve(Sigma, c);   %This gives cumbersome results
%!hhe = solve(Sigma, he);
%simplify(hhe)

%received solution:
%Sigma 
%Sigma = (he*(C0 - C1*(c - he)))/2 + ((2*c + h - 2*he)*(4*C0 - 2*C1*c +
%C1*h + 2*C1*he))/8 - ((c - h/2)^2*(C0 - C1*(c - he)))/(2*he);
 
%cc(he)
%cc1 = root(4*C1*z^3 - 8*C1*he*z^2 - 4*C1*h*z^2 - 4*C0*z^2 + 4*C1*h*he*z + 4*C1*he^2*z + C1*h^2*z + 8*C0*he*z + 4*C0*h*z + 4*C0*h*he - 4*C0*he^2 - C0*h^2, z, 1);

%!!!!!coefficients of equation (polynom) for finding c  
% A3 = 4*C1;                           %at c^3
% A2 = -8*C1*he - 4*C1*h - 4*C0;       %at c^2
% A1 = 4*C1*h*he +4*C1*he^2 + C1*h^2 + 8*C0*he + 4*C0*h;   %at c
% A0 = 4*C0*h*he -4*C0*he^2 - C0*h^2;                      %at c^0

%he(c)
%!!! hhe = (2*C0*c + C0*h - 2*C1*c^2 - 2*(C0*c*h*(2*C0 - 2*C1*c + C1*h))^(1/2) + C1*c*h)/(2*(C0 - C1*c))  % this is needed root!!!!
% (2*C0*c + C0*h - 2*C1*c^2 + 2*(C0*c*h*(2*C0 - 2*C1*c + C1*h))^(1/2) + C1*c*h)/(2*(C0 - C1*c))


%Bending moment
%!M = int(k*y^2, y, -(h/2-c), he) + int(Say*y, y, he, h/2+c);
% Mb = simplify(M - c*Sigma)  %this expression is more cumbersome
% Mb2 = simplify(M + c*Sigma)  %this expression is more cumbersome
% Mb3 = simplify(M - h*Sigma/2)  %this expression is more cumbersome
% Mb4 = simplify(M + h*Sigma/2)  %this expression is more cumbersome
% Mb5 = simplify(M - h*Sigma)  %this expression is more cumbersome
% Mb6 = simplify(M + h*Sigma)  %this expression is more cumbersome
% Mb7 = simplify(M + he*Sigma)  %this expression is more cumbersome
% Mb8 = simplify(M - he*Sigma)  %this expression is more cumbersome

%Solution for M
%M = (C0/2 - (C1*c)/2)*(c + h/2)^2 - he^2*(C0/2 - (C1*c)/2) - (C1*he^3)/3 + (C1*(c + h/2)^3)/3 + (he^2*(C0 - C1*(c - he)))/3 - ((c - h/2)^3*(C0 - C1*(c - he)))/(3*he);
%!hhe2 = solve(M-Mgiven, he);
%cc2 = solve(M-Mgiven, c);  %this gives the 4th order equation for c
% hheb2 = solve(Mb2-Mgiven, he)
% hheb3 = solve(Mb3-Mgiven, he)
% hheb4 = solve(Mb4-Mgiven, he)
% hheb5 = solve(Mb5-Mgiven, he)
% hheb6 = solve(Mb6-Mgiven, he)
% hheb7 = solve(Mb7-Mgiven, he)
% hheb8 = solve(Mb8-Mgiven, he)

%Solution for hh
%hh = root(4*C1*c*z^3 - 4*C0*z^3 + 12*C0*c*h*z + 12*C1*c^2*h*z -
%3*C1*c*h^2*z + 2*C1*h^3*z + 3*C0*h^2*z - 12*C1*c^3*z + 12*C0*c^2*z - 24*Mgiven*z - 12*C1*c^3*h + 12*C0*c^2*h - 6*C0*c*h^2 + 6*C1*c^2*h^2 - C1*c*h^3 + C0*h^3 + 8*C1*c^4 - 8*C0*c^3, z, 1);

%!!!!!coefficients of equation (polynom) for finding he
%!B3 = 4*C1*c - 4*C0;   %at he^3
%!B2 = 0;               %at he^2
%!B1 = 12*C0*c*h + 12*C1*c^2*h - 3*C1*c*h^2 + 2*C1*h^3 + 3*C0*h^2 - 12*C1*c^3 + 12*C0*c^2 - 24*Mgiven;    %at he
%= B1 = 12*C0*(h/2 + c)^2 + C1*(- 12*c^3 + 12*c^2*h - 3*c*h^2 + 2*h^3) - 24*Mgiven;  %!!! the same!!!
%!B0 = - 12*C1*c^3*h + 12*C0*c^2*h - 6*C0*c*h^2 + 6*C1*c^2*h^2 - C1*c*h^3 + C0*h^3 + 8*C1*c^4 - 8*C0*c^3; %at he^0
%= B0 = 8*(C0-C1*c)*(h/2 - c)^3;  %!!! the same!!!

%The real solution of the above cubic equation, hope, is (hhe2=hhe3)
% Q = (h/2-c)^6 - (C0*(h/2+c)^2 - C1(c^3-c^2*h+c*h^2/4-h^3/6)-2*Mgiven)^3/(C0-C1*c)^3;
% QQ = sqrt(Q);
% hhe3 = ((h/2 - c)^3 + QQ)^(1/3) + ((h/2 - c)^3 - QQ)^(1/3);

%We can equate both solutions for hhe and thus find c numerically !!!!!!
%solve(hhe-hhe3, c)
%after that to find he from hhe.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stage2B: compressed plastic zone 
%!kc = (C0 - C1*(he+c))/he;  %withour multiplaer d

%!Sigma = int(kc*y, y, -he, h/2+c) - int(Say, y, -h/2+c, -he)
%cc = solve(Sigma, c)   %This gives cumbersome results
%!hhe = solve(Sigma, he)
%simplify(hhe)

%received solution:
%Sigma = ((C0 - C1*(c + he))*(c + h/2)^2)/(2*he) - (he*(C0 - C1*(c + he)))/2 -
%((2*c - h + 2*he)*(2*C1*c - 4*C0 + C1*h + 2*C1*he))/8;
 
%cc1 = root(4*C1*z^3 + 8*C1*he*z^2 + 4*C1*h*z^2 - 4*C0*z^2 + 4*C1*h*he*z +
%4*C1*he^2*z + C1*h^2*z - 8*C0*he*z - 4*C0*h*z + 4*C0*h*he - 4*C0*he^2 - C0*h^2, z, 1);

%!!! hhe = (C0*h - 2*C0*c + 2*C1*c^2 - 2*(C0*c*h*(2*C1*c - 2*C0 + C1*h))^(1/2) + C1*c*h)/(2*(C0 - C1*c)) % this is needed root!!!!
%hhe = (C0*h - 2*C0*c + 2*C1*c^2 + 2*(C0*c*h*(2*C1*c - 2*C0 + C1*h))^(1/2) + C1*c*h)/(2*(C0 - C1*c))

%!!!!!coefficients of equation (polynom) for finding c  
% A3 = 4*C1;                           %at c^3
% A2 = -8*C1*he - 4*C1*h - 4*C0;       %at c^2
% A1 = 4*C1*h*he +4*C1*he^2 + C1*h^2 + 8*C0*he + 4*C0*h;   %at c
% A0 = 4*C0*h*he -4*C0*he^2 - C0*h^2;                      %at c^0

%Bending moment
%!M = d*int(kc*y^2, y, -he, h/2+c) - d*int(Say*y, y, -h/2+c, -he)
% Mb = simplify(M - c*Sigma)  %this expression is more cumbersome
% Mb2 = simplify(M + c*Sigma)  %this expression is more cumbersome
% Mb3 = simplify(M - h*Sigma/2)  %this expression is more cumbersome
% Mb4 = simplify(M + h*Sigma/2)  %this expression is more cumbersome
% Mb5 = simplify(M - h*Sigma)  %this expression is more cumbersome
% Mb6 = simplify(M + h*Sigma)  %this expression is more cumbersome
% Mb7 = simplify(M + he*Sigma)  %this expression is more cumbersome
% Mb8 = simplify(M - he*Sigma)  %this expression is more cumbersome

%Solution for M
%M = (C0/2 - (C1*c)/2)*(c + h/2)^2 - he^2*(C0/2 - (C1*c)/2) - (C1*he^3)/3 + (C1*(c + h/2)^3)/3 + (he^2*(C0 - C1*(c - he)))/3 - ((c - h/2)^3*(C0 - C1*(c - he)))/(3*he);
%!hhe2 = solve(M-Mgiven, he)
% cc2 = solve(M-Mgiven, c);  %this gives the 4th order equation for c
% hheb2 = solve(Mb2-Mgiven, he)
% hheb3 = solve(Mb3-Mgiven, he)
% hheb4 = solve(Mb4-Mgiven, he)
% hheb5 = solve(Mb5-Mgiven, he)
% hheb6 = solve(Mb6-Mgiven, he)
% hheb7 = solve(Mb7-Mgiven, he)
% hheb8 = solve(Mb8-Mgiven, he)

%Solution for hh
%hh =  root(4*C1*c*d*z^3 - 4*C0*d*z^3 - 12*C1*c^2*d*h*z - 3*C1*c*d*h^2*z -
%12*C0*c*d*h*z - 2*C1*d*h^3*z + 3*C0*d*h^2*z - 12*C1*c^3*d*z + 12*C0*c^2*d*z - 24*Mgiven*z - 6*C1*c^2*d*h^2 - 12*C1*c^3*d*h + 12*C0*c^2*d*h + 6*C0*c*d*h^2 - C1*c*d*h^3 - 8*C1*c^4*d + 8*C0*c^3*d + C0*d*h^3, z, 1);

%!!!!!coefficients of equation (polynom) for finding he
%!B3 = 4*C1*c*d - 4*C0*d;   %at he^3  or B3 = -4*d*(C0-C1*c)
%!B2 = 0;
%!B1 = - 12*C1*c^2*d*h - 3*C1*c*d*h^2 - 12*C0*c*d*h - 2*C1*d*h^3 + 3*C0*d*h^2 - 12*C1*c^3*d + 12*C0*c^2*d - 24*Mgiven;
%or B1 = 12*d*C0*(h/2-c)^2 - 12*d*C1*(c^3+c^2*h+c*h^2/4+h^3/6) - 24*Mgiven;
%!B0 = - 6*C1*c^2*d*h^2 - 12*C1*c^3*d*h + 12*C0*c^2*d*h + 6*C0*c*d*h^2 - C1*c*d*h^3 - 8*C1*c^4*d + 8*C0*c^3*d + C0*d*h^3;
%or B0 = 8*d*(C0-C1*c)*(h/2+c)^2;


%The real solution of the above cubic equation, hope, is (hhe2=hhe3)
% Q = (h/2+c)^6 - (C0*(h/2-c)^2 - C1(c^3+c^2*h+c*h^2/4+h^3/6))^3/(C0-C1*c)^3;
% QQ = sqrt(Q);
% hhe3 = ((h/2 + c)^3 + QQ)^(1/3) + ((h/2 + c)^3 - QQ)^(1/3);

%We can equate both solutions for hhe and thus find c numerically !!!!!!
%solve(hhe-hhe3, c)
%after that to find he from hhe.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%stage 3: two plastic zones
%!hem = d*he*(C0-C1*c)/(C0-C1*c + C1*he*(d+1));   %he^-

%!Sigma3 = d*int(k*y, y, -(h/2-c), he) + int(Say, y, he, h/2+c);
%cc = solve(Sigma, c);   %This gives cumbersome results
%!hhe = solve(Sigma, he);
%simplify(hhe)

