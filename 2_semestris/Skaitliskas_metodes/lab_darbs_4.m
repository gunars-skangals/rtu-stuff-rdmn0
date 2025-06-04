%% 1.piemçrs 
% Ņûtona interpolâcijas polinoms 
clear all  clc ,close all, format compact
ynodes=[2.3,3.4,5.1,6.3,7.5,8.2,7.4]; m = length(ynodes); 
xnodes=[1:m];
coef=ynodes;
for k = 2:m
    coef(k:m) = (coef(k:m) - coef(k-1:m-1))./(xnodes(k:m) - xnodes(1:m+1-k));
end

syms x
pol= coef(m); % polinoma konstuçðana
for k = m-1:-1:1
pol=pol*(x-xnodes(k))+coef(k);
end
polyn(x)=collect(pol)
polyn(x)=expand(pol)
coefpol=sym2poly(polyn(x)) % polinoma koeficienti

disp('Atbilde:')
disp('0.0043x^6-0.1121x^5+1.1451x^4-5.8563x^3+15.6006x^2-18.6817x+10.2')

x_pr = xnodes(1):0.01:xnodes(m);
plot(x_pr,polyn(x_pr),'r-',xnodes,ynodes,'o','LineWidth',3)
title('Newton''s interpolating polynomial')
legend('polynomial','nodes') 
grid on, ylim([1 10])

%disp([num2str(koefpol(1)),'x^6+',num2str(koefpol(2)),'x^5+',num2str(koefpol(3))...
%    'x^4+',num2str(koefpol(4)),'x^3+',num2str(koefpol(5)),'x^2+',num2str(koefpol(6))...
%    'x+',num2str(koefpol(7))])
X_main = sym('x_%d',[1 4])


%% 2.piemçrs.  
% Òûtona interpolâcijas polinoms 
clear all, clc ,close all, format compact
syms x
y=@(x)log(1+x.^2);
xnodes=[0:3];
ynodes=y(xnodes); coef=ynodes;
for k = 2:4
    coef(k:4) = (coef(k:4) - coef(k-1:3))./(xnodes(k:4) - xnodes(1:5-k));
end

pol= coef(4); % polinoma konstuçðana
for k = 3:-1:1
  pol=pol*(x-xnodes(k))+coef(k);
end
polyn(x)=collect(pol)
coefpol=sym2poly(polyn) % polinoma koeficienti

disp('Atbilde:')
disp(' -0.0744x^3+0.3347x^2+0.4328x')

x_pr = 0:0.01:3;
plot(x_pr,polyn(x_pr),'r-',x_pr,y(x_pr),'k--',xnodes,ynodes,'o','LineWidth',3)
title('Newton''s interpolating polynomial')
legend('polynomial','y=ln(1+x^2)','nodes') 
ylim([-0.5 3]), grid

% turpinâjums
figure  % Uzzîmçt interpolâcijas kïûdas grafiku
plot(x_pr,y(x_pr)-polyn(x_pr),'g:','LineWidth',3)
title('Interpolation error'), grid

figure % x vçrtîbas ârpus intervâla [0,3], kur ir doti
       % interpolâcijas mezgli
x1=0:0.01:10;
plot(x1,polyn(x1),'r-',x1,y(x1),'k--','LineWidth',3)
legend('polynomial','y=log(1+x^2)')
ylim([0, 5]), grid

%% 3.piemçrs.  
% Òûtona interpolâcijas polinoms 
clear all, clc ,close all, format compact
syms x, f=@(x)sqrt(2+x.^7);
xnodes=[1:5]; ynodes=f(xnodes);
coef=ynodes;

for k = 2:5
    coef(k:5) = (coef(k:5) - coef(k-1:4))./...
                (xnodes(k:5) - xnodes(1:6-k));
end

pol= coef(5); % polinoma konstuçðana
for k = 4:-1:1
pol=pol*(x-xnodes(k))+coef(k);
end
polyn(x)=collect(pol)
coefpol=sym2poly(polyn) % polinoma koeficienti

int_res=double(polyn(2.5)), f_res=f(2.5)
int_error=abs(int_res-f_res)
disp('Atbilde:')
disp(['1) koef_x^2 = ', num2str(coefpol(3))])
disp(['2) interpolâcijas rezultâts punktâ(2.5) = ', num2str(int_res)])
disp(['3) interpolâcijas kïûda  = ', num2str(int_error)])
disp('4) interpolâcijas kïûda punktâ(x=3) = 0')
disp('tâ kâ punkts x=3 ir interpolâcijas mezgla punkts')

%% 4.piemçrs.  
% Ermita interpolâcija( interp1 )
clear all, close all , clc, format compact
xnodes=[0 pi/8 pi/6 pi/4 pi/3 pi/2] ;
ynodes=[0 0.7486 0.9286 0.8849 -0.1010 0.67];
f=@(x)sin(x.^3+2*x);
x=linspace(0,pi/2,200);
y=interp1(xnodes,ynodes,x,'pchip');

plot(x,y,'r-',xnodes,ynodes,'o',x,f(x),'k-','LineWidth',3)
legend('int.polynomial','nodes','y=sin(x^3+2x)')
title('Hermite interpolation( interp1 )'), grid
 
figure  % Uzzîmçt interpolâcijas kïûdas grafiku 
plot(x,f(x)-y,'g--','LineWidth',3)
%plot(x,f(x)-y,'LineWidth',3,'Color','g','LineStyle','--')
title('Interpolation error'), grid

%% 4.a.uzdevums.  Ermita interpolâcija( interp1 )
clear all, close all ,clc
xnodes=linspace(0,pi/2,20);
f=@(x1)sin(x1.^3+2*x1);
ynodes=f(xnodes);

x=linspace(0,pi/2,200);
y=interp1(xnodes,ynodes,x,'pchip');
plot(x,y,'r-',xnodes,ynodes,'o',x,f(x),'k-','LineWidth',3)
legend('int.polynomial','nodes','y=sin(x^3+2x)')
title('Hermite interpolation( interp1 )'),grid
 
figure  % klûda
plot(x,f(x)-y,'g--','LineWidth',3)
title('Interpolation error'), grid

%% Runge piemçrs
clear all, close all, clc, format compact
% xnodes=linspace(-1,1,7); % vienmçrîgais reþìis ar 7 punktiem intervâlâ [-1,1]
xnodes=linspace(-1,1,21); % vienmçrîgais reþìis ar 21 punktiem intervâlâ [-1,1]
y=@(x)1./(1+25*x.^2); ynodes=y(xnodes);
coef=ynodes; m=length(xnodes);
syms x
pol=0; % Lagranþa interpolâcijas polinoms
for k = 1:m
   w=1;
     for j = [1:k-1 k+1:m]
         w = (x-xnodes(j))/(xnodes(k)-xnodes(j))*w;
     end
   pol=pol + w*ynodes(k);
end
polyn(x)=collect(pol);
coefpol=sym2poly(polyn) % polinoma koeficienti
x_pr = -1:0.01:1;
plot(x_pr,polyn(x_pr),'r-',x_pr,y(x_pr),'k--',xnodes,ynodes,'o','LineWidth',3)
title('Lagrange interpolating polynomial n=20')
legend('polynomial','y=1/(1+25*x^2)','nodes') 
ylim([-1 2.1]),grid

% turpinâjums
figure  % Uzzîmçt interpolâcijas kïûdas grafiku
plot(x_pr,y(x_pr)-polyn(x_pr),'g:','LineWidth',4)
title('Interpolation error'), grid
%ylim([-1 20])

%% Bernðteina piemçrs (1916.g.)
clear all, close all , clc, format compact
% xnodes=linspace(-1,1,7)
xnodes=linspace(-1,1,11);
y=@(x)abs(x);

ynodes=y(xnodes); coef=ynodes;
m = length(xnodes)
syms x
pol=0; % polinoma konstuçðana
for k = 1:m
 w= 1; 
 for j = [1:k-1 k+1:m]
   w = (x-xnodes(j))/(xnodes(k)-xnodes(j))*w;
 end
 pol = pol + w*ynodes(k);
end 
pol
polyn(x)=collect(pol);
coefpol=sym2poly(polyn) % polinoma koeficienti

x_pr = -1:0.01:1;
plot(x_pr,polyn(x_pr),'r-',x_pr,y(x_pr),'k--',xnodes,ynodes,'o','LineWidth',3)
title('Lagrange interpolating polynomial n=10')
legend('polynomial','y=abs(x)','nodes') 
ylim([0 1.3]), grid

% turpinâjums
figure  % Uzzîmçt interpolâcijas kïûdas grafiku
plot(x_pr,y(x_pr)-polyn(x_pr),'g:','LineWidth',4)
title('Interpolation error'), grid


%% 5.piemçrs
% Interpolâcija ar kubiskiem splainiem.
clear all, close all , format compact,clc
f=@(x)abs(x)
xnodes=linspace(-1,1,11); ynodes=f(xnodes);

x=linspace(-1,1,200);
y=interp1(xnodes,ynodes,x,'spline');
% y_inter=interp1(xmezg,ymezg,x,'pchip')
% y_inter=spline(xmezg,ymezg,x)
plot(x,y,'r-',x,f(x),'k--',xnodes,ynodes,'o','LineWidth',3)
legend('spline function','y=abs(x)','nodes')
title('Spline interpolation')
ylim([0 1.3]),grid
% turpinâjums
figure  % klûda
plot(x,f(x)-y,'g--','LineWidth',4)
title('Interpolation error'), grid

%% 6.uzdevums. Interpolâcija ar kubiskiem splainiem.
clear all, clc ,close all, format, format compact
f=@(x)sqrt(2+x.^7);
xnodes=[1:5]; ynodes=f(xnodes);
disp([xnodes' ynodes'])
x0=2.5;
int_res=interp1(xnodes,ynodes,x0,'spline'), f_res=f(x0)
int_error=abs(int_res-f_res)
disp('Atbilde:')
disp(['1) interpolâcijas rezultâts punktâ(x=2.5) = ', num2str(int_res)])
disp(['2) interpolâcijas kïûda  = ', num2str(int_error)])
disp('3) interpolâcijas kïûda punktâ(x=3) = 0, ')
disp('tâ kâ punkts x=3 ir interpolâcijas mezgla punkts')

%% Uzcevumi pastâvîgai risinâðanai
%% 1.uzdevums.
clear all  clc ,close all, format compact
xnodes = ([ 1:11 ] - 10) / 10;
fn = @(x) erf(x);
ynodes= fn(xnodes);
m = length(ynodes);

coef=ynodes;
for k = 2:m
    coef(k:m) = (coef(k:m) - coef(k-1:m-1))./(xnodes(k:m) - xnodes(1:m+1-k));
end

syms x
pol= coef(m); % polinoma konstuçðana
for k = m-1:-1:1
pol=pol*(x-xnodes(k))+coef(k);
end
polyn(x)=collect(pol)
polyn(x)=expand(pol)
coefpol=sym2poly(polyn(x)) % polinoma koeficienti

disp('Atbilde:')

x_pr = xnodes(1):0.01:xnodes(m);
plot(x_pr,polyn(x_pr),'r-',xnodes,ynodes,'o','LineWidth',3)
title('Newton''s interpolating polynomial')
legend('polynomial','nodes') 
grid on
figure
plot(x_pr,polyn(x_pr) - fn(x_pr),'g:','LineWidth',4)
grid on
%% 2.uzdevums.
clear all  clc ,close all, format compact
xnodes = 3:1:8;
fn = @(x) x.^3 .* sqrt(cos(x) + x);
ynodes= fn(xnodes);
m = length(ynodes);

coef=ynodes;
for k = 2:m
    coef(k:m) = (coef(k:m) - coef(k-1:m-1))./(xnodes(k:m) - xnodes(1:m+1-k));
end

syms x
pol= coef(m); % polinoma konstuçðana
for k = m-1:-1:1
pol=pol*(x-xnodes(k))+coef(k);
end
polyn(x)=collect(pol)
polyn(x)=expand(pol)
coefpol=sym2poly(polyn(x)) % polinoma koeficienti

disp('Atbilde:')

x_pr = xnodes(1):0.01:xnodes(m);
plot(x_pr,polyn(x_pr),'r-',xnodes,ynodes,'o','LineWidth',3)
title('Newton''s interpolating polynomial')
legend('polynomial','nodes') 
grid on
figure
plot(x_pr,polyn(x_pr) - fn(x_pr),'g:','LineWidth',4)
grid on

disp(['interpolācijas vērtība pie x = 6.5  ' num2str(double(polyn(6.5)))])
disp(['koeficiens pie 3. pakāpes ' num2str(coefpol(length(coefpol) - 3))])
%% 3.uzdevums.
clear all,  clc ,close all, format compact, clearvars
xnodes = 3:1:8;
f = @(x) x.^3 .* sqrt(cos(x) + x);
ynodes = f(xnodes);
m = length(ynodes);

x=linspace(3,8,200);
y=interp1(xnodes,ynodes,x,'spline');
% y_inter=interp1(xmezg,ymezg,x,'pchip')
% y_inter=spline(xmezg,ymezg,x)
plot(x,y,'r-',x,f(x),'k--',xnodes,ynodes,'o','LineWidth',3)
legend('spline function','y=abs(x)','nodes')
title('Spline interpolation'),grid
% turpinâjums
figure  % klûda
plot(x,f(x)-y,'g--','LineWidth',4)
title('Interpolation error'), grid

disp(['interpolācijas vērtība pie x = 6.5  ' num2str(double(y(6.5)))])
disp(['interpolācijas kļūda pie x = 6.5  ' num2str(double(y(6.5) - f(6.5)))])
disp(['interpolācijas vērtība pie x = 6.5  ' num2str(double(y(6.5)))])
%% 4.uzdevums. 
clear all  clc ,close all, format compact
xnodes = [ 2.1 2.6 3.1 3.6 4.1 4.6 5.1 5.6 ];
ynodes= [ 2.15 3.35 3.87 4.15 3.74 3.30 2.71 2.33 ];
m = length(ynodes);

coef=ynodes;
for k = 2:m
    coef(k:m) = (coef(k:m) - coef(k-1:m-1))./(xnodes(k:m) - xnodes(1:m+1-k));
end

syms x
pol= coef(m); % polinoma konstuçðana
for k = m-1:-1:1
pol=pol*(x-xnodes(k))+coef(k);
end
polyn(x)=expand(collect(pol));

x_pr = xnodes(1):0.01:xnodes(m);

iterp = interp1(xnodes,ynodes,x_pr,'spline');
plot(x_pr,polyn(x_pr),'r-',xnodes,ynodes,'o', x_pr, iterp, 'b--', 'LineWidth',3)
legend('polynomial','nodes', 'spline') 
grid on

%% 5.uzdevums. 

clear all  clc ,close all, format compact
xnodes = [ -1 -0.96 -0.86 -0.79 0.22 0.5 0.93 ];
ynodes=  [ -1 -0.151 0.894 0.986 0.895 0.5 -0.306 ];
m = length(ynodes);


coef=ynodes;
for k = 2:m
    coef(k:m) = (coef(k:m) - coef(k-1:m-1))./(xnodes(k:m) - xnodes(1:m+1-k));
end

syms x
pol= coef(m); % polinoma konstuçðana
for k = m-1:-1:1
pol=pol*(x-xnodes(k))+coef(k);
end
polyn(x)=expand(collect(pol));

x_pr = xnodes(1):0.01:xnodes(m);

iterp = interp1(xnodes,ynodes,x_pr,'spline');
plot(x_pr,polyn(x_pr),'r-',xnodes,ynodes,'o', x_pr, iterp, 'b--', 'LineWidth',3)
legend('polynomial','nodes', 'spline') 
grid on

%% 6.uzdevums.

clear all  clc ,close all, format compact
xnodes = -3:0.5:3;
ynodes=  [ -1 -1 -1 -0.98 -0.79 -0.46 0 0.46 0.79 0.98 1 1 1 ];
m = length(ynodes);

x_pr = xnodes(1):0.01:xnodes(m);

interp_spline = interp1(xnodes,ynodes,x_pr,'spline');
inter_erm = interp1(xnodes,ynodes,x_pr,'pchip');

plot(xnodes,ynodes,'o', x_pr, interp_spline, 'r', x_pr, inter_erm, 'b--', 'LineWidth',3)
legend('nodes', 'spline') 
grid on

%% 7.uzdevums. 

