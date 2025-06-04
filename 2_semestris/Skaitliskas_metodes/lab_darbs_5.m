%% Aproksim�cija  & �Skaitlisk� integr��ana"


%% 1.piem�rs. 
% Aproksim�cija- polyfit
clear all, clc, close all,format compact
xpoints=[1:0.5:5]; ypoints=[5.5,4.7,4.9,6.4,8.1,9.3,11.4,13.5,16.0];
plot(xpoints,ypoints,'or','LineWidth',3),xlim([0.7 5.3])
ylim([3 17]),grid on, title('Data points')

% turpin�jums - aproksim�cija ar 1.k�rtas polinomu
p=polyfit(xpoints,ypoints,1)
xapprox = 1:0.01:5; yapprox=polyval(p,xapprox);
figure
plot(xpoints,ypoints,'or',xapprox,yapprox,'k','LineWidth',3)
legend('data points','polynomial'), xlim([0.7 5.3]),ylim([3 19])
title('Approximation by a polynomial of degree 1'),grid on
disp('Atbilde:')
disp([' 1.k.polinoms: ', num2str(p(1)),'x+', num2str(p(2))])

% turpin�jums - aproksim�cija ar 2.k�rtas polinomu
%{
...
%}

disp([' 2.k.polinoms: ', num2str(p(1)),'x^2', num2str(p(2)),'x+', num2str(p(3))])

% turpin�jums - aproksim�cija ar 3.k�rtas polinomu 
%{
...
%}
disp([' 3.k.polinoms: ', num2str(p(1)),'x^3+', num2str(p(2)),'x^2', num2str(p(3)),'x+', num2str(p(4))])

%% 2.piem�rs.
% Aproksim�cija_fit
clear all, close all, clc
xpoints=[0,1,1.5,2:6,6.5,7]'; 
ypoints=[8.2,7.0,6.3,5.0,4.7,7,8,8.5,8.8,9]';
plot(xpoints,ypoints,'or','LineWidth',3)
xlim([-0.7 9]),ylim([4 10]),grid on
title('Data points')

% turpin�jums - aproksim�cija ar 3.k�rtas polinomu
[p,reg]=fit(xpoints,ypoints,'poly3')
figure
plot(p,xpoints,ypoints)
xlim([-0.7 9]),ylim([4 10]),grid on
title('Approximation by a polynomial of degree 3')

% turpin�jums. Aproksim�cija ar trig.funkciju.
ft=fittype({'1','sin(x)','cos(x)'});
[p,reg]=fit(xpoints,ypoints,ft)
figure
plot(p,xpoints,ypoints)
xlim([-0.7 9]),ylim([4 10]),grid on
title('Approximation by trig.function y=a+b*sinx+c*cosx ')

% turpin�jums. Aproksim�cija ar trig.funkciju.
%{
...
%}

disp('Atbilde :')
disp('Lab�k� aproksim�cija ir 2.trig.funkcija ')
disp([num2str(p.a),num2str(p.b),'sinx+',num2str(p.c),'cosx+',...
    num2str(p.d),'sin2x',num2str(p.e),'cos2x'])

%% 3.piem�rs.
% Komanda int un integral.
clear all, clc,format compact
syms x
fun=@(x)sin(x)./x.*exp(-x.^2);
def_int=int(fun(x),0,2)
num_int=integral(fun,0,2)

%% 4.piem�rs. 
% Komanda int un integral-sin(1000x)
clear all, clc, close all, format compact, format long
syms x,fun=@(x)sin(10^3.*x);
def_int=double(int(fun(x),0,5))
num_int=integral(fun,0,5)
fplot(fun,[0, 5],'r')
xlim([-0.5 5.5]),ylim([-1.3 1.3]),grid on
title('y=sin(10^3x)')
format

%% 5.piem�rs. 
% Noteikt� integr��a aproksim�cija ar 4.k�rtas polinomu 
clear all, clc, close all, format compact
syms t, y=@(t)sqrt(1+t.^2).*cos(t);
i=0;
for x=0:0.5:5
    i=i+1;
    f(i)=integral(y,0,x);
end
xpoints=[0:0.5:5]'; ypoints=f';
[p,reg]=fit(xpoints,ypoints,'poly4')
plot(p,xpoints,ypoints)
xlim([-0.6 6]),grid on
title('Approximation by a polynomial of degree 4')

% turpin�jums
% Uzz�m�t funk. f(x) un aproksim�jo�� polinoma grafikus interv�l� 
% [ 0,5 ] ar soli 0.1  (2)
i=0;
for x=0:0.1:5
    i=i+1;
    f(i)=integral(y,0,x);
end
figure
x_pr=[0:0.1:5];
plot(x_pr,p(x_pr),'r-',x_pr,f,'b--','LineWidth',3)
legend('polynomial','f(x)'), xlim([-0.6 6]),grid on
title('Graphs of the function and approximating polynomial')

% turpin�jums
% aproksim�jo�� polinoma koeficients pie otr�s pak�pes main�g� (3)  
coef_x_2=p.p3
% polinoma v�rt�ba punkt� x0=1.2  (4)
x0=1.2; pol_value=p(x0)
% funkcijas v�rt�ba punkt� x0=1.2  (5)
f_value=integral(y,0,x0)
% aproksim�cijas k��du punkt� x1=3.2 (6)
app_error=abs(p(3.2)-integral(y,0,3.2))
disp('Atbilde:')
disp([' 3) koeficients pie otr�s pak�pes main�g� = ',num2str(p.p3)])
disp([' 4) polinoma v�rt�ba punkt�(1.2) = ',num2str(pol_value)])
disp([' 5) funkcijas v�rt�ba punkt�(1.2) = ',num2str(f_value)])
disp([' 6) aproksim�cijas k��da punkt�(3.2) = ',num2str(app_error)])

%% 6.piem�rs.
% Interpol�t integr�li ar tre��s k�rtas ��tona
% interpol�cijas polinomu ar mezgliem punktos 
clc, clear all, close all ,format compact
xnodes=[0.2,0.6,1.2,2.0];
syms t, y=@(t)sin(t)./t;
for i=1:4
    ynodes(i)=integral(y,0,xnodes(i));
end
disp('Funkcijas v�rt�bas interpol�cijas mezglos:')
disp(ynodes)
coef=ynodes;
for k = 2:4
    coef(k:4) = (coef(k:4) - coef(k-1:3))./(xnodes(k:4) - xnodes(1:5-k));
end

syms x, pol= coef(4); % polinoma konstu��ana
for k = 3:-1:1
pol=pol*(x-xnodes(k))+coef(k);
end
polyn(x)=collect(pol);
coefpol=sym2poly(polyn) % polinoma koeficienti

% turpin�jums
disp('Atbilde :')
disp([' 3.k.polinoms: ', num2str(coefpol(1)),'x^3', num2str(coefpol(2)),'x^2+', num2str(coefpol(3)),'x', num2str(coefpol(4))])

% Uzz�m�t funk. f(x) un interpol�cijas polinoma grafikus interv�l� 
% [ 0.2,2 ] ar soli 0.1  (3)
i=0;
for x=0.2:0.1:2
    i=i+1;
    f(i)=integral(y,0,x);
end
x_pr=[0.2:0.1:2];
plot(x_pr,double(polyn(x_pr)),'r-',x_pr,f,'g--', xnodes,ynodes,'o','LineWidth',3)
title('Graphs of the function and interpolating polynomial')
legend('polynomial','f(x)','nodes'), ylim([0.2 2]), xlim([0 2.2])
grid on

% turpin�jums
figure  % Uzz�m�t interpol�cijas k��das grafiku (4)
plot(x_pr,f-double(polyn(x_pr)),'g--','LineWidth',3)
title('Interpolation error'), xlim([0 2.2])
grid on

%% Uzdevumi patst�v�gai risin��anai