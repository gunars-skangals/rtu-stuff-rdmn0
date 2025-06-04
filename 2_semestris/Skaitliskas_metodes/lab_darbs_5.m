%% Aproksimâcija  & «Skaitliskâ integrçðana"


%% 1.piemçrs. 
% Aproksimâcija- polyfit
clear all, clc, close all,format compact
xpoints=[1:0.5:5]; ypoints=[5.5,4.7,4.9,6.4,8.1,9.3,11.4,13.5,16.0];
plot(xpoints,ypoints,'or','LineWidth',3),xlim([0.7 5.3])
ylim([3 17]),grid on, title('Data points')

% turpinâjums - aproksimâcija ar 1.kârtas polinomu
p=polyfit(xpoints,ypoints,1)
xapprox = 1:0.01:5; yapprox=polyval(p,xapprox);
figure
plot(xpoints,ypoints,'or',xapprox,yapprox,'k','LineWidth',3)
legend('data points','polynomial'), xlim([0.7 5.3]),ylim([3 19])
title('Approximation by a polynomial of degree 1'),grid on
disp('Atbilde:')
disp([' 1.k.polinoms: ', num2str(p(1)),'x+', num2str(p(2))])

% turpinâjums - aproksimâcija ar 2.kârtas polinomu
%{
...
%}

disp([' 2.k.polinoms: ', num2str(p(1)),'x^2', num2str(p(2)),'x+', num2str(p(3))])

% turpinâjums - aproksimâcija ar 3.kârtas polinomu 
%{
...
%}
disp([' 3.k.polinoms: ', num2str(p(1)),'x^3+', num2str(p(2)),'x^2', num2str(p(3)),'x+', num2str(p(4))])

%% 2.piemçrs.
% Aproksimâcija_fit
clear all, close all, clc
xpoints=[0,1,1.5,2:6,6.5,7]'; 
ypoints=[8.2,7.0,6.3,5.0,4.7,7,8,8.5,8.8,9]';
plot(xpoints,ypoints,'or','LineWidth',3)
xlim([-0.7 9]),ylim([4 10]),grid on
title('Data points')

% turpinâjums - aproksimâcija ar 3.kârtas polinomu
[p,reg]=fit(xpoints,ypoints,'poly3')
figure
plot(p,xpoints,ypoints)
xlim([-0.7 9]),ylim([4 10]),grid on
title('Approximation by a polynomial of degree 3')

% turpinâjums. Aproksimâcija ar trig.funkciju.
ft=fittype({'1','sin(x)','cos(x)'});
[p,reg]=fit(xpoints,ypoints,ft)
figure
plot(p,xpoints,ypoints)
xlim([-0.7 9]),ylim([4 10]),grid on
title('Approximation by trig.function y=a+b*sinx+c*cosx ')

% turpinâjums. Aproksimâcija ar trig.funkciju.
%{
...
%}

disp('Atbilde :')
disp('Labâkâ aproksimâcija ir 2.trig.funkcija ')
disp([num2str(p.a),num2str(p.b),'sinx+',num2str(p.c),'cosx+',...
    num2str(p.d),'sin2x',num2str(p.e),'cos2x'])

%% 3.piemçrs.
% Komanda int un integral.
clear all, clc,format compact
syms x
fun=@(x)sin(x)./x.*exp(-x.^2);
def_int=int(fun(x),0,2)
num_int=integral(fun,0,2)

%% 4.piemçrs. 
% Komanda int un integral-sin(1000x)
clear all, clc, close all, format compact, format long
syms x,fun=@(x)sin(10^3.*x);
def_int=double(int(fun(x),0,5))
num_int=integral(fun,0,5)
fplot(fun,[0, 5],'r')
xlim([-0.5 5.5]),ylim([-1.3 1.3]),grid on
title('y=sin(10^3x)')
format

%% 5.piemçrs. 
% Noteiktâ integrâïa aproksimâcija ar 4.kârtas polinomu 
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

% turpinâjums
% Uzzîmçt funk. f(x) un aproksimçjoðâ polinoma grafikus intervâlâ 
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

% turpinâjums
% aproksimçjoðâ polinoma koeficients pie otrâs pakâpes mainîgâ (3)  
coef_x_2=p.p3
% polinoma vçrtîba punktâ x0=1.2  (4)
x0=1.2; pol_value=p(x0)
% funkcijas vçrtîba punktâ x0=1.2  (5)
f_value=integral(y,0,x0)
% aproksimâcijas kïûdu punktâ x1=3.2 (6)
app_error=abs(p(3.2)-integral(y,0,3.2))
disp('Atbilde:')
disp([' 3) koeficients pie otrâs pakâpes mainîgâ = ',num2str(p.p3)])
disp([' 4) polinoma vçrtîba punktâ(1.2) = ',num2str(pol_value)])
disp([' 5) funkcijas vçrtîba punktâ(1.2) = ',num2str(f_value)])
disp([' 6) aproksimâcijas kïûda punktâ(3.2) = ',num2str(app_error)])

%% 6.piemçrs.
% Interpolçt integrâli ar treðâs kârtas Òûtona
% interpolâcijas polinomu ar mezgliem punktos 
clc, clear all, close all ,format compact
xnodes=[0.2,0.6,1.2,2.0];
syms t, y=@(t)sin(t)./t;
for i=1:4
    ynodes(i)=integral(y,0,xnodes(i));
end
disp('Funkcijas vçrtîbas interpolâcijas mezglos:')
disp(ynodes)
coef=ynodes;
for k = 2:4
    coef(k:4) = (coef(k:4) - coef(k-1:3))./(xnodes(k:4) - xnodes(1:5-k));
end

syms x, pol= coef(4); % polinoma konstuçðana
for k = 3:-1:1
pol=pol*(x-xnodes(k))+coef(k);
end
polyn(x)=collect(pol);
coefpol=sym2poly(polyn) % polinoma koeficienti

% turpinâjums
disp('Atbilde :')
disp([' 3.k.polinoms: ', num2str(coefpol(1)),'x^3', num2str(coefpol(2)),'x^2+', num2str(coefpol(3)),'x', num2str(coefpol(4))])

% Uzzîmçt funk. f(x) un interpolâcijas polinoma grafikus intervâlâ 
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

% turpinâjums
figure  % Uzzîmçt interpolâcijas kïûdas grafiku (4)
plot(x_pr,f-double(polyn(x_pr)),'g--','LineWidth',3)
title('Interpolation error'), xlim([0 2.2])
grid on

%% Uzdevumi patstâvîgai risinâðanai