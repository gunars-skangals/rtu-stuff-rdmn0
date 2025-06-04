%% 1.uzdevums
clc, clear all, format longG
integ = @(a,x) exp(a .* x) ./ (x + 2);
for a=1:11
    a_adj = (a - 1) / 10;
    fun= @(x) integ(a_adj, x);
    int_val(a) = integral(fun, 0, 1) - 0.5;   
end
a_pr=0:0.1:1;
plot(a_pr,int_val,'r','LineWidth',3)
xlabel('a'),ylabel('f(a)'),grid on

a_xn(1, :) = [0.4 0.5];
for i=1:2
    fun = @(x) integ(a_xn(i),x);
    f(i) = integral(fun, 0, 1) - 0.5;
end

k = 1;
delta = 1;
while (delta > 0.0001)

    k1=k+1;
    i=k+2;

    a_xn(i) = a_xn(k1) - f(k1) * (a_xn(k1) - a_xn(k)) / (f(k1) - f(k));
    fun = @(x) integ(a_xn(i), x);
    f(i) = integral(fun, 0, 1) - 0.5; 
    
    M_pr(k,1)=a_xn(i);
    M_pr(k,2)=f(i);

    k = k + 1;

    delta = abs(a_xn(i) - a_xn(k1));
end

disp('Atbilde: ')
disp([' a = ',num2str(M_pr(k - 1, 1), 8)])
%% 2.uzdevums
clc, clear all,format longG
f1 = @(x1,x2) x1 .^ 3 + x2 .^ 3 - 1;
f2 = @(x1,x2) sin(x1 + 0.4) - x2 - 0.2;

fimplicit(f1,[-5 5 -5 5],'r','LineWidth',3)
hold on
fimplicit(f2,[-5 5 -5 5],'b','LineWidth',3)
hold off
grid on,xlabel('x1'),ylabel('x2')

%% 
syms x1 x2

xapp= [ 0.8 0.7 ];
xapp_pr=xapp;
xpr=[x1 x2];

% funkcijas
fun=[ x1 ^ 3 + x2 ^ 3 - 1, sin(x1 + 0.4) - x2 - 0.2 ];
% parciālie atvasinājumi
fun_pr=[diff(fun(1),xpr(1)) diff(fun(1),xpr(2))
        diff(fun(2),xpr(1)) diff(fun(2),xpr(2))];

epsi=10^(-4);
k=0; 
sol_norm=1;
while sol_norm > epsi
    for i=1:2
        B(i,1)=-double(subs(fun(i),xpr,xapp));
        for j=1:2      
            A(i,j)=double(subs(fun_pr(i,j),xpr,xapp));
        end
    end
    xdelta=A\B; xapp=xapp+xdelta';
    c=double(subs(fun,xpr,xapp));
    sol_norm=norm(c);
    k=k+1;
    M_pr(k,1:2)=xapp(1:2); 
    M_pr(k,3)=sol_norm; 
end

% turpinâjums
disp(['saknes tuvinājumi: ', num2str(xapp_pr(1:2))])
disp('            x1                       x2                  kļūdas norma')
disp(M_pr)

% turpinâjums
disp(['Atbilde :'])
disp([' x1 = ',num2str(M_pr(k,1),5),',  x2 = ',num2str(M_pr(k,2),5)])
format

%% 3.uzdevums

K = 200;
N = 100;
x0 = 0.9;
rs = [ 0.5 1.9 2.3 2.6 3 ];

% x(j) = x_prev * exp(r * (1 - x_prev) / K) - izskatās kļūda formulā
for i = 1:length(rs)
    r = rs(i);
    
    x = ones(1, N) * x0;

    for j = 2:N
        x_prev = x(j-1);
        x(j) = x_prev * exp(r * (1 - x_prev) / K);
    end
    
    figure
    plot(x)
    title([ 'x pie r = ' num2str(r) ])
end
% secinājums: jo lielāks r, jo izteiktāks (bet joprojām neliels) izliekums

% ja būtu x_prev / K, nevis (1 - x_prev) / K, kā Ricker model formulā
for i = 1:length(rs)
    r = rs(i);
    
    x = ones(1, N) * x0;

    for j = 2:N
        x_prev = x(j-1);
        x(j) = x_prev * exp(r * (1 - x_prev) / K);
    end
    
    figure
    plot(x)
    title([ 'x pie r = ' num2str(r) ])
end

% secinājumi:
% sākotnējā r izvēle ļoti krasi maina x uzvedību:
% pie r = 0.5 - funkcija sasniedz K vērtību, vizuāli līdzīga loģistiskajai
% līknei
% pie r = 1.9 - svārstības ap K vērtību, kas norimst
% pie r = 2.3 - nerimstošanas svārstības starp diviem stāvokļiem, ~ K +- 125
% pie r = 2.6 nerimstošas svārstības starp 4 stāvokļiem, K +- 125, K +- 175
% pie r = 3 nerimstošas svārstības starp vairākiem stāvokļiem, virspusēji
% izskatās, ka tie atkārtojas, līdzīgi kā pie iepriekšējām r vērtībām

%% 4.uzdevums
clc, clear all
% 1.punkts
% V_t+1 = 0.998 * V_t + 0.45 * S1_t
% S1_t+1 = 0.002 * V_t + 0.1 * S1_t
% S2_t+1 = 0.45 * S1_t + 1 * S2_t
A = [ 0.998 0.45 0; 0.002 0.1 0; 0 0.45 1 ];

% 2.punkts
n0 = [ 100 0 0 ]';
n_t = n0;
for i=1:100
    n_t = A * n_t;
    if i == 10 || i == 50 || i == 100
        disp(['Pēc ' num2str(i) ' gadiem:'])
        disp(n_t)
    end
end
% 3.punkts
[P, D] = eig(A);
c = P \ n0;

for i = 1:3
    disp (['c_' num2str(i) ' = ' num2str(c(i))])
    disp (['V_' num2str(i) ' = '])
    disp (P(:, i))
    disp (['mu_' num2str(i) ' = ' num2str(D(i, i))])
end
% ja t tiecas uz bezgalību, (mu_1 ^ t) tiecas
% uz 1, pārējie tiecas uz 0 (mu_2, mu_3 < 1), attiecīgi sadalījumu nosaka
% pirmais saskaitāmais, sadalījums tiecas uz c_1 * V_1 - visi cilvēki stāvoklī S2

%% 4.punkts

% V_t+1 = 0.998 * V_t + 0.65 * S1_t
% S1_t+1 = 0.002 * V_t + 0.1 * S2_t
% S2_t+1 = 0.35 * S1_t + 0.9 * S2_t
A = [ 0.998 0.65 0; 0.002 0 0.1; 0 0.35 0.9 ];

n0 = [ 100 0 0 ]';
n_t = n0;
for i=1:100
    n_t = A * n_t;
    if i == 10 || i == 50 || i == 100
        disp(['Pēc ' num2str(i) ' gadiem:'])
        disp(n_t)
    end
end
% 3.punkts
[P, D] = eig(A);
c = P \ n0;

for i = 1:3
    disp (['c_' num2str(i) ' = ' num2str(c(i))])
    disp (['V_' num2str(i) ' = '])
    disp (P(:, i))
    disp (['mu_' num2str(i) ' = ' num2str(D(i, i))])
end
% ja t tiecas uz bezgalību: 
% mu_1 < 1, 1. saskaitamais tiecas uz 0
% mu_2 = 1, mu_2 ^ t būs 1
% mu_3 < 1, 3.saskaitāmais tiecas uz 0
% tātad sadalījums tiecas uz c_2 * V_2
disp('c_2 * V_2 = ')
disp(c(2) * P(:, 2))

%% 5.uzdevums 
clc, clear all
% a) ja vilku nav:
% V = 0
% dT/ dt = 2*T * (1 - 0.0001 * T) - 0.01 * T * 0
% tātad T turpinās pieaugt līdz (1 - 0.0001 * T) kļūs 0 vai negatīvs -
% sasniegs 10000 trušu un pieaugums apstāsies

% b)
syms T V t

T_0 = 1000;
V_0 = 200;

dT = 2 * T * (1 - 0.0001 * T) - 0.01 * T * V;
dV = -0.5 * V + 0.0001 * T * V;
dT_b = subs(subs(dT, T, T_0), V, V_0);
dV_b = subs(subs(dV, T, T_0), V, V_0);

[X,Y] = meshgrid(800:50:1200,0:50:400);

U_mesh = 2 * X * (1 - 0.0001 * X) - 0.01 * X * Y;
V_mesh = -0.5 * X + 0.0001 * X * Y;

quiver(X, Y, U_mesh, V_mesh, 1)
title('Fāzes diagramma ap sākuma stāvokli')
% secinājums - sākuma nosacījuma tuvumā abas populācijas samazinās

% c)
n = 50;
Tv = ones(1, n) * T_0;
Vv = ones(1, n) * V_0;

for i = 2:n
   T_prev = Tv(1, i - 1);
   V_prev = Vv(1, i - 1);
   deltaT = subs(subs(dT, T, T_prev), V, V_prev);
   deltaV = subs(subs(dV, T, T_prev), V, V_prev);
   Tv(i) = T_prev + deltaT;
   Vv(i) = V_prev + deltaV;
end

figure
plot(Vv)
grid on,xlabel('t'),ylabel('V')
title('Vilku skaits')
figure
plot(Tv)
grid on,xlabel('t'),ylabel('T')
title('Trušu skaits')

% d)
% vilku skaits samazinās un trušu skaits eksplozīvi pieaug līdz palikuši
% ~27 vilki un 8500 truši - pēc tam populācijas līdzsvarojas pie 100
% vilkiem un 5000 trušiem

