%%
clc
clearvars
syms x
y = x^2 * sin(3*x);
atv = diff(y)
atv3 = diff(y, 3)
rez = simplify(atv3)
%%
clc
clearvars
syms t
x = 2 * cos(t) ^ 3;
y = 2 * sin(t) ^ 3;
atv = diff(y, t) / diff(x, t)
atv2 = diff(atv, t) / diff(x, t)
av2_rez = simplify(atv2)
%% 
clc
clearvars
syms x y
z = log(x ^3 + 3 * y ^ 2);
zx = diff(z,x);
zy = diff(z,y);
zxxx = diff(z, x, 3)
zxyy = diff(diff(z, x), y, 2)
%% 
clc
clearvars
syms x y
F = x*exp(x*y) + y^2 * tan(x) - 2;
yx = - diff(F,x) / diff(F, y)
%% 1.uzd a
clc
clearvars
syms x
y = x * asin(sqrt(x / (x + 1)))- sqrt(x) + atan(sqrt(x));
yatv = diff(y)
rez = simplify(yatv, 'IgnoreAnalyticConstraints', true)
%% 1.uzd b
clc
clearvars
syms x
y = sqrt(x^2 + 2) / x^2 - 1/sqrt(2) * log((sqrt(2) + sqrt(x^2 + 2)/x));
yatv = diff(y, x)
rez = simplify(yatv, 'IgnoreAnalyticConstraints', true)
%% 2.uzd a
clc
clearvars
syms x
y = exp(x) * (cos(2*x) - 3 * sin(2*x));
atv3 = diff(y, 3);
rez = simplify(atv3, 'IgnoreAnalyticConstraints', true)

%% 3.uzd a
clc
clearvars
syms t
x = cos(t) ^ 2;
y = tan(t) ^ 2;
atv1 = simplify( diff(y,t) / diff(x, t))
atv2 = simplify(diff(atv1, t) / diff(x, t))

%% 4.uzd a
clc
clearvars
syms x y
z = 1/(x * sqrt(y))*atan(exp(3*x)) + (tan(2*x))^(2*y) - 1/sin(3*y);
zx = simplify(diff(z,x))
zy = simplify(diff(z,y))
%% 5.uzd a
clc
clearvars
syms x y
z = exp(x * y ^ 4);
zxx = diff(z, x, 2)
zyy = diff(z, y, 2)
zyx = diff(diff(z, y), x)
%% 6.uzd a
clc
clearvars
syms x y
z = 1 / (1- x * y);
atv = diff(diff(z, x, 2), y)
%% 7.uzd a
clc
clearvars
syms x y
F = cos(x * y / 4) + exp(x) * y - sqrt(3 + x^2) - y ^ 2;
atv = - diff(F, x) / diff(F, y)
%% 8.uzd a
clc
clearvars
syms x
f = (x - 3) / 2 * sqrt(6 * x - x^2 - 8) + asin(sqrt(x / 2 - 1));
g = ((1 + x)* atan(sqrt(x)) - sqrt(x)) / x;
rez = simplify(diff(f) - diff(g))
%% 9.uzdevums
clc
clearvars
syms x
y = x ^ 2 * cos (x / 2);
yatv = diff(y);
hold off
fplot(y, 'b--', [-10, 10])
hold on
fplot(yatv, 'r:', [-10, 10])
