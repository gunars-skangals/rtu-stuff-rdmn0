%%
syms t
att = laplace(exp(-3*t)*sin(2*t)^2)

%% 
clc
clearvars
syms x p
y(x) = x^2 * sin(3*x);
att = laplace(y, x, p)
%% ilaplace
clc
clearvars
syms s 
X = (2*s + 1) / (s^3 - 4*s^2);
orig = ilaplace(X)
%% ilaplace
clc
clearvars
syms s x
X = (2*s + 1) / (s^3 - 4*s^2);
orig = ilaplace(X, s, x)
%% diferenciālvienādojuma risināšanas shēma
clc
clearvars
syms t x(t) s Xs
x1=diff(x,t);
vien = diff(x, t, 2) + 3 * diff(x, t) + 2 * x == 1 - t + 2 * t^2;
lapl = laplace(vien);
lapl1 = subs(lapl, [ laplace(x(t), t, s), x(0), x1(0)], [ Xs, 0, 1 ])
att = solve(lapl1, Xs);
orig = ilaplace(att)
%% 1.uzd a
clc
clearvars
syms t
laplace(7 * cos(t)^2 - 4 * sin(t) ^ 2)
%% 1.uzd b
clc
clearvars
syms t
laplace(exp(2* t) * (4 * exp(-5 * t) + 6 * sin(3*t) - 9*cos(6*t)))
%% 2.uzd a
clc
clearvars
syms s
ilaplace((8 * s ^ 2 - 19 * s - 73)/ (s + 3) / (s - 1) / (s - 4))
%% 2.uzd b
clc
clearvars
syms s
ilaplace((s - 3 * s ^ 3 + 15 * s ^ 2 - 3)/s^3/(s - 3))
%% 3.uzd a
clc
clearvars
syms t x(t) s Xs
x1=diff(x,t);
vien = diff(x, t, 2) -4 * diff(x, t) -5 * x == 9*exp(2*t);
lapl = laplace(vien);
lapl1 = subs(lapl, [ laplace(x(t), t, s), x(0), x1(0)], [ Xs, 5, 4 ]);
att = solve(lapl1, Xs);
orig = ilaplace(att)
fplot(orig)
%% 3.uzd b
clc
clearvars
syms t x(t) s Xs
x1=diff(x,t);
vien = diff(x, t, 2) + 3 * diff(x, t) == 18 *t - 6;
lapl = laplace(vien);
lapl1 = subs(lapl, [ laplace(x(t), t, s), x(0), x1(0)], [ Xs, -4, -7 ]);
att = solve(lapl1, Xs);
orig = ilaplace(att)
fplot(orig, [0, 5])