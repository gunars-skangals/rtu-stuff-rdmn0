%% 1.uzdevums
clearvars
syms x
y(x) = 2 * acot(x^3 + 3);
z(x) = 2 * acot(x)^3 + 3;
fplot(y(x), 'b--', [-10, 10], 'LineWidth', 2.5)
hold on
fplot(z(x), 'r:', [-10, 10], 'LineWidth', 2)
legend('y = 2arcctg(x^3 + 3)', 'z = 2arcctg(x)^3 + 3')
xlabel x-ass
ylabel y-ass
grid on
title 'Divu funkciju grafiki'
axis ([-10, 10, -6, 12])
hold off
%% 2.uzdevums
clearvars
syms x
y(x) = cos(2^x);
fplot(y(x), 'r-x', [-2*pi, 2*pi], 'LineWidth', 2)
xlabel x-ass
ylabel y-ass
title 'Funkcijas 𝑦 = cos(2^x) grafiks'
%% 3. uzdevums
clc
clearvars
syms x
y1(x) = cos(x^3 + 3);
y2(x) = sin(x^3) + 3;
fplot(y1(x), 'r-', [-pi, pi], 'LineWidth', 3)
hold on
fplot(y2(x), 'g:', [-pi, pi], 'LineWidth', 3)
legend('y1 = cos(x^3 + 3)', 'y2 = sin(x^3) + 3')
xlabel x-ass
ylabel y-ass
grid on
title 'Trigonometriskās funkcijas'
hold off
%% 5.uzdevums
clc
clearvars
x = -2*pi:0.01:2*pi;
y1 = @(x) sin(x.^2 + 1);
y2 = @(x) cos(2*x).^2;
zim = plot(x, y1(x), 'b--', x, y2(x), 'r-', 'LineWidth', 2)
set(zim(2), 'LineWidth', 3);
legend('y = sin(x^2 + 1)', 'y2 = cos(2*x)^2')
xlabel x-ass
ylabel y-ass
grid on
title 'Divu funkciju grafiki'
%% 7.uzdevums
clc
clearvars
syms t
x1(t) = 2* cos(t)^2 + cos(t);
y1(t) = sin(2* t) + sin(t);
x2(t) = -cos(t) - 2 * cos (0.5*t);
y2(t) = -sin(t) + 2 * sin (0.5*t);
fplot(x1, y1, [-pi, pi], 'g', 'LineWidth', 3)
hold on 
fplot(x2, y2, [-pi, pi], 'b', 'LineWidth', 4)
xlabel x-ass
ylabel y-ass
grid on
legend ('x1(t)y1(t)', 'x2(t)y2(t)')
title 'Divu parametriski dotu funkciju grafiki'
hold off
%% 8.uzdevums
clc
clearvars
syms t
x(t) = 2 * (t^2 - 1)/(1 + t^2);
y(t) = 2*t*(t^2 - 1)/(1 + t^3);
f(t) = sin(t)*cos(2*t);
fplot(x, y, [0, 2*pi], 'r-x', 'LineWidth', 1.5)
hold on 
fplot(f, [0, 2*pi], 'g-o', 'LineWidth', 2)
xlabel x-ass
ylabel y-ass
grid on
legend ('x1(t)y1(t)', 'y=sin(x)*cos(2*x)')
title 'Divu funkciju grafiki'
axis ([-2,7,-1,2])
hold off
%% 9. uzdevums
clc
clearvars
x = -9:0.01:9;
y = @(x) acot(2*x+5).^3;
f = @(x) 2*acot(x.^3) + 5;
plot(x, y(x), 'r:', x, f(x), 'g-', 'LineWidth', 3)
xlabel x-ass
ylabel y-ass
grid on
legend ('y(x)', 'f(x)')
title 'Divu funkciju grafiki'
axis ([-10,10,-4,10])
%% 10.uzdevums
clc
clearvars
x = 3:0.01:11;
y = @(x) cos(x+1) .* log (x.^2 - 1);
f = @(x) (sin(x) + 1) .* (log(x.^2) - 1);
plot(x, y(x), 'k-', x, f(x), 'b--', 'LineWidth', 5)
xlabel x-ass
ylabel y-ass
grid on
legend ('y(x)', 'f(x)')
title 'Divu funkciju grafiki'
axis ([3,11,-6,8])
