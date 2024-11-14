%% fplot3
clc
clearvars
syms t
x = t * cos(2*t);
y = t * sin(2*t);
z = 0.5 * t;
fplot3(x, y, z, [0, 8*pi], "r--", ...
    'LineWidth', 5)
%% 
clc
clearvars
syms t
t = 0:0.1:8*pi;
x = t .* cos(2.*t); % uzdodam kaa handle, taapeec punkts pirms nelinearam darrbiibaam
y = t .* sin(2.*t);
z = 0.5 .* t;
plot3(x,y,z)
view([1,2,3])
%% virsmas - dekarta
clc
clearvars
syms x y
z(x, y) = sin(x.^2 + y.^2);
% fmesh, fsurf, mesh ,surf
% fsurf(z, [-pi, pi])
fsurf(z, [-pi, pi, -pi, pi])
%% virsmas 
clc
clearvars
[x,y] = meshgrid(-pi:0.1:pi, -pi:0.1:pi);
z = sin(x.^2 + y.^2);
mesh(x, y, z)
%% virsmas - parametriski
clc
clearvars
syms u v
x(u, v) = 2 * (1 - exp(u/(4*pi))) * cos(u) * cos(v / 2)^2;
y(u, v) = 2 * (-1 + exp(u/(4*pi))) * sin(u) * cos(v / 2)^2;
z(u, v) = 1 - exp(u / (3 * pi)) - sin(v) + exp(u/(4*pi))*sin(v);
fsurf(x, y, z, [0, 8*pi, -pi, pi])
%%
clc
clearvars
[u,v] = meshgrid(0:0.1:3, 0:0.1:pi);
x = u .* cos(v);
y = u .* sin(v);
z = 2 .* u;
surf(x, y, z)
%% liimenjliinijas - contour
clc
clearvars
syms x y
%fsurf(x*y, [-5, 5])
fcontour(x*y, [-5, 5]) % nav opciju att'elo'sanai
%%
clc
clearvars
[x, y] = meshgrid(-5:0.1:5, -5:0.1:5);
contour(x.*y, "--", "LineWidth", 3)
%% liimenjliinijas - iekraasotas
clc
clearvars
[x, y] = meshgrid(-5:0.1:5, -5:0.1:5);
contourf(x.*y)
%% noraadit veertiibas
clc
clearvars
[x, y] = meshgrid(-5:0.1:5, -5:0.1:5);
contour(x.*y, [-4,-3,-2,-1,1,2,3])
%% noraadit skaitu
clc
clearvars
[x, y] = meshgrid(-5:0.1:5, -5:0.1:5);
contourf(x.*y, 30)
%% paraadiit veertiibas
clc
clearvars
[x, y] = meshgrid(-5:0.1:5, -5:0.1:5);
contour(x, y, x.*y, "showText", "on")
%% 1.uzdevums
clc
clearvars
syms x y
z(x, y) = x .* exp(-x.^2 - y.^2);
% fmesh, fsurf, mesh ,surf
% fsurf(z, [-pi, pi])
fsurf(z, [-2, 2])
view([-3,2,1])
xlabel x
ylabel y
zlabel z
grid on
title 'Necaurspīdīga virsma'
%% 2.uzd - skip
%% 3.uzdevums
clc
clearvars
[x,y] = meshgrid(-6:0.2:6, -5:0.1:5);
f = log(1 + sqrt(x .^2 + 2.* y  .^ 2)) + cos(sqrt(1 + x.^4)) + y.^3 .* exp(x .* y);
contour(f, [-2,0,1,3,5,7], "LineWidth", 2)
title 'Līmeņlīnijas ar vērtībām -2, 0, 1, 3, 5, 7'
%% 4.uzdevums
clc
clearvars
syms t
a = 0.3;
x = cos(t) / (sqrt(1+ a^2*t^2));
y = sin(t) / (sqrt(1+ a^2*t^2));
z = a*t / (sqrt(1+ a^2*t^2));
fplot3(x, y, z, [-12*pi, 12*pi], "r", 'LineWidth', 2)
xlabel x-ass
ylabel y-ass
zlabel z-ass
grid on
title '3D līnija'
%% 5.uzdevums
clc
clearvars
[x,y] = meshgrid(-10:0.1:10, -10:0.1:10);
z = sin(sqrt(x.^2 + y.^2 + 4))./sqrt(x.^2 + y.^2 + 4);
mesh(z)
title 'Caurspīdīga virsma'
xlabel x-ass
ylabel y-ass
zlabel z-ass
%% 6.uzdevums
clc
clearvars
[t,u] = meshgrid(0:0.1:2*pi, -pi/2:0.1:pi/2);
x = 5 .* cos(t) .* cos(u);
y = 3 .* sin(t) .* cos(u);
z = 7 .* sin(u);
surf(x,y,z)
view([1,-1,3])
title 'Necaurspīdīga virsma'
xlabel x
ylabel y
zlabel z

%% 6.uzd vēlereiz
clc
clearvars
syms t u
x = 5 .* cos(t) .* cos(u);
y = 3 .* sin(t) .* cos(u);
z = 7 .* sin(u);
fsurf(x,y,z, [0, 2*pi, -pi/2, pi/2])
view([1,-1,3])
title 'Necaurspīdīga virsma'
xlabel x
ylabel y
zlabel z
%% 7.uzd 
clc
clearvars
[x, y] = meshgrid(-2:0.1:4, -3:0.1:-1);
contour(x.*y, )
