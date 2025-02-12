% fplot - funkciju var uzdot ka'simbolisku, kā hanflde vai ka'izteiksmi
% plot
syms x
y(x) = sin(x)^2;
fplot(y, [-pi,pi] ) % default intervāls -5, 5

%% fplot ar handle
clearvars
syms x
y = @(x) sin(x)^2;
fplot(y, [-pi,pi])

%% fplot ar handle - vektorizēts
clc
clearvars
syms x
y = @(x) sin(x).^2;
fplot(y, [-pi,pi])
%% 
clearvars
syms x
fplot(cos(x^2), [-pi,pi])
figure
y(x) = sin(x)^2;
fplot(y, [-pi,pi] ) % default intervāls -5, 5
% fplot atļauj tikai 1 līniju 1 koordina'tu sistēmā
%% 
clearvars
syms x
y1 = sin(x^2);y2 = cos(x^2)
fplot(y1,[-pi,pi])
hold on % hold off lai izslēgtu
fplot(y2,[-pi,pi])
hold off
%% plot
clearvars
syms x
x1 = -pi:0.01:pi; % 0.1 arī standarta
y(x) = sin(x)^2;
plot(x1, y(x1))
%%
clc
x1=linspace(-pi,pi,100);
y = @(x) cos(x).^2
plot(x1, y(x1))
%%
clc
clearvars
x1 = -pi:0.01:pi; % 0.1 arī standarta
y1 = @(x) sin(x).^2;
y2 = @(x) cos(x.^2);
plot(x1, y1(x1), x1, y2(x1))
%% 
clearvars
syms x
y(x) = sin(x)^2;
fplot(y, [-pi,pi], 'r:o','LineWidth', 2) % melnai krāsai lieto burtu k
title 'y=sin(x^2)'
xlabel 'x'
ylabel 'y'
%%
clc
clearvars
x1 = -pi:0.01:pi; % 0.1 arī standarta
y1 = @(x) sin(x).^2;
y2 = @(x) cos(x.^2);
plot(x1, y1(x1), 'g--', x1, y2(x1), 'm-.', 'LineWidth', 3)
% dažādiem platumiem - zīmē atsevišķi un lieto set komandu
%% parametriski
clc
clearvars
syms t
x(t) = cos(t)^3
y(t) = sin(t)^3
fplot(x, y, [0,2*pi])
%% parametriski
clc
clearvars
t = 0:0.1:2*pi;
x = cos(t).^3;
y = sin(t).^3;
plot(x, y, 'k--*', 'LineWidth', 2.5)