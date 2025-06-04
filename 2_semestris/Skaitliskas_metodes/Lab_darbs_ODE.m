% Parasto diferenci�lvien�?dojumu skaitlisk�s risin��anas metodes "
%              165. - 174. lpp.
%  Diferenci�lvien�dojumu sist�mu risin��nas metodes
%              174. - 189. lpp.

%% 1.piem�rs. ( 4.piem�rs(gr�mat�) 172.lpp. ) 
% Ko�� probl�ma.( M - fail)
% Pirm�s k�rtas diferenci�lvien�dojumu sist�ma 
clear all, clc, close all, format compact
x_int = [0 2];   % interv�ls
y0 = [1];        % s�kuma nosac�jums
% uzdevuma skaitlisk� risin��ana ar komandu ode45
sol = ode45(@st_fun_prob1,x_int,y0)

% turpin�jums
sol_x=sol.x'
sol_y=sol.y'

% turpinajums
t = [0:0.01:2];       % x v�rt�bu vektors
y = deval(sol,t);     % y v�rt�bu vektors
plot(t,y,'r','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on

%% 2.piem�rs.( 5.piem�rs(gr�mat�) 173.lpp. ) 
% Pirm�s kartas diferenci�lvien�dojumu sist�ma ( function handle )
clear all, clc, close all, format compact
% defin�sim diferenci�lvien�dojuma labo pusi:  
fun_prob2  = @(x,y) [x+y];  % function handle
x_int = [0 1];    % interv�ls
y0 = 0;           % sakuma nosac�jums
% uzdevuma skaitlisk� risin��ana ar komandu ode45
sol = ode45(fun_prob2 ,x_int,y0) 

% turpin�?jums
y_02 = deval(sol,0.2)    % y vert�ba punkt� 0.2
t = [0:0.01:1];          % x v�rt�bu vektors
y = deval(sol,t);        % y vērtību vektors
plot(t,y,'r','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on
disp('Atbilde:')
disp([' funkcijas v�rt�ba punkt�(0.2) =  ',num2str(y_02)])

%% 3.piem�rs( 6.piem�rs(gr�mat�) 177.lpp. )
% Pirm�s k�rtas diferenci�lvien�dojumu sist�ma( M-fails )
clear all, clc, close all,format compact 
x_int = [0,12];        % interv�ls
y0 = [0;1;1];          % s�kuma nosac�jums
% uzdevuma skaitlisk� atrisin��ana ar komandu ode45
sol = ode45(@st_fun_prob3,x_int,y0)

% turpin�jums
sol_x=sol.x'
sol_y=sol.y'

% turpinajums
t = [0:0.01:12];       % x v�rt�bu vektors
y = deval(sol,t);      % y1, y2 un y3 v�rt�bu vektors
plot(t,y(1,:),'r', t,y(2,:),'g', t,y(3,:),'b','LineWidth',3)
xlim([-0.5 12.5]), ylim([-1 1.2])
legend('y1(x)','y2(x)','y3(x)'),xlabel('x'),ylabel('y')
title('Graphs of the functions y1(x),y2(x) un y3(x)'),grid on

%% 4.piem�rs. ( 7.piem�rs(gr�mat�) 178.lpp. )
% Pirm�s kartas diferenci�lvien�dojumu sist�ma
clear all, clc, close all,format compact
% defin�sim diferenci�lvien�dojumu sist�mas labo pusi:  
dxy_dt = @(t,y) [-y(2)-y(1).^2; 2*y(1)-y(2)] % function handle
t_int = [0,10];       % interv�ls
y0 = [1;1];           % sakuma nosac�jums
% uzdevuma skaitlisk� atrisin��ana ar komandu ode45
sol = ode45(dxy_dt,t_int,y0)

% turpinajums
sol_t=sol.x'
sol_xy=sol.y'

% turpinajums
t = [0:0.01:10];        % t vert�bu vektors
xy = deval(sol,t);      % x un y vert�bu vektors
plot(t,xy(1,:),'r',t,xy(2,:),'g','LineWidth',3)
legend('x(t)','y(t)'), xlabel('t'), ylabel('x(t),y(t)')
title('Graphs of the functions x(t),y(t)'),grid on

% turpin�jums
xy_2_7=deval(sol,2.7)
disp('Atbilde:')
disp(' funkciju vert�bas punkt�?(2.7)')
disp(['  x(2.7) = ',num2str(xy_2_7(1))])
disp(['  y(2.7) = ',num2str(xy_2_7(2))])

%% 5.piem�rs.( 9.piem�rs(gr�mata)  183.lpp. )
% Augst�ku k�rtu diferencialvienadojumi 
clear all, clc, close all,format compact
x_int = [1 4];             % intervals
y0 = [5/4; 3/4; 11/2];     % sakuma nosac�jumi
% numerical solution using ode45
sol = ode45(@st_fun_prob5,x_int,y0)

% turpin�jums
sol_x = sol.x', sol_y = sol.y'

% turpinajums
t = [1:0.01:4];         % x v�rt�bu vektors
y = deval(sol,t);       % y, y' un y'' v�rt�bu vektors
plot(t,y(1,:),'r',t,y(2,:),'g',t,y(3,:),'b','LineWidth',3)
legend('y(x)','y''(x)','y''''(x)'),xlabel('x'),ylabel('y(x), y''(x), y''''(x)')
title('Graphs of  the functions y(x), y''(x) and y''''(x)'),grid on

% turpin�jums
y_value=deval(sol,1.87)

disp(' Atbilde:')
disp('  funkciju v�rt�bas punkt�?(1.87)')
disp(['  y(1.87) = ',num2str(y_value(1))])
disp(['  y''(1.87) = ',num2str(y_value(2))])
disp(['  y''''(1.87) = ',num2str(y_value(3))])

%% Uzdevumi patst�v�gai risin��anai

%% 1. uzdevums
% y' +y * cos x = sin x * cos x, y(0) = 1
% y' = sin x * cos x - y * cos x

clear all, clc, close all, format compact
x_int = [0 5];   % interv�ls
y0 = [1];        % s�kuma nosac�jums
% uzdevuma skaitlisk� risin��ana ar komandu ode45

fun = @(x, y) sin(x) .* cos(x) - y .* cos(x);
sol = ode45(fun, x_int, y0)

% turpin�jums
sol_x=sol.x';
sol_y=sol.y';

% turpinajums
t = [0:0.01:5];       % x v�rt�bu vektors
y = deval(sol,t);     % y v�rt�bu vektors
plot(t,y,'r','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on

y2 = deval(sol, [2])


%% 2.uzdevums
% y'  - x^2 - y ^ 2 = 0; y (0) = 1
clear all, clc, close all, format compact
x_int = [0 0.8];   % interv�ls
y0 = [1];        % s�kuma nosac�jums
% uzdevuma skaitlisk� risin��ana ar komandu ode45

fun = @(x, y) x .^ 2 + y ^ 2;
sol = ode45(fun, x_int, y0)

% turpin�jums
sol_x=sol.x';
sol_y=sol.y';

% turpinajums
t = [0:0.01:0.8];       % x v�rt�bu vektors
y = deval(sol,t);     % y v�rt�bu vektors

asda = zeros(1, length(t))
for i = 1:length(t)
    asda(i) = fun(t(i), y(i));
end

plot(t,y,'r','LineWidth',3)
hold on
plot(t, asda, 'b')
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on

y2 = deval(sol, [0.5])

%% 7.3
clear all, clc, close all, format compact
x_int = [0 2];   % interv�ls
y0 = [ 1; -1 ];        % s�kuma nosac�jums

fun = @(t,y) [y(2); 10 .* y(2) + 11 .* y(1)];

sol = ode45(fun, x_int, y0)

% turpin�jums
sol_x=sol.x';
sol_y=sol.y';

% turpinajums
t = [0:0.01:2];       % x v�rt�bu vektors
y = deval(sol,t);     % y v�rt�bu vektors

plot(t,y,'r','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on

y2 = deval(sol, [0.7])

%% 7.5
clear all, clc, close all, format compact
x_int = [0 1.5];   % interv�ls
y0 = [ 0; 1 ];        % s�kuma nosac�jums
% y '' + y / (x ^ 2 -4) = 0

fun = @(x, y) [ y(2) ; - y(1) ./ (x .^ 2 - 4) ];

sol = ode45(fun, x_int, y0)

% turpin�jums
sol_x=sol.x';
sol_y=sol.y';

% turpinajums
t = [0:0.01:1.5];       % x v�rt�bu vektors
y = deval(sol,t);     % y v�rt�bu vektors

plot(t,y(1, :),'r','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on

y2 = deval(sol, [0.7])

%% 7.6 
clear all, clc, close all, format compact
x_int = [1 3];   % interv�ls
y0 = [ 1; 2 ];        % s�kuma nosac�jums

fun = @(x, y) [ y(2); (- y(2) - x .* y(1)) ./ x ];

sol = ode45(fun, x_int, y0);

% turpin�jums
sol_x=sol.x';
sol_y=sol.y';

% turpinajums
t = [x_int(1):0.01:x_int(2)];       % x v�rt�bu vektors
y = deval(sol,t);     % y v�rt�bu vektors

plot(t,y(1, :),'r','LineWidth',3)
hold on
plot(t,y(2, :),'g','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on

y2 = deval(sol, [2.3])