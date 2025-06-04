% Parasto diferenciâlvienÄ?dojumu skaitliskâs risinâðanas metodes "
%              165. - 174. lpp.
%  Diferenciâlvienâdojumu sistçmu risinâðnas metodes
%              174. - 189. lpp.

%% 1.piemçrs. ( 4.piemçrs(grâmatâ) 172.lpp. ) 
% Koðî problçma.( M - fail)
% Pirmâs kârtas diferenciâlvienâdojumu sistçma 
clear all, clc, close all, format compact
x_int = [0 2];   % intervâls
y0 = [1];        % sâkuma nosacîjums
% uzdevuma skaitliskâ risinâðana ar komandu ode45
sol = ode45(@st_fun_prob1,x_int,y0)

% turpinâjums
sol_x=sol.x'
sol_y=sol.y'

% turpinajums
t = [0:0.01:2];       % x vçrtîbu vektors
y = deval(sol,t);     % y vçrtîbu vektors
plot(t,y,'r','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on

%% 2.piemçrs.( 5.piemçrs(grâmatâ) 173.lpp. ) 
% Pirmâs kartas diferenciâlvienâdojumu sistçma ( function handle )
clear all, clc, close all, format compact
% definçsim diferenciâlvienâdojuma labo pusi:  
fun_prob2  = @(x,y) [x+y];  % function handle
x_int = [0 1];    % intervâls
y0 = 0;           % sakuma nosacîjums
% uzdevuma skaitliskâ risinâðana ar komandu ode45
sol = ode45(fun_prob2 ,x_int,y0) 

% turpinÄ?jums
y_02 = deval(sol,0.2)    % y vertîba punktâ 0.2
t = [0:0.01:1];          % x vçrtîbu vektors
y = deval(sol,t);        % y vÄ“rtÄ«bu vektors
plot(t,y,'r','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on
disp('Atbilde:')
disp([' funkcijas vçrtîba punktâ(0.2) =  ',num2str(y_02)])

%% 3.piemçrs( 6.piemçrs(grâmatâ) 177.lpp. )
% Pirmâs kârtas diferenciâlvienâdojumu sistçma( M-fails )
clear all, clc, close all,format compact 
x_int = [0,12];        % intervâls
y0 = [0;1;1];          % sâkuma nosacîjums
% uzdevuma skaitliskâ atrisinâðana ar komandu ode45
sol = ode45(@st_fun_prob3,x_int,y0)

% turpinâjums
sol_x=sol.x'
sol_y=sol.y'

% turpinajums
t = [0:0.01:12];       % x vçrtîbu vektors
y = deval(sol,t);      % y1, y2 un y3 vçrtîbu vektors
plot(t,y(1,:),'r', t,y(2,:),'g', t,y(3,:),'b','LineWidth',3)
xlim([-0.5 12.5]), ylim([-1 1.2])
legend('y1(x)','y2(x)','y3(x)'),xlabel('x'),ylabel('y')
title('Graphs of the functions y1(x),y2(x) un y3(x)'),grid on

%% 4.piemçrs. ( 7.piemçrs(grâmatâ) 178.lpp. )
% Pirmâs kartas diferenciâlvienâdojumu sistçma
clear all, clc, close all,format compact
% definçsim diferenciâlvienâdojumu sistçmas labo pusi:  
dxy_dt = @(t,y) [-y(2)-y(1).^2; 2*y(1)-y(2)] % function handle
t_int = [0,10];       % intervâls
y0 = [1;1];           % sakuma nosacîjums
% uzdevuma skaitliskâ atrisinâðana ar komandu ode45
sol = ode45(dxy_dt,t_int,y0)

% turpinajums
sol_t=sol.x'
sol_xy=sol.y'

% turpinajums
t = [0:0.01:10];        % t vertîbu vektors
xy = deval(sol,t);      % x un y vertîbu vektors
plot(t,xy(1,:),'r',t,xy(2,:),'g','LineWidth',3)
legend('x(t)','y(t)'), xlabel('t'), ylabel('x(t),y(t)')
title('Graphs of the functions x(t),y(t)'),grid on

% turpinâjums
xy_2_7=deval(sol,2.7)
disp('Atbilde:')
disp(' funkciju vertîbas punktÄ?(2.7)')
disp(['  x(2.7) = ',num2str(xy_2_7(1))])
disp(['  y(2.7) = ',num2str(xy_2_7(2))])

%% 5.piemçrs.( 9.piemçrs(grâmata)  183.lpp. )
% Augstâku kârtu diferencialvienadojumi 
clear all, clc, close all,format compact
x_int = [1 4];             % intervals
y0 = [5/4; 3/4; 11/2];     % sakuma nosacîjumi
% numerical solution using ode45
sol = ode45(@st_fun_prob5,x_int,y0)

% turpinâjums
sol_x = sol.x', sol_y = sol.y'

% turpinajums
t = [1:0.01:4];         % x vçrtîbu vektors
y = deval(sol,t);       % y, y' un y'' vçrtîbu vektors
plot(t,y(1,:),'r',t,y(2,:),'g',t,y(3,:),'b','LineWidth',3)
legend('y(x)','y''(x)','y''''(x)'),xlabel('x'),ylabel('y(x), y''(x), y''''(x)')
title('Graphs of  the functions y(x), y''(x) and y''''(x)'),grid on

% turpinâjums
y_value=deval(sol,1.87)

disp(' Atbilde:')
disp('  funkciju vçrtîbas punktÄ?(1.87)')
disp(['  y(1.87) = ',num2str(y_value(1))])
disp(['  y''(1.87) = ',num2str(y_value(2))])
disp(['  y''''(1.87) = ',num2str(y_value(3))])

%% Uzdevumi patstâvîgai risinâðanai

%% 1. uzdevums
% y' +y * cos x = sin x * cos x, y(0) = 1
% y' = sin x * cos x - y * cos x

clear all, clc, close all, format compact
x_int = [0 5];   % intervâls
y0 = [1];        % sâkuma nosacîjums
% uzdevuma skaitliskâ risinâðana ar komandu ode45

fun = @(x, y) sin(x) .* cos(x) - y .* cos(x);
sol = ode45(fun, x_int, y0)

% turpinâjums
sol_x=sol.x';
sol_y=sol.y';

% turpinajums
t = [0:0.01:5];       % x vçrtîbu vektors
y = deval(sol,t);     % y vçrtîbu vektors
plot(t,y,'r','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on

y2 = deval(sol, [2])


%% 2.uzdevums
% y'  - x^2 - y ^ 2 = 0; y (0) = 1
clear all, clc, close all, format compact
x_int = [0 0.8];   % intervâls
y0 = [1];        % sâkuma nosacîjums
% uzdevuma skaitliskâ risinâðana ar komandu ode45

fun = @(x, y) x .^ 2 + y ^ 2;
sol = ode45(fun, x_int, y0)

% turpinâjums
sol_x=sol.x';
sol_y=sol.y';

% turpinajums
t = [0:0.01:0.8];       % x vçrtîbu vektors
y = deval(sol,t);     % y vçrtîbu vektors

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
x_int = [0 2];   % intervâls
y0 = [ 1; -1 ];        % sâkuma nosacîjums

fun = @(t,y) [y(2); 10 .* y(2) + 11 .* y(1)];

sol = ode45(fun, x_int, y0)

% turpinâjums
sol_x=sol.x';
sol_y=sol.y';

% turpinajums
t = [0:0.01:2];       % x vçrtîbu vektors
y = deval(sol,t);     % y vçrtîbu vektors

plot(t,y,'r','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on

y2 = deval(sol, [0.7])

%% 7.5
clear all, clc, close all, format compact
x_int = [0 1.5];   % intervâls
y0 = [ 0; 1 ];        % sâkuma nosacîjums
% y '' + y / (x ^ 2 -4) = 0

fun = @(x, y) [ y(2) ; - y(1) ./ (x .^ 2 - 4) ];

sol = ode45(fun, x_int, y0)

% turpinâjums
sol_x=sol.x';
sol_y=sol.y';

% turpinajums
t = [0:0.01:1.5];       % x vçrtîbu vektors
y = deval(sol,t);     % y vçrtîbu vektors

plot(t,y(1, :),'r','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on

y2 = deval(sol, [0.7])

%% 7.6 
clear all, clc, close all, format compact
x_int = [1 3];   % intervâls
y0 = [ 1; 2 ];        % sâkuma nosacîjums

fun = @(x, y) [ y(2); (- y(2) - x .* y(1)) ./ x ];

sol = ode45(fun, x_int, y0);

% turpinâjums
sol_x=sol.x';
sol_y=sol.y';

% turpinajums
t = [x_int(1):0.01:x_int(2)];       % x vçrtîbu vektors
y = deval(sol,t);     % y vçrtîbu vektors

plot(t,y(1, :),'r','LineWidth',3)
hold on
plot(t,y(2, :),'g','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on

y2 = deval(sol, [2.3])