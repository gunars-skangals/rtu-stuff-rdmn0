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
x = [0:0.01:2];       % x vçrtîbu vektors
y = deval(sol,x);     % y vçrtîbu vektors
plot(x,y,'r','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on

%% 2.piemçrs.( 7.5.piemçrs(grâmatâ) 201.lpp. ) 
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
x = [0:0.01:1];          % x vçrtîbu vektors
y = deval(sol,x);        % y vÄ“rtÄ«bu vektors
plot(x,y,'r','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on
disp('Atbilde:')
disp([' funkcijas vçrtîba punktâ(0.2) =  ',num2str(y_02)])

%% 3.piemçrs( 7.6.piemçrs(grâmatâ) 204.lpp. )
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
x = [0:0.01:12];       % x vçrtîbu vektors
y = deval(sol,x);      % y1, y2 un y3 vçrtîbu vektors
plot(x,y(1,:),'r', x,y(2,:),'g', x,y(3,:),'b','LineWidth',3)
%xlim([-0.5 12.5]), ylim([-1 1.2])
legend('y1(x)','y2(x)','y3(x)'),xlabel('x'),ylabel('y')
title('Graphs of the functions y1(x),y2(x) un y3(x)'),grid on

%% 4.piemçrs. ( 7.7.piemçrs(grâmatâ) 206.lpp. )
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

%% 5.piemçrs.( 7.9.piemçrs(grâmata)  211.lpp. )
% Augstâku kârtu diferencialvienadojumi 
clear all, clc, close all,format compact
x_int = [1 4];             % intervals
y0 = [5/4; 3/4; 11/2];     % sakuma nosacîjumi
% numerical solution using ode45
sol = ode45(@st_fun_prob5,x_int,y0)

% turpinâjums
sol_x = sol.x', sol_y = sol.y'

% turpinajums
x = [1:0.01:4];         % x vçrtîbu vektors
y = deval(sol,x);       % y, y' un y'' vçrtîbu vektors
plot(x,y(1,:),'r',x,y(2,:),'g',x,y(3,:),'b','LineWidth',3)
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
%% 7.1. uzdevums
% Pirmâs kartas diferenciâlvienâdojumu sistçma ( function handle )
clear all, clc, close all, format compact
% definçsim diferenciâlvienâdojuma labo pusi:  
fun_prob2  = @(x,y)-y.*cos(x)+sin(x).*cos(x);  % function handle
x_int = [0 5];    % intervâls
y0 = 1;           % sakuma nosacîjums
% uzdevuma skaitliskâ risinâðana ar komandu ode45
sol = ode45(fun_prob2 ,x_int,y0) 

% turpinÄ?jums
y_02 = deval(sol,2)    % y vertîba punktâ 0.2
x = [0:0.01:5];          % x vçrtîbu vektors
y = deval(sol,x);        % y vÄ“rtÄ«bu vektors
plot(x,y,'r','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the function y(x)'),grid on
disp('Atbilde:')
disp([' funkcijas vçrtîba punktâ(2) =  ',num2str(y_02)])
%% 7.2. uzdevums
% Pirmâs kartas diferenciâlvienâdojumu sistçma ( function handle )
clear all, clc, close all, format compact
% definçsim diferenciâlvienâdojuma labo pusi:  
fun_prob2  = @(x,y)x.^2+y.^2;  % function handle
x_int = [0 0.8];    % intervâls
y0 = 1;           % sakuma nosacîjums
% uzdevuma skaitliskâ risinâðana ar komandu ode45
sol = ode45(fun_prob2 ,x_int,y0) 

% turpinÄ?jums
y_02 = deval(sol,0.5)    % y vertîba punktâ 0.2
x = [0:0.01:0.8];          % x vçrtîbu vektors
y = deval(sol,x);        % y vÄ“rtÄ«bu vektors
y1=x.^2+y.^2;
plot(x,y,'r',x,y1,'b','LineWidth',3)
xlabel('x-axis'), ylabel('y-axis')
title('Graph of the functions '),grid on
disp('Atbilde:')
disp([' funkcijas vçrtîba punktâ(2) =  ',num2str(y_02)])
%% 7.3. uzdevums
% Pirmâs kartas diferenciâlvienâdojumu sistçma
clear all, clc, close all,format compact
% definçsim diferenciâlvienâdojumu sistçmas labo pusi:  
dxy_dt = @(t,y) [y(2); 10*y(2)+11*y(1)] % function handle
t_int = [0,2];       % intervâls
y0 = [1;-1];           % sakuma nosacîjums
% uzdevuma skaitliskâ atrisinâðana ar komandu ode45
sol = ode45(dxy_dt,t_int,y0)

% turpinajums
t = [0:0.01:2];        % t vertîbu vektors
xy = deval(sol,t);      % x un y vertîbu vektors
plot(t,xy(1,:),'r',t,xy(2,:),'g','LineWidth',3)
legend('x(t)','y(t)'), xlabel('t'), ylabel('x(t),y(t)')
title('Graphs of the functions x(t),y(t)'),grid on

% turpinâjums
xy_2_7=deval(sol,0.7)
disp('Atbilde:')
disp(' funkciju vertîbas punktÄ?(2.7)')
disp(['  x(0.7) = ',num2str(xy_2_7(1))])
disp(['  y(0.7) = ',num2str(xy_2_7(2))])
%% 7.4. uzdevums
% Pirmâs kartas diferenciâlvienâdojumu sistçma
clear all, clc, close all,format compact
% definçsim diferenciâlvienâdojumu sistçmas labo pusi:  
dxy_dt = @(t,y) [2*y(1)-0.01*y(1).*y(2); -y(2)+0.01*y(1).*y(2)] % function handle
t_int = [0,1];       % intervâls
y0 = [1;1];           % sakuma nosacîjums
% uzdevuma skaitliskâ atrisinâðana ar komandu ode45
sol = ode45(dxy_dt,t_int,y0)

% turpinajums
t = [0:0.01:1];        % t vertîbu vektors
xy = deval(sol,t);      % x un y vertîbu vektors
plot(t,xy(1,:),'r',t,xy(2,:),'g','LineWidth',3)
legend('x(t)','y(t)'), xlabel('t'), ylabel('x(t),y(t)')
title('Graphs of the functions x(t),y(t)'),grid on

% turpinâjums
xy_2_7=deval(sol,0.3)
disp('Atbilde:')
disp(' funkciju vertîbas punktÄ?(0.3)')
disp(['  x(0.3) = ',num2str(xy_2_7(1))])
disp(['  y(0.3) = ',num2str(xy_2_7(2))])
%% 7.6. uzdevums
% Pirmâs kartas diferenciâlvienâdojumu sistçma
clear all, clc, close all,format compact
% definçsim diferenciâlvienâdojumu sistçmas labo pusi:  
dxy_dx = @(x,y) [y(2); -y(2)./x-y(1)] % function handle
x_int = [1,3];       % intervâls
y0 = [1;2];           % sakuma nosacîjums
% uzdevuma skaitliskâ atrisinâðana ar komandu ode45
sol = ode45(dxy_dx,x_int,y0)

% turpinajums
x = [1:0.01:3];        % t vertîbu vektors
xy = deval(sol,x);      % x un y vertîbu vektors
plot(x,xy(1,:),'r',x,xy(2,:),'g','LineWidth',3)
legend('x(t)','y(t)'), xlabel('t'), ylabel('x(t),y(t)')
title('Graphs of the functions x(t),y(t)'),grid on

% turpinâjums
xy_2_3=deval(sol,2.3)
disp('Atbilde:')
disp(' funkciju vertîbas punktÄ?(2.3)')
disp(['  x(2.3) = ',num2str(xy_2_3(1))])
disp(['  y(2.3) = ',num2str(xy_2_3(2))])
%% 7.8. uzdevums
% Pirmâs kartas diferenciâlvienâdojumu sistçma
clear all, clc, close all,format compact
% definçsim diferenciâlvienâdojumu sistçmas labo pusi:  
dxy_dx = @(x,y) [y(2); 3/2*y(1).^2] % function handle
x_int = [-2,1];       % intervâls
y0 = [1;-1];           % sakuma nosacîjums
% uzdevuma skaitliskâ atrisinâðana ar komandu ode45
sol = ode45(dxy_dx,x_int,y0)

% turpinajums
x = [-2:0.01:1];        % t vertîbu vektors
xy = deval(sol,x);      % x un y vertîbu vektors
plot(x,xy(1,:),'r',x,xy(2,:),'g','LineWidth',3)
legend('x(t)','y(t)'), xlabel('t'), ylabel('x(t),y(t)')
title('Graphs of the functions x(t),y(t)'),grid on

% turpinâjums
xy_2_3=deval(sol,-1.25)
disp('Atbilde:')
disp(' funkciju vertîbas punktÄ?(2.3)')
disp(['  x(-1.25) = ',num2str(xy_2_3(1))])
disp(['  y(-1.25) = ',num2str(xy_2_3(2))])
%% 7.8. uzdevums
% Pirmâs kartas diferenciâlvienâdojumu sistçma
clear all, clc, close all,format compact
% definçsim diferenciâlvienâdojumu sistçmas labo pusi:  
dxy_dx = @(x,y) [y(2); y(2)./x-y(1)./x.^2+8*x] % function handle
x_int = [1,3];       % intervâls
y0 = [3;4];           % sakuma nosacîjums
% uzdevuma skaitliskâ atrisinâðana ar komandu ode45
sol = ode45(dxy_dx,x_int,y0)

% turpinajums
x = [1:0.01:3];        % t vertîbu vektors
xy = deval(sol,x);      % x un y vertîbu vektors
plot(x,xy(1,:),'r',x,xy(2,:),'g','LineWidth',3)
legend('x(t)','y(t)'), xlabel('t'), ylabel('x(t),y(t)')
title('Graphs of the functions x(t),y(t)'),grid on

% turpinâjums
xy_1_86=deval(sol,1.86)
xy_2_68=deval(sol,2.68)
disp('Atbilde:')
disp(' funkciju vertîbas')
disp(['  y(1.86) = ',num2str(xy_1_86(1))])
disp(['  yprim(2.68) = ',num2str(xy_2_68(2))])