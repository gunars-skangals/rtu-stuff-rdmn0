%% Nelineâro vienâdojumu risinâðana "     

%% 1.piemçrs. 
% Atrisinât vienâdojumu ar Òûtona metodi 
clc, clear all, format compact, close all
format longG
f = @(x)x.^3-2.*x-5;
x_pr = -10:0.01:10;
plot(x_pr,f(x_pr),'r','LineWidth',3)
title(['Graph of the function on',...
       ' the interval [-10,10]'])
grid on

% turpinâjums
figure
x_pr = -3:0.01:3;
plot(x_pr,f(x_pr),'r','LineWidth',3)
title('Graph of the function on the interval [-3,3]')
grid on

% turpinâjums
x_app =[2]; % saknes tuvinâjums 
iter=6;     % max iterâciju skaits
syms x, fpr(x)=diff(f(x),x);
xn=x_app; M =zeros(iter,2);
for i=1:iter
   xn=xn-f(xn)/double(fpr(xn));
   M(i,1)=xn; M(i,2)=f(xn);
end
disp(['saknes tuvinâjums:  ',num2str(x_app(:))])
disp('              x                 f(x)')
disp(M)

% turpinâjums
disp('Atbilde:')
disp([' x = ',num2str(M(iter,1),6)]),format

%% 2.piemçrs. 
% Atrisinât vienâdojumu ar Òûtona metodi
clc, clear all,close all,format compact,format longG
f = @(x)x.^3-7*x.^2+11.01*x-5;
x_pr=-10:0.01:10;
plot(x_pr,f(x_pr),'r','LineWidth',3)
title('Graph of the function on the interval [-10,10]')
grid on

% turpinâjums
figure
x_pr = -1:0.01:6;
plot(x_pr,f(x_pr),'r','LineWidth',3)
title(['Graph of the function on the interval [-1,6]'])
grid on

% turpinâjums
figure
x_pr =0.8:0.01:1.2;
plot(x_pr,f(x_pr),'r','LineWidth',3)
title(['Graph of the function on the interval [0.8, 1.2]'])
grid on, grid minor

% turpinâjums
x_app =[0.95 1.05 5]; % sakòu tuvinâjumi 
iter=6;
syms x, fpr(x)=diff(f(x),x);
for j=1:3
    M =zeros(iter,2);
    xn=x_app(j);
    for i=1:iter
       xn=xn-f(xn)/double(fpr(xn));
       M(i,1)=xn; M(i,2)=f(xn);
    end
    disp(['saknes tuvinâjums:  ',num2str(x_app(j))])
    disp('              x                 f(x)')
    disp(M)
    xapp_val(j)=M(iter,1);
end  
disp('Atbilde:')
disp([' x1 = ',num2str(xapp_val(1),5),', x2 = ',num2str(xapp_val(2),6),...
      ',  x3 = ',num2str(xapp_val(3),6)])
format
 
%% 3.piemçrs.
% Atrisinât vienâdojumu ar hordu metodi .
clc, clear all,close all,format compact, format longG
f = @(x)x.^3-2*x-5;
x_pr = -10:0.01:10;
plot(x_pr,f(x_pr),'r','LineWidth',3)
title(['Graph of the function on the interval [-10,10]'])
grid on

% turpinâjums
figure
x_pr = -3:0.01:3;
plot(x_pr,f(x_pr),'r','LineWidth',3)
title(['Graph of the function on the interval [-3,3]'])
grid on

% turpinâjums
iter = 6;
xn = zeros(1,iter+2);
xn = [2 2.2]; 
for k=1:iter
    k_1=k+1; k_2=k+2;
    xn(k_2)=xn(k_1)-f(xn(k_1))*(xn(k_1)-xn(k))/(f(xn(k_1))-f(xn(k)));
end
M_pr(:,1)=xn(3:iter+2); M_pr(:,2)=f(xn(3:iter+2));
disp(['saknes tuvinâjumu intervâls: ', num2str(xn(1:2))])
disp('              x                 f(x)')
disp(M_pr)

% turpinâjums
disp('Atbilde: ')
disp([' x = ',num2str(M_pr(iter,1),5)])
format

%% 4.piemçrs.
% Izmantojot hordu metodi, atrast minimâlo pozitîvo skaitli a,
% kas apmierina vienâdojumu
clc, clear all,close all, format compact, format longG
integ=@(a,x)(2*x+3).^(1/3).*log(a.*x);
for a=1:10
    fun=@(x)integ(a,x);
    int_val(a) = integral(fun,2,3)-5;   
end
a_pr=1:10;
plot(a_pr,int_val,'r','LineWidth',3)
title(['Graph of the function f(a) on the interval [1 10]'])
xlabel('a'),ylabel('f(a)'),grid on

% turpinâjums
iter=6;
a_xn =zeros(1,iter+2);
a_xn=[4.7 5.1]; 
for i=1:2
    fun=@(x)integ(a_xn(i),x); f(i)=integral(fun,2,3)-5;
end

% turpinâjums
for k=1:iter
    k1=k+1; i=k+2;
    a_xn(i)=a_xn(k1)-f(k1)*(a_xn(k1)-a_xn(k))/(f(k1)-f(k));
    fun = @(x)integ(a_xn(i),x);  f(i)=integral(fun,2,3)-5; 
    M_pr(k,1)=a_xn(i);M_pr(k,2)=f(i);
end
disp(['skaitïa a  tuvinâjumu intervâls: ', num2str(a_xn(1:2))])
disp('            a                  f(x)')
disp(M_pr)

% turpinâjums
disp('Atbilde: ')
disp([' a = ',num2str(M_pr(iter,1),6)])
format

%% 5.piemçrs. 
% Izmantojot hordu metodi, atrast minimâlo pozitîvo skaitli a,
% kas apmierina vienâdojumu
clc, clear all,close all, format compact, format longG
integ=@(a,x)(2+a.*x.^4).^(1/2).*sin(x.^2);
for a=1:10
    fun=@(x)integ(a,x);
    int_val(a)=integral(fun,2,4)-5;   
end
a_pr=1:10;
plot(a_pr,int_val,'r','LineWidth',3)
title(['Graph of the function f(a) on the interval [1 10]'])
xlabel('a'),ylabel('f(a)'),grid on

%% turpinâjums
for a=11:20
    fun=@(x)integ(a,x);
    int_val(a)=integral(fun,2,4)-5;   
end
figure
a_pr = 11:20;
plot(a_pr,int_val(11:20),'r','LineWidth',3)
title(['Graph of the function f(a) on the interval [11 20]'])
xlabel('a'),ylabel('f(a)'),grid on

%% turpinâjums
iter=6;
a_xn =zeros(1,iter+2);
a_xn=[14 14.5]; 
for i=1:2
    fun=@(x)integ(a_xn(i),x); f(i)=integral(fun,2,4)-5;
end
for k=1:iter
    k1=k+1; i=k+2;
    a_xn(i)=a_xn(k1)-f(k1)*(a_xn(k1)-a_xn(k))/(f(k1)-f(k));
    fun=@(x)integ(a_xn(i),x);  f(i)=integral(fun,2,4)-5; 
    M_pr(k,1)=a_xn(i); M_pr(k,2)=f(i);
end
disp(['skaitïa a  tuvinâjumu intervâls: ', num2str(a_xn(1:2))])
disp('            a                  f(x)')
disp(M_pr)

% turpinâjums
disp('Atbilde: ')
disp([' a = ',num2str(M_pr(iter-1,1),8)])
format

%% Tirgus lîdsvars

%% Uzdevumi patstâvîgai risinâðanai
%% 1.UZDEVUMS
clc, clear all, format compact, close all
format longG
f = @(x) x / 4 + 3 - sin (x + 0.2);
x_pr = -10:0.01:10;
plot(x_pr,f(x_pr),'r','LineWidth',3)
title(['Graph of the function on',...
       ' the interval [-10,10]'])
grid on

%% turpinâjums
figure
x_pr = -17:0.01:-8;
plot(x_pr,f(x_pr),'r','LineWidth',3)
title('Graph of the function on the interval [-3,3]')
grid on

%% turpinâjums - sakne 1
x_app =[-15]; % saknes tuvinâjums 1
iter=6;     % max iterâciju skaits
syms x, fpr(x)=diff(f(x),x);
xn=x_app; M =zeros(iter,2);
for i=1:iter
   xn=xn-f(xn)/double(fpr(xn));
   M(i,1)=xn; M(i,2)=f(xn);
end
disp(['saknes tuvinâjums:  ',num2str(x_app(:))])
disp('              x                 f(x)')
disp(M)

% turpinâjums
disp('Atbilde:')
disp([' x1 = ',num2str(M(iter,1),6)]),format
%% turpinâjums - sakne 2
x_app =[-13]; % saknes tuvinâjums 1
iter=6;     % max iterâciju skaits
syms x, fpr(x)=diff(f(x),x);
xn=x_app; M =zeros(iter,2);
for i=1:iter
   xn=xn-f(xn)/double(fpr(xn));
   M(i,1)=xn; M(i,2)=f(xn);
end
disp(['saknes tuvinâjums:  ',num2str(x_app(:))])
disp('              x                 f(x)')
disp(M)

% turpinâjums
disp('Atbilde:')
disp([' x2 = ',num2str(M(iter,1),6)]),format
%% turpinâjums - sakne 3
x_app =[-10]; % saknes tuvinâjums 1
iter=6;     % max iterâciju skaits
syms x, fpr(x)=diff(f(x),x);
xn=x_app; M =zeros(iter,2);
for i=1:iter
   xn=xn-f(xn)/double(fpr(xn));
   M(i,1)=xn; M(i,2)=f(xn);
end
disp(['saknes tuvinâjums:  ',num2str(x_app(:))])
disp('              x                 f(x)')
disp(M)

% turpinâjums
disp('Atbilde:')
disp([' x3 = ',num2str(M(iter,1),6)]),format
%% 2. UZDEVUMS
clc, clear all, format compact, close all
format longG
f = @(x) x .^ 2 - 4 * x + 4 - 2 * cos(2*x);
x_pr = -10:0.01:10;
plot(x_pr,f(x_pr),'r','LineWidth',3)
title(['Graph of the function on',...
       ' the interval [-10,10]'])
grid on
%% 
x_pr = -2:0.01:4;
plot(x_pr,f(x_pr),'r','LineWidth',3)
title(['Graph of the function on',...
       ' the interval [-10,10]'])
grid on
%% turpinâjums - sakne 1
x_app =[2.3]; % saknes tuvinâjums 1
iter=6;     % max iterâciju skaits
syms x, fpr(x)=diff(f(x),x);
xn=x_app; M =zeros(iter,2);
for i=1:iter
   xn=xn-f(xn)/double(fpr(xn));
   M(i,1)=xn; M(i,2)=f(xn);
end
disp(['saknes tuvinâjums:  ',num2str(x_app(:))])
disp('              x                 f(x)')
x1 = M(iter,1)
%% turpinâjums - sakne 2
x_app =[3.3]; % saknes tuvinâjums 1
iter=6;     % max iterâciju skaits
syms x, fpr(x)=diff(f(x),x);
xn=x_app; M =zeros(iter,2);
for i=1:iter
   xn=xn-f(xn)/double(fpr(xn));
   M(i,1)=xn; M(i,2)=f(xn);
end
disp(['saknes tuvinâjums:  ',num2str(x_app(:))])
disp('              x                 f(x)')
x2 = M(iter,1)

%% 
disp(['Atbilde: x1 = ' num2str(x1, 6) ', x2 = ' num2str(x2, 6)])
%% 3. UZDEVUMS
clc, clear all,close all, format compact, format longG
f = @(x) sin(x) - 2 * log(x);

a_pr=-5:0.1:3;
plot(a_pr,f(a_pr),'r','LineWidth',3)
title(['Graph of the function f(a) on the interval [1 10]'])
xlabel('a'),ylabel('f(a)'),grid on

%%
iter = 6;
xn = zeros(1,iter+2);
xn = [ 1.5 1.9 ];
for k=1:iter
    k_1=k+1; k_2=k+2;
    xn(k_2)=xn(k_1)-f(xn(k_1))*(xn(k_1)-xn(k))/(f(xn(k_1))-f(xn(k)));
end
M_pr(:,1)=xn(3:iter+2); M_pr(:,2)=f(xn(3:iter+2));
disp(['saknes tuvinâjumu intervâls: ', num2str(xn(1:2))])
disp('              x                 f(x)')
disp(M_pr)
x1 = M_pr(iter, 1)
%% skatāmies tikai x > 0
% iter = 6;
% xn = zeros(1,iter+2);
% xn = [ -0.9 -0.5 ];
% for k=1:iter
%     k_1=k+1; k_2=k+2;
%     xn(k_2)=xn(k_1)-f(xn(k_1))*(xn(k_1)-xn(k))/(f(xn(k_1))-f(xn(k)));
% end
% M_pr(:,1)=xn(3:iter+2); M_pr(:,2)=f(xn(3:iter+2));
% disp(['saknes tuvinâjumu intervâls: ', num2str(xn(1:2))])
% disp('              x                 f(x)')
% disp(M_pr)
% x1 = M_pr(iter, 1)
x_pr = -2:0.01:4;
plot(x_pr,f(x_pr),'r','LineWidth',3)
title(['Graph of the function on',...
       ' the interval [-10,10]'])
grid on
%% 
disp(['Atbilde: x1 = ' num2str(x1, 6)])
%% 4. UZDEVUMS
clc, clear all,close all, format compact, format longG
integ = @(a,x) sqrt(a .* sin(x) + x .^ 3);
for a=1:10
    fun=@(x) integ(a,x);
    int_val(a)=integral(fun,0,1)-1;   
end
a_pr=1:10;
plot(a_pr,int_val,'r','LineWidth',3)
title(['Graph of the function f(a) on the interval [1 10]'])
xlabel('a'),ylabel('f(a)'),grid on

%% turpinâjums
iter=6;
a_xn =zeros(1,iter+2);
a_xn=[0.5 1.5]; 
for i=1:2
    fun=@(x)integ(a_xn(i),x); f(i)=integral(fun,0,1)-1;
end
for k=1:iter
    k1=k+1; i=k+2;
    a_xn(i)=a_xn(k1)-f(k1)*(a_xn(k1)-a_xn(k))/(f(k1)-f(k));
    fun=@(x)integ(a_xn(i),x); 
    f(i)=integral(fun,0,1)-1; 
    M_pr(k,1)=a_xn(i); M_pr(k,2)=f(i);
end
disp(['skaitïa a  tuvinâjumu intervâls: ', num2str(a_xn(1:2))])
disp('            a                  f(x)')
disp(M_pr)

% turpinâjums
disp('Atbilde: ')
disp([' a = ',num2str(M_pr(iter-1,1),8)])
format
%% 5.uzdevums
clc, clear all,close all, format compact, format longG
integ = @(a,x) sqrt(a .* x .^ 3 + 1);

vv = zeros(10, 2);

for i=1:10
    a = i - 1;
    fun=@(x) integ(a,x);
    vv(i, 1) = a;
    vv(i, 2) = integral(fun,1,2) - 2;   
end

plot(vv(:,1),vv(:, 2),'r','LineWidth',3)
title(['Graph of the function f(a) on the interval [1 10]'])
xlabel('a'),ylabel('f(a)'),grid on

%% turpinâjums
iter=6;
a_xn =zeros(1,iter+2);
a_xn=[0.5 1]; 
for i=1:2
    fun=@(x)integ(a_xn(i),x);
    f(i)=integral(fun,1,2) - 2;
end
for k=1:iter
    k1=k+1; i=k+2;
    a_xn(i)=a_xn(k1)-f(k1)*(a_xn(k1)-a_xn(k))/(f(k1)-f(k));
    fun=@(x)integ(a_xn(i),x);
    f(i)=integral(fun,1,2) - 2;
    M_pr(k,1)=a_xn(i); M_pr(k,2)=f(i);
end
disp(['skaitïa a  tuvinâjumu intervâls: ', num2str(a_xn(1:2))])
disp('            a                  f(x)')
disp(M_pr)

% turpinâjums
disp('Atbilde: ')
disp([' a = ',num2str(M_pr(iter-1,1),8)])
format
%% 6.uzdevums
clc, clear all, format compact, close all
format longG
f = @(x) sin(x) .* log(x + 5) - sqrt(x + 1);
x_pr = -1:0.01:10;
plot(x_pr,f(x_pr),'r','LineWidth',3)
title(['Graph of the function on',...
       ' the interval [-10,10]'])
grid on
%% 
x_pr = 0:0.01:3;
plot(x_pr,f(x_pr),'r','LineWidth',3)
title(['Graph of the function on',...
       ' the interval [-10,10]'])
grid on
%% turpinâjums - sakne 1
x_app =[ 0.88 ]; % saknes tuvinâjums 1
iter=6;     % max iterâciju skaits
syms x, fpr(x)=diff(f(x),x);
xn=x_app; M =zeros(iter,2);
for i=1:iter
   xn=xn-f(xn)/double(fpr(xn));
   M(i,1)=xn; M(i,2)=f(xn);
end
disp(['saknes tuvinâjums:  ',num2str(x_app(:))])
disp('              x                 f(x)')
x1 = M(iter,1)
%% turpinâjums - sakne 2
x_app =[2.04]; % saknes tuvinâjums 1
iter=6;     % max iterâciju skaits
syms x, fpr(x)=diff(f(x),x);
xn=x_app; M =zeros(iter,2);
for i=1:iter
   xn=xn-f(xn)/double(fpr(xn));
   M(i,1)=xn; M(i,2)=f(xn);
end
disp(['saknes tuvinâjums:  ',num2str(x_app(:))])
disp('              x                 f(x)')
x2 = M(iter,1)

%% 
disp(['Atbilde: x1 = ' num2str(x1, 6) ', x2 = ' num2str(x2, 6)])
%% 7.UZDEVUMS
clc, clear all,close all, format compact, format longG
f = @(x) 4 .* cos(x .^2) - sqrt(1 + x .^ 3);

a_pr=2.3:0.1:2.6;
plot(a_pr,f(a_pr),'r','LineWidth',3)
title(['Graph of the function f(a) on the interval [1 10]'])
xlabel('a'),ylabel('f(a)'),grid on

%%
iter = 6;
xn = zeros(1,iter+2);
xn = [ 1 1.2 ];
for k=1:iter
    k_1=k+1; k_2=k+2;
    xn(k_2)=xn(k_1)-f(xn(k_1))*(xn(k_1)-xn(k))/(f(xn(k_1))-f(xn(k)));
end
M_pr(:,1)=xn(3:iter+2); M_pr(:,2)=f(xn(3:iter+2));
disp(['saknes tuvinâjumu intervâls: ', num2str(xn(1:2))])
disp('              x                 f(x)')
disp(M_pr)
x1 = M_pr(iter, 1)
%% 8.UZDEVUMS
clc, clear all,close all, format compact, format longG
integ = @(a,x) log(1 + a .* x + x .^ 3);
for a=1:10
    fun=@(x) integ(a,x);
    int_val(a)=integral(fun,0,1) - 1;   
end
a_pr=1:10;
plot(a_pr,int_val,'r','LineWidth',3)
title(['Graph of the function f(a) on the interval [1 10]'])
xlabel('a'),ylabel('f(a)'),grid on

%% turpinâjums
iter=6;
a_xn =zeros(1,iter+2);
a_xn=[ 3.2 3.9 ]; 
for i=1:2
    fun=@(x) integ(a_xn(i),x);
    f(i)=integral(fun,0,1) - 1; 
end
for k=1:iter
    k1=k+1; i=k+2;
    a_xn(i)=a_xn(k1)-f(k1)*(a_xn(k1)-a_xn(k))/(f(k1)-f(k));

    fun=@(x) integ(a_xn(i),x);
    f(i)=integral(fun,0,1) - 1; 

    M_pr(k,1)=a_xn(i);
    M_pr(k,2)=f(i);
end
disp(['skaitïa a  tuvinâjumu intervâls: ', num2str(a_xn(1:2))])
disp('            a                  f(x)')
disp(M_pr)

% turpinâjums
disp('Atbilde: ')
disp([' a = ',num2str(M_pr(iter-1,1),8)])
format
%% 9.UZDEVUMS
clc, clear all,close all, format compact, format longG
f = @(x) sin(3 .* x) - x .^ 2 ./ 5 + x - 0.2;

a_pr=-5:0.1:10;
plot(a_pr,f(a_pr),'r','LineWidth',3)
xlabel('a'),ylabel('f(a)'),grid on

%%
a_pr=1:0.01:2;
plot(a_pr,f(a_pr),'r','LineWidth',3)
xlabel('a'),ylabel('f(a)'),grid on
%%
a_pr=3:0.01:4;
plot(a_pr,f(a_pr),'r','LineWidth',3)
xlabel('a'),ylabel('f(a)'),grid on
%% x0
a_pr=-0.5:0.01:0.5;
plot(a_pr,f(a_pr),'r','LineWidth',3)
xlabel('a'),ylabel('f(a)'),grid on
%% x5
a_pr=4.5:0.01:6;
plot(a_pr,f(a_pr),'r','LineWidth',3)
xlabel('a'),ylabel('f(a)'),grid on
%% sakne 0
iter = 6;
xn = zeros(1,iter+2);
xn = [ 0 0.1 ];
for k=1:iter
    k_1=k+1; k_2=k+2;
    xn(k_2)=xn(k_1)-f(xn(k_1))*(xn(k_1)-xn(k))/(f(xn(k_1))-f(xn(k)));
end
M_pr(:,1)=xn(3:iter+2); M_pr(:,2)=f(xn(3:iter+2));
disp(['saknes tuvinâjumu intervâls: ', num2str(xn(1:2))])
disp('              x                 f(x)')
disp(M_pr)
x0 = M_pr(iter, 1)
%% sakne 1
iter = 6;
xn = zeros(1,iter+2);
xn = [ 1.3 1.4 ];
for k=1:iter
    k_1=k+1; k_2=k+2;
    xn(k_2)=xn(k_1)-f(xn(k_1))*(xn(k_1)-xn(k))/(f(xn(k_1))-f(xn(k)));
end
M_pr(:,1)=xn(3:iter+2); M_pr(:,2)=f(xn(3:iter+2));
disp(['saknes tuvinâjumu intervâls: ', num2str(xn(1:2))])
disp('              x                 f(x)')
disp(M_pr)
x1 = M_pr(iter, 1)
%% sakne 2
iter = 6;
xn = zeros(1,iter+2);
xn = [ 1.65 1.75 ];
for k=1:iter
    k_1=k+1; k_2=k+2;
    xn(k_2)=xn(k_1)-f(xn(k_1))*(xn(k_1)-xn(k))/(f(xn(k_1))-f(xn(k)));
end
M_pr(:,1)=xn(3:iter+2); M_pr(:,2)=f(xn(3:iter+2));
disp(['saknes tuvinâjumu intervâls: ', num2str(xn(1:2))])
disp('              x                 f(x)')
disp(M_pr)
x2 = M_pr(iter, 1)
%% sakne 3
iter = 6;
xn = zeros(1,iter+2);
xn = [ 3.4 3.6 ];
for k=1:iter
    k_1=k+1; k_2=k+2;
    xn(k_2)=xn(k_1)-f(xn(k_1))*(xn(k_1)-xn(k))/(f(xn(k_1))-f(xn(k)));
end
M_pr(:,1)=xn(3:iter+2); M_pr(:,2)=f(xn(3:iter+2));
disp(['saknes tuvinâjumu intervâls: ', num2str(xn(1:2))])
disp('              x                 f(x)')
disp(M_pr)
x3 = M_pr(iter, 1)

%% sakne 4
iter = 6;
xn = zeros(1,iter+2);
xn = [ 3.9 4 ];
for k=1:iter
    k_1=k+1; k_2=k+2;
    xn(k_2)=xn(k_1)-f(xn(k_1))*(xn(k_1)-xn(k))/(f(xn(k_1))-f(xn(k)));
end
M_pr(:,1)=xn(3:iter+2); M_pr(:,2)=f(xn(3:iter+2));
disp(['saknes tuvinâjumu intervâls: ', num2str(xn(1:2))])
disp('              x                 f(x)')
disp(M_pr)
x4 = M_pr(iter, 1)
%% sakne 5
iter = 6;
xn = zeros(1,iter+2);
xn = [ 5 5.5 ];
for k=1:iter
    k_1=k+1; k_2=k+2;
    xn(k_2)=xn(k_1)-f(xn(k_1))*(xn(k_1)-xn(k))/(f(xn(k_1))-f(xn(k)));
end
M_pr(:,1)=xn(3:iter+2); M_pr(:,2)=f(xn(3:iter+2));
disp(['saknes tuvinâjumu intervâls: ', num2str(xn(1:2))])
disp('              x                 f(x)')
disp(M_pr)
x5 = M_pr(iter, 1)
%% 
disp(['Atbilde: '])
disp(['x0 = ' num2str(x0, 6)])
disp(['x1 = ' num2str(x1, 6)])
disp(['x2 = ' num2str(x2, 6)])
disp(['x3 = ' num2str(x3, 6)])
disp(['x4 = ' num2str(x4, 6)])
disp(['x5 = ' num2str(x5, 6)])