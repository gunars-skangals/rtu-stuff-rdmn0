%% Nelineâro vienâdojumu sistçmas

%% Tirgus lîdsvars

%% Vienâdojumu sistçmas
%% 1.piemçrs. 
% Nelineâru vienâdojumu sistçmas. 
clc, clear all, close all, format compact,format longG
f1=@(x1,x2)x1.^3+x2.^3-1;
f2=@(x1,x2)sin(x1+0.1)-x2-0.3;
fimplicit(f1,[0 1.5 0 1.5],'r','LineWidth',3)
hold on
fimplicit(f2,[0 1.5 0 1.5],'b','LineWidth',3)
hold off
grid on,xlabel('x1'),ylabel('x2')
legend('x1^3+x2^3-1','sin(x1+0.1)-x2-0.3')
title('Intersection of two curves')
%%
clear all, close all
zim1 = ezplot('x1^3+x2^3-1',[0 1.5 0 1.5])
hold on
zim2 = ezplot('sin(x1+0.1)-x2-0.3',[0,1.5,0,1.5])
set(zim1,'Color','r','LineWidth',3)
set(zim2,'Color','b','LineWidth',3)
grid
f1=@(x1,x2)x1.^3+x2.^3-1;
f2=@(x1,x2)sin(x1+0.1)-x2-0.3;
%% turpinâjums
syms x1 x2 
xapp=[0.9 0.6]; xapp_pr=xapp; xpr=[x1 x2];
fun=[x1^3+x2^3-1,sin(x1+0.1)-x2-0.3];            % funkcijas
fun_pr=[diff(fun(1),xpr(1)) diff(fun(1),xpr(2))
        diff(fun(2),xpr(1)) diff(fun(2),xpr(2))] % funkciju atvasinâjumi

% turpinâjums
epsi=10^(-5);  k=0; 
sol_norm=1;
while sol_norm > epsi
    for i=1:2
        B(i,1)=-double(subs(fun(i),xpr,xapp));
        for j=1:2      
            A(i,j)=double(subs(fun_pr(i,j),xpr,xapp));
        end
    end
    xdelta=A\B; xapp=xapp+xdelta';
    c=double(subs(fun,xpr,xapp)); sol_norm=norm(c);
    k=k+1;
    M_pr(k,1:2)=xapp(1:2); 
    M_pr(k,3)=sol_norm; 
end

% turpinâjums
disp(['saknes tuvinâjumi: ', num2str(xapp_pr(1:2))])
disp('            x1                       x2                  kïûdas norma')
disp(M_pr)

% turpinâjums
disp('Atbilde: ')
disp([' x1 = ',num2str(M_pr(k,1),5),',  x2 = ',num2str(M_pr(k,2),5)])
format

%% patstavigai risināšanai
%% 1. uzdevums 

% Nelineâru vienâdojumu sistçmas. 
clc, clear all, close all, format compact,format longG
f1 = @(x1,x2) x1 .^ 2 - x2.^ 2 - 1;
f2 = @(x1,x2) x1 .^ 3 + 2 .* x2 .^ 2 - 5;

fimplicit(f1,[0 5 -5 5],'r','LineWidth',3)
hold on
fimplicit(f2,[0 5 -5 5],'b','LineWidth',3)
hold off
grid on,xlabel('x1'),ylabel('x2')
title('Intersection of two curves')

%% tuvinājumi
% x1 = 1.4, x2 = 1
% x1 = 1.4, x2 = -1
%% turpinâjums
syms x1 x2
tuvinajumi = [ 1.4 1; 1.4 -1 ];
for tv = 1:2
    xapp=tuvinajumi(tv, 1:2);
    xapp_pr=xapp;
    xpr=[x1 x2];

    fun=[ x1 ^ 2 - x2 ^ 2 - 1 , x1 ^ 3 + 2 * x2 ^ 2 - 5 ]; % funkcijas
    fun_pr=[diff(fun(1),xpr(1)) diff(fun(1),xpr(2))
            diff(fun(2),xpr(1)) diff(fun(2),xpr(2))] % funkciju atvasinâjumi
    
    % turpinâjums
    epsi=10^(-5);  k=0; 
    sol_norm=1;
    while sol_norm > epsi
        for i=1:2
            B(i,1)=-double(subs(fun(i),xpr,xapp));
            for j=1:2      
                A(i,j)=double(subs(fun_pr(i,j),xpr,xapp));
            end
        end
        xdelta=A\B; xapp=xapp+xdelta';
        c=double(subs(fun,xpr,xapp)); sol_norm=norm(c);
        k=k+1;
        M_pr(k,1:2)=xapp(1:2); 
        M_pr(k,3)=sol_norm; 
    end
    
    % turpinâjums
    disp(['saknes tuvinâjumi: ', num2str(xapp_pr(1:2))])
    disp('            x1                       x2                  kïûdas norma')
    disp(M_pr)
    
    % turpinâjums
    disp(['Atbilde (' num2str(tv) '):'])
    disp([' x1 = ',num2str(M_pr(k,1),5),',  x2 = ',num2str(M_pr(k,2),5)])
    format

end

%% 2. uzdevums 

% Nelineâru vienâdojumu sistçmas. 
clc, clear all, close all, format compact,format longG
f1 = @(x1,x2) x2 .^ 2 - x1.^ 3 ./ (4 - x1);
f2 = @(x1,x2) x1 .^ 2 - x2 .^ 2 - 0.5;

fimplicit(f1,[0 4 0 4],'r','LineWidth',3)
hold on
fimplicit(f2,[0 4 0 4],'b','LineWidth',3)
hold off
grid on,xlabel('x1'),ylabel('x2')
title('Intersection of two curves')

%% tuvinājumi
% x1 = 0.81, x2 = 0.41
% x1 = 1.84, x2 = 1.7
%% turpinâjums
syms x1 x2
tuvinajumi = [ 0.81 0.41; 1.84 1.7 ];
for tv = 1:2
    xapp=tuvinajumi(tv, 1:2);
    xapp_pr=xapp;
    xpr=[x1 x2];

    fun=[ x2 ^ 2 - x1 ^ 3 / (4 - x1), x1 ^ 2 - x2 ^ 2 - 0.5]; % funkcijas
    fun_pr=[diff(fun(1),xpr(1)) diff(fun(1),xpr(2))
            diff(fun(2),xpr(1)) diff(fun(2),xpr(2))]; % funkciju atvasinâjumi
    
    % turpinâjums
    epsi=10^(-5);  k=0; 
    sol_norm=1;
    while sol_norm > epsi
        for i=1:2
            B(i,1)=-double(subs(fun(i),xpr,xapp));
            for j=1:2      
                A(i,j)=double(subs(fun_pr(i,j),xpr,xapp));
            end
        end
        xdelta=A\B; xapp=xapp+xdelta';
        c=double(subs(fun,xpr,xapp)); sol_norm=norm(c);
        k=k+1;
        M_pr(k,1:2)=xapp(1:2); 
        M_pr(k,3)=sol_norm; 
    end
    
    % turpinâjums
    disp(['saknes tuvinâjumi: ', num2str(xapp_pr(1:2))])
    disp('            x1                       x2                  kïûdas norma')
    disp(M_pr)
    
    % turpinâjums
    disp(['Atbilde (' num2str(tv) '):'])
    disp([' x1 = ',num2str(M_pr(k,1),5),',  x2 = ',num2str(M_pr(k,2),5)])
    format

end

%% 3.uzdevums
% Nelineâru vienâdojumu sistçmas. 
clc, clear all, close all, format compact,format longG
f1 = @(x1,x2) x1 .^ 3 + x2 .^ 3 - 3 .* x1 .* x2;
f2 = @(x1,x2) x1 .^ 2 + 3 .* x2 .^ 2 - 16;

fimplicit(f1,[-5 5 -5 5],'r','LineWidth',3)
hold on
fimplicit(f2,[-5 5 -5 5],'b','LineWidth',3)
hold off
grid on,xlabel('x1'),ylabel('x2')
title('Intersection of two curves')

%% turpinâjums
syms x1 x2
tuvinajumi = [ 
    -2.7 1.7; 
    1.3  -2.2
    ];
for tv = 1:2
    xapp=tuvinajumi(tv, 1:2);
    xapp_pr=xapp;
    xpr=[x1 x2];

    fun=[ x1 ^ 3 + x2 ^ 3 - 3 * x1 * x2, x1 ^ 2 + 3 * x2 ^ 2 - 16 ]; % funkcijas
    fun_pr=[diff(fun(1),xpr(1)) diff(fun(1),xpr(2))
            diff(fun(2),xpr(1)) diff(fun(2),xpr(2))]; % funkciju atvasinâjumi
    
    % turpinâjums
    epsi=10^(-5);  k=0; 
    sol_norm=1;
    while sol_norm > epsi
        for i=1:2
            B(i,1)=-double(subs(fun(i),xpr,xapp));
            for j=1:2      
                A(i,j)=double(subs(fun_pr(i,j),xpr,xapp));
            end
        end
        xdelta=A\B; xapp=xapp+xdelta';
        c=double(subs(fun,xpr,xapp)); sol_norm=norm(c);
        k=k+1;
        M_pr(k,1:2)=xapp(1:2); 
        M_pr(k,3)=sol_norm; 
    end
    
    % turpinâjums
    disp(['saknes tuvinâjumi: ', num2str(xapp_pr(1:2))])
    disp('            x1                       x2                  kïûdas norma')
    disp(M_pr)
    
    % turpinâjums
    disp(['Atbilde (' num2str(tv) '):'])
    disp([' x1 = ',num2str(M_pr(k,1),5),',  x2 = ',num2str(M_pr(k,2),5)])
    format

end