%% uzd 1
clear all, clc

f1 = @(x, y) 4 .* x .^ 2 - 2 .* x .* y + 3 .* y .^ 2 + 10 .* x + 12 .* y - 8;
f2 = @(x, y) (x .^ 2) ./ 2 + (y .^ 2) ./ 3 - 1;

fimplicit(f1,[-5 5 -5 5],'r','LineWidth',3)
hold on
fimplicit(f2,[-5 5 -5 5],'b','LineWidth',3)
hold off

%% 

%% turpinâjums
syms x1 x2
tuvinajumi = [ 
    -1.25 0.75; 
    1 -1.1
    ];
for tv = 1:2
    xapp=tuvinajumi(tv, 1:2);
    xapp_pr=xapp;
    xpr=[x1 x2];

    fun=[ 4 .* x1 .^ 2 - 2 .* x1 .* x2 + 3 .* x2 .^ 2 + 10 .* x1 + 12 .* x2 - 8, (x1 .^ 2) ./ 2 + (x2 .^ 2) ./ 3 - 1 ]; % funkcijas
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
    % disp(['saknes tuvinâjumi: ', num2str(xapp_pr(1:2))])
    % disp('            x1                       x2                  kïûdas norma')
    % disp(M_pr)
    
    % turpinâjums
    disp(['Atbilde (' num2str(tv) '):'])
    disp([' x1 = ',num2str(M_pr(k,1),5),',  x2 = ',num2str(M_pr(k,2),5)])
    format

end

%%
fun=[ 4 .* x1 .^ 2 - 2 .* x1 .* x2 + 3 .* x2 .^ 2 + 10 .* x1 + 12 .* x2 - 8, (x1 .^ 2) ./ 2 + (x2 .^ 2) ./ 3 - 1 ]; % funkcijas
jakobi=[diff(fun(1),xpr(1)) diff(fun(1),xpr(2))
        diff(fun(2),xpr(1)) diff(fun(2),xpr(2))];

B1 = double(subs(jakobi, [x1 x2], [ -1.2374 0.83866 ]));
eig(B1)

B2 = double(subs(jakobi, [x1 x2], [ 1.0716 -1.1303 ]));
eig(B2)

% for i=1:2
%     B(i,1) = -double(subs(fun(i),xpr,xapp));
%     for j=1:2      
%         A(i,j)=double(subs(fun_pr(i,j),xpr,xapp));
%     end
% end
%% 2.uzdevums
clear all, clc

f1 = @(x, y) 4 .* x .^ 2 + ;
f2 = @(x, y) (x .^ 2) ./ 2 + (y .^ 2) ./ 3 - 1;

fimplicit(f1,[-5 5 -5 5],'r','LineWidth',3)
hold on
fimplicit(f2,[-5 5 -5 5],'b','LineWidth',3)
hold off