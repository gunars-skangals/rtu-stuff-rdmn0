clear all
x=[3 6 9 12 15 18 24 27 30];
% x ir saraþotas produkcijas apjoms
ps=[6.24 7.41 8.5 10.07 12.39 14.61 20.2 24.83 29.94];
% ps is vienas vienîbas cena (piedâvâjums)
pd=[39.48 38.65 34.1 29.43 24.41 20.39 14.21 10.36 4.97];
% pd ir vienas vienîbas cena (pieprâsîjums)
plot(x,ps,'o')
hold on
plot(x,pd,'o')
%%
[p1,reg]=fit(x',ps','poly2')
hold on
x1=[3:0.1:30];
plot(x1,p1(x1))
[p2,reg]=fit(x',pd','poly2')
hold on
x1=[3:0.1:30];
plot(x1,p2(x1))

xapp=[21 17];
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