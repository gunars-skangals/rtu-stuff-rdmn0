%% 1.piem�rs. Faktoriz�cijas metode.
clear all, clc, format compact
a=[0 2 1 4]; b=[-5 -6 -9 -6].*(-1); c=[1 3 4 0];
d=[1;7;-2;9]; n=length(a);
ksi=zeros(1,n+1); eta=zeros(1,n+1);
for i=1:n
    ksi(i+1)=c(i)/(b(i)-a(i)*ksi(i));
    eta(i+1)=(a(i)*eta(i)-d(i))/(b(i)-a(i)*ksi(i));
end
X(n)=eta(n+1); % apr��in�t x
for i=n-1:-1:1
    X(i)=ksi(i+1)*X(i+1)+eta(i+1);
end
X
disp('Atbilde:'),disp(' Atrisin�jums='), disp(X)
% p�rbaude: Faktoriz�cijas metode ir stabila 
for i=1:n
    if abs(b(i)) < ( abs(a(i))+abs(c(i)) )
        disp(' Faktoriz�cijas metode nav stabila'), return
    end
end
disp(' Faktoriz�cijas metode ir stabila')        

%% 2.piem�rs.Vektoru un matricu normas
clear all, clc, format compact
u=[1;4;6;8];
A=[3,7,-1,2;8,4,3,5;-6,9,0,4;5,8,11,3];
norm_v1=norm(u,1), norm_v2=norm(u,2)
norm_vInf=norm(u,inf),  norm_m2=norm(A,2)% Ctrl+Enter
disp('Atbilde:')
disp([' vektora norma(p=1) = ' , num2str(norm_v1) ])
disp([' vektora norma(p=2) = ' , num2str(norm_v2) ])
disp([' vektora norma(inf) = ' , num2str(norm_vInf) ])
disp([' matricas norma(p=2) = ' , num2str(norm_m2) ])

%% 3.piem�rs.Jakobi metode.
clear all, clc, format compact
A=[-4,1,2;3,-7,3;1,2,5];  B=[2;29;17]; n=length(B);
x_app=zeros(n,1); % s�kuma tuvin�jums 
itermax=20;       % max iter�ciju skaits
epsi=0.001;       % apr��inu precizit�te
iter=0;           % iter�ciju skaits
solnorm=1;        % k��das norma 
prnorm=zeros(1,2);% iter�ciju skaits un k��das norma
k=1;
while solnorm > epsi && iter < itermax
    k=k+1;  iter=iter+1;
    for i=1:3
        res_sum=0;
        for j=1:3
            if j~=i
               res_sum=res_sum+x_app(j,k-1)*A(i,j);
            end
        end
        x_app(i,k)=(B(i,1)-res_sum)/A(i,i);
    end
    solnorm=norm((x_app(:,k)-x_app(:,k-1)),2);
    prnorm(iter,:)=[iter,solnorm];
end
disp('Iter. numurs  K��da'),disp(prnorm)
x_approx=x_app(:,k)
x_sol=linsolve(A,B)
[row,col]=size(prnorm);
disp('Atbilde:')
disp([' iter.skaits = ' , num2str(prnorm(row,1)) ])
disp([' k��da = ',num2str(prnorm(row,2)) ])
disp([' x_tuvin�jumi : {'...
        num2str(x_approx(:)') '}'])
fun_prob3(A);  % p�rbaude: Jakobi metode konver��

%% 4.piem�rs. Jakobi metode.
clear all, clc,format, format compact, n=20;
A=zeros(n,n);B=ones(n,1);
for i=1:n-1
    A(i+1,i)=1;A(i,i)=-3;A(i,i+1)=1;
end
A(n,n)=-3;
x_app=zeros(n,1); k=1;itermax=8;
solnorm=1; prnorm=zeros(1,2);
for iter = 1:itermax
    k=k+1;
    for i=1:n
        rez_sum=0;
        for j=1:n
            if j~=i
               rez_sum=rez_sum+x_app(j,k-1)*A(i,j);
            end
        end
        x_app(i,k)=(B(i,1)-rez_sum)/A(i,i);
    end
    solnorm=norm((x_app(:,k)-x_app(:,k-1)),2);
    prnorm(iter,:)=[iter,solnorm];
end
disp('Iter. numurs     K��da')
disp(prnorm)
x_approx=x_app(:,k)
epsi3_norm=norm(x_app(:,4)-x_app(:,3)) 
% atrisin�juma norma 8.iter�cij�
X_app8_norm=norm(x_app(:,9))
disp('Atbilde:')
disp([' iter.skaits = ' num2str(prnorm(8,1)),...
      ' -> k��da = ',num2str(prnorm(8,2))  ])
disp([' iter.skaits = ' num2str(prnorm(3,1)),...
      ' -> k��da = ',num2str(epsi3_norm)  ])
disp([' X_app8_norm = ',num2str(X_app8_norm) ])
disp([' x_tuvin�jumi : {'...
        num2str(x_approx(:)') '}'])
fun_prob3(A);  % p�rbaude: Jakobi metode konver��, matrica ir pozit�vi defin�ta 
       
%% Uzdevumi patst�v�gai risin��anai
%% 1.uzdevums.
%% 2.uzdevums
%% 3.uzdevums
%% 4.uzdevums
%% 5.uzdevums
%% 6.uzdevums
%% 7.uzdevums
%% 8.uzdevums
%% 9.uzdevums


%% p�rbaude: Jakobi metode konver��
function fun_prob3(A_mat)
 [row,col]=size(A_mat);
 for i=1:row
    sum=0;
    for j=1:col
        if i~=j
            sum=sum+abs(A_mat(i,j));
        end
    end
    if abs(A_mat(i,i)) <= sum
      disp(' Jakobi metode diver��'), 
      disp([' rindas numurs ' num2str(i) ':' '  --> ' num2str(A_mat(i,i)) ...
          ' < ' num2str(sum) ])
      return
    end
 end
 disp(' Jakobi metode konver��')
end