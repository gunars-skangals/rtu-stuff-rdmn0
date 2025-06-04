%% 1.piem�rs. Vienk�r�� iter�cijas metode
clear all,format ,clc, format compact
A=[9 2 5;2 18 6;5 6 27];b=[1;6;2];
ni=fun_prob1(A); % p�rbaude vai matrica ir pozitiv� defin�ta
if ni==2 
    disp('Koeficientu matrica nav pozit�vi defin�ta'), return
end
check=isequal(A,A');
if check==0 
    disp('Koeficientu matrica nav simetrisk�'), return
end
disp('Koeficientu matrica ir simetriska un pozit�vi defin�ta ')
tau=0.01;  n=length(b);
x_app=zeros(n,1);
epsi=10^(-3);itermax=12;
k_iter=0; resid=b-A*x_app;
for k_iter=1:itermax
%while norm(resid) > epsi & k_iter < itermax
    x_app=x_app+tau*resid;  resid=b-A*x_app;
end
k_iter, x_app, x_sol=linsolve(A,b)
x_app(1),x_app(2),x_app(3)
norma_pec_12_iter=norm(x_app)
nesaistes_norma=norm(b-A*x_app)
disp('Atbilde:')
disp([' iter. skaits = ' , num2str(k_iter) ])
disp([' x_tuvin�jumi : {', num2str(x_app(:)') '}'])
disp(' vai ' )
disp([' x_tuvin�jumi :  x1=', num2str(x_app(1)) ,', x2=' num2str(x_app(2)),...
      ', x3=' num2str(x_app(3))])
x_str=num2str(x_app)
disp([' x_tuvin�jumi :  ' x_str(1,1:6) ' ' x_str(2,1:6) '  '  x_str(3,1:6) ])

fprintf('\nThe result is: [') 
fprintf(' %d ', x_app) 
fprintf(']\n') 

%% Piem�rs. Vienk�r�� iter�cijas metode
clear all,format longG,clc, format compact
A=[5,3,7;3,12,-4;7,-4,52];b=[32;15;155];
ni=fun_prob1(A); % p�rbaude vai matrica ir pozitiv� defin�ta
if ni==2 
    disp('Koeficientu matrica nav pozit�vi defin�ta'), return
end
check=isequal(A,A');
if check==0 
    disp('Koeficientu matrica nav simetrisk�'), return
end
disp('Koeficientu matrica ir simetriska un pozit�vi defin�ta ')
tau=0.04;  x_app=zeros(3,1);
epsi=10^(-3);itermax=1000;
k_iter=0; resid=b-A*x_app;
while norm(resid) > epsi & k_iter < itermax
    x_app=x_app+tau*resid; resid=b-A*x_app;
    k_iter=k_iter+1;
end
tau,k_iter, x_app, x_sol=linsolve(A,b)

disp('Atbilde:')
disp(' tau=0.01,0.02,0.03 - metode konver��, bet tau=0.04 - metode diver��')
format

%% 2.uzdevums.Matricas �pa�v�rt�bas un �pa�vektori
clear all,clc, format compact
A=[1,-2,3,-4;3,7,-8,2;13,5,-9,0;6,4,10,-3];
[V,D]=eig(A)
lambda=eig(A) 
lambda=diag(D)

%% 3.& 4.uzdevums.Apr��in�t tau_max un k��das samazin��anas koeficientu 1.piem�ram
% Apr��in�t tau_max 1.piem�ram
clear all, clc, format compact
A=[5,3,7;3,12,-4;7,-4,52];
eig_val=eig(A),tau_max=2/max(eig_val)
disp('Atbilde:')
disp([' �pa�v�rt�bas = {' num2str(eig_val(:)') '}'])
disp([' tau_max =' num2str(tau_max)])

%% turpin�jums
% Apr��in�t k��das samazin��anas koeficientu
E_mat=eye(3); tau=0.01; 
F=E_mat-tau*A;
ro=norm(F,2)
disp('Atbilde:')
disp([' ro = '  num2str(ro) '-�oti tuvs vieniniekam, '])
disp(' t�p�c iter�ciju skaits ir tik liels, n = 251')

%% 5.piem�rs. Vienk�r�� iter�cijas metode.Pie�emt optim�lo parametra tau v�rt�bu    
clear all, clc, format compact
A=[5,3,7;3,12,-4;7,-4,52];b=[32;15;155];
ni=fun_prob1(A); % p�rbaude vai matrica ir pozitiv� defin�ta
if ni==2 
    disp('Koeficientu matrica nav pozit�vi defin�ta'), return
end
check=isequal(A,A');
if check==0 
    disp('Koeficientu matrica nav simetrisk�'), return
end
disp('Koeficientu matrica ir simetriska un pozit�vi defin�ta ')
x_app=zeros(3,1); epsi=10^(-3); itermax=1000;
lambda=eig(A); tau_opt=2/(max(lambda)+min(lambda))
k_iter=0; resid=b-A*x_app;
while norm(resid) > epsi && k_iter < itermax
    x_app=x_app+tau_opt*resid; resid=b-A*x_app;
    k_iter=k_iter+1;
end
k_iter, x_app
disp('Atbilde:')
disp([' iter. skaits = ' num2str(k_iter) ])
disp([' x_tuvinajumi : {' num2str(x_app(:).') '}'])
disp([' optim�lo parametra tau v�rt. = ' num2str(tau_opt)])

%% 6.piem�rs. Minim�l�s nesaistes metode
clear all, clc, format compact
A=[4,-2,3;-2,6,-1;3,-1,12];B=[5;3;14];
ni=fun_prob1(A); % p�rbaude vai matrica ir pozitiv� defin�ta
if ni==2 
    disp('Koeficientu matrica nav pozit�vi defin�ta'), return
end
check=isequal(A,A');
if check==0 
    disp('koeficientu matrica nav simetrisk�'), return
end
disp('A  ir simetriska un pozit�vi defin�ta ')

k_iter=0; epsi=10^(-3);itermax=300;n=length(B);
x_app=zeros(n,1);
r=A*x_app-B;norm_r=norm(r);
while norm_r > epsi & k_iter < itermax
    k_iter=k_iter+1
    tau=((A*r)'*r)/norm(A*r)^2
    x_app=x_app-(tau*r')'
    r=A*x_app-B; norm_r=norm(r)
end
x_sol=linsolve(A,B)

disp('Atbilde:')
disp([' iter. skaits = ' num2str(k_iter),...
    ', nesaistes norma = ' num2str(norm_r)])
disp([' x_tuvin�jumi : {' num2str(x_app(:).') '}'])
x_minres=minres(A,B)

%% 7.piem�rs. Matricas �pa�v�rt�bu apr��in��ana. 
clear all, clc, format compact
A=[4,1,3,2;2,5,6,1;7,4,8,3;2,1,5,9];n=length(A(1,:));
x_app=ones(n,1); % sak.tuvin�jums �pa�vektoram 
e_mas(:,1)=x_app/norm(x_app); % norm�tais vektors
x_app(:,2)=A*e_mas;
e_mas(:,2)=x_app(:,2)/norm(x_app(:,2));
k=2; 
epsi=10^(-3);
iter_max=20;
k_iter=0;
while norm(e_mas(:,k)-e_mas(:,k-1)) > epsi & k <= iter_max
    x_app(:,k+1)=A*e_mas(:,k); 
    e_mas(:,k+1)=x_app(:,k+1)/norm(x_app(:,k+1));
    k_iter=k_iter+1;
    lambda=dot( x_app(:,k+1)',e_mas(:,k));
    x_pr=x_app(:,k+1);
    e_pr=e_mas(:,k+1); 
    k=k+1;
end
k_iter,lambda,x_pr,e_pr

% turpin�jums
eig_val_max=eigs(A,1)  % p�c modu�a liel�k� �pa�v�rt�ba
disp('Atbilde:')
disp([' iter�ciju skaits = ' num2str(k_iter) ])
disp([' liel�k� �pa�v�rt�ba = ' num2str(lambda) '( ar prec.=10^(-3 ))'])

%% Uzdevumi patst�v�gai risin��anai
%% 1. uzdevums
clear all,format ,clc, format compact
A=[9 2 5;2 18 6;5 6 27];b=[1;6;2];
ni=fun_prob1(A); % p�rbaude vai matrica ir pozitiv� defin�ta
if ni==2 
    disp('Koeficientu matrica nav pozit�vi defin�ta'), return
end
check=isequal(A,A');
if check==0 
    disp('Koeficientu matrica nav simetrisk�'), return
end
disp('Koeficientu matrica ir simetriska un pozit�vi defin�ta ')
tau=0.01;  n=length(b);
x_app=zeros(n,1);
epsi=10^(-3);itermax=12;
k_iter=0; resid=b-A*x_app;
for k_iter=1:itermax
%while norm(resid) > epsi & k_iter < itermax
    x_app=x_app+tau*resid;  resid=b-A*x_app;
end
k_iter, x_app, x_sol=linsolve(A,b)
x_app(1),x_app(2),x_app(3)
norma_pec_12_iter=norm(x_app)
nesaistes_norma=norm(b-A*x_app)
disp('Atbilde:')
disp([' iter. skaits = ' , num2str(k_iter) ])
disp([' x_tuvin�jumi : {', num2str(x_app(:)') '}'])
disp(' vai ' )
disp([' x_tuvin�jumi :  x1=', num2str(x_app(1)) ,', x2=' num2str(x_app(2)),...
      ', x3=' num2str(x_app(3))])
x_str=num2str(x_app)
disp([' x_tuvin�jumi :  ' x_str(1,1:6) ' ' x_str(2,1:6) '  '  x_str(3,1:6) ])

fprintf('\nThe result is: [') 
fprintf(' %d ', x_app) 
fprintf(']\n') 

%% 2.uzdevums. 
clear all,format ,clc, format compact
A=[9 2 5;2 18 6;5 6 27];b=[1;6;2];
ni=fun_prob1(A); % p�rbaude vai matrica ir pozitiv� defin�ta
if ni==2 
    disp('Koeficientu matrica nav pozit�vi defin�ta'), return
end
check=isequal(A,A');
if check==0 
    disp('Koeficientu matrica nav simetrisk�'), return
end
disp('Koeficientu matrica ir simetriska un pozit�vi defin�ta ')
lambda_max=eigs(A,1);
tau=2/lambda_max;  n=length(b);
x_app=zeros(n,1);
epsi=10^(-3);itermax=15;
k_iter=0; resid=b-A*x_app;
for k_iter=1:itermax
%while norm(resid) > epsi & k_iter < itermax
    x_app=x_app+tau*resid;  resid=b-A*x_app;
end
k_iter, x_app, x_sol=linsolve(A,b)
x_app(1),x_app(2),x_app(3)
norma_pec_15_iter=norm(x_app)
nesaistes_norma=norm(b-A*x_app)
disp('Atbilde:')
disp([' iter. skaits = ' , num2str(k_iter) ])
disp([' x_tuvin�jumi : {', num2str(x_app(:)') '}'])
disp(' vai ' )
disp([' x_tuvin�jumi :  x1=', num2str(x_app(1)) ,', x2=' num2str(x_app(2)),...
      ', x3=' num2str(x_app(3))])
x_str=num2str(x_app)
disp([' x_tuvin�jumi :  ' x_str(1,1:6) ' ' x_str(2,1:6) '  '  x_str(3,1:6) ])

fprintf('\nThe result is: [') 
fprintf(' %d ', x_app) 
fprintf(']\n') 

%% 3.uzdevums. 
clear all,format ,clc, format compact
A=[9 2 5;2 18 6;5 6 27];b=[1;6;2];
ni=fun_prob1(A); % p�rbaude vai matrica ir pozitiv� defin�ta
if ni==2 
    disp('Koeficientu matrica nav pozit�vi defin�ta'), return
end
check=isequal(A,A');
if check==0 
    disp('Koeficientu matrica nav simetrisk�'), return
end
disp('Koeficientu matrica ir simetriska un pozit�vi defin�ta ')
lambda=eig(A);
tau=2/(max(lambda)+min(lambda));  n=length(b);
x_app=zeros(n,1);
epsi=10^(-3);itermax=10;
k_iter=0; resid=b-A*x_app;
for k_iter=1:itermax
%while norm(resid) > epsi & k_iter < itermax
    x_app=x_app+tau*resid;  resid=b-A*x_app;
end
k_iter, x_app, x_sol=linsolve(A,b)
x_app(1),x_app(2),x_app(3)
norma_pec_15_iter=norm(x_app)
nesaistes_norma=norm(b-A*x_app)
disp('Atbilde:')
disp([' iter. skaits = ' , num2str(k_iter) ])
disp([' x_tuvin�jumi : {', num2str(x_app(:)') '}'])
disp(' vai ' )
disp([' x_tuvin�jumi :  x1=', num2str(x_app(1)) ,', x2=' num2str(x_app(2)),...
      ', x3=' num2str(x_app(3))])
x_str=num2str(x_app)
disp([' x_tuvin�jumi :  ' x_str(1,1:6) ' ' x_str(2,1:6) '  '  x_str(3,1:6) ])

fprintf('\nThe result is: [') 
fprintf(' %d ', x_app) 
fprintf(']\n') 
%% 4.uzdevums.
clear all, clc, format compact
A=[1 4 7 -1 5;8 2 -3 4 10;...
    6 12 -2 7 4; 3 -4 9 0 11;...
    5 6 -4 2 8];n=length(A(1,:));
x_app=ones(n,1); % sak.tuvin�jums �pa�vektoram 
e_mas(:,1)=x_app/norm(x_app); % norm�tais vektors
x_app(:,2)=A*e_mas;
e_mas(:,2)=x_app(:,2)/norm(x_app(:,2));
k=2; 
epsi=10^(-3);
iter_max=20;
k_iter=0;
while norm(e_mas(:,k)-e_mas(:,k-1)) > epsi & k <= iter_max
    x_app(:,k+1)=A*e_mas(:,k); 
    e_mas(:,k+1)=x_app(:,k+1)/norm(x_app(:,k+1));
    k_iter=k_iter+1;
    lambda=dot( x_app(:,k+1)',e_mas(:,k));
    x_pr=x_app(:,k+1);
    e_pr=e_mas(:,k+1); 
    k=k+1;
end
k_iter,lambda,x_pr,e_pr

% turpin�jums
eig_val_max=eigs(A,1)  % p�c modu�a liel�k� �pa�v�rt�ba
disp('Atbilde:')
disp([' iter�ciju skaits = ' num2str(k_iter) ])
disp([' liel�k� �pa�v�rt�ba = ' num2str(lambda) '( ar prec.=10^(-3 ))'])

%% 5.uzdevums.
clear all, clc, format compact
A=[16 3 4 2;3 12 2 -1;4 2 8 -1;2 -1 -1 2];B=[25;16;13;2];
ni=fun_prob1(A); % p�rbaude vai matrica ir pozitiv� defin�ta
if ni==2 
    disp('Koeficientu matrica nav pozit�vi defin�ta'), return
end
check=isequal(A,A');
if check==0 
    disp('koeficientu matrica nav simetrisk�'), return
end
disp('A  ir simetriska un pozit�vi defin�ta ')

k_iter=0; epsi=10^(-3);itermax=300;n=length(B);
x_app=zeros(n,1);
r=A*x_app-B;norm_r=norm(r);
while norm_r > epsi & k_iter < itermax
    k_iter=k_iter+1
    tau=((A*r)'*r)/norm(A*r)^2
    x_app=x_app-(tau*r')'
    r=A*x_app-B; norm_r=norm(r)
end
x_sol=linsolve(A,B)

disp('Atbilde:')
disp([' iter. skaits = ' num2str(k_iter),...
    ', nesaistes norma = ' num2str(norm_r)])
disp([' x_tuvin�jumi : {' num2str(x_app(:).') '}'])
%% 6.uzdevums. 
clear all, clc, format compact
A=[7 2.5 0.5;2.5 14 1.5;0.5 1.5 21];B=[9;-4;22];
ni=fun_prob1(A); % p�rbaude vai matrica ir pozitiv� defin�ta
if ni==2 
    disp('Koeficientu matrica nav pozit�vi defin�ta'), return
end
check=isequal(A,A');
if check==0 
    disp('koeficientu matrica nav simetrisk�'), return
end
disp('A  ir simetriska un pozit�vi defin�ta ')

k_iter=0; epsi=10^(-4);itermax=300;n=length(B);
x_app=[25;-10;0];
r=A*x_app-B;norm_r=norm(r);
while norm_r > epsi & k_iter < itermax
    k_iter=k_iter+1
    tau=((A*r)'*r)/norm(A*r)^2
    x_app=x_app-(tau*r')'
    r=A*x_app-B; norm_r=norm(r)
end
x_sol=linsolve(A,B)

disp('Atbilde:')
disp([' iter. skaits = ' num2str(k_iter),...
    ', nesaistes norma = ' num2str(norm_r)])
disp([' x_tuvin�jumi : {' num2str(x_app(:).') '}'])
%% 7.uzdevums.
clear all, clc, format compact
A=[6 3 5;3 12 6;5 6 18];B=[-12;8;35];
ni=fun_prob1(A); % p�rbaude vai matrica ir pozitiv� defin�ta
if ni==2 
    disp('Koeficientu matrica nav pozit�vi defin�ta'), return
end
check=isequal(A,A');
if check==0 
    disp('koeficientu matrica nav simetrisk�'), return
end
disp('A  ir simetriska un pozit�vi defin�ta ')

k_iter=0; epsi=10^(-3);itermax=300;n=length(B);
x_app=[-15;0;22];
r=A*x_app-B;norm_r=norm(r);
while norm_r > epsi & k_iter < itermax
    k_iter=k_iter+1
    tau=((A*r)'*r)/norm(A*r)^2
    x_app=x_app-(tau*r')'
    r=A*x_app-B; norm_r=norm(r)
end
x_sol=linsolve(A,B)

disp('Atbilde:')
disp([' iter. skaits = ' num2str(k_iter),...
    ', nesaistes norma = ' num2str(norm_r)])
disp([' x_tuvin�jumi : {' num2str(x_app(:).') '}'])
%% 8.uzdevums. 
clear all,format ,clc, format compact
A=[16 3 4 2;3 12 2 -1;4 2 8 -1;2 -1 -1 2];b=[25;16;13;2];
ni=fun_prob1(A); % p�rbaude vai matrica ir pozitiv� defin�ta
if ni==2 
    disp('Koeficientu matrica nav pozit�vi defin�ta'), return
end
check=isequal(A,A');
if check==0 
    disp('Koeficientu matrica nav simetrisk�'), return
end
disp('Koeficientu matrica ir simetriska un pozit�vi defin�ta ')
tau=0.02;  n=length(b);
x_app=zeros(n,1);
epsi=10^(-3);itermax=15;
k_iter=0; resid=b-A*x_app;
for k_iter=1:itermax
%while norm(resid) > epsi & k_iter < itermax
    x_app=x_app+tau*resid;  resid=b-A*x_app;
end
k_iter, x_app, x_sol=linsolve(A,b)
x_app(1),x_app(2),x_app(3)
norma_pec_15_iter=norm(x_app)
nesaistes_norma=norm(b-A*x_app)
lambda=eig(A);
t_max=2/max(lambda)
t_opt=2/(max(lambda)+min(lambda))
disp('Atbilde:')
disp([' iter. skaits = ' , num2str(k_iter) ])
disp([' x_tuvin�jumi : {', num2str(x_app(:)') '}'])
disp(' vai ' )
disp([' x_tuvin�jumi :  x1=', num2str(x_app(1)) ,', x2=' num2str(x_app(2)),...
      ', x3=' num2str(x_app(3))])
x_str=num2str(x_app)
disp([' x_tuvin�jumi :  ' x_str(1,1:6) ' ' x_str(2,1:6) '  '  x_str(3,1:6) ])

fprintf('\nThe result is: [') 
fprintf(' %d ', x_app) 
fprintf(']\n') 
%% 9.uzdevums.
function ni = fun_prob1(A_mat)
  ni=1; 
  [row,col]=size(A_mat);
   for i=1:row
      if det(A_mat(1:i,1:i))>0
      else ni=2; return
      end
   end
end