%% 1.piemçrs. Vienkârðâ iterâcijas metode
clear all,format ,clc, format compact
A=[5,3,7;3,12,-4;7,-4,52];b=[32;15;155];
ni=fun_prob1(A); % pârbaude vai matrica ir pozitivâ definçta
if ni==2 
    disp('Koeficientu matrica nav pozitîvi definçta'), return
end
check=isequal(A,A');
if check==0 
    disp('Koeficientu matrica nav simetriskâ'), return
end
disp('Koeficientu matrica ir simetriska un pozitîvi definçta ')
tau=0.01;  x_app=zeros(3,1);
epsi=10^(-3);itermax=1000;
k_iter=0; resid=b-A*x_app;
while norm(resid) > epsi & k_iter < itermax
    x_app=x_app+tau*resid;  resid=b-A*x_app;
    k_iter=k_iter+1;
end
k_iter, x_app, x_sol=linsolve(A,b)
x_app(1),x_app(2),x_app(3)

disp('Atbilde:')
disp([' iter. skaits = ' , num2str(k_iter) ])
disp([' x_tuvinâjumi : {', num2str(x_app(:)') '}'])
disp(' vai ' )
disp([' x_tuvinâjumi :  x1=', num2str(x_app(1)) ,', x2=' num2str(x_app(2)),...
      ', x3=' num2str(x_app(3))])
x_str=num2str(x_app)
disp([' x_tuvinâjumi :  ' x_str(1,1:6) ' ' x_str(2,1:6) '  '  x_str(3,1:6) ])

fprintf('\nThe result is: [') 
fprintf(' %d ', x_app) 
fprintf(']\n') 

%% Piemçrs. Vienkârðâ iterâcijas metode
clear all,format longG,clc, format compact
A=[5,3,7;3,12,-4;7,-4,52];b=[32;15;155];
ni=fun_prob1(A); % pârbaude vai matrica ir pozitivâ definçta
if ni==2 
    disp('Koeficientu matrica nav pozitîvi definçta'), return
end
check=isequal(A,A');
if check==0 
    disp('Koeficientu matrica nav simetriskâ'), return
end
disp('Koeficientu matrica ir simetriska un pozitîvi definçta ')
tau=0.04;  x_app=zeros(3,1);
epsi=10^(-3);itermax=1000;
k_iter=0; resid=b-A*x_app;
while norm(resid) > epsi & k_iter < itermax
    x_app=x_app+tau*resid; resid=b-A*x_app;
    k_iter=k_iter+1;
end
tau,k_iter, x_app, x_sol=linsolve(A,b)

disp('Atbilde:')
disp(' tau=0.01,0.02,0.03 - metode konverìç, bet tau=0.04 - metode diverìç')
format

%% 2.uzdevums.Matricas îpaðvçrtîbas un îpaðvektori
clear all,clc, format compact
A=[1,-2,3,-4;3,7,-8,2;13,5,-9,0;6,4,10,-3];
[V,D]=eig(A)
lambda=eig(A) 
lambda=diag(D)

%% 3.& 4.uzdevums.Aprçíinât tau_max un kïûdas samazinâðanas koeficientu 1.piemçram
% Aprçíinât tau_max 1.piemçram
clear all, clc, format compact
A=[5,3,7;3,12,-4;7,-4,52];
eig_val=eig(A),tau_max=2/max(eig_val)
disp('Atbilde:')
disp([' îpaðvçrtîbas = {' num2str(eig_val(:)') '}'])
disp([' tau_max =' num2str(tau_max)])

% turpinâjums
% Aprçíinât kïûdas samazinâðanas koeficientu
E_mat=eye(3); tau=0.01; 
F=E_mat-tau*A;
ro=norm(F,2)
disp('Atbilde:')
disp([' ro = '  num2str(ro) '-ïoti tuvs vieniniekam, '])
disp(' tâpçc iterâciju skaits ir tik liels, n = 251')

%% 5.piemçrs. Vienkârðâ iterâcijas metode.Pieòemt optimâlo parametra tau vçrtîbu    
clear all, clc, format compact
A=[5,3,7;3,12,-4;7,-4,52];b=[32;15;155];
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
disp([' optimâlo parametra tau vçrt. = ' num2str(tau_opt)])

%% 6.piemçrs. Minimâlâs nesaistes metode
clear all, clc, format compact
A=[4,-2,3;-2,6,-1;3,-1,12];B=[5;3;14];
ni=fun_prob1(A); % pârbaude vai matrica ir pozitivâ definçta
if ni==2 
    disp('Koeficientu matrica nav pozitîvi definçta'), return
end
check=isequal(A,A');
if check==0 
    disp('koeficientu matrica nav simetriskâ'), return
end
disp('A  ir simetriska un pozitîvi definçta ')

k_iter=0; epsi=10^(-3);itermax=300;x_app=zeros(3,1);
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
disp([' x_tuvinâjumi : {' num2str(x_app(:).') '}'])


%% 7.piemçrs. Matricas îpaðvçrtîbu aprçíinâðana. 
clear all, clc, format compact
A=[4,1,3,2;2,5,6,1;7,4,8,3;2,1,5,9];n=length(A(1,:));
x_app=ones(n,1); % sak.tuvinâjums îpaðvektoram 
e_mas(:,1)=x_app/norm(x_app); % normçtais vektors
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

% turpinâjums
eig_val_max=eigs(A,1)  % pçc moduïa lielâkâ îpaðvçrtîba
disp('Atbilde:')
disp([' iterâciju skaits = ' num2str(k_iter) ])
disp([' lielâkâ îpaðvçrtîba = ' num2str(lambda) '( ar prec.=10^(-3 ))'])

%% Uzdevumi patstâvîgai risinâðanai
%% 1. uzdevums
%% 2.uzdevums. 
%% 3.uzdevums. 
%% 4.uzdevums.
%% 5.uzdevums.
%% 6.uzdevums. 
%% 7.uzdevums.
%% 8.uzdevums. 
%% 9.uzdevums.
