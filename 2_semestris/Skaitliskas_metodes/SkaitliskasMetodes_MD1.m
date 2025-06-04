%% 1.uzdevums
clear all, clc
A=[ 1 -2 -1 3; 2 -4 -2 6; 2 1 0 1 ];
B=[ 5; 10; 20 ];
Ap = [ A B ];
rank_A = rank(A);
rank_Ap = rank(Ap);
if rank_A ~= rank_Ap
    disp('Sistēma ir nesaderīga'), return
end;
if rank_Ap < length(A)
    disp('Sistēma ir nenoteikta')
end;

sol = rref(Ap)

disp('Atbilde:')
disp(' x3, x4 - jebkuri reāli skaitļi')
disp([' x1 = ' num2str(sol(1, 5)) ' + ' num2str(-sol(1, 3)) '*x_3'  ' +  ' num2str(-sol(1, 4)) '*x_4' ])
disp([' x2 = ' num2str(sol(2, 5)) ' + ' num2str(-sol(2, 3)) '*x_3'  ' +  ' num2str(-sol(2, 4)) '*x_4' ])

%% 2.uzdevums
clear all, clc
N = 20;
A = diag(ones(N, 1) * -5) + diag(ones(N - 1, 1),-1) + diag(ones(N - 1, 1), 1);
B = ones(N,1) * 3;

for i = 1:N
    if det(A(1:i, 1:i)) == 0
        disp([ i '. galvenais minors ir vienāds ar 0!' ]);
        return;
    end;
end;

[L U P] = lu(A);
Y=L\(P*B);
X=U\Y;
X_labels = ('x_' + string(1:20))';
sol = table(X_labels, X);

disp('Atbilde:')
disp(sol)
%% 3.uzdevums
clear all, clc
A = [ 5 5 0.5; 5 5.5 1.5; 0.5 1.5 50.5 ];
b = [ 18; 32; 7 ];
x_app = [ -25; 29; -0.5 ];
eps = 0.001;

[row,col]=size(A);
for i=1:row
  if det(A(1:i,1:i))>0
  else 
      disp('Matrica nav pozitīvi definēta!')
      return
  end
end

if not(isequal(A,A'))
    disp('Koeficientu matrica nav simetriskâ'), return
end

disp('Koeficientu matrica ir simetriska un pozitīvi definēta ')

eig_val=eig(A); tau = 2 / (max(eig_val) + min(eig_val));

k_iter = 0;
resid = b - A * x_app;
while norm(resid) > eps
    x_app = x_app + tau * resid;
    resid = b - A * x_app;
    k_iter = k_iter + 1;
end

X_labels = ('x_' + string(1:3))';
sol = table(X_labels, x_app);

disp('Atbilde:')
disp(sol)
disp(['Iterāciju skaits: ' num2str(k_iter)])
%% 4.uzdevums 
clc, clear all
A = [ 14 11 14; 11 9 11; 14 11 15 ];
B = [ 3; 23; 4 ];
n = length(A)
x_app = zeros(n, 1);
iter_max = 40;

epsi=0.001;       % aprçíinu precizitâte
iter=0;           % iterâciju skaits
solnorm=1;        % kïûdas norma 

k=1;
while solnorm > epsi && iter < iter_max
    k = k + 1;  iter = iter + 1;
    for i = 1:n
        res_sum = 0;
        for j = 1:n
            if j ~= i
               res_sum = res_sum + x_app(j, k-1) * A(i, j);
            end
        end
        x_app(i, k) = (B(i, 1) - res_sum) / A(i, i);
    end
    solnorm=norm((x_app(:, k) - x_app(:, k-1)), 2);
end

X_labels = ('x_' + string(1:n))';
sol = table(X_labels, x_app(:, k));

disp('Atbilde:')
disp(sol)

disp(' Metode diverģē - neizpildas arī konverģences nosacījums par dominējošu diagonāli')
disp(' Vienkāršo iterāciju metode būtu atbilstošāka')

[row,col]=size(A);
for i=1:row
    sum=0;
    for j=1:col
        if i~=j
            sum=sum+abs(A(i,j));
        end
    end
    if abs(A(i,i)) <= sum
        disp(' Jakobi metode diverģē'), 
        disp([' rindas numurs ' num2str(i) ':' '  --> ' num2str(A(i,i)) ' < ' num2str(sum) ])
    return
    end
end
disp(' Jakobi metode konverģē')

%% 4. uzdevums ar vienkāršo iterāciju metodi
clear all, clc
A = [ 14 11 14; 11 9 11; 14 11 15 ];
b = [ 3; 23; 4 ];
x_app = [ 0;0;0 ];
eps = 0.001;

[row,col]=size(A);
for i=1:row
  if det(A(1:i,1:i))>0
  else 
      disp('Matrica nav pozitīvi definēta!')
      return
  end
end

if not(isequal(A,A'))
    disp('Koeficientu matrica nav simetriskâ'), return
end

disp('Koeficientu matrica ir simetriska un pozitīvi definēta ')

eig_val=eig(A); tau = 2 / (max(eig_val) + min(eig_val));

k_iter = 0;
resid = b - A * x_app;
while norm(resid) > eps
    x_app = x_app + tau * resid;
    resid = b - A * x_app;
    k_iter = k_iter + 1;
end

X_labels = ('x_' + string(1:3))';
sol = table(X_labels, x_app);

disp('Atbilde:')
disp(sol)
disp(['Iterāciju skaits: ' num2str(k_iter)])

%% 5.uzdevums
clc, clear all
x_i = [ 0 0.6 1.2 1.8 2.4 3 3.6 4.2 ];
y_i = [ 1 2.26 2.75 2.33 1.03 -0.91 -3.12 -5.15 ];
m = length(y_i); 

coef=y_i;
for k = 2:m
    coef(k:m) = (coef(k:m) - coef(k-1:m-1)) ./ (x_i(k:m) - x_i(1:m+1-k));
end

syms x
pol = coef(m);
for k = m-1:-1:1
    pol = pol * (x - x_i(k)) + coef(k);
end
polyn(x) = expand(collect(pol));

x_pr = x_i(1):0.01:x_i(m);
spline=interp1(x_i, y_i, x_pr, 'spline');


plot(x_pr, polyn(x_pr), 'r-', x_pr, spline, 'b--', x_i, y_i, 'o', 'LineWidth',2)
legend('polinoms', 'splaini', 'mezgli') 
grid on, ylim([-6 3])