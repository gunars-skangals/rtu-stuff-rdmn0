%% asda
A = [ 1,2,3; 6,4,-2; 3,7,9 ]

B = [ 4 7 1; -3 8 0; 2 5 -1 ]

%% īpašas matricas
O = zeros(3, 3)
eye(3) % vienības matrica
V = ones (2, 5) % vieninieku matrica

%% elementi
a23 = A(2, 3)
rinda2 = A(2, :)
kolonna3 = A(:, 3)
d = diag(A)
[r,k] = size(A)
%% rindas matica
C = [ 3 6 9 2 7 ]
length(C)
%% darbības - determinants
determinants = det(A)
D = 5 * A - 3 * B
P = A * B
inversa = inv(A)
inversa2 = inv(sym(A))
inversa3 = A ^ (-1)
%transponet
A_t = A'
r = rank(A) % lineāri neatkarīgu rindu skaits

%% 
clc
A = [ 1 2 -3; 3 -1 4; 5 7 -1 ];
B = [ 1; -4; 2 ];
x = linsolve(A, B)
AP = [ A B ];
AP_1 = rref(AP)
%% 
clc
clearvars
A = [ 1 2 -3; 3 -1 4; 4 1 1 ];
B = [ 1; -4; -3 ];
% x = linsolve(A, B) % pārāk zems rank
AP = [ A B ];
AP_1 = rref(sym(AP))
%% pārveidot - pāra'rindām peiskaita 5, nepaŗa - piereizina
clc
clearvars
A = [ 1 2 -3 4; 3 -1 4 5; 4 1 1 6; 2 3 4 5 ];
[r, k] = size(A);
for i=1:r
   for j=1:k
     if (rem(i, 2) == 0)
         A(i, j) = A(i, j) + 5;
     else
         A(i, j) = A(i, j) * 5;
     end;        
   end;
end;
A
%% funkcijas vērtību tabula
clc
clearvars
syms x
fn(x) = x^2 + 5 *x;
x_V = 0:0.5:10
tabula(:, 1) = x_V;
tabula(:, 2) = fn(x_V);
disp ("x_vērtības   fn_vērtības")
disp(tabula)
%% uzdevums 1
clc
clearvars
A = [ 1 -2; 3 -4 ]
A^4 + 5 * eye(2)
%% uzdevums 2
A = [ 1 3 4; 2 6 5; 1 7 9 ];
B = [ 2 -1 3; 0 4 -2; 4 2 8];
A * B' + A ^ (-1) * B
%% uzdevums 3
syms x
f(x) = sqrt(cos(x)^2 + x^(3 / 4));
f_0 = double(f(2))
x_V = 1:2:21;
tabula(:, 1) = x_V;
tabula(:, 2) = f(x_V);
disp("x_vērtības   f_vērtības")
disp(tabula)
%% uzdevums 4
clearvars
syms x
f(x) = 2* x ^ 2 + sin(x);
g(x) = cos(x) ^ 2 + 3 * x;
y(x) = f(x) + g(x);
y_0 = double(y(1 / 3))
x_V = 1:2:11;
tabula(:, 1) = x_V;
tabula(:, 2) = y(x_V);
disp("x_vērtības   y_vērtības")
disp(tabula)
%% uzdevums 5
clc
A = [ 1 3 5 7; 2 4 6 8; 13 14 15 16; 9 10 11 12 ];
[r, k] = size(A);
B = zeros(r, k);
for i=1:r
    for j=1:k
        B(i, j) = A(i, j) + j ^ 2 + 2 * j;
    end;
end;
B
B2 = zeros(r, k);
for i=1:r
    for j=1:k
        B2(i, j) = A(i, j) ^ 2 + 3 * i ^ 3;
    end;
end;
B2
%% uzdevums 6
clc
A = [ 3 9; 5 4; 7 2 ];
[r, k] = size(A);
A2 = zeros(r, k * 5);
for i = 1:r
    for j = 1:(k*5)
        if (j < 3)
            A2(i, j) = A(i, j);
        else 
            A2(i, j) = 2 * A2(i, j - 2) - 3 * A2(i, j - 1);
        end
    end
end
disp(A2)
%% uzdevums 7
A = zeros(6, 8);
for i=1:6
    A(i, i) = -1;
    A(i, i+1) = 2;
    A(i, i+2) = 5;
end;
A
%% uzdevums 8
A = [ 5 6 2 3; 11 3 7 2; 8 1 12 9; 0 4 1 5 ];
[r, k] = size(A);
B = ones(r, k);
for i=1:4
    if (rem(i, 2) == 0)
        B(i, i) = A(i, i)^2 + A(i, i - 1);
    else
        B(i, i) = A(i, i)^2 + A(i, i + 1);
    end;
end
B
%% uzdevums 9
A = zeros(10, 10)
for i=1:10
    A(i, i) = 5;    
    if ( i + 1 <= 10)
        A(i, i + 1) = 2;
        A(i + 1, i) = 1;
    end
end
A
%% uzdevums 10
A = [ 3 -3 2 -4 5; 72 -2 11 7 -8; 9 22 12 -2 5; -1 14 -6 10 -2 ];
[r k] = size(A);
B = zeros(1, r);
for i = 1:r
    for j = 1:k
      elem = A(i,j);
      if rem(i, 2) == 0
          if (elem < 0)
              B(i) = B(i) + elem;
          end
      else
          if (elem > 0)
              B(i) = B(i) + elem;
          end
      end
    end
end
disp(B)
rez = 1;
for i = 1:4
    rez = rez * B(i);
end
rez

