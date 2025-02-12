% vienāojumu risināšana
% x^2 -6x+8=0
syms x
atr = solve(x^2 - 6*x + 8)

%%
clc
clearvars
syms x
atr = solve(x^2 - 6*x + 8 == 0);
x1 = atr(1), x2 = atr(2)
%% 
clc
clearvars
syms x
atr = solve(x^3 == 1, 'Real', true)
%%
clc
clearvars
syms x y
[x , y] = solve(x^2 + y^2 == 5, x*y == 2);
atr = [x, y]
%%
clc
clearvars
syms a b c x
solve (a * x^2 + b * x + c, a)
%%
% for i=vērtību vektors
%   komandas
A = [3 6 1 -8 32 104 93];
sum = 0;
n = length(A);
for i = 1:n
    sum = sum + A(i);
end
sum

%%
% for i=vērtību vektors
%   komandas
A = [3 6 1 -8 32 104 93];
sum = 0;
for j = A
    sum = sum + j;
end
sum
%%
clc
clearvars
A = [ 1 2 3 4; 9 8 7 6; 3 5 7 9 ];
[r, k] = size(A);
result = 1;
for i=1:r
   for j=1:k
       result = result * A(i, j);
   end
end
result
%% nosacījuma operators
% if nosacījumi
%     komandas
% end
% if nosacījumi
%    komandas1
% elseif
%    komandas2    
% else
%    komandas3
% end
% a == b - vienādība
% a ~= b - nav vienāds
% < > normāli
% a <= b, a >= b
% nosacījumi1 && nosacījumi2 - UN
% nosacījumi || nosacījumi2  - VAI
% pāra skaitlis - rem (n , 2)  == 0
% nepāra skaitlis - rem (n , 2)  == 1
%% 
clc
clearvars
A = [ 1 2 3 4; 9 8 7 6; 3 5 7 9; 2 4 6 8 ];
[r, k] = size(A);
% - nepāra rindās pieskaitīt 1
% pāra rindās sareizina't ar 2
for i=1:r
   for j=1:k
       if rem(i, 2) == 0
           A(i , j) = A(i , j) * 2;
       else
           A(i , j) = A(i , j) + 1;
       end
   end
end
A
%% 
clc
clearvars
% while nosacījums
% komandas
% end
n = 1;
summa = 0;
while (1 / n ^ 2 > 10 ^(-3))
    summa = summa + 1 / n ^ 2
    n = n + 1;
end
%% 1. uzdevums a
clc
clearvars
syms x
solve(x^2 + 3 * x == 4*x + 2)
%% 1. uzdevums b
clc
clearvars
syms x y z
[xa, ya, za] = solve(2 * (x + y) == x * y * z, 6 * (y + z) == 5 * x * y * z, 3 * (x + z) == 2 *x *y *z);
[xa, ya, za]
%% 5.uzdevums
clc
clearvars
A = [2, 3, 5, -1, 4, 7, 10, -15, 0, 25, 62]
sum = 0
reiz = 1
n = length(A);
for i=1:n
    if rem(A(i), 2) == 0
        sum = sum + A(i);
    else
        reiz = reiz * A(i);
    end
end
sum
reiz
%% 6. uzdevuma
clc
clearvars
summa = 0;
for i = 1:40
    if (i <= 10) || (i >= 30)
       summa = summa + i^2 + 2 * i ;
    else
        summa = summa + i *3 + 1;
    end
end
summa
%% 7.uzdevums
clc
clearvars
syms x
s3 = 0;
s6 = 0;
s9 = 0;
for n=1:9
    elem = 4 * (-1)^n / (n ^2 * sym(pi) ^ 2) * cos(n * pi * x);
    if n <= 3
        s3 = s3 + elem;
    end
    if n <= 6
        s6 = s6 + elem;
    end
    
    s9 = s9 + elem;
end
s3
s6
s9
fplot(s3)
hold on
fplot(s6)
fplot(s9)
%axis ([0, 5, -0.5, 1])
hold off
%% 8.uzdevums
clc
clearvars
k = 1;
% ak = 4 / (k ^ 3 + 3);
summa = 0;
% k = k + 1;
while (4 / (k ^ 3 + 3) > 0.0001)
    summa = summa + 4 / (k ^ 3 + 3);
    k = k + 1;
end 
k - 1
summa
%% 8.uzdevums vēlreiz
clc
clearvars
k = 0;
summa = 0;
while true
    k = k + 1;
    ak = 4 / (k ^ 3 + 3);
    
    if ak < 0.0001
        break
    end
    
    summa = summa + ak;
end
k - 1
summa
