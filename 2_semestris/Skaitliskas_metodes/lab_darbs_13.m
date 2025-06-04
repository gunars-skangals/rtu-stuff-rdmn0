%% 10. lab darbs
%% 1. uzdevums (zāļu koncentrācija)
% Vienādojums ir C(n+1)=0.3*C(n)+0.2
% Pieņemsim, ka C(1)=0.2
c=0;
m=100;
for n=1:m
    c=0.3*c+0.2;
end
c
%% 3. uzdevums (vaļu populācijas dinamika)
x=[0 0 1]';
A=[0 0 0.3;0.072 0.8 0.88;0 0.19 0];
k=10000;
for i=1:k
    x=A*x;
end
x
%%
% clc
% A = [ 1 0.5; 1 1.5];
% [P, D]=eig(A)
% Pinv = 
% C = 
%% 1.uzdevums
% Cn = Cn-1 * 0.3 + 0.2
syms t
b = 0.2; a = 0.3;
x0 = 0; % ????
x_nekust = b / (1-a)
C(t) = x_nekust + (a ^ t) * (x0 - x_nekust);
double(C(4))

%% 3.uzdevums
clc
A = [0 0 0.3 ;0.072 0.8 0.88; 0 0.19 0]

[D, V] = eig(A)

%% 4.UZDEVUMS
clc
A = [ 0.998 0.45 0; 0.002 0.1 0; 0 0.45 1 ];

x=[100 0 0]';
k=10000;
for i=1:k
    x=A*x;
end
x


%%
clc
A = [ 0.998 0.65 0; 0.002 0 0.1; 0 0.35 0.9 ];

x=[100 0 0]';
k=1000000;
for i=1:k
    x=A*x;
end
x


% A = [ 1 0.5; 1 1.5];
% [P, D]=eig(A)
% Pinv = 
% C = 
%% 4.uzdevums
clc
A = [ 0 2 4 ; 0.25 0 0 ; 0 0.5 0 ]

x=[100 0 0]';
k=10000;
for i=1:k
    x=A*x;
end
x

%% 5
k = 100;
R = 3.8;
Xt = ones(1, k) / 2;
for i=2:k
    X_prev = Xt(i - 1);
    Xt(i) = R * X_prev * (1 - X_prev);
end
plot(Xt)

