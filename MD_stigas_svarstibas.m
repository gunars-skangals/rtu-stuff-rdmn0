%%
syms s n1 
c = 13;
L = 6;
f(s) = 1 / 55 * s * (6 - s); % f(x)
g(s) = 1 / 99 * sin(7 *pi * s / 6); % g(x)
% u(x, 0) = 
% ut(x, 0) = 
% u (0, t) = 0
% u (6, t) = 0
% koeficientu aprēķināšana
bn(n1) = 2 / L * int(f(s) * sin(n1 * pi * s / L), s , 0, L);
cn(n1) = 2 / (c * n1 * pi) * int (g(s) * sin(n1 * pi * s / L), s, 0, L);
eps = 10 ^ (-3); % precizitāte
y1 = linspace(0, L);
z1 = 1 ./ 55 .* y1 .* (6 - y1); % f(x) - sākuma stāvoklis
hold off
plot(y1, z1)
hold on
t=0.5
% aprekjinam 5 t vertibam
for k=1:5
    A = zeros(101); x = - L / 100; % intervalu (0, L) sadalām 100 daļās
    for m=1:101
        x = x + L / 100;
        S=0; n=1; u=1;
        % koeficientu reekjinaashana
        while (u > eps)
            Bn = double(subs(bn(n1), n1, n));
            Cn = double(subs(cn(n1), n1, n));
            u = (Bn * cos(c * n * pi * t / L) + Cn * sin(c * n * pi * t / L)) * sin(n * pi * x / L);
            S = S + u;
            n = n + 1;
        end;
        
        A(m) = S;
    end;
    p = 1:101;
    y = p * L / 100;
    plot(y, A)
    t = t + 0.5;    
end;
hold off
