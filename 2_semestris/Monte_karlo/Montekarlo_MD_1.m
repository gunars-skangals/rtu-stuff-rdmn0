%% 1.spēlmaņu uzdevums
clc, clear all
q = 0.495; %  pirmā spēlētaja uzvaras varbūtība
S_A = 5; % pirmā spēlētaja sākuma summa

ns = [ 10 100 1000 10000 ]; % spēļu skaiti
n = ns(end); % simulējam vienreiz, ar maksimālo spēļu skaitu
SK_A = S_A;
results = zeros(n, 1);
for i = 1:n
   if rand > q
        SK_A = SK_A - 1;
   else
        SK_A = SK_A + 1;
   end;
   
   results(i) = SK_A;
end;

% veidojam atsevišķus grafikus
figure
for j = 1:length(ns)
    n = ns(j);
    
    subplot(2, 2, j);
    plot(results(1:n));
    title(sprintf('Spēles gaita pie n = %d', n));
    xlabel('N')
    ylabel('S_A')
    grid on
end;

%% 1.uzdevums, otrā daļa
trials = 1000;
M = ones(trials, 1);

for i = 1:trials
    SK_A = S_A;

    while (SK_A >= 0)
       if rand > q
            SK_A = SK_A - 1;
       else
            SK_A = SK_A + 1;
       end;
       M(i) = M(i) + 1;
    end;
end;

figure
histogram(M)
title('Pirmās negatīvās summas spēles numurs - pilnā histogramma')

figure
histogram(M, [ 50 100 150 200 250 300 350 400 ])
title('Pirmās negatīvās summas spēles numurs')

%% 2. Eiropas rulete
clc, clear all
Ns = [ 100 100000 ];
trials = 10;
for j = 1:length(Ns)
    N = Ns(j);
    zaud_skaits = 0;
    figure
    hold on
    for t = 1:trials
        results = zeros(N, 1);    
        S = 0;
        
        for i=1:N            
            if rand * 37 > 1
                S = S - 1;
                zaud_skaits = zaud_skaits + 1;
            else
                S = S + 36;
            end;
            results(i) = S;
        end;
        
        plot(results)
    end;

    xlabel('N')
    title(sprintf('Spēļu gaita pie N = %d', N))
    hold off
    disp([ 'Zaudējuma varbūtība pie N = ' num2str(N) ': ' num2str(zaud_skaits / (N * trials)) ]);
end;
disp([ 'Teorētiskā zaudējuma varbūtība: ' num2str(36/37)]);

%% 3. uzdevums - Bifona adata
clc
clearvars

ns = [ 100 10000 1000000 ];
conv_n_idx = 2;
deltas = zeros(ns(conv_n_idx), 1);
d = 2; l = 1.5;
for j = 1:length(ns)
    n = ns(j);
    h = 0;
    pi_est = 0;
    for i = 1:n
        x = rand * d / 2; % random x - attālums no vidusspunkta līdz tuvākajai līnijai x = U[0, d/2]
        alpha = rand * pi / 2; % random alpha - šaurais leņķis starp adatu un līniju alpha = U[0, pi/2]
    
        if (x <= l / 2 * cos(alpha)) % vai adata trāpa uz līnijas
            h = h + 1;
        end

        if (h > 0)
            pi_est = 2 * l / d * i / h;
        end

        if (j == conv_n_idx)
            deltas(i) = pi_est;
        end
    end
    
    disp([ 'Pie N = ' num2str(n) ' novērtētais Pi = ' num2str(pi_est) ])
end
hold on
plot(deltas)
plot(ones(ns(conv_n_idx), 1)* pi, '--g')
eps = 0.05;
plot(ones(ns(conv_n_idx), 1)* pi - eps, '--r')
plot(ones(ns(conv_n_idx), 1)* pi + eps, '--r')
hold off
title('Konverģences līkne')

%% 4. uzdevums - noteiktais integrālis
clc, clear all
syms x
fn2 = @(x) (x >= -7 & x <= -3).* (2 .* sqrt(x + 7)) + ...
     (x > -3 & x <= -1).*(0.5 .* (x + 1) .^ 2 + 2) + ...
     (x >= -1 & x < -0.75).*((x + 2) .* 8 - 6) + ...
     (x >= -0.75 & x < -0.5).*((-1.75-x) .* 4 + 8) + ...
     (x >= -0.5 & x < 0.5) .* 3 + ...
     (x >= 0.5 & x < 0.75) .* ((-1.75+x) .* 4 + 8) + ...
     (x >= 0.75 & x < 1) .* ((2 - x) .* 8 - 6) + ...
     (x >= 1 & x < 3) .* (0.5 .* (x - 1) .^ 2 + 2) + ...
     (x >= 3 & x <= 7) .* (2 .* sqrt(7-x));
fplot(fn2, [-8, 8])
ylim([-1 6])
grid on

ns = [ 100, 10000, 1000000 ];
for j = 1:length(ns)
    h = 0;
    n = ns(j);
    for i = 1:n
        x = rand * 14 - 7; % x ~ U[-7,7]
        y = rand * 4; % y ~ U[0, 4]
    
        if fn2(x) >= y
            h = h + 1;
        end;
    end;
    
    est_area = 14 * 4 * h / n; % -7 līdz 7 = 14
    disp( [ 'pie N = ' num2str(n) ' novērtētais laukums = ' num2str(est_area) ])
end;
% integrāļā precīzas vērtības aprēķins
syms fn(x)
fn(x) = piecewise( ... 
    x >= -7 & x <= -3, 2 * sqrt(x + 7), ...
    x > -3 & x <= -1, 0.5 * (x + 1) ^ 2 + 2, ...
    x >= -1 & x < -0.75, (x + 2) * 8 - 6, ...
    x >= -0.75 & x < -0.5, (-1.75-x) * 4 + 8, ...
    x >= -0.5 & x < 0.5, 3, ...
    x >= 0.5 & x < 0.75, (-1.75+x) * 4 + 8, ...
    x >= 0.75 & x < 1, (2 - x) * 8 - 6, ...
    x >= 1 & x < 3, 0.5 * (x - 1) ^ 2 + 2, ...
    x >= 3 & x <= 7, 2 * sqrt(7-x), ...
    0);
exact_area = double(int(fn, -7, 7));
disp([ 'Precīza laukuma vērtība = ' num2str(exact_area) ])

%% 5. uzdevums - loka šaušana
clc, clear all
x_mean = 0; x_sigma = 60; y_mean = 0; y_sigma = 30;
xi = [0 -46.6 -4 33 33 73 33 ]; yi = [-38 0 55 55 30 30 -38 ];
pgon = polyshape(xi, yi);
plot(pgon)
grid on
Ns = [ 1000 100000 ];
for j = 1:length(Ns)
    N = Ns(j); attempts = zeros(N, 1);
    for i = 1:N
        hits = 0;
        for k = 1:10
            x = normrnd(x_mean, x_sigma); y = normrnd(y_mean, y_sigma);
            if (pgon.isinterior(x, y)) 
                hits = hits + 1;
            end
        end
        attempts(i) = hits;
    end
    figure, histogram(attempts)
    title(sprintf('N = %d', N))
    disp([ 'N = ' num2str(N) ])
    disp([ 'Varbūtība - tieši 5 trāpījumi = ' num2str(numel(attempts(attempts == 5)) / N) ])
    disp([ 'Varbūtība - vismaz 5 trāpījumi = ' num2str(numel(attempts(attempts >= 5)) / N) ])
end

%% 6. uzdevums - papīra ražošana
clc, clear all

t_max = 900; c_base = 2500; c_workday = 1500; c_sunday = 250;
l_mean = 24.4; l_sigma = sqrt(4.5);

trials = 100;
min_c_avg = +inf;
best_M = 0;
best_f = 0;

for M = 1:6 % nedēļu skaits līdz preventīvai nomaiņai
    for f = 0:0.05:1 % cik liela daļa no nedēļas skaitās "pagājusi nedēļa" priekš preventīvas nomaiņas
        cs = zeros(trials, 1);
        for i = 1:trials
            t = 0;
            c = 0;
            t_p = t + M * 7;
            while (t < t_max)
                l = normrnd(l_mean, l_sigma); % lentas dzīves ilgums
                t_B = t + l;% sagaidāmās avārijas laiks        
                if (t_B < t_p)
                    % avārija
                    t = t_B;
                    c = c + c_workday;
                    N = floor(t); % eksperimenta diena
                    D = mod (N, 7); % nedēļas diena
                    f_fakt = D / 6; % cik liela daļa no nedēļas ir pagājusi
                    if (f_fakt < f)
                        t_p = N + M * 7 - D;
                    else
                        t_p = N + M * 7  - D + 7;
                    end
                    % disp([ 'Avārija @ t = ' num2str(t) ', c = ' num2str(c) ])
                else
                    % preventīvā nomaiņa
                    t = t_p;
                    c = c + c_sunday;
                    t_p = t + M * 7;
                    % disp([ 'Preventīva nomaiņa @ t = ' num2str(t) ', c = ' num2str(c) ])
                end
                
                c = c + c_base;
            end
        
            cs(i) = c / t;
        end

        avg_c = mean(cs);

        % disp([ 'Vidējā cena = ' num2str(avg_c) ' pie M = ' num2str(M) ' un f = ' num2str(f) ' (t = '  num2str(t) ')' ])

        if (avg_c < min_c_avg)
            min_c_avg = avg_c;
            best_M = M;
            best_f = f;
        end
    end
end

disp([ 'Labākā vidējā cena = ' num2str(min_c_avg) ' pie M = ' num2str(best_M) ' un f = ' num2str(best_f)  ])