%% 1.spēlmaņu uzdevums
clc
q = 0.51 %  pirmā spēlētaja uzvaras varbūtība
S_A = 100 % pirmā spēlētaja sākuma summa
n = 10000 % spēļu skaits
trials = 500 % meģinajumu skaits
end_S = zeros(1, trials);
for t = 1:trials
   results = zeros(1, trials);
   SK_A = S_A;
   
   for i = 1:n
       if rand > q
            SK_A = SK_A - 1;
       else
            SK_A = SK_A + 1;
       end;
       
       results(i) = SK_A;
   end;
   
   plot(results);   
   end_S(t) = SK_A;   
end;

histogram(end_S)

% historgramma !! 
%% 2. Eiropas rulete
clc
clearvars

N = 100000;
trials = 3;

for t = 1:trials
    results = zeros(1, N);    
    S = 0;
    
    for i=1:N
        
        if rand * 37 > 1
            S = S - 1;
        else
            S = S + 36;
        end;
        results(i) = S;
    end;

    plot(results);
end;

