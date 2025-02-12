%% 1.spēlmaņu uzdevums
clc
q = 0.495; %  pirmā spēlētaja uzvaras varbūtība
S_A = 5; % pirmā spēlētaja sākuma summa
n = 10; % spēļu skaits
trials = 1000; % procesa realizācijas
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
   
   % plot(results);   
   end_S(t) = SK_A;   
end;

histogram(end_S) % histogramma ar beigu summām pie dažādiem N

% otrā daļa - pirmās negatīvās naudas summas spēļu numuri
m = zeros(1, trials);
for t = 1:trials
   results = zeros(1, trials);
   SK_A = S_A;
   
   i = 0;   
   while SK_A > 0
       if rand > q
            SK_A = SK_A - 1;
       else
            SK_A = SK_A + 1;
       end;
       
       i = i + 1;
   end;
   m(t) = i;
end;

histogram(m) % TODO: better bins
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

%% 
% 2*sqrt(-abs(abs(x)-1)*abs(3-abs(x))/((abs(x)-1)*(3-abs(x))))(1+abs(abs(x)-3)/(abs(x)-3))sqrt(1-(x/7)^2)+(5+0.97(abs(x-.5)+abs(x+.5))-3(abs(x-.75)+abs(x+.75)))(1+abs(1-abs(x))/(1-abs(x)))
% (2.71052+(1.5-.5abs(x))-1.35526sqrt(4-(abs(x)-1)^2))sqrt(abs(abs(x)-1)/(abs(x)-1))+0.9 

% -3sqrt(1-(x/7)^2)sqrt(abs(abs(x)-4)/(abs(x)-4))
% abs(x/2)-0.0913722(x^2)-3+sqrt(1-(abs(abs(x)-2)-1)^2)

