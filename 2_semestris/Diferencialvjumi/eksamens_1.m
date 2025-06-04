%% 1.uzdevums
clc, clear all

% dotie
S0 = 66;
K = 70;
r = 0.03;
Tm = 10;
T = Tm / 12;
sigma = 0.7;
Sb = 45;

Smax = S0 * 10; % Smax būtiski ietekmē opcijas cenu - sigma ir diezgan liels, attiecīgi iespējams sasniegt diezgan augstu cenu
dS = 0.5
dT = 0.001
EuCallImpl2(S0, K, r, T, sigma, Smax, dS, dT) % 15.6469
% izskatās, ka šāds cenas apreķins konverģē uz ~15.65
% tālāka dS,dT samazināšana rezultātu vairs būtiski neietekmē
Call = blsprice(S0, K, r, T, sigma) % 15.7150

%% 2.uzdevums
clc, clear all

% dotie
S0 = 66;
K = 70;
r = 0.03;
Tm = 10;
T = Tm / 12;
sigma = 0.7;
Sb = 45;

Smax = S0 * 5; % opcija ir bezvērtīga ja S > K, augstākiem Smax vairs nav būtiskas nozīmes
dS = 0.5;
dT = 0.001;
ckPut = KnockOutPutOption(S0, K, r, T, sigma, Sb, Smax, dS, dT);
ckPut % 0.8360

[Call Put] = blsprice(S0, K, r, T, sigma);
Put % 17.9867

% Parastās Eiropas put opcijas cena ir būtiski lielāka
% To varētu skaidrot ar augsto sigma vērtību - pie šādiem nosacījumiem būtu ļoti
% iespējams sasniegt barjeras cenu, kas knockout opciju padarītu nevērtīgu

%% 3.uzdevums
clc, clear all
% dotie
S0 = 66;
K = 70;
r = 0.03;
Tm = 10;
T = Tm / 12;
sigma = 0.7;
Sb = 45;

Smax = S0 * 5;
dS = 0.5;
dT = 0.001;
omega = 1.2;
tol = 0.001;

putPrice = AmPutCK(S0, K, r, T, sigma, Smax, Sb, dS, dT, omega, tol);
putPrice % 16.5347

% Cena būtiski lielāka nekā Eiropas veida knock-out put opcijai un mazliet mazāka
% nekā parastai Eiropas put opcijai. Opcija šobrīd ir in-the-money
% Izskaidrojams ar to, ka šādu opciju var izpildīt pirms norādītā termiņa
% beigām - kad pamata vērtspapīra cena tuvojas barjerai, opciju var
% paspēt izpildīt un tā neiet zudumā kā tas būtu ar Eiropas veida opciju.
% Arī pretējā virzienā - ja vērtspapīra cena tuvojas strike price tuvu opcijas 
% derīguma termiņa beigam, opciju var izpildīt un iegūt kādu peļņu, kamēr
% Eiropas opcijas gadījumā svarīgs ir tikai stāvoklis opcijas termiņa
% beigās un augstas sigma vērtības apstākļos ļoti iespējams, ka tad opcija
% būs nevērtīga.

%% izmainīts EuPutImpl2
function  price = EuCallImpl2 (S0 ,K, r, T, sigma, Smax, dS , dt)
    % set up grid and adjust increments if necessary
    M = round(Smax/dS);
    % dS = Smax/M;
    N = round(T/dt) ;
    dt = T/N;
    matval = zeros(M+1,N+1);
    vetS = linspace(0,Smax,M+1)';
    veti = 0:M;
    vetj = 0:N;
    % set up boundary conditions
    matval(:,N+1) = max(vetS - K,0); % termiņa beigās vērtība ir (S - K)
    matval(1, :) = 0; % ja cena ir 0, opcijas vērtība ir 0
    matval(M+1,:) = (Smax - K)*exp(-r*dt*(N-vetj)); % pie max cenas: diskontē (Smax - K) atbilstoši laika momentam
    % set up coefficients
    a = 0.5*dt*(-sigma^2*veti + r).*veti;
    b = 1+ dt*(sigma^2*veti.^2 + r);
    c =- 0.5*dt*(sigma^2*veti + r).*veti;
    coeff = diag(a(3:M) ,-1) + diag(b(2:M)) + diag(c(2:M-1),1) ;
    [L,U] = lu(coeff);
    % solve the sequence of linear systems
    aux = zeros(M-1,1);
    for j=N:-1:1
        aux(1) = - a(2) * matval(1,j); % other term from BC is zero
        matval(2:M,j) = U\(L \ (matval(2:M,j+1) + aux));
    end
    % return price, possibly by linear interpolation outside the grid 
    price = interp1(vetS, matval(:,1), S0);
end;
%% DOPutCK.m ar komentāriem priekš sevis
function  price = KnockOutPutOption (S0 ,K, r, T, sigma,Sb, Smax, dS , dt)
    % set up grid and adjust increments if necessary
    M = round((Smax-Sb)/dS);
    dS = (Smax-Sb)/M;
    N = round(T/dt) ;
    dt = T/N;
    matval = zeros(M+1,N+1);
    vetS = linspace(Sb,Smax,M+1)'; % opcija pārtrauc pastāvēt, ja cena nokrītas līdz Sb, tādēļ tas ir minimums
    veti = vetS/dS;
    % vetj = 0:N;
    % set up boundary conditions
    matval(:,N+1) = max(K-vetS,0);
    matval(1, :)= 0 ; % pie min cenas (Sb) opcijas vērtība ir 0
    matval(M+1,:) = 0; % pie max cenas (Smax) opcijas vērtība ir 0
    % set up coefficients
    aalpha = 0.25*dt*(sigma^2*(veti.^2) - r*veti);
    bbeta = - dt*0.5*(sigma^2*veti.^2 + r);
    ggamma = 0.25*dt*(sigma^2*(veti.^2) + r*veti);
    M1 = -diag(aalpha(3:M) ,-1) + diag(1-bbeta(2:M)) - diag(ggamma(2:M-1),1) ;
    [L,U] = lu(M1);
    M2=diag(aalpha(3:M) ,-1) + diag(1+bbeta(2:M)) + diag(ggamma(2:M-1),1) ;
    % solve the sequence of linear systems
    
    for j=N:-1:1
        matval(2:M,j) = U\(L \ (M2*matval(2:M,j+1)));
    end
    % return price, possibly by linear interpolation outside the grid
    price = interp1(vetS, matval(:,1), S0);
end
%% AmPutCK.m
function price = AmPutCK(S0,K,r,T,sigma,Smax,Sb,dS,dt,omega,tol)    
    M = round((Smax-Sb)/dS);
    dS = (Smax-Sb)/M;
    N = round(T/dt) ;
    dt = T/N;
    vetS = linspace(Sb,Smax,M+1)'; % opcija pārtrauc pastāvēt, ja cena nokrītas līdz Sb, tādēļ tas ir minimums
    veti = vetS/dS;

    % set up boundary conditions
    payoff = max(K-vetS (2 : M),0);
    pastval = payoff; % values for the last layer

    % set up coefficients
    aalpha = 0.25*dt*(sigma^2*(veti.^2) - r*veti);
    bbeta = - dt*0.5*(sigma^2*veti.^2 + r);
    ggamma = 0.25*dt*(sigma^2*(veti.^2) + r*veti);
    M2=diag(aalpha(3:M) ,-1) + diag(1+bbeta(2:M)) + diag(ggamma(2:M-1),1) ;
    aux = zeros (M-1,1) ;
    newval=zeros(1,M-1);
    for j=N: -1 : 1
        aux(1) = 0; % robežvērtības ir 0
    % set up right hand side and initialize
        rhs = M2*pastval(:) + aux;
        oldval = pastval;
        error = realmax;
        while tol < error
            
            newval(1) = max ( payoff (1),...
                oldval(1) + omega/(1-bbeta(2)) * (...
                rhs(1) - (1-bbeta(2))*oldval(1) + ggamma(2)*oldval(2))) ;
            
            for k=2 : M-2
                newval(k) = max ( payoff (k),...
                oldval(k) + omega/(1-bbeta(k+1)) * (...
                rhs(k) + aalpha(k+1)*newval(k-1) - ...
                (1-bbeta(k+1))*oldval(k) + ggamma(k+1)*oldval(k+1))) ;
            end
            newval(M-1) = max( payoff (M-1), ...
                oldval(M-1) + omega/(1-bbeta(M)) * (...
                rhs(M-1) + aalpha(M)*newval(M-2) - ...
                (1-bbeta(M) ) *oldval (M-1))) ;
            % newval=newval';
            error = norm(newval - oldval);
            oldval = newval;
        end
        pastval = newval;
    end
    newval = [0, newval , 0] ; % add missing values - robežvērtība ir 0
    % return price, possibly by linear interpolation outside the grid
    price = interp1(vetS, newval, S0);
end