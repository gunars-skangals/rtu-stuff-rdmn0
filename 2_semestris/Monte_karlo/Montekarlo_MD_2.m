%% 1. gadījuma skaitļu ģenerēšana
clear all, clc
n = 10000;

% a
A = zeros(n, 1);
for i=1:n
    A(i) = custom_rng();
end;
plot(A)
title('Ģenerēšanas rezultāti kā secība')
figure
histogram(A, 10)
title('Ģenerēšanas rezultātu histogramma')

% b
Bx = zeros(n, 1);
By = zeros(n, 1);
for i=1:n
    Bx(i) = custom_rng();
    By(i) = custom_rng();
end;
figure
scatter(Bx, By)
title('Ģenerēšanas rezultāti plaknē')

% c
Cx = zeros(n, 1);
Cy = zeros(n, 1);
Cz = zeros(n, 1);
for i=1:n
    Cx(i) = custom_rng();
    Cy(i) = custom_rng();
    Cz(i) = custom_rng();
end;
figure
scatter3(Cx, Cy, Cz)
title('Ģenerēšanas rezultāti 3 dimensijās')

%% d., e. - periods / aperiods
clear all, clc

P_length = 2000;
P = zeros(P_length, 1);
match_start = 0;
match_offset = 0;
match_count = 0;

tic
max_i = 2 ^ 33;
for i = 1:max_i
    rand = custom_rng();

    if (i < P_length)
        P(i) = rand;
    end

    if (match_start > 0)
        k = i - match_start + match_offset;
        if (P(k) == rand)
            match_count = match_count + 1;
            if (match_count > 100)
                disp(['sakritība sākot ar elementu: ' num2str(match_start)])
                disp(['aperiods: ' num2str(match_offset)])
                disp(['periods: ' num2str(match_start - match_offset)])
                break
            end
        else
            % atkārtojumu sērija pārtraukta
            match_count = 0;
            match_start = 0;
            match_offset = 0;
        end

    end

    if (i > P_length & match_start == 0)
        for j = 1:length(P)
            if (rand == P(j))
                match_start = i;
                match_offset = j;
                break;
            end
        end
    end
end
toc
% sakritība sākot ar elementu: 2147483399
% aperiods: 1
% periods: 2147483398
% Elapsed time is 4896.515391 seconds.
%% f. vienmērīguma analīze
clc, clear all
n = 100000;
sample = zeros(n, 1);
for i=1:n
    sample(i) = custom_rng();
end;

nbins = 100; % number of bin
edges = linspace(0,1, nbins+1); % edges of the bins
E = n/nbins*ones(nbins,1); % expected value (equal for uniform dist)
[h,p,stats] = chi2gof(sample, 'Expected',E,'Edges',edges);
h
p
% h = 0, p = 0.9817
% ļoti liela varbūtība, ka sadalījums ir vienmērīgs
%% g. korelācijas īpašības
clc, clear all
n = 1000;
sample = zeros(n, 1);
for i=1:n
    sample(i) = custom_rng();
end;
[c,lags] = xcorr(sample);
stem(lags,c)
title('Autokorelācija')
%% g. korelācijas īpašības - iebūvētais rand
clc, clear all
n = 1000;
sample = zeros(n, 1);
for i=1:n
    sample(i) = rand;
end;
[c,lags] = xcorr(sample);
stem(lags,c)
title('Autokorelācija (iebūvētais rand)')
%% h. uzdevuma risināšana
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
ns = [ 100, 10000, 1000000 ];
for j = 1:length(ns)
    h = 0;
    n = ns(j);
    for i = 1:n
        x = custom_rng() * 14 - 7; % x ~ U[-7,7]
        y = custom_rng() * 4; % y ~ U[0, 4]
    
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
% iebūvētais rand
% pie N = 100 novērtētais laukums = 39.76
% pie N = 10000 novērtētais laukums = 37.9176
% pie N = 1000000 novērtētais laukums = 38.2412
% Precīza laukuma vērtība = 38.25
% custom_rand:
% pie N = 100 novērtētais laukums = 38.64
% pie N = 10000 novērtētais laukums = 38.276
% pie N = 1000000 novērtētais laukums = 38.2725
%% II - normāla gadījuma lieluma ģenerēšana
clear all, clc
syms x
n = 10000;
polya_sample = zeros(n, 1);
builtin_sample = randn(n, 1);
polya_inv = @(x) (1 + (x >= 0.5) .* -2) .* sqrt(-pi ./ 2 .* log(1 - (x .* 2 - 1) .^ 2));
for i=1:n
    polya_sample(i) = polya_inv(rand);
end;
disp([ 'Kolmogorova-Smirnova tests' ])
[h, p] = kstest(polya_sample); % Kolmogorova-Smirnova tests, noklusējuma parametri - salīdzina ar normālsadalījumu
disp([ 'Polya: h = ' num2str(h) ', p = ' num2str(p) ])
[h, p] = kstest(builtin_sample);
disp([ 'Iebūvētais: h = ' num2str(h) ', p = ' num2str(p) ])
disp([ 'Shapiro-Wilks tests' ])
[h, p, W] = swtest(polya_sample);
disp([ 'Polya: h = ' num2str(h) ', p = ' num2str(p) ', W = ' num2str(W) ])
[h, p, W] = swtest(builtin_sample);
disp([ 'Iebūvētais: h = ' num2str(h) ', p = ' num2str(p) ', W = ' num2str(W) ])
%% III - Eiropas option līguma cena
clear all, clc

strike_price = 100;
T = 1;
s_0 = 100;
volatility = 0.3;
risk_free_rate = 0.06;
n = 100000;
sample = zeros(n, 1);
estimates = zeros(n - 1, 1);
for i = 1:n
    % akcijas cena 
    S_T = s_0 .* exp((risk_free_rate - 0.5 .* (volatility .^ 2)) .* T + volatility .* randn .* sqrt(T));
    % diskontēta starpība starp vienošanās cenu un akcijas cenu
    optPrice = exp (- risk_free_rate * T) * (S_T - strike_price); 
    if (optPrice < 0)
        optPrice = 0;
    end
    sample(i) = optPrice;

    estimates(i) = mean(sample(1:i)); % novērtējums - vidējā vērtība no līdzšinējām
end

% iebūvētā funkcija šim aprēķinam
[Call,Put] = blsprice(s_0,strike_price,risk_free_rate,T,volatility);
Call

plot(estimates)
hold on
eps = 0.10;
convLine = ones(length(estimates), 1) * Call;
plot(convLine - eps, '--r')
plot(convLine + eps, '--r')
% plot(convLine, '--g')
title('Konverģence')
hold off
mean(sample)
figure
histogram(sample, [0 eps:5:100 100])
title('Histogramma')

%%
function result = custom_rng
    persistent current_x
    if isempty(current_x)
       current_x = 1;
    end
    % L'Ecuyer, 14. slaids
    m = 2147483399;
    a = 40692;

    current_x = rem(a * current_x, m);
    result = current_x / m;  % [0,1)
end

% https://github.com/antagomir/scripts/blob/master/matlab/swtest.m
function [H, pValue, W] = swtest(x, alpha, tail)
    %SWTEST Shapiro-Wilk parametric hypothesis test of composite normality.
    %   [H, pValue, SWstatistic] = SWTEST(X, ALPHA, TAIL) performs
    %   the Shapiro-Wilk test to determine if the null hypothesis of
    %   composite normality is a reasonable assumption regarding the
    %   population distribution of a random sample X. The desired significance 
    %   level, ALPHA, is an optional scalar input (default = 0.05).
    %   TAIL indicates the type of test (default = 1).
    %
    %   The Shapiro-Wilk hypotheses are: 
    %   Null Hypothesis:        X is normal with unspecified mean and variance.
    %      For TAIL =  0 (2-sided test), alternative: X is not normal.
    %      For TAIL =  1 (1-sided test), alternative: X is upper the normal.
    %      For TAIL = -1 (1-sided test), alternative: X is lower the normal.
    %
    %   This is an omnibus test, and is generally considered relatively
    %   powerful against a variety of alternatives.
    %   Shapiro-Wilk test is better than the Shapiro-Francia test for
    %   Platykurtic sample. Conversely, Shapiro-Francia test is better than the
    %   Shapiro-Wilk test for Leptokurtic samples.
    %
    %   When the series 'X' is Leptokurtic, SWTEST performs the Shapiro-Francia
    %   test, else (series 'X' is Platykurtic) SWTEST performs the
    %   Shapiro-Wilk test.
    % 
    %    [H, pValue, SWstatistic] = SWTEST(X, ALPHA, TAIL)
    %
    % Inputs:
    %   X - a vector of deviates from an unknown distribution. The observation
    %     number must exceed 3 and less than 5000.
    %
    % Optional inputs:
    %   ALPHA - The significance level for the test (default = 0.05).
    %
    %   TAIL  - The type of the test (default = 1).
    %  
    % Outputs:
    %  SWstatistic - The test statistic (non normalized).
    %
    %   pValue - is the p-value, or the probability of observing the given
    %     result by chance given that the null hypothesis is true. Small values
    %     of pValue cast doubt on the validity of the null hypothesis.
    %
    %     H = 0 => Do not reject the null hypothesis at significance level ALPHA.
    %     H = 1 => Reject the null hypothesis at significance level ALPHA.
    %
    
    %
    % References: Royston P. "Algorithm AS R94", Applied Statistics (1995) Vol. 44, No. 4.
    %   AS R94 -- calculates Shapiro-Wilk normality test and P-value
    %   for sample sizes 3 <= n <= 5000. Handles censored or uncensored data.
    %   Corrects AS 181, which was found to be inaccurate for n > 50.
    %
    
    %
    % Ensure the sample data is a VECTOR.
    %
    
    if numel(x) == length(x)
        x  =  x(:);               % Ensure a column vector.
    else
        error(' Input sample ''X'' must be a vector.');
    end
    
    %
    % Remove missing observations indicated by NaN's and check sample size.
    %
    
    x  =  x(~isnan(x));
    
    if length(x) < 3
       error(' Sample vector ''X'' must have at least 3 valid observations.');
    end
    
    if length(x) > 5000
        warning('Shapiro-Wilk test might be inaccurate due to large sample size ( > 5000).');
    end
    
    %
    % Ensure the significance level, ALPHA, is a 
    % scalar, and set default if necessary.
    %
    
    if (nargin >= 2) && ~isempty(alpha)
       if numel(alpha) > 1
          error(' Significance level ''Alpha'' must be a scalar.');
       end
       if (alpha <= 0 || alpha >= 1)
          error(' Significance level ''Alpha'' must be between 0 and 1.'); 
       end
    else
       alpha  =  0.05;
    end
    
    %
    % Ensure the type-of-test indicator, TAIL, is a scalar integer from 
    % the allowable set {-1 , 0 , 1}, and set default if necessary.
    %
    
    if (nargin >= 3) && ~isempty(tail)
       if numel(tail) > 1
          error('Type-of-test indicator ''Tail'' must be a scalar.');
       end
       if (tail ~= -1) && (tail ~= 0) && (tail ~= 1)
          error('Type-of-test indicator ''Tail'' must be -1, 0, or 1.');
       end
    else
       tail  =  1;
    end
    
    % First, calculate the a's for weights as a function of the m's
    % See Royston (1995) for details in the approximation.
    
    x       =   sort(x); % Sort the vector X in ascending order.
    n       =   length(x);
    mtilde  =   norminv(((1:n)' - 3/8) / (n + 0.25));
    weights =   zeros(n,1); % Preallocate the weights.
    
    if kurtosis(x) > 3
        
        % The Shapiro-Francia test is better for leptokurtic samples.
        
        weights =   1/sqrt(mtilde'*mtilde) * mtilde;
    
        %
        % The Shapiro-Francia statistic W is calculated to avoid excessive rounding
        % errors for W close to 1 (a potential problem in very large samples).
        %
    
        W   =   (weights' * x) ^2 / ((x - mean(x))' * (x - mean(x)));
    
        nu      =   log(n);
        u1      =   log(nu) - nu;
        u2      =   log(nu) + 2/nu;
        mu      =   -1.2725 + (1.0521 * u1);
        sigma   =   1.0308 - (0.26758 * u2);
    
        newSFstatistic  =   log(1 - W);
    
        %
        % Compute the normalized Shapiro-Francia statistic and its p-value.
        %
    
        NormalSFstatistic =   (newSFstatistic - mu) / sigma;
        
        % the next p-value is for the tail = 1 test.
        pValue   =   1 - normcdf(NormalSFstatistic, 0, 1);
        
    else
        
        % The Shapiro-Wilk test is better for platykurtic samples.
    
        c    =   1/sqrt(mtilde'*mtilde) * mtilde;
        u    =   1/sqrt(n);
    
        PolyCoef_1   =   [-2.706056 , 4.434685 , -2.071190 , -0.147981 , 0.221157 , c(n)];
        PolyCoef_2   =   [-3.582633 , 5.682633 , -1.752461 , -0.293762 , 0.042981 , c(n-1)];
    
        PolyCoef_3   =   [-0.0006714 , 0.0250540 , -0.39978 , 0.54400];
        PolyCoef_4   =   [-0.0020322 , 0.0627670 , -0.77857 , 1.38220];
        PolyCoef_5   =   [0.00389150 , -0.083751 , -0.31082 , -1.5861];
        PolyCoef_6   =   [0.00303020 , -0.082676 , -0.48030];
    
        PolyCoef_7   =   [0.459 , -2.273];
    
        weights(n)   =   polyval(PolyCoef_1 , u);
        weights(1)   =   -weights(n);
    
        % Special attention when n=3 (this is a special case).
        if n == 3
            weights(1)  =   0.707106781;
            weights(n)  =   -weights(1);
        end
    
        if n >= 6
            weights(n-1) =   polyval(PolyCoef_2 , u);
            weights(2)   =   -weights(n-1);
        
            count  =   3;
            phi    =   (mtilde'*mtilde - 2 * mtilde(n)^2 - 2 * mtilde(n-1)^2) / ...
                    (1 - 2 * weights(n)^2 - 2 * weights(n-1)^2);
        else
            count  =   2;
            phi    =   (mtilde'*mtilde - 2 * mtilde(n)^2) / ...
                    (1 - 2 * weights(n)^2);
        end
    
        %
        % The vector 'WEIGHTS' obtained next corresponds to the same coefficients
        % listed by Shapiro-Wilk in their original test for small samples.
        %
    
        weights(count : n-count+1)  =  mtilde(count : n-count+1) / sqrt(phi);
    
        %
        % The Shapiro-Wilk statistic W is calculated to avoid excessive rounding
        % errors for W close to 1 (a potential problem in very large samples).
        %
    
        W   =   (weights' * x) ^2 / ((x - mean(x))' * (x - mean(x)));
    
        %
        % Calculate the significance level for W (exact for n=3).
        %
    
        newn    =   log(n);
    
        if (n > 3) && (n <= 11)
        
            mu      =   polyval(PolyCoef_3 , n);
            sigma   =   exp(polyval(PolyCoef_4 , n));    
            gam     =   polyval(PolyCoef_7 , n);
        
            newSWstatistic  =   -log(gam-log(1-W));
        
        elseif n >= 12
        
            mu      =   polyval(PolyCoef_5 , newn);
            sigma   =   exp(polyval(PolyCoef_6 , newn));
        
            newSWstatistic  =   log(1 - W);
        
        elseif n == 3
            mu      =   0;
            sigma   =   1;
            newSWstatistic  =   0;
        end
    
        %
        % Compute the normalized Shapiro-Wilk statistic and its p-value.
        %
    
        NormalSWstatistic       =   (newSWstatistic - mu) / sigma;
        
        % The next p-value is for the tail = 1 test.
        pValue       =   1 - normcdf(NormalSWstatistic, 0, 1);
    
        % Special attention when n=3 (this is a special case).
        if n == 3
            pValue  =   1.909859 * (asin(sqrt(W)) - 1.047198);
            NormalSWstatistic =   norminv(pValue, 0, 1);
        end
        
    end
    
    % The p-value just found is for the tail = 1 test.
    if tail == 0
        pValue = 2 * min(pValue, 1-pValue);
    elseif tail == -1
        pValue = 1 - pValue;
    end
    
    %
    % To maintain consistency with existing Statistics Toolbox hypothesis
    % tests, returning 'H = 0' implies that we 'Do not reject the null 
    % hypothesis at the significance level of alpha' and 'H = 1' implies 
    % that we 'Reject the null hypothesis at significance level of alpha.'
    %
    
    H  = (alpha >= pValue);
end