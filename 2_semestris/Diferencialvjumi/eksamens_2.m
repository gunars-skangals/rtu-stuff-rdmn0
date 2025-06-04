clc
m = 1.5;
% mean-reverting ātrums
k1 = m * 0.1;
k2 = m * 0.3;
% ilgtermiņa vidējie
theta1 = m * 0.025;
theta2 = m * 0.01;
% svārstīgums
sigma1 = m * 0.01;
sigma2 = m * 0.015;
% korelācija starp svārstībām
rho = 0.4;

y10 = m * 0.02;
y20 = m * 0.01;

T = 10;
N = 120;

t = linspace(0, T, N + 1);
A = twoFactorVasicek(k1, theta1, sigma1, y10, k2, theta2, sigma2, y20, T, N, rho);
plot(t, A, 'LineWidth', 1.5);
xlabel('Laiks');
ylabel('Procentu likme');
legend('a', 'b', 'c')
title('Divu faktoru Vasicek modelis');
grid on

%% 2.punkts
fwd_rate_fn = @(r1, t1, r2, t2) (r2 .* t2 - r1 .* t1) ./ (t2 - t1);

rate_3m = interp1(t, A, 0.25); % likme 3 mēnešiem
disp([ 'Procentu likme pēc 3 mēnešiem: ' num2str(rate_3m)])
rate_9m = interp1(t, A, 0.75); % likme 9 mēnešiem
disp([ 'Procentu likme pēc 9 mēnešiem: ' num2str(rate_9m)])
rate_15m = interp1(t, A, 1.25); % likme 15 mēnešiem
disp([ 'Procentu likme pēc 15 mēnešiem: ' num2str(rate_15m)])

fwd_3_9 = fwd_rate_fn(rate_3m, 0.25, rate_9m, 0.75);
disp([ 'Procentu likme 3-9 mēnešu periodam: ' num2str(fwd_3_9)])

fwd_9_15 = fwd_rate_fn(rate_9m, 0.75, rate_15m, 1.25);
disp([ 'Procentu likme 9-15 mēnešu periodam: ' num2str(fwd_9_15)])

%% 
function P = twoFactorVasicek(k1, theta1, sigma1, b1_0, ...
                             k2, theta2, sigma2, b2_0, ...
                             T, N, rho)
    dt = T / N; % laika solis
    
    P = zeros(N+1,1);

    b1 = @(tau) 1 / k1 * (1 - exp( -k1 * tau));
    b2 = @(tau) 1 / k2 * (1 - exp( -k2 * tau));
    
    % pieņem lambda = 0, attiecīgi theta_dash_i = theta_i

    gamma_i = @(k, theta_dash, sigma) k .^ 2 .* theta_dash - sigma .^ 2 ./ 2;

    gamma1 = gamma_i(k1, theta1, sigma1);
    gamma2 = gamma_i(k2, theta2, sigma2);

    sigma12 = sigma1*sigma2*rho;

    a = @(t, B1, B2) gamma1 .* (B1 - t) ./ (k1 .^ 2) ...
          - (sigma1 .^ 2 .* B1 .^ 2) ./ 4 ./ k1 ...
          + gamma2 .* (B2 - t) ./ (k2 .^ 2) ...
          - sigma2 .^ 2 .* B2 .^ 2 ./ 4 ./ k2 ...
          + sigma12 ./ k1 ./ k2 * (t - B1 - B2 + 1 ./ (k1 + k2) .* (1 - exp(-(k1 + k2) .* t)));

    t = 0;

    for i = 0:N
        
        t = t + dt;

        B1 = b1(t);
        B2 = b2(t);

        P(i + 1) = (-a(t, B1, B2) + B1 * b1_0 + B2 * b2_0) / t;
    end
    
end

