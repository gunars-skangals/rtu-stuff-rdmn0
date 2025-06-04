clc
% mean-reverting ātrums
k1 = 0.1;
k2 = 0.3;
% ilgtermiņa vidējie
theta1 = 0.03;
theta2 = 0.02;
% svārstīgums
sigma1 = 0.01;
sigma2 = 0.015;
% korelācija starp svārstībām
rho = 0.4;

T = 30;
N = 120;

% 1)
y_i0 = [ 0.02 0; 0.02 0.01; 0.02 0.02 ];
lineStyles = [ "r" "g" "b" ];
t = linspace(0, T, N + 1);
for i = 1:3
   A = twoFactorVasicek(k1, theta1, sigma1, y_i0(i, 1), k2, theta2, sigma2, y_i0(i, 2), T, N, rho);
   plot(t, A, lineStyles(i), 'LineWidth', 1.5);
   hold on
end
xlabel('Laiks');
ylabel('Procentu likme');
legend('a', 'b', 'c')
title('Divu faktoru Vasicek modelis - 1');
grid on
hold off

% 2)
y_i0 = [ 0 0; 0 0.01; 0 0.02 ];
lineStyles = [ "r" "g" "b" ];
t = linspace(0, T, N + 1);
figure
for i = 1:3
   A = twoFactorVasicek(k1, theta1, sigma1, y_i0(i, 1), k2, theta2, sigma2, y_i0(i, 2), T, N, rho);
   plot(t, A, lineStyles(i), 'LineWidth', 1.5);
   hold on
end
xlabel('Laiks');
ylabel('Procentu likme');
legend('a', 'b', 'c')
title('Divu faktoru Vasicek modelis - 2');
grid on
hold off

% 3)
y_i0 = [ 0.04 0; 0.04 0.01; 0.04 0.02 ];
lineStyles = [ "r" "g" "b" ];
t = linspace(0, T, N + 1);
figure
for i = 1:3
   A = twoFactorVasicek(k1, theta1, sigma1, y_i0(i, 1), k2, theta2, sigma2, y_i0(i, 2), T, N, rho);
   plot(t, A, lineStyles(i), 'LineWidth', 1.5);
   hold on
end
xlabel('Laiks');
ylabel('Procentu likme');
legend('a', 'b', 'c')
title('Divu faktoru Vasicek modelis - 3');
grid on
hold off

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

