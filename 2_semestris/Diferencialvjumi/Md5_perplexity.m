% https://www.mathworks.com/help/finance/hwv.html
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

T = 5;
N = 100;

twoFactorVasicek(k1, theta1, sigma1, 0.02, k2, theta2, sigma2, 0, T, N);
%%
function r = twoFactorVasicek(kappa_x, theta_x, sigma_x, x0, ...
                             kappa_y, theta_y, sigma_y, y0, ...
                             T, N)
    dt = T / N;
    x = zeros(N+1,1);
    y = zeros(N+1,1);
    r = zeros(N+1,1);
    x(1) = x0;
    y(1) = y0;
    for i = 1:N
        % Factor x
        mu_x = theta_x + (x(i) - theta_x)*exp(-kappa_x*dt);
        var_x = (sigma_x^2)/(2*kappa_x)*(1 - exp(-2*kappa_x*dt));
        x(i+1) = mu_x + sqrt(var_x)*randn;
        
        % Factor y
        mu_y = theta_y + (y(i) - theta_y)*exp(-kappa_y*dt);
        var_y = (sigma_y^2)/(2*kappa_y)*(1 - exp(-2*kappa_y*dt));
        y(i+1) = mu_y + sqrt(var_y)*randn;
        
        % Short rate is sum of factors
        r(i+1) = x(i+1) + y(i+1);
    end
    
    % Plot the short rate path
    t = linspace(0, T, N+1);
    plot(t, r, 'LineWidth', 1.5);
    xlabel('Time');
    ylabel('Short Rate r(t)');
    title('Two-Factor Vasicek Model Simulation');
    grid on;
end