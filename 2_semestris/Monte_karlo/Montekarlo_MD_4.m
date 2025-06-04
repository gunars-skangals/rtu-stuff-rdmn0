clc, clear all
% https://www.kaggle.com/datasets/marianadeem755/stock-market-data
data = readtable('/home/g/Documents/15 Years Stock Data of NVDA AAPL MSFT GOOGL and AMZN.csv');
deltas = diff(log(data.Close_AAPL));
hold on
histogram(deltas,'Normalization','pdf')
[f1,x1] = ksdensity(deltas);
plot(x1,f1, 'LineWidth', 3)
title('Ieejas datu histogramma un varbūtības blīvuma funkcija')
hold off
% izvēlas trokšņa std.deviation - mazāku par datu kopā novēroto
sigma = std(deltas) * 0.75; 
n = 25000; %simulāciju skaits
x=zeros(n,1);
x(1) = 0; %sākuma punkts
for i = 1:n-1
    x_c = normrnd(x(i), sigma);
    alpha = interp1(x1,f1,x_c) / interp1(x1,f1,x(i));
    if rand < min(1,alpha)
        x(i+1) = x_c; % ņem jauno vērtību
    else
        x(i+1) = x(i); % noraida jauno vērtību
    end
end
figure;
hold on
plot(x1,f1, 'LineWidth', 2)
histogram(x,'Normalization','pdf');
title('Simulētie dati un sākotnējā VBF')
xlim([-0.2 0.2]) % gara, negtiva aste
hold off