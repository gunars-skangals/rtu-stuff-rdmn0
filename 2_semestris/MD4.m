%% 1.uzdevums
% https://finance.yahoo.com/quote/BRK-B/options/?date=1750377600
% skatīts 24.03.2025 19:32
S0 = 523.28;
K = 520;
r = 0.041820;
T = wrkdydif(today, 'june 20 2025' ) / 252;
sigma = 0.2193;

[OptCall, OptPut] = blsprice(S0, K, r, T, sigma)
% aprēķināta cena: 27.7416
% yahoo cena: 23.67 

[CallVal, PutVal] = blsdelta(S0, K, r, T, sigma, 0)
% Delta = 0.5828
GammaVal = blsgamma(S0, K, r, T, sigma, 0)
% Gamma = 0.0067
VegaVal = blsvega(S0, K, r, T, sigma, 0)
% Vega = 103.7327

N=1000;% number call sold
Cal_position= N*OptCall
%Delta neutral portfolio
N*CallVal*S0
D_neutral=-Cal_position+N*CallVal*S0

disp([' Delta-neitrāls portfelis: ' num2str(N) ' call opcijas, ' num2str(D_neutral) ' akcijas'])
%  Delta-neitrāls portfelis: 1000 call opcijas, 277212.277 akcijas