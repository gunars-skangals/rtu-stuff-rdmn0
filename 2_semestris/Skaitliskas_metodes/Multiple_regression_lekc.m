%% Transporta izdevumi
clear all
a=readtable('/home/g/source/uni/faili/transportation.xlsx');
%% Divfaktoru regresija
% X = km 
%Y = izdevumi
X1=a.Km;
X2=a.Days;
Y=a.Cost;
x=[ones(size(X1)) X1 X2];
b=regress(Y,x)

%% Grafiks
scatter3(X1,X2,Y,'filled')
hold on
x1fit = min(X1):100:max(X1);
x2fit = min(X2):2:max(X2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('Km')
ylabel('Days')
zlabel('Cost')
view(50,10)
hold off
%% HATCO
clear all, clc
a=readtable('/home/g/source/uni/faili/hatcodata.xls');
% X1 = delivery speed 
% X2 = price level
% X3 = price flexibility
% X4 = manufacturer image
% X5 = overall service
% X6 = salesforce image
% X7 = product quality
x1=a.X1;
x2=a.X2;
x3=a.X3;
x4=a.X4;
x5=a.X5;
x6=a.X6;
x7=a.X7;
y=a.X9;
x=[x1 x2 x3 x4 x5 x6 x7];
b=fitlm(x,y)

%% HATCO (stepwise model)
c=stepwiselm(x,y,'constant','Upper','linear')
% Linear stepwise regression without interaction terms
x8=a.X8;
x=[x1 x2 x3 x4 x5 x6 x7 x8];
b1=stepwiselm(x,y,'constant','Upper','linear')