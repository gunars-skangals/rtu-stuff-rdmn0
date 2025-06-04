%% 1. uzdevums
clear all
a=readtable('/home/g/source/uni/faili/electr.xlsx');
test=a.Test;
output=a.Output;
[p,reg]=fit(test,output,'poly1')
%% Linear model
plot(p,test,output)
b=table(test,output);
x=[test];
c=fitlm(x,output)
%% Prediction interval
x1=[480:10:650];
p1=predint(p,x1,0.95);
% plot predicition interval
hold on
plot(x1,p1,'m--')
%% Histogram of residuals
figure
plotResiduals(c)
% plot histogram of residuals
%% Normal probability plot
figure
plotResiduals(c,'probability')
% normal probability plot for residuals
%% 2. uzdevums (Piegâde)
clear all
a=readtable('/home/g/source/uni/faili/piegade.xlsx');
test=a.Skaits;
sk = sqrt(a.Skaits);

test = sk;
output=a.Laiks;
[p,reg]=fit(test,output,'poly1')
%% Linear model

plot(p,test,output)
b=table(test,output);
x=[test];
c=fitlm(x,output)
%% Prediction interval
x1=[sqrt(50):sqrt(10):sqrt(300)];
p1=predint(p,x1,0.95);
% plot predicition interval
hold on
plot(x1,p1,'m--')
%% Histogram of residuals
figure
plotResiduals(c)
% plot histogram of residuals
%% Normal probability plot
figure
plotResiduals(c,'probability')
%% 3. uzdevums (Reklâma)
clear all
a=readtable('/home/g/source/uni/faili/reklama.xlsx');
test=a.Izdevumi;
output=a.Apjoms;
[p,reg]=fit(test,output,'poly1')
%% Linear model
plot(p,test,output)
b=table(test,output);
x=[test];
c=fitlm(x,output)
%% Prediction interval
x1=[20:10:100];
p1=predint(p,x1,0.95);
% plot predicition interval
hold on
plot(x1,p1,'m--')
%% Histogram of residuals
figure
plotResiduals(c)
% plot histogram of residuals
%% Normal probability plot
figure
plotResiduals(c,'probability')
%% 4. uzdevums (Izdevumi)
clear all
a=readtable('/home/g/source/uni/faili/transportation.xlsx');
test=a.Km;
output=a.Cost;
[p,reg]=fit(test,output,'poly1')
%% Linear model
plot(p,test,output)
b=table(test,output);
x=[test];
c=fitlm(x,output)
%% Prediction interval
min = floor(min(test) / 10) * 10;
max = ceil(max(test) / 10) * 10;
step = (max - min) / 100;
x1=[min:step:max];
p1=predint(p,x1,0.95);
% plot predicition interval
hold on
plot(x1,p1,'m--')
%% Histogram of residuals
figure
plotResiduals(c)
% plot histogram of residuals
%% Normal probability plot
figure
plotResiduals(c,'probability')
%% 5. uzdevums (Latvijas iedzîvotâju skaits)
clear all
a=readtable('/home/g/source/uni/faili/population.xlsx');
test=a.X;
output=a.Y;
[p,reg]=fit(test,output,'poly1')
%% Linear model
plot(p,test,output)
b=table(test,output);
x=[test];
c=fitlm(x,output)
%% Prediction interval
min_test = min(test);
max_test = max(test);
step = (max_test - min_test) / 50;
x1=[1:1:25];
p1=predint(p,x1,0.95);
% plot predicition interval
hold on
plot(x1,p1,'m--')
%% Histogram of residuals

figure
plotResiduals(c)
% plot histogram of residuals
%% Normal probability plot
figure
plotResiduals(c,'probability')
%% HATCO
clear all
a=readtable('/home/g/source/uni/faili/hatcodata.xls');
test=a.X5;
output=a.X9;
[p,reg]=fit(test,output,'poly1')
%% Linear model
plot(p,test,output)
b=table(test,output);
x=[test];
c=fitlm(x,output)
%% Prediction interval
min_test = min(test);
max_test = max(test);
step = (max_test - min_test) / 50;
x1=[min_test:step:max_test];
p1=predint(p,x1,0.95);
% plot predicition interval
hold on
plot(x1,p1,'m--')
%% Histogram of residuals

figure
plotResiduals(c)
% plot histogram of residuals
%% Normal probability plot
figure
plotResiduals(c,'probability')