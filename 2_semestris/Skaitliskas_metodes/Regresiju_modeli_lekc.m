%% Nonlinear model
clear all, clc
a=readtable('/home/g/source/uni/faili/Cobb_Douglas.xlsx');
x1=a.lnlabor;
x2=a.lncapital;
y=a.lnoutput;
x=[x1 x2];
b=fitlm(x,y)

%% Linear model
x3=a.labor;
x4=a.capital;
y1=a.output;
xx=[x3 x4];
b1=fitlm(xx,y1)
%% GDP growth
% Log-lin modelis
clear all, clc
a=readtable('/home/g/source/uni/faili/GDP_growth.xlsx');
x1=a.time;
y=a.lnrgdp;
x=[x1];
b=fitlm(x,y)
%% Lin. modelis
x1=a.time;
y1=a.rgdp;
x=[x1];
b1=fitlm(x,y1)
%% CEO 
clear all, clc
a=readtable('/home/g/source/uni/faili/ceosal2.xls');
%%
x1=a.ceoten;
y=log(a.salary);
x=[x1];

b=fitlm(x,y)
%% Lin. modelis
b1=fitlm(x1,a.salary)

scatter(a.ceoten, a.salary)