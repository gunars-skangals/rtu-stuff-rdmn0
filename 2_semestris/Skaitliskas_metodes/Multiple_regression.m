%% 1. uzdevums
% Referendums Kvebek� 1995.g.
clear all
a=readtable('/home/g/source/uni/faili/Quebec_ref.xlsx');
% 1. regresija
% X = �Par� balsu �patsvars 
%Y = neder�go bi�etenu �patsvars
X=a.PctYes;
Y=a.PctReject;
x=[X];
c=fitlm(x,Y)
plotResiduals(c,'fitted')
%% 1. regresija (bez Nr. 27)
clear all
a=readtable('/home/g/source/uni/faili/Quebec_m27.xlsx');
% X = �Par� balsu �patsvars 
%Y = neder�go bi�etenu �patsvars
X=a.PctYes;
Y=a.PctReject;
x=[X];
c=fitlm(x,Y)
plotResiduals(c,'fitted')
%% 2. regresija (bez Nr. 27)
% X = balsot�ju �patsvars, kuriem dzimt� valoda ir ang�u
% Y = neder�go bi�etenu �patsvars
X=a.PctAnglophone;
Y=a.PctReject;
x=[X];
c=fitlm(x,Y)
plotResiduals(c,'fitted')
%% 3. regresija (bez Nr. 27)
%X = balsot�ju �patsvars, kuriem dzimt� valoda 
%nav ne ang�u, ne fran�u
%Y = neder�go bi�etenu �patsvars
X=a.PctAllophone;
Y=a.PctReject;
x=[X];
c=fitlm(x,Y)
plotResiduals(c,'fitted')

%%
clear all
a=readtable('/home/g/source/uni/faili/hourly_wages.xlsx');
wage = a.wage;
isFemale = a.female;
isNonWhite = a.nonwhite;
isUnion = a.union;
education = a.education;
experience = a.exper;
age = a.age;
%%
x = [experience education isFemale isNonWhite isUnion];
% x=[isFemale isNonWhite isUnion education experience age];
b=fitlm(x,wage, 'Intercept',true)

% c = stepwiselm(x,wage,'constant','Upper','linear')
%%
