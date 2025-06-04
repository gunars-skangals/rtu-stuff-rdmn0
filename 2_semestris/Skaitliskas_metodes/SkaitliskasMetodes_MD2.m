%% 1.uzdevums
clear all, clc

xnodes = [ 1.1 2.6 4.1 5.6 7.1 8.6 10.10 11.6 ]';
ynodes = [ -0.86 2.7 3.22 1.68 3.07 7.02 8.48 7.06 ]';

plot(xnodes, ynodes, '*')

[p,reg]=fit(xnodes,ynodes,'poly3');
figure
plot(p,xnodes,ynodes)
title('Approximation by a polynomial of degree 3')

[p,reg]=fit(xnodes,ynodes,'poly5');
figure
plot(p,xnodes,ynodes)
title('Approximation by a polynomial of degree 5')

[p,reg]=fit(xnodes,ynodes,'poly6');
figure
plot(p,xnodes,ynodes)
title('Approximation by a polynomial of degree 6')

disp('Labākā aproksimācija tika iegūta ar 6. pakāpes polinomu:')
p
%% 2.uzdevums
clear all, clc
a=readtable('/home/g/source/uni/faili/Apartments.xlsx'); % taka specifiski Linux sistēmai
test=a.Floor;
output=a.Price;
plot(test,output, '*')
%% Linear model
x=[test];
c=fitlm(x,output)
plot(c) % g) punkts redzams grafikā 
p1 = c.Coefficients(1, 1).Estimate;
p2 = c.Coefficients(2, 1).Estimate;
disp([ 'a) Cena = ' num2str(p1) ' + ' num2str(p2) ' * Floor'])
disp([ 'b) Koeficientu interpretācija: Dzīvokļa pamata cena ir ' num2str(p1) ', plus piemaksa par katru stāvu ' num2str(p2) ' apmērā'])
disp([ 'c) Abu koeficientu p vērtības ir būtiski mazākas par 0.05 - ' num2str(c.Coefficients(1, 4).pValue) ' un ' num2str(c.Coefficients(2, 4).pValue) ', tātad modelis ir statistiski nozīmīgs' ])
disp([ 'd) Modeļa determinācijas koeficients ir ' num2str(c.Rsquared.Ordinary) ' - modelis izskaidro apmēram ceturto daļu no cenu variācijas' ])
%% e)
plot(c.Residuals.Raw)
title ('Kļūdas')
figure
plot(c.Residuals.Standardized)
title ('Standartizētas kļūdas')
figure
plotResiduals(c, 'probability')
disp([ 'f) Modeļa kļūdas ir ļoti tuvas normālam sadalījumam' ])
%% g)
[p,reg]=fit(test,output,'poly1') % cita fit metode, bet iegūtās koeficientu vērtības ir tādas pašas
f1 = 1:30;
p1=predint(p,f1,0.95); % 
plot(test,output,'*')
indexes = [1:length(test)]';
labels = cellstr(num2str(indexes));
text(test, output, labels, 'Vert','bottom', 'Horiz','left', 'FontSize',7) 
hold on
plot(f1,p1,'m--')
hold off
disp([ 'h) 3., 21. un 50. novērojumi ir izlēcēji - atrodas ārpus prognozes zonas' ])

%% 3.uzdevums
clear all, clc
a=readtable('/home/g/source/uni/faili/Cobb_Douglas_model.xlsx'); % taka specifiski Linux sistēmai
% ln (Q) = ln (alpha * L ^ beta * K ^ gamma)
% ln (Q) = ln (alpha) + beta * ln(L) + gamma * ln(K)
lnoutput = log(a.output);
lnlabor = log(a.labor);
lncapital = log(a.capital);
inputs = [lnlabor lncapital];
b=fitlm(inputs, lnoutput)
disp(' ')
disp('Atbilde:')
disp(['Koeficientu p-vērtības ir ļoti zemas (< 10^-4), modelis ir statistiski nozīmīgs'])
disp(['Modeļa determinācijas koeficients ir ļoti augsts: ' num2str(b.Rsquared.Ordinary) ', modelis izskaidro vairāk 95% no variācijas (logaritmiskajā skalā) un ir ļoti precīzs'])