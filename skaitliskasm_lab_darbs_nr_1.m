%% Lineâru vienâdojumu sistçmas 
% Sadaïa "Tieðâs metodes" 


%% 1.piemçrs (Gausa metode)
clear all, clc, format compact
A=[2,2,-1,1;4,3,-1,2;8,5,-3,4;3,3,-2,2];B=[4;6;12;6];
[row,col]=size(A), Aaug=[A B];
A_rank=rank(A), Aaug_rank=rank(Aaug)
sol=rref(Aaug)

% turpinâjums
X_name=['x1';'x2';'x3';'x4'];
X_value=sol(:,col+1);
disp('Atbilde:')
disp(' LVS ir saderîga un noteikta')
disp(' ( viens vienîgs atrisinâjums )')
solution = table(X_name,X_value)

%% 2.piemçrs (Gausa metode)


%% 3.piemçrs (Gausa metode)
clear all, clc, format compact
A=[9,-3,5,6;6,-2,3,1;3,-1,3,14];B=[4;5;-8];
Aaug=[A B];
A_rank=rank(A), Aaug_rank=rank(Aaug)
sol=sym(rref(Aaug))

% turpinâjums
syms x2 x4, X_gen=sol(:,5)-sol(:,4).*x4-sol(:,2).*x2

% turpinâjums
disp('Atbilde:')
disp(' LVS ir saderîga un nenoteikta(bezg.daudz atrisinâjumu )')
disp(' x1=(13+13x4+x2)/3, x3=-7-9x4, x2 un x4-jebkuri reâli skaitïi')

%% 4.piemçrs (Gausa metode)    
clear all, clc, format compact
A=[1 -3 2;2 1 -4;5 -8 2];B=[1;5;8];
Aaug=[A B];
A_rank=rank(A), Aaug_rank=rank(Aaug)
sol=sym(rref(Aaug))

% bezgalīgi daudz atrisinājumu, z patvaļīgi

syms z, X_gen=sol(:,4)-sol(:,3).*z
subs(X_gen(2), z = 3)


%% 5.piemçrs (a)     
% Salîdzinâjums: rref un linsolve
clear all,clc,format compact
A=[4,2,-1;1,2,1;0,1,-1];B=[0;1;-3];
Aaug=[A B];
sol=sym(rref(Aaug))
sol_1=linsolve(A,B)

%% 5.piemçrs (b)    
% Salîdzinâjums: rref un linsolve
clear all,  clc, format compact
A=[ 1 5 4; 2 10 8; 3 15 12 ]; B=[ 1; 3; 5 ];
Aaug=[A B];
sol=sym(rref(Aaug)) % nav atrisinajumu
sol_1=linsolve(A,B)

%% 5.piemçrs (c)     
% Salîdzinâjums: rref un linsolve
clear all,  clc, format compact
A=[ 1 -3 2; 1 9 6; 1 3 4 ]; B=[ -1; 3; 1 ];
Aaug=[A B];
sol=sym(rref(Aaug)) % bezgalīgi daudz atrisinajumu
sol_1=linsolve(A,B)

%% 6.piemçrs. (LU metode)
clc, clear all, format compact
A=[2,1,1;1,-4,3;3,2,2];B=[7;2;13];
ni = parb(A);% pârbaude: galvenie minori ir vienâdi ar nulli
if ni==2 
  disp('Galvenie minori ir vienâdi ar nulli'), return
end
disp('Galvenie minori nav vienâdi ar nulli')
[L,U,P]=lu(A)

% turpinâjums
Amat=P*L*U % A matrica
Y=linsolve(L,P*B) % vai Y=L\(P*B)  
X=U\Y

% turpinâjums
X_name =['x1';'x2';'x3'];
disp('Atbilde:')
solution = table(X_name,X)

%% 7.piemçrs. Hoïecka metode
clear all, format compact, clc
A=[1,2,6;2,7,3;6,3,64]; B=[7;2;13];
check = isequal(A,A'); % pârbaude: vai matrica ir simetriskâ
if check ==0 % check = 1(TRUE) vai check = 0(FALSE)
disp('Koeficientu matrica nav simetriskâ'), return
end
ni=Hmetode(A); % pârbaude: vai matrica ir pozitîvi definçta
if ni==2
disp('Koeficientu matrica nav pozitîvi definçta'), return
end
disp(' Koeficientu matrica ir simetriskâ un pozitîvi definçta ')
% turpinâjums
L=chol(A,'lower')
Amat=L*L' % matrica A
Y=L\B % vai Y=linsolve(L,B)
X=L'\Y % vai X=linsolve(L’,Y)
% Ctrl+Enter
%% 8.piemçrs. Atstaroðanas metode 
clear all,format compact, clc
A=[1,2,6;2,7,3;6,3,64]; B=[7;2;13];
if det(A) == 0
disp(' Koeficientu matricas determinats ir = 0 '),return
end
disp(' Koeficientu matricas determinats nav vienâds ar 0 ')
[Q,R]=qr(A)
qtb = (Q') * B
X = rref([ R qtb ])
% Ctrl+Enter
%% Uzdevumi patstâvîgai risinâðanai
%% 1.uzd
clear all
m = [ 3 -1 4 -3; 2 3 1 5; 1 -5 -3 3 ]
rref(m)
%% 2.uzd
clear all
m = [ 3 -1 4 -1; 2 3 1 1; 7 5 6 1 ]
sol = sym(rref(m)) % bezgalīgi daudz, izsakām pēc z
%sol=sym(rref(Aaug))

% bezgalīgi daudz atrisinājumu, z patvaļīgi
syms z, X_gen=sol(:,4)-sol(:,3).*z

%% 3.uzd
clear all
A = [ 3 -3 3 6; 3 6 -3 7; 3 -9 3 1 ];
b = [ 6; 1; 7 ];
Aaug = [ A b ];
rankA = rank(A);
rankAaug = rank(Aaug);
if (rankA == rankAaug)
    disp('Ranki sakrīt')
else
    disp('Ranki nesakrīt')
    return
end

sol = sym(rref(Aaug))

syms x4, X_gen=sol(:,5)-sol(:,4).*x4
x1 = subs(X_gen(1), x4 = -2)

%% 4.uzd
clear all, clc
A = [1 2 1 -3; 2 1 -2 9; 3 3 -1 6; 4 5 0 3]; b = [ 5; -2; 3; 8];
Aaug = [ A b ];
rankA = rank(A);
rankAaug = rank(Aaug);
if (rankA == rankAaug)
    disp('Ranki sakrīt')
else
    disp('Ranki nesakrīt')
    return
end

sol = sym(rref(Aaug)) % bezgalīgi daudz atrisinājumu, izvēlas patvaļīgus x3, x4
syms x3 x4, X_gen=sol(:,5)-sol(:,4).*x4 - sol(:,3).*x3

X_part = subs(subs(X_gen, x3 = -1), x4 = 3)
x1 = X_part(2)

%% 5.uzdevums
clear all, clc
A = [ 4 4 5 5; 2 0 3 -1; 1 1 -5 0; 0 3 2 0 ]; b = [ 0; 10; -10; 1 ];
Aaug = [ A b ];
rankA = rank(A);
rankAaug = rank(Aaug);
if (rankA == rankAaug)
    disp('Ranki sakrīt')
else
    disp('Ranki nesakrīt')
    return
end

rref(Aaug) % Gausa

ni = parb(A);% pârbaude: galvenie minori ir vienâdi ar nulli
if ni==2 
  disp('Galvenie minori ir vienâdi ar nulli'), return
end
disp('Galvenie minori nav vienâdi ar nulli')
[L,U,P] = lu(A);
% turpinâjums
Amat=P*L*U % A matrica
Y=linsolve(L,P*b) % vai Y=L\(P*B)  
X=U\Y
% turpinâjums
X_name =['x1';'x2';'x3';'x4'];
disp('Atbilde:')
solution = table(X_name,X) % LU

%% 5 uzd - holecka
check = isequal(A,A'); % pârbaude: vai matrica ir simetriskâ
if check ==0 % check = 1(TRUE) vai check = 0(FALSE)
disp('Koeficientu matrica nav simetriskâ'), return
end
ni=Hmetode(A); % pârbaude: vai matrica ir pozitîvi definçta
if ni==2
disp('Koeficientu matrica nav pozitîvi definçta'), return
end
disp(' Koeficientu matrica ir simetriskâ un pozitîvi definçta ')
% turpinâjums
L=chol(A,'lower')
Amat=L*L' % matrica A
Y=L\b % vai Y=linsolve(L,B)
X=L'\Y % vai X=linsolve(L’,Y)

%% 5. uzd - atstarošanas
if det(A) == 0
disp(' Koeficientu matricas determinats ir = 0 '),return
end
disp(' Koeficientu matricas determinats nav vienâds ar 0 ')
[Q,R]=qr(A)
qtb = (Q') * b
X = rref([ R qtb ])


%% Ârçjâ junkcija. Holecka metode
function ni = Hmetode(A_mat)
ni=1; % pârbaude: vai matrica ir pozitîvi definçta
[row,col]=size(A_mat);
for i=1:row
if det(A_mat(1:i,1:i))>0
else ni=2; break
end
end
end
%% %% Ârçja funkcija. LU metode
function ni = parb(A_mat)
ni=1; % pârbaude: galvenie minori nav vienâdi ar nulli
[row,col]=size(A_mat);
for i=1:row
if det(A_mat(1:i,1:i))~=0
else ni=2; break
end
end
end