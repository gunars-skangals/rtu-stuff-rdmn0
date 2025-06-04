%% 1. gadījuma skaitļu ģeneracija
clear all, clc
n = 100000;

% a
A = zeros(n, 1);
for i=1:n
    A(i) = custom_rng();
end;
plot(A)
title('Ģenerēšanas rezultāti kā secība')
figure
histogram(A, 10)
title('Ģenerēšanas rezultātu histogramma')

% b
Bx = zeros(n, 1);
By = zeros(n, 1);
for i=1:n
    Bx(i) = custom_rng();
    By(i) = custom_rng();
end;
figure
scatter(Bx, By)
title('Ģenerēšanas rezultāti plaknē')

% c
Cx = zeros(n, 1);
Cy = zeros(n, 1);
Cz = zeros(n, 1);
for i=1:n
    Cx(i) = custom_rng();
    Cy(i) = custom_rng();
    Cz(i) = custom_rng();
end;
figure
scatter3(Cx, Cy, Cz)
title('Ģenerēšanas rezultāti 3 dimensijās')

%% d
clear all, clc
c = 2^32;
P_length = 2000;
P = zeros(P_length, 1);
match_start = 0;
match_offset = 0;
match_count = 0;

for i = 1:(c*2)
    rand = custom_rng();

    if (i < P_length)
        P(i) = rand;
    end

    if (match_start > 0)
        k = i - match_start + match_offset;
        if (P(k) == rand)
            match_count = match_count + 1;
            if (match_count > 100)
                disp(['100 secigi atkārtojumi'])
                disp(['sakritība sākot ar elementu: ' num2str(match_start)])
                disp(['aperiods: ' num2str(match_offset)])
                disp(['periods: ' num2str(match_start - match_offset)])
                return
            end
        end
    end

    if (i > P_length & match_start == 0)
        for j = 1:length(P)
            if (rand == P(j))
                match_start = i;
                match_offset = j;
                break;
            end
        end
    end
end


%% periods
clc
asda = 0;
for j = 1:(length(P) - c)
    for i = (j + 1):length(P)
        if (P(i) == P(j)) 
            disp(['Atkartota ' num2str(j) '. vertiba  @ index:  ' num2str(i) ' periods ' num2str(i - j)])
            asda = asda + 1; 
            break;
        end;

        if (asda > 100)
            break
        end
    end
end
%%
for j = 1:(length(P) - c)
    for i = c:length(P)
        if (P(i) == P(j)) 
            disp(['Atkartota ' num2str(j) ' vertiba @ index:  ' num2str(i)])
        end;
    end
end
% 
% min_guess = 5001;
% max_guess = 10000;
% 
% while (min_guess < max_guess)
% 
%     current_guess = ceil((max_guess + min_guess) / 2);
%     start_idx = current_guess;
%     disp(['Guess for aperiod: ' num2str(current_guess)])
% 
%     repeat_start = 0;
%     repeat_length = 0;
%     tic
%     for i = 1:(length * 2)
%         rand = custom_rng();
% 
%         if (i > start_idx & i <= end_idx)
%             P(i - start_idx) = rand;
%         end;
% 
%         if (i > end_idx)
%             if repeat_length == 0            
%                 if rand == P(repeat_length + 1)
%                     % atkārtojums sākas
%                     repeat_start = i;
%                     repeat_length = 1;
%                 end
%             else
%                 if repeat_length >= length
%                     disp(['max length repeat ' num2str(repeat_start)])
%                     break
%                 end
% 
%                 if rand == P(repeat_length + 1)
%                     % atkārtojums turpinās
%                     repeat_length = repeat_length + 1;
%                 else
%                     % atkārtojums beidzas
%                     repeat_length = 0;
%                     repeat_start = 0;
%                 end
%             end;
%         end;
%     end;
%     % repeat_start
%     % repeat_length
% 
%     if (repeat_length == length)
%         % atrada atkārtojumu
%         max_guess = current_guess;
%     else
%         min_guess = current_guess;
%     end;
% 
%     disp([ 'min guess: ' num2str(min_guess) ' max guess ' num2str(max_guess) ])
% 
%     toc
% end
% disp([ 'aperiod length : ' num2str(min_guess) ])


% 536 870 912 - 2^29
% 536 880 913 - perioda garums m + 1
% start idx = 10000
% max length repeat 536880913
% 
% repeat_start =
% 
%    536880913
% 
% 
% repeat_length =
% 
%    536870912
% 
% Elapsed time is 203.085706 seconds.
% 5000
% max length repeat 536875913
% 
% repeat_start =
% 
%    536875913
% 
% 
% repeat_length =
% 
%    536870912
% 
% Elapsed time is 212.420947 seconds.
% 2500
% repeat_start =
% 
%    536873413
% 
% 
% repeat_length =
% 
%    536868412
% 
% Elapsed time is 203.985967 seconds.
% 1250
% repeat_start =
% 
%    536872163
% 
% 
% repeat_length =
% 
%    536869662
% 
% Elapsed time is 212.113977 seconds.
% 625
% repeat_start =
% 
%    536871538
% 
% 
% repeat_length =
% 
%    536870287
% 
% Elapsed time is 202.347594 seconds.
% 310
% repeat_start =
% 
%    536871538
% 
% 
% repeat_length =
% 
%    536870287
% 
% Elapsed time is 202.347594 seconds.
%% 2. normālā sadalījuma generēšana

%% 3. Call option cenas modelis
2^29
num2str(2^32)
num2str(rem(5 ^ 96810627, 2^32))
% 536880913 - 536870912

%%
function result = custom_rng

    persistent current_x
    if isempty(current_x)
       current_x = 1;
    end

    max_value = 4294967296; % 2^32
    % 536870912; %
    a = 1732073221; % no lekcijas meteriāla, priekš 32-bit skaitļa

    current_x = rem(a * current_x, max_value);
    result = current_x / max_value;  % [0,1)
end
