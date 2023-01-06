clc;
clf;
clear;
close all; 
%.......... test signal...............
 % Sine Wave generation
N=128;
T1=5;
 n_n=0:1:N-1;
 x1=sin(2*pi*n_n/T1);
 subplot(4,1,1);
 plot(n_n,x1);
 
 % Impulse Train generation
 T2=7;
n_n=0:1:N-1;
x2(1)=0;
for i=2:N
    if rem(i,T2)==0
        x2(i)=1;
    else
        x2(i)=0;
    end
end
 subplot(4,1,2);
stem(n_n,x2);

% addition of white gaussian noise
rng(0);
x3=rand(1,N);
subplot(4,1,3);
plot(n_n,x3);

x = (x1 + x2 +x3)-mean(x1+x2+x3) ; %substract dc component

xl=length(x);
subplot(4,1,4);
plot(n_n,x);


% ...... Calculation of RAMANUJAN SUMS COEFFICIENT ......
% using Euler totient function & mobieus function ......
% create user defined eulertotient & mobieus function .........


for q = 1:128
    for n_n = 1:128
        d1(q,n_n) = gcd(q,n_n);                            % ....find greatest common divisor of(q,n).. 
         mob(q,n_n) = mobiusfunction(q/d1(q,n_n));         %... call mobiusfunction & compute µ (q/((q,n))..
    eul_num(q) = eulertotient(q);                      % ... call eulertotient function & compute ø(q).....
    eul_den(q,n_n) = eulertotient(q/d1(q,n_n));            %... call eulertotient function &compute ø (q/((q,n))) .....
    cq(q,n_n) = (mob(q,n_n) * eul_num(q))./eul_den(q,n_n);   % find ramanujan sums coefficient ........ ( µ (q/((q,n) * ø(q)) / ø (q/((q,n)))
    end   
end

% ............... RFT ...................

for q=1:N
    s1=zeros(1,1);   
for n_n=1:N
       s1=s1+(x(n_n).*cq(q,n_n)); %S(j,ii));
end
       RFT1(q)=(1/(N*eulertotient(q)))*s1;
end
RFT = abs(RFT1);
figure;
stem(1:128,RFT);
axis([0 20 0 0.4]);
title('RFT'); xlabel('period length[aa]');ylabel('|R|');