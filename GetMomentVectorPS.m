function [v] = GetMomentVectorPS(seq)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n=length(seq);

%Binary indicator sequences of A, C, G, T (each value will be either 1 or 0
%; 1 if that nucleotide appears in that position and 0 otherwise)
uA=zeros(1, n);
uC=zeros(1, n);
uG=zeros(1, n);
uT=zeros(1, n);

%Sequences' length (of A, C, G, T respectively)
nA=0;
nC=0;
nG=0;
nT=0;

%Get the binary indicator sequences and their lengths
for i=1:n
    nu=seq(i);
   switch upper(nu)
      case {'A'}
           uA(i)=1;
           uC(i)=0;
           uG(i)=0;
           uT(i)=0;
           nA=nA+1;
      case {'C'}
           uC(i)=1;
           uA(i)=0;
           uG(i)=0;
           uT(i)=0;
           nC=nC+1;     
      case {'G'}
           uG(i)=1;
           uA(i)=0;
           uC(i)=0;
           uT(i)=0;
           nG=nG+1;
      case {'T'}
           uT(i)=1;
           uA(i)=0;
           uC(i)=0;
           uG(i)=0;
           nT=nT+1;
   end
end

%Ramanujan Fourier Transforms

%.......... test signal...............
 % Sine Wave generation
N=n;
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


for q = 1:n
    for n_n = 1:n
        d1(q,n_n) = gcd(q,n_n);                            % ....find greatest common divisor of(q,n).. 
         mob(q,n_n) = mobiusfunction(q/d1(q,n_n));         %... call mobiusfunction & compute µ (q/((q,n))..
    eul_num(q) = eulertotient(q);                      % ... call eulertotient function & compute ø(q).....
    eul_den(q,n_n) = eulertotient(q/d1(q,n_n));            %... call eulertotient function &compute ø (q/((q,n))) .....
    cq(q,n_n) = (mob(q,n_n) * eul_num(q))./eul_den(q,n_n);   % find ramanujan sums coefficient ........ ( µ (q/((q,n) * ø(q)) / ø (q/((q,n)))
    end   
end

% ............... RFT ...................

for q=1:uA
    s1=zeros(1,1);   
for n_n=1:uA
       s1=s1+(x(n_n).*cq(q,n_n)); %S(j,ii));
end
       RFT1(q)=(1/(uA*eulertotient(q)))*s1;
end
RFT1 = UT;
 

%
for q=1:uC
    s1=zeros(1,1);   
for n_n=1:uC
       s1=s1+(x(n_n).*cq(q,n_n)); %S(j,ii));
end
       RFT2(q)=(1/(uC*eulertotient(q)))*s1;
end
RFT2 = UA;

%
for q=1:uG
    s1=zeros(1,1);   
for n_n=1:uG
       s1=s1+(x(n_n).*cq(q,n_n)); %S(j,ii));
end
       RFT3(q)=(1/(uG*eulertotient(q)))*s1;
end
RFT3 = UG;

%
for q=1:uT
    s1=zeros(1,1);   
for n_n=1:uT
       s1=s1+(x(n_n).*cq(q,n_n)); %S(j,ii));
end
       RFT4(q)=(1/(uT*eulertotient(q)))*s1;
end
RFT4 = UT;





%Exclude the first term
%Discrete Fourier Transforms


%Exclude the first term
UA(1)=[];
UC(1)=[];
UG(1)=[];
UT(1)=[];

% Get the first half of the transform (since it's symmetric)
m=floor((n-1)/2);
UA1=UA(1:m);
UC1=UC(1:m);
UG1=UG(1:m);
UT1=UT(1:m);


%Power spectrums
PSA=abs(UA).^2;     
PSC=abs(UC).^2;     
PSG=abs(UG).^2;     
PST=abs(UT).^2;     


%Normalized moments
MA=zeros(1,3);   
MC=zeros(1,3);
MG=zeros(1,3);
MT=zeros(1,3);

%Moment vector
for j=1:3
   MA(j)=(nA*(n-nA))*sum(PSA.^j)/(nA^(j)*(n-nA)^(j)); 
   MC(j)=(nC*(n-nC))*sum(PSC.^j)/(nC^(j)*(n-nC)^(j)); 
   MG(j)=(nG*(n-nG))*sum(PSG.^j)/(nG^(j)*(n-nG)^(j)); 
   MT(j)=(nT*(n-nT))*sum(PST.^j)/(nT^(j)*(n-nT)^(j)); 
end


v=[MA(1), MC(1), MG(1), MT(1), MA(2), MC(2), MG(2), MT(2), MA(3), MC(3), MG(3), MT(3)];


end

