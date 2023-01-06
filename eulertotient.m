function eutot = eulertotient(n)
%n=input('the value of eulertotient function');
 eutot = 0;
 for i=1:n
     if gcd(i,n) == 1;
         eutot = eutot + 1;
     end
 end

 %disp('eulertotient=');
  %  disp(eutot)
end