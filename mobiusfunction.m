function mobius = mobiusfunction(n)
%n=input('the value of mobeius function');
if n==1
    mobius=1;
else
    f=factor(n);
    l=length(f);
    if l == 1
        mobius = -1;
        %disp('mobeius=');
    %disp(mobius)
        %return
    else
    for i=1:l-1;
        if f(1)==f(i+1)
            mobius=0;
            break
                else
            mobius=(-1)^l;
        end
    end
    end
end
    %disp('mobeius=');
    %disp(mobius)
end
            
        