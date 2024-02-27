function [xOut] = valInNewInterval(x, a1,b1,a2,b2)
%Takes a value x in the interval [a1, b1] and returns its corresponding
%point in the interval [a2, b2].

    %slope
    m=(b2-a2)./(b1-a1);
    %y-intercept
    c=a2-m.*a1;
    xOut=m.*x+c;

end