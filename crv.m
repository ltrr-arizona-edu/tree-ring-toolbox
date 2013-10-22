t=[1:100]';
k=.1;
k=k(ones(length(t),1),:);
a=.1;
b=.05;
b=b(ones(length(t),1),:);

u= -a * exp(-b .* t) + k;

plot(t,u)
title('-0.1 * exp (-.05 * t) + .1')
pause

u= a * log(b .*t) +k;
plot(t,u)
title('0.1 *  log(.05 * t) +0.1')
