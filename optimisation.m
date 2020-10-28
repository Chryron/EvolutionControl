D = 0:0.01:1;
mu = 5;
f = (2-exp(D)).*sqrt(mu.*D).*(1-exp(-D));
plot(D,f)
xlabel("D"); ylabel("f(D)");
