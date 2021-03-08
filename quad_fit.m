function val = quad_fit(samples, k)
    s0 = (-3*samples(1)+12*samples(2)+17*samples(3)+12*samples(4)-3*samples(5))/35;
    s1 = (-2*samples(1)-samples(2)+samples(4)+2*samples(5))/10;
    s2 = (2*samples(1)-samples(2)-2*samples(3)-samples(4)+2*samples(5))/14;
    val = s0 + k*s1 + s2*k^2;
end