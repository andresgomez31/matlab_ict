function yapprox = RK2 (m, x0, xf, y0, step, a2)
    x = x0:step:xf;
    y = zeros(size(x));
    y(1) = y0;
    
    if a2 ~= 0
        q = 1/(2*a2);
    else
        q = 0;
    end
     
    a1 = 1 - a2;

    for i = 1:length(x) - 1
        k1 = m(x(i), y(i));
        k2 = m(x(i) + q*step, y(i) + q*step*k1);

        y(i + 1) = y(i) + (a1*k1 + a2*k2)*step;
    end

    yapprox = y;
end