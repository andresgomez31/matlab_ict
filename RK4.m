function yapprox = RK4(m, x0, xf, y0, step)
    x = x0:step:xf;
    y = zeros(size(x));
    y(1) = y0;

    for i = 1:length(x)-1
        k1 = m(x(i), y(i));
        k2 = m(x(i) + step/2, y(i) + step*k1/2);
        k3 = m(x(i) + step/2, y(i) + step*k2/2);
        k4 = m(x(i) + step, y(i) + step*k3);

        y(i+1) = y(i) + (step/6)*(k1 + 2*k2 + 2*k3 + k4);
    end

    yapprox = y;
end
