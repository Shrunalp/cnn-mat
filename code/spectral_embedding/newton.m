function [x] = newton(y)
    MAXIT = 100;
    SMALL1 = 1.e-14;
    SMALL2 = 1.e-12;
    
    % Initial search interval
    xl = y - pi / 4.;
    xh = y + pi / 4.;
    
    % Initial guess
    x = y;
    
    % Initialize step lengths
    f  = x - sin(2. * x) / 2. - y;
    df = max(1. - cos(2. * x), SMALL1);
    
    dx_old = pi / 4.;
    dx = -f / df;
    
    for a = 1:MAXIT      
        if (xh - x) < dx || dx < (xl - x) || abs(dx) > abs(dx_old / 2.)
            dx = (xh - xl) / 2.;
            x = xl;
        end
        x = x + dx;
        if abs(dx) < SMALL2
            return;
        end
        
        f  = x - sin(2. * x) / 2. - y;
        df = max(1. - cos(2. * x), SMALL1);
        
        dx_old = dx;
        dx = -f / df;
        
        if dx > 0.
            xl = x;
        else
            xh = x;
        end
    end
    warning('newton: failed to converge');
end