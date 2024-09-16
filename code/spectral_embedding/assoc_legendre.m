function [Plm] = assoc_legendre(l, m, x)
% assoc_legendre calculates the associated Legendre function of degree l and 
% order m. More efficient than legendre when a single associated Legendre is 
% required. Assumes that all arguments are single values.
% 
% Inputs:
%   l - degree of the associated Legendre function.
%   m - order of the associated Legendre function.
%   x - argument of the associated Legendre function.
%
% Outputs:
%   Plm - value of the specified associated Legendre function. Follows the 
%         Condon-Shortley phase convention.
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.

    % Start with explicit formula for P^m_m(x) and use recursion to increase
    % the degree. Recursion relation on degree is known to be stable in this
    % direction.
    if abs(m) < eps
        Pm = 1.;
    else
        Pm  = (-1)^m * prod(m:(2 * m - 1)) / 2.^(m - 1) * (1. - x * x)^(m / 2.);
    end
    if abs(l - m) < eps
        Plm = Pm;
        return;
    end
    
    % Given current degree l1, the three functions appearing in the recursion
    % relation are Pm = P^m_{l1 - 1}, Po = P^m_l1, and Pp = P^m_{l1 + 1}.
    Po = x * (2 * m + 1) * Pm;
    l1 = m + 1;
    while abs(l - l1) > eps
        Pp = ((2 * l1 + 1) * x * Po - (l1 + m) * Pm) / (l1 - m + 1);
        Pm = Po;
        Po = Pp;
        l1 = l1 + 1;
    end
    
    Plm = Po;
end
