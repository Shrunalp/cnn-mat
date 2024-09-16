function [Z] = get_hsh(n, l, m, w, th, ph)
% get_hsh calculates the value of the hyperspherical harmonic Z^n_{lm} for the
% specified rotation angle and spherical angles of the rotation axis. Assumes 
% that all arguments are single values.
% 
% Inputs:
%   n  - upper index of the hyperspherical harmonic.
%   l  - degree of the spherical harmonic appearing in the definition.
%   m  - order of the spherical harmonic appearing in the definition.
%   w  - rotation angle in the interval [0, 2 \pi].
%   th - polar angle of rotation axis in the interval [0, \pi].
%   ph - aximuthal angle of rotation axis in the interval [0, 2 \pi].
%
% Outputs:
%   Z  - value of the specified hyperspherical harmonic. Follows the
%        conventions established in J. K. Mason, Acta Cryst A 65, 259 (2009).
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.

    n_2 = round(n / 2);
    abs_m = abs(m);

    % The expected usage pattern involves repeated function calls for the same
    % values of n. Clebsch-Gordan coefficients should be memoized. The usage
    % pattern below only requires that C^{l, 0}_{n / 2, m, n / 2, -m}.
    persistent C;
    if isempty(C)
        C = cell(0);
    end
    if length(C) < n_2 + 1 || isempty(C{n_2 + 1})
        C{n_2 + 1} = zeros(n + 1, n + 1);
        for a = 1:(n + 1)
            C{n_2 + 1}(a, :) = clebsch_gordan(n_2, n_2, a - 1, 0)';
        end
    end
    
    % Assuming that this is called to construct the symmetrized versions, there
    % is no consistent pattern in the values of l or m. assoc_legendre is much
    % more efficient than returning the associated Legendre functions for all
    % values of m.
    p = assoc_legendre(l, abs_m, cos(th));
    p = p / realsqrt(prod((l - abs_m + 1):(l + abs_m)));
    if m < 0
        p = (-1)^abs_m * p;
    end
    Y = realsqrt((2 * l + 1) / (4 * pi)) * p * exp(1i * m * ph);
    
    % With the spherical harmonic is calculated, the only remaining part is to
    % find the generalized character of order l for the irrep of rank n / 2.
    % The relevant formula for the hyperspherical harmonics is implicit in Eq.
    % 8 on page 81 of D. A. Varshalovich et al, Quantum Theory of Angular
    % Momentum, 1988. The generalized character is calculated by the 
    % trigonometric expansion in Eq. 2 on page 106 of the same reference, with
    % the Clebsh-Gordan coefficient transformed into the memoized form.
    C1 = C{n_2 + 1}(l + 1, :);
    m1 = -n_2:n_2;
    Z = sum((-1).^m1 .* C1 .* exp(-1i * m1 * w));
    Z = (-1)^(n_2) * realsqrt((2. * (n + 1)) / (pi * (2 * l + 1))) * Z * Y;
end
