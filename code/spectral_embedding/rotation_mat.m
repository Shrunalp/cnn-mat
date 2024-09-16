function [U] = rotation_mat(j, w, th, ph)
% rotation_mat converts a rotation angle and spherical angles of the rotation
% axis into the equivalent irrep of SO(3).
% 
% Inputs:
%   j  - specifies dimension of the irrep.
%   w  - rotation angle in the interval [0, 2 \pi].
%   th - polar angle of rotation axis in the interval [0, \pi].
%   ph - aximuthal angle of rotation axis in the interval [0, 2 \pi].
%
% Outputs:
%   U  - (2 j + 1)-dimensional representation of SO(3), using the conventions
%        established in Eq. 6 on page 81 of D. A. Varshalovich et al, Quantum
%        Theory of Angular Momentum, 1988. Rows and columns ordered in
%        increasing values of m' and m.
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    tmp = tan(w / 2.) * cos(th);
    tmp = (1. - 1i * tmp) / realsqrt(1 + tmp^2.);
    r_base = 1i * exp(-1i * ph) * tmp;
    c_base = -1i * exp(1i * ph) * tmp;
    
    % Require w to be in [-\pi, \pi]
    w = mod(w + pi, 2. * pi) - pi;
    xi = 2. * asin(sin(w / 2.) * sin(th));
    U = wigner_little_d(j, xi);
    
    m = -j:j;
    n = 2 * j + 1;
    for a = 1:n
        U(a, :) = U(a, :) * r_base^m(a);
        U(:, a) = U(:, a) * c_base^m(a);
    end
end
