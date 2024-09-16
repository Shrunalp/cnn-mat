function [C, m1, m2] = clebsch_gordan(j1, j2, j, m)
% clebsch_gordan returns the Clebsch-Gordan coefficients for the specified
% total angular momenta and coupled z angular momentum component. This form
% is particularly convenient from a computational standpoint. Follows the
% approach of W. Straub in viXra:1403.0263.
% 
% Inputs:
%   j1 - first uncoupled total angular momentum.
%   j2 - second uncoupled total angular momentum.
%   j  - coupled total angular momentum.
%   m  - coupled z angular momentum component.
%
% Outputs:
%   C  - all of the nonzero Clebsch-Gordan coefficients for the specified
%        inputs. Ordered in increasing values of m1.
%   m1 - first uncoupled z angular momentum components.
%   m2 - second uncoupled z angular momentum components.
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    assert(isscalar(j) && isscalar(m) && isscalar(j1) && isscalar(j2), 'All inputs must be scalars.')

    if j1 < 0 || mod(j1, 0.5) ~= 0 || ...
       j2 < 0 || mod(j2, 0.5) ~= 0 || ...
       j  < 0 || mod(j,  0.5) ~= 0 || ...
       mod(j + m, 1) ~= 0
        % Input doesn't make sense
        error(sprintf('You probably made a mistake. \n 1. j1, j2, j should be non-negative integers or half-integers. \n 2. j + m should be an integer.')); %#ok<SPERR>
    end
    
    if j < abs(j1 - j2) || j > j1 + j2 || abs(m) > j
        % Nothing to do
        m1 = [];
        m2 = [];
        C  = [];
        return;
    end

    m11 = (m - j1 - j2 + abs(j1 - j2 + m)) / 2;
    m1n = (m + j1 + j2 - abs(j1 - j2 - m)) / 2;
    
    m1 = (m11:m1n)';
    m2 = m - m1;
    
    j_const = j1 * (j1 + 1) + j2 * (j2 + 1) - j * (j + 1);
    
    n = m1n - m11 + 1;
    % A is a tridiagonal symmetric matrix
    A = zeros(n, n);
    for a = 1:n
        A(a, a) = j_const + 2 * m1(a) * m2(a);
    end
    for a = 1:(n - 1)
        tmp = realsqrt(j1 * (j1 + 1) - m1(a) * m1(a + 1)) * realsqrt(j2 * (j2 + 1) - m2(a) * m2(a + 1));
        A(a, a + 1) = tmp;
        A(a + 1, a) = tmp;
    end
    
    % A determines C up to sign and normalization
    C = null(A);
    C = sign(C(n)) * C / sqrt(C' * C);
end
