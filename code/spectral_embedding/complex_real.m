function [U] = complex_real(n)
% complex_real constructs a sparse matrix that performs the similarity
% transformation to convert from the complex harmonics to the real
% harmonics. Rows ordered by increasing l and by increasing m. Columns
% ordered by increasing l, by increasing m, and by c before s.
% 
% Inputs:
%   n - maximum value of l.
% 
% Outputs:
%   U - sparse matrix that performs the transformation. 
% 
% Copyright 2023 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.

sqrt2 = realsqrt(2.);

U = sparse((n + 1)^2, (n + 1)^2);
for L = 0:n
    U((L + 1)^2 - L, (L + 1)^2 - 2 * L) = (1i)^L;
    for M = 1:L
        U((L + 1)^2 - L - M, (L + 1)^2 - 2 * L + 2 * M - 1) = (-1)^M * (1i)^L / sqrt2;
        U((L + 1)^2 - L + M, (L + 1)^2 - 2 * L + 2 * M - 1) = (1i)^L / sqrt2;
        U((L + 1)^2 - L - M, (L + 1)^2 - 2 * L + 2 * M) = (-1)^(M - 1) * (1i)^(L - 1) / sqrt2;
        U((L + 1)^2 - L + M, (L + 1)^2 - 2 * L + 2 * M) = (1i)^(L - 1) / sqrt2;
    end
end
