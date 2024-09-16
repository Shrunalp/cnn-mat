function [d, err] = wigner_little_d(j, theta)
% wigner_little_d constructs a Wigner little d matrix given a total angular 
% momenum and an angle. This corresponds to the irrep of SO(3) for a rotation
% about the y axis.
% 
% Inputs:
%   j     - specifies dimension of the matrix.
%   theta - rotation angle in the interval [0, 2 \pi].
%
% Outputs:
%   d     - Wigner little d matrix, with the rows and columns ordered in
%           increasing values of m' and m. Follows the approach of X. M. Feng
%           et al in 10.1103/PhysRevE.92.043307.
%   err   - simple estimate of the numerical error in the calculation of the
%           entries. Usually not calculated or returned.
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    m = -j:j;
    n = 2 * j + 1;
    
    X = sqrt((j + m) .* (j - m + 1)) / (2 * 1i);
    % Jy is a tridiagonal Hermitian matrix
    Jy = zeros(n, n);
    for a = 2:n
        b = n - a + 1;
        Jy(a - 1, a) = -X(a);
        Jy(b + 1, b) =  X(a);
    end
    % Requires that eigenvectors be ordered with increasing eigenvalues
    [V, ~] = eig(Jy);
    
    W = V;
    for a = 1:n
        W(:, a) = W(:, a) * exp(-1i * m(a) * theta);
    end
    
    d = W * V';
    
    if nargout > 1
        err = max(abs(imag(d(:))));
    end
    d = real(d);
end
