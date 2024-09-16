function [A] = clean(A, TOL)
% clean attempts to construct an orthonormal column space for A with the 
% minimum number of nonzero entries. Assumes that A has full column rank.
% 
% Inputs:
%   A   - matrix specifying the column space.
%   TOL - threshold below which an entry is considered insignificant.
%
% Outputs:
%   A   - matrix whose columns form an orthonormal column space for A.
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    [m, n] = size(A);
    
    % Initial processing
    for a = 1:n
        A(:, a) = A(:, a) / sqrt(A(:, a)' * A(:, a));
        A = clean_col(A, a, TOL);
    end
    
    % Forward sweep
    r1 = 0;
    for a = 1:n
        % Find column containing most significant entry
        p = 0.;
        while p < TOL && r1 < m
            r1 = r1 + 1;
            [p, b] = max(abs(A(r1, a:n)));
        end
        b = (a - 1) + b;
        
        % Should not happen
        if r1 == m
            disp('A is not full column rank');
            break;
        end
        
        % Swap column into leading position
        swap = A(:, a);
        A(:, a) = A(:, b);
        A(:, b) = swap;
        
        % Make leading entry real
        A(:, a) = A(:, a) * (conj(A(r1, a)) / abs(A(r1, a)));

        % Cancel leading entries, breaks orthogonality and normality
        for b = a + find(abs(A(r1, (a + 1):n)) > TOL)
            A(:, b) = A(:, b) - (A(r1, b) / A(r1, a)) * A(:, a);
            A = clean_col(A, b, TOL);
        end
        
    end
    
    % Backward sweep to restore orthonormality
    [Q, ~] = qr(fliplr(A), 0);
    A = fliplr(Q);

    for a = 1:n
        % Strip out small values introduced by orthonormalization
        A = clean_col(A, a, TOL);
        
        % Make leading entry positive
        A(:, a) = sign(A(find(A(:, a), 1), a)) * A(:, a);
    end
end

function [A] = clean_col(A, a, TOL)
% clean_col strips out any real and imaginary components of the entries in
% A(:, a) with magnitudes below the threshold.
% 
% Inputs:
%   A   - matrix containing the column.
%   a   - column index.
%   TOL - threshold below which a component is considered insignificant.
%
% Outputs:
%   A   - matrix containing the column.
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    I = find(A(:, a));
    
    Re = real(A(I, a));
    Re(abs(Re) < TOL) = 0.;
    Im = imag(A(I, a));
    Im(abs(Im) < TOL) = 0.;
    
    A(I, a) = Re + 1i * Im;
end