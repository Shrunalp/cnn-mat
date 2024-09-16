function [N] = sp_null(A, opt, TOL)
% sp_null finds the left and right null spaces of the sparse matrix A. The 
% algorithm is not well-established in the literature, but performs the
% decomposition A = L U Q where L and Q are invertible matrices, U and A are
% the same size, and U has the block structure U = [U_{11}, 0; 0, 0]. If U is
% rank r, then since A Q^{-1} = L U the n - r columns of Q^{-1} corresponding
% to pivotless columns of U are a basis for the right null space. Similarly,
% since L^{-1} A = U Q the m - r rows of L^{-1} corresponding to pivotless
% rows of U are a basis for the left null space.
% 
% Inputs:
%   A   - a sparse matrix with m rows and n columns.
%   opt - 0 or 1 specifies that the left or right null space is desired.
%   TOL - threshold for a pivot to be considered significant.
%
% Outputs:
%   N   - the left or right null space of the matrix A depending on the
%         value of opt. The sparsity is increased at the expense of the columns
%         not being orthonormal.
% 
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    switch opt
        case 0
            I = find(A);
            A(I(abs(A(I)) < TOL)) = 0.;
            
            [L, U, ~] = luq(A, TOL);
            r = nnz(abs(spdiags(U, 0)) > TOL);
            inv_L = inv(L);
            N = inv_L((r + 1):size(A, 1), :);
            
            I = find(N);
            N(I(abs(N(I)) < TOL)) = 0.;
        case 1
            N = sp_null(A.', 0, TOL).';
        otherwise
            disp('Invalid option. Consult documentation.');
    end
end

function [L, U, Q] = luq(A, TOL)
% luq is the workhorse of sp_null, and actually performs the decomposition
% A = L * U * Q. This is constructed using the A = P' * L * U decomposition
% given by the lu function in MATLAB.
% 
% Inputs:
%   A   - a sparse matrix with m rows and n columns.
%   TOL - threshold for a pivot to be considered significant.
%
% Outputs:
%   L   - an invertible matrix with the same number of rows as A.
%   U   - a matrix with the block structure U = [U_{11}, 0; 0, 0] where U_{11}
%         is a square matrix indicating the rank of A.
%   Q   - an invertible matrix with the same number of columns as A.
% 
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    if ~issparse(A)
        A = sparse(A);
    end
            
    % This should be as stable as possible
    [L, U, P] = lu(A, 1.);
    
    % Convert decomposition to A = L * U * Q where L is m x m, U is m x n,
    % and Q is n x n.
    [m, n] = size(A);
    if n < m
        p = m - n;
        L = [L, [sparse(n, p); speye(p)]];
        U = [U; sparse(p, n)];
    end
    L = P' * L;
    Q = speye(n);
    
    % Find the row and col indices of significant pivots
    p1 = find(abs(spdiags(U, 0)) > TOL)';
    p = numel(p1);
    
    % Permute rows and columns to give U = [U_{11}, U_{12}; U_{21}, U_{22}]
    % with all of the significant pivots in U_{11}. 
    r2 = 1:m;
    r2(p1) = [];
    r_perm = [p1, r2];
    L = L(:, r_perm);
    U = U(r_perm, :);
        
    c2 = 1:n;
    c2(p1) = [];
    c_perm = [p1, c2];
    U = U(:, c_perm);
    Q = Q(c_perm, :);
    
    % Update indices to reflect reordering
    p1 = 1:p;
    r2 = (p + 1):m;
    c2 = (p + 1):n;
    
    % The U_{21} block is converted to zeros by U = [I, 0; -X, I] * U and 
    % L = L * [I, 0; X, I], preserving the product A = L * U * Q.
    X = U(r2, p1) / U(p1, p1);
    L(:, p1) = L(:, p1) + L(:, r2) * X;
    U(r2, p1) = 0.;
    U(r2, c2) = U(r2, c2) - X * U(p1, c2);
    
    % The U_{1} block is converted to zeros by U = U * [I, -X; 0, I] and 
    % Q = [I, X; 0, I] * Q, preserving the product A = L * U * Q. U_{22} is
    % unchanged as a consequence of the previous step.
    X = U(p1, p1) \ U(p1, c2);
    U(p1, c2) = 0.;
    Q(p1, :) = Q(p1, :) + X * Q(c2, :);
    
    % At this point U = [U_{11}, 0; 0, U_{22}] but U_{22} could still have
    % nonzero values off of the diagonal. Permute rows and columns to
    % change block structure to U = [U_{11}, 0, 0; 0, U_{22}, 0; 0, 0, 0].
    r3 = p + find(max(abs(U(r2, c2)), [], 2) > TOL)';
    c3 = p + find(max(abs(U(r2, c2)), [], 1) > TOL);
    if ~isempty(r3) || ~isempty(c3)
        r4 = r2;
        r4(r3 - p) = [];
        r_perm = [r3, r4];
        L(:, (p + 1):m) = L(:, r_perm);
        U((p + 1):m, :) = U(r_perm, :);

        c4 = c2;
        c4(c3 - p) = [];
        c_perm = [c3, c4];
        U(:, (p + 1):n) = U(:, c_perm);
        Q((p + 1):n, :) = Q(c_perm, :);
        
        % Update indices to reflect reordering
        r3 = p + (1:numel(r3));
        c3 = p + (1:numel(c3));
        
        % Perform decomposition of U_{22} = L_r * U_r * Q_r, resulting in
        % A = [L_1, L_2 * L_r, L_3] * ...
        %     [U_{11}, 0, 0; 0, U_r, 0; 0, 0, 0] * ...
        %     [Q_1; Q_r * Q_2; Q_3]
        [Lr, Ur, Qr] = luq(U(r3, c3), TOL);
        U(r3, c3) = Ur;
        L(:, r3) = L(:, r3) * Lr;
        Q(c3, :) = Qr * Q(c3, :);
        
        % Find the number of significant pivots
        q = nnz(abs(spdiags(U, 0)) > TOL);
        U = [U(1:q, 1:q), sparse(q, n - q); sparse(m - q, q), sparse(m - q, n - q)];
    else
        U(r2, c2) = 0.;
    end
end
