function [X, L, M] = make_symm_pq(n, TOL)
% make_symm_pq finds sets of coefficients that make the hyperspherical
% harmonic expansion invariant to the specified point symmetry group
% generators. These can then be interpreted as defining a compact basis for an
% orientation distribution obeying the required symmetry, namely, the
% symmetrized hyperspherical harmonics.
% 
% The difference with make_symm is in the way the eigenspaces are combined. The
% more formal procedure followed here constructs projection matrices and finds
% the subspace belonging to the intersection. This is slower, but has the
% advantage of being supported by well-established theorems.
% 
% Inputs:
%   n   - upper index of the hyperspherical harmonics. The calculation is
%         formally independent for each value of n.
%   lx  - rotation angle and spherical angles of the rotation axis for one of 
%         the generators left multiplying the argument of the hyperspherical
%         harmonics. That is, a symmetry of the sample. The angles [w, th, ph]
%         should be in the usual intervals. x can be a, b, c or d.
%   rx  - rotation angle and spherical angles of the rotation axis for one of 
%         the generators right multiplying the argument of the hyperspherical
%         harmonics. That is, a symmetry of the crystal. The angles [w, th, ph]
%         should be in the usual intervals. x can be a, b, c or d.
%   TOL - tolerance for the desired eigenvalue of a matrix. For all practical
%         cases the desired eigenvalues are well-separated from the rest of the
%         spectrum, so this does not need to be excessively small.
%
% Outputs:
%   X   - matrix with orthonormal columns that forms a basis for the space of 
%         expansion coefficients satisfying the requested symmetries. That is,
%         a column gives the expansion coefficients required to define one of
%         the symmetrized hyperspherical harmonics. The expansion uses complex
%         hyperspherical harmonics and orders the coefficients
%         lexicographically in (l, m).
%   L   - values of l for the rows of X.
%   M   - values of m for the rows of X.
% 
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.

%#ok<*SPRIX>
%#ok<*UNRCH>

    % TOL = 1e-12;

    % n = 8;
    l = n / 2;
        
    % Symmetry generators. The crystal and sample point group symmetries do not
    % require more than two generators each.
    la = [0., 0., 0.];
    ra = [pi / 2, 0., 0.];
    
    lb = [0., 0., 0.];
    rb = [pi / 2, pi / 2, 0.];
    
    lc = [0., 0., 0.];
    rc = [0., 0., 0.];
    
    ld = [0., 0., 0.];
    rd = [0., 0., 0.];
    
    % The equation for the irreps of SO(4) below allows the Clebsch-Gordan
    % coefficients to be collected into a unitary matrix. Performing a similarity
    % transformation with this matrix has the effect of converting from the 
    % uncoupled to the coupled basis.
    CG = sparse((n + 1)^2, (n + 1)^2);
    j = 0:n;
    for a = 1:(n + 1)
        m = -j(a):j(a);
        for b = 1:(2 * j(a) + 1)
            row = (j(a) + 1)^2 - j(a) + m(b);
            [C, m1, m2] = clebsch_gordan(l, l, j(a), m(b));
            for c = 1:length(C)
                col = (m1(c) + l) * (2 * l + 1) + (m2(c) + l) + 1;
                CG(row, col) = C(c);
            end
        end
    end
    
    CR = complex_real(n);

    X = speye((n + 1)^2);
    P = X * X';
    
    % Construct simultaneous eigenvectors of eigenvalue one of the irreps of SO(4).
    % Observe that the eigenspace of eigenvalue one of a matrix A is equivalent to
    % the right nullspace of the matrix (A - I). The irrep of SO(4) is constructed
    % using the following equation (a transform of Eq. 2 of the manuscript);
    % 
    % R^{a a}_{c \gamma d \delta}(g_l, g_r) = \sum_{\alpha' \beta' \alpha \beta}
    % C^{c \gamma}_{a \alpha' a \beta'} U^a_{\alpha' \alpha}(g_r^{-1})
    % U^a_{\beta' \beta}(g_l) C^{d \delta}_{a \alpha b \beta}
    % 
    % That is, the desired irrep of SO(4) is given by a similarity transformation
    % of the Kroneker product of irreps of SO(3). The effect of right multiplying
    % a row vector of hyperspherical harmonics with this irrep is to carry every
    % rotation g \rightarrow g_l g g_r. Since the eigenspectrum of a matrix is
    % unchanged by a unitary transformation, the calculation can be performed in
    % the uncoupled basis and transformed to the coupled basis afterwards.
    Ul = rotation_mat(l,  la(1), la(2), la(3));
    Ur = rotation_mat(l, -ra(1), ra(2), ra(3));
    Ra = kron(Ur, Ul);
    Y = sp_orth(sp_null(Ra - eye(size(Ra)), 1, TOL));
    Q = Y * Y';
    X = sp_orth(sp_null(P * Q - speye(size(Q)), 1, TOL));
    P = X * X';
    
    Ul = rotation_mat(l,  lb(1), lb(2), lb(3));
    Ur = rotation_mat(l, -rb(1), rb(2), rb(3));
    Rb = kron(Ur, Ul);
    Y = sp_orth(sp_null(Rb - eye(size(Rb)), 1, TOL));
    Q = Y * Y';
    X = sp_orth(sp_null(P * Q - speye(size(Q)), 1, TOL));
    P = X * X';
    
    Ul = rotation_mat(l,  lc(1), lc(2), lc(3));
    Ur = rotation_mat(l, -rc(1), rc(2), rc(3));
    Rc = kron(Ur, Ul);
    Y = sp_orth(sp_null(Rc - eye(size(Rc)), 1, TOL));
    Q = Y * Y';
    X = sp_orth(sp_null(P * Q - speye(size(Q)), 1, TOL));
    P = X * X';
    
    Ul = rotation_mat(l,  ld(1), ld(2), ld(3));
    Ur = rotation_mat(l, -rd(1), rd(2), rd(3));
    Rd = kron(Ur, Ul);
    Y = sp_orth(sp_null(Rd - eye(size(Rd)), 1, TOL));
    Q = Y * Y';
    X = sp_null(P * Q - speye(size(Q)), 1, TOL);
    
    % Convert to the coupled basis
    X = CG * X;
    X = full(CR * clean(CR' * X, TOL));
    
    err_a = max(max(abs((CG * Ra * CG' - eye(size(Ra))) * X)));
    err_b = max(max(abs((CG * Rb * CG' - eye(size(Rb))) * X)));
    err_c = max(max(abs((CG * Rc * CG' - eye(size(Rc))) * X)));
    err_d = max(max(abs((CG * Rd * CG' - eye(size(Rd))) * X)));
    
    if ~isempty(X)
        disp(['Numerical error: ', num2str(max([err_a, err_b, err_c, err_d]))]);
    else
        disp('No symmetrized harmonics found.');
    end
    
    % Assign indices
    L = zeros((n + 1)^2, 1);
    M = zeros((n + 1)^2, 1);
    for l = 0:n
        for m = -l:l
            ind = (l + 1)^2 - l + m;
            L(ind) = l;
            M(ind) = m;
        end
    end
end

function [Q] = sp_orth(A)
% sp_orth finds an orthogonal basis for the column space of A, but is 
% actually just a wrapper for qr.
% 
% Inputs:
%   A - matrix specifying the column space.
%
% Outputs:
%   Q - matrix whose columns form an orthonormal column space for A.
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    [Q, ~] = qr(A, 0);
end
