function [symm] = collect_symm(n_dim)
% collect_symm finds at least n_dim sets of coefficients that make the
% hyperspherical harmonic expansion invariant to the cubic crystal point
% group generators. generators. These can  be interpreted as defining a
% compact basis for an orientation distribution obeying the required
% symmetry, namely, the symmetrized hyperspherical harmonics.
% 
% Inputs:
%   n_dim - number of dimensions required for the spectral embedding.
% 
% Outputs:
%   summ  - a cell array containing three columns containing the outputs X,
%           L and M of the make_symm_pq function. The value of N starts at
%           zero and increases with the row.
% 
% Copyright 2023 Jeremy Mason
% 
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    n_found = 0;

    n = 0;
    symm = cell(1, 3);
    while n_found < n_dim
        [X, L, M] = make_symm_pq(n, 1.e-12);
        symm{n + 1, 1} = X;
        symm{n + 1, 2} = L;
        symm{n + 1, 3} = M;

        n_found = n_found + size(X, 2);
        n = n + 2;
    end
end