function [coords] = get_spectral(q, n_dim, symm, DREAM3D_order)
% get_spectral finds the first n_dim coordinates of the quaternion q in the
% spectral embedding.
% 
% Inputs:
%   q      - four components of the normalized quaternion q in the form of a
%            row vector
%   n_dim  - number of coordinates of the quaternion q in the spectral
%            embedding.
%   symm   - a cell array containing three columns containing the outputs X,
%            L and M of the make_symm_pq function. The value of N starts at
%            zero and increases with the row.
%   DREAM3D_order - boolean that indicates the quaternion uses DREAM3D
%            vector-scalar ordering.
% 
% Outputs:
%   coords - first n_dim coordinates of the quaternion q in the spectral
%            embedding in the form of a row vector
% 
% Copyright 2023 Jeremy Mason
% 
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.

    if DREAM3D_order
        % DREAM3D uses vector-scalar ordering
        w = 2. * acos(q(4));
        if w > 10. * eps
            th = acos(q(3) / sin(w / 2.));
            ph = atan2(q(2), q(1));
        else
            th = 0.;
            ph = 0.;
        end
    else
        % normal scalar-vector ordering
        w = 2. * acos(q(1));
        if w > 10. * eps
            th = acos(q(4) / sin(w / 2.));
            ph = atan2(q(3), q(2));
        else
            th = 0.;
            ph = 0.;
        end
    end
    
    coords = zeros(1, n_dim);

    ind = 0;
    % start at n = 1 since n = 0 gives a constant coordinate
    for n = 1:(size(symm, 1) - 1)
        coeffs = symm{n + 1, 1};
        L = symm{n + 1, 2};
        M = symm{n + 1, 3};
        for c = 1:size(coeffs, 2)
            val = 0.;
            for r = 1:size(coeffs, 1)
                if abs(coeffs(r, c)) > eps 
                    val = val + coeffs(r, c) * get_hsh(n, L(r), M(r), w, th, ph);
                end
            end

            ind = ind + 1;
            coords(ind) = val;
            if abs(ind - n_dim) < eps
                break;
            end
        end
        if abs(ind - n_dim) < eps
            break;
        end
    end
    
    if max(abs(imag(coords))) < 1.e-14
        coords = real(coords);
    end
end