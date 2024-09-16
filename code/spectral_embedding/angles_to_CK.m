function [CK] = angles_to_CK(w, th, ph)
% angles_to_CK converts a rotation angle and spherical angles of the rotation
% axis into the equivalent Cayley-Klein matrix.
% 
% Inputs:
%   w  - rotation angle in the interval [0, 2 \pi].
%   th - polar angle of rotation axis in the interval [0, \pi].
%   ph - aximuthal angle of rotation axis in the interval [0, 2 \pi].
%
% Outputs:
%   CK - 2-dimensional representation of SU(2), using the conventions
%        established in Section 14-2 of Simon L. Altmann, Rotations, 
%        Quaternions and Double Groups, 1986, but with the rows and columns
%        ordered in increasing values of m' and m.
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    a = cos(w / 2.) - 1i * sin(w / 2.) * cos(th);
    b = -1i * sin(w / 2.) * sin(th) * exp(-1i * ph);
    
    % Reorder elements to be consistent with current convention
    CK = [conj(a), -conj(b); b, a];
end
