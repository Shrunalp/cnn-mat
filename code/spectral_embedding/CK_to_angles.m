function [w, th, ph] = CK_to_angles(CK)
% CK_to_angles converts a rotation expressed by a Cayley-Klein matrix to the
% equivalent rotation angle and spherical angles of the rotation axis.
% 
% Inputs:
%   CK - 2-dimensional representation of SU(2), using the conventions
%        established in Section 14-2 of Simon L. Altmann, Rotations, 
%        Quaternions and Double Groups, 1986, but with the rows and columns
%        ordered in increasing values of m' and m.
%
% Outputs:
%   w  - rotation angle in the interval [0, 2 \pi].
%   th - polar angle of rotation axis in the interval [0, \pi].
%   ph - aximuthal angle of rotation axis in the interval [-\pi, \pi].
%
% Copyright 2019 Jeremy Mason
%
% Licensed under the Apache License, Version 2.0, <LICENSE-APACHE or
% http://apache.org/licenses/LICENSE-2.0> or the MIT license <LICENSE-MIT or
% http://opensource.org/licenses/MIT>, at your option. This file may not be
% copied, modified, or distributed except according to those terms.
    a = CK(2, 2);
    b = CK(2, 1);
    
    w  = 2. * acos(max(min(real(a), 1.), -1.));
    if w > eps
        tmp = -sin(w / 2.);
        th = acos(max(min(imag(a) / tmp, 1.), -1.));
        if th > eps
            tmp = tmp * sin(th);
            ph = atan2(real(b) / tmp, imag(b) / tmp);
        else
            ph = 0.;
        end
    else
        th = 0.;
        ph = 0.;
    end
end
