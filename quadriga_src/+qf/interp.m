function zi = interp( x, y, z, xc, yc, use_double )
%INTERP 2D linear interpolation optimized for performace
%
% Calling object:
%   None (static method)
%
% Description:
%   This function implements a 2D linear interpolation which is highly optimized for fast
%   execution. All calculations are done in single-precision floating point (30% faster than double
%   precision, but less accurate), and multiple data sets can be interpolated simultaneously. One-
%   dimensional linear interpolation can be done by using zi = interp( x, 0, z, xc )
%
% Input:
%   x
%   Vector of x sample points; size [ 1, nx ] or [ nx, 1 ]
%
%   y
%   Vector of y sample points; size [ 1, ny ] or [ ny, 1 ]
%
%   z
%   The input data matrix; size [ ny, nx, ne ]; the 3rd dimension allows for interpolations of
%   multiple data-sets; for one-dimensional interpolation, the size must be [ 1, nx, ne ]
%
%   xc
%   Vector of x sample points after interpolation; size [ 1, nxi ] or [ nxi , 1 ]
%
%   yc
%   Vector of y sample points after interpolation; size [ 1, nyi ] or [ nyi, 1 ]
%
%   use_double
%   If set to true, double precision is used instead of single precision.
%
% Output:
%   zi
%   The interpolated data; size [ nyi, nxi, ne ]
%
%
% QuaDRiGa Copyright (C) 2011-2023
% Fraunhofer-Gesellschaft zur Foerderung der angewandten Forschung e.V. acting on behalf of its
% Fraunhofer Heinrich Hertz Institute, Einsteinufer 37, 10587 Berlin, Germany
% All rights reserved.
%
% e-mail: quadriga@hhi.fraunhofer.de
%
% This file is part of QuaDRiGa.
%
% The Quadriga software is provided by Fraunhofer on behalf of the copyright holders and
% contributors "AS IS" and WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES, including but not limited to
% the implied warranties of merchantability and fitness for a particular purpose.
%
% You can redistribute it and/or modify QuaDRiGa under the terms of the Software License for 
% The QuaDRiGa Channel Model. You should have received a copy of the Software License for The
% QuaDRiGa Channel Model along with QuaDRiGa. If not, see <http://quadriga-channel-model.de/>.

if ~exist( 'yc','var' )
    yc = [];
end

if exist( 'use_double','var' )
    if use_double
        zi = quadriga_lib.interp( x, y, double(z), xc, yc );
    else
        zi = quadriga_lib.interp( x, y, single(z), xc, yc );
    end
else % Auto detect based on type of "z"
    zi = quadriga_lib.interp( x, y, z, xc, yc );
end

end
