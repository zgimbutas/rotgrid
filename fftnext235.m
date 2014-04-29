function next235=fftnext235(n)
%FFTNEXT235 Returns the next multiple of 2, 3, and 5.
%
%  next235=fftnext235(n) returns the next integer that 
%                     is a multiple of 2, 3, and 5.
%
next235 = n;
mex_id_ = 'fftnext235(i int[x], io int[x])';
[next235] = rotgrid_r2014a(mex_id_, n, next235, 1, 1);

