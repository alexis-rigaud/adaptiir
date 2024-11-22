function r = lat2cor(N,nu,sth,cth)
    % r = lat2cor(N,nu,sth,cth)
    % returns to column vector r the first N terms of the autocorrelation
    % sequence of the transfer function parametrized in normalized lattice form.
    % Filter order (M, say) is taken as length(sth).
    % nu is a vector of M+1 elements containing the tap coefficients.
    % sth is a vector of M elements containing the sines of the rotation angles,
    % a.k.a. reflection coefficients.
    % cth contains the corresponding cosines, and may be omitted from the
    % input list.
    %
    % Because of indexing conventions, the zeroth lag r_0 will be in position r(1),
    % the first lag r_1 will be in position r(2), and so on.
    %
    % Requires file `lattice.m'.

    if nargin==3,
     cth = sqrt(1.0 - sth.^2);
    end
    r = zeros(N,1);
    M = length(sth);
    r(1) = sum(nu.^2);
    xst = zeros(M+1,1);
    xst(1:M) = nu(1:M);
    for k=[2:N]
     xst = lattice(xst,nu,sth,cth);
     r(k) = xst(M+1);
     xst(M+1) = 0.0;
    end
    % end of lat2cor().
    %
    % Last modified: February 1997, Phil Regalia