function h = lat2imp(N,nu,sth,cth)
    % h = lat2imp(N,nu,sth,cth)
    % returns to column vector h the first N terms of the impulse
    % response of a transfer function parametrized in normalized lattice form.
    % Filter order (M, say) is taken as length(sth).
    % nu is a vector of length M+1 containing the tap coefficients.
    % sth is a vector of length M containing the sines of the rotation angles,
    % a.k.a. reflection coefficients.
    % cth contains the corresponding cosines, and may be omitted from
    % the input list.
    %
    % Because of indexing conventions, sample h_0 will be in position h(1),
    % sample h_1 will be in position h(2), and so forth.
    %
    % Requires file `lattice.m'.

    if nargin==3,
     cth = sqrt(1.0 - sth.^2);
    end
    h = zeros(N,1);
    M = length(sth);
    xst = zeros(M+1,1);
    xst(M+1) = 1.0;
    for k=[1:N]
     xst = lattice(xst,nu,sth,cth);
     h(k) = xst(M+1);
     xst(M+1) = 0.0;
    end
    % end of lat2imp().
    %
    % Last modified: March 1997, Phil Regalia