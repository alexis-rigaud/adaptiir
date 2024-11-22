function [a,b] = lat2dir(nu,sth)
    % [a,b] = lat2dir(nu,sth)
    % converts lattice parameters to direct form parameters.
    % Filter order M is taken as length(sth).
    % nu is a vector of M+1 elements containing the tap coefficients,
    % while sth is a vector of length M containing sines of rotation angles
    % (a.k.a. reflection coefficients).
    % On output,
    % a and b contain, respectively, the denominator and numerator
    % coefficients, according to
    % A(z) = a(1) + a(2)*z + a(3)*z^2 + ... + a(M+1)*z^M
    % B(z) = b(1) + b(2)*z + b(3)*z^2 + ... + b(M+1)*z^M
    % a(1) = 1.0 should result on output.
    %
    % If size(nu) = [M+1 k] on input,
    % then each column is associated to a different transfer function,
    % each having the same poles.
    % On output, b is a matrix whose columns correspond to the numerator
    % polynomials of the different transfer functions.
    %
    % Requires file `lat2tdl.m'.

    M = length(sth);
    Tdl = lat2tdl(sth);
    a = Tdl(M+1:-1:1,M+1);
    ksize = size(nu);
    if ksize(1) == (M+1)
     b = Tdl * nu;
    elseif ksize(2) == (M+1)
     b = Tdl*nu';
    end
    % end of lat2dir().
    %
    % Last modified: February 1997, Phil Regalia