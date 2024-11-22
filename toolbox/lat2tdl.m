function Tdl = lat2tdl(sth)
    % Tdl = lat2tdl(sth)
    % calculates the transformation matrix Tdl given the sines of
    % the rotation angles of the lattice filter, stored in the
    % vector sth.
    % Used by files `lat2dir.m' and `hsval.m'.

    M = length(sth);
    cth = sqrt(1.0 - sth.^2);
    Tdl = zeros(M+1,M+1);
    Tdl(1,1) = prod(cth);
    for k=2:M+1
      Tdl(1:k,k) = ([0;Tdl(1:k-1,k-1)] + ...
        sth(k-1) * [Tdl(k-1:-1:1,k-1);0])/cth(k-1);
    end
    % end of lat2tdl().
    %
    % Last modified: February 1997, Phil Regalia