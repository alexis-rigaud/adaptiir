function [sval,Wroot] = hsval(nu,sth)
    % [sval,Wroot] = hsval(nu,sth)
    % calculates the Hankel singular values of the rational
    % function described by its normalized lattice parameters
    % nu and sth.  Filter order (M, say) is taken as length(sth).
    % On input:
    %  nu = vector of M+1 elements containing the tap coefficients;
    %  sth = vector of M elements containing sines of rotation
    %        angles, a.k.a. the reflection coefficients.
    % On output:
    %  Wroot = M by M square root of observability gramian W, such that
    %          W = Wroot' * Wroot;
    %  sval = M vector containing Hankel singular values = svd(Wroot).
    % Controllability gramian of normalized lattice is the identity matrix,
    % which need not be computed.
    %
    % sval = hsval(nu,sth) returns a vector containing the Hankel singular values.
    % See also `hsval2' in this package, as well as `dbalreal' in control toolbox.
    %
    % Requires file `lat2tdl.m'.

    % Form Q matrix:
    %
    M = length(sth); %McMillan degree
    cth = sqrt(1.0 - sth.^2);
    Q = eye(M+1);
    for k=1:M
     Q(:,k:k+1) = Q(:,k:k+1) * [-sth(k) cth(k); cth(k) sth(k)];
    end
    %
    % Get matrix Tdl
    %
    Tdl = lat2tdl(sth);
    %
    % form A matrix and c vector
    %
    Amat = Q(1:M,1:M);
    [k,n] = size(nu);
    if k < n,
     cvec = nu * Q;
    else
     cvec = nu' * Q;
    end
    cvec = cvec(1:M);
    X = zeros(M,M);
    X(:,M) = cvec';
    for k = 1:M-1
     X(:,M-k) = Amat' * X(:,M+1-k) + Tdl(M+1-k,M+1) * cvec';
    end
    Wroot = Tdl(1:M,1:M)\X'; %square root of observability gramian.
    sval = svd(Wroot); % Hankel singular values.
    % end of hsval().
    %
    % Last modified: February 1997, Phil Regalia