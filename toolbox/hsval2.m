function [sval,Wroot] = hsval2(nu,sth)
    % [sval,Wroot] = hsval2(nu,sth)
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
    % sval = hsval2(nu,sth) returns a vector containing the Hankel singular values.
    % See also hsval in this package, as well as dbalreal in control toolbox.

    % Form Q matrix = orthogonal Hessenberg matrix:
    %
    M = length(sth); %McMillan degree
    cth = sqrt(1.0 - sth.^2);
    Q = eye(M+1);
    for k=1:M
     Q(:,k:k+1) = Q(:,k:k+1) * [-sth(k) cth(k); cth(k) sth(k)];
    end
    %
    % Form Cholesky factor of Toeplitz matrix:
    %
    Toe = zeros(M+1,M+1);
    Toe(1,1) = 1.0/prod(cth);
    for k=2:M+1
     Toe(:,k) = Q * Toe(:,k-1);
    end                %Toe contains Cholesky factor of Toeplitz matrix.
    den = zeros(M+1,1);
    den(M+1) = 1.0;
    den = Toe\den; %characteristic polynomial
    %
    % Form c vector; Amat is Q(1:M,1:M).
    %
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
     X(:,M-k) = (Q(1:M,1:M))' * X(:,M+1-k) + den(M+1-k) * cvec';
    end
    Wroot = Toe(1:M,1:M) * X'; %square root of observability gramian.
    sval = svd(Wroot); % Hankel singular values.
    % end of hsval2().
    %
    % Last modified: March 1997, Phil Regalia