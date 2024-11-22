function [pole,zero] = lat2pz(nu,sth)
    % [pole,zero] = lat2pz(nu,sth)
    % converts the lattice description to poles and zeros of the
    % corresponding transfer function.
    % Filter order (M, say) is taken as length(sth).
    % nu is a vector of M+1 entries, containing the tap coefficients.
    % sth is an M vector containing the sines of the rotation angles, a.k.a.
    % reflection coefficients.
    % On output, pole and zero are M vectors containing, respectively,
    % the poles and zeros of the transfer function, expressed in the
    % z^(-1) plane, so that a stable system will yield each entry of the
    % `pole' vector with modulus less than one.
    % If transfer function has a pure delay factor of order k,
    % then the transfer function has k zeros at infinity in the z^(-1) plane,
    % and the first k entries of the vector `zero' should contain `Inf'.

    M = length(sth);
    cth = sqrt(1.0 - sth.^2);
    % build Q matrix
    Q = eye(M+1);
    for k = 1:M
     Q(:,k:k+1) = Q(:,k:k+1) * [-sth(k) cth(k);cth(k) sth(k)];
    end
    % compute poles
    Amat = Q(1:M,1:M);
    pole = eig(Amat);
    bvec = Q(1:M,M+1);
    ksize = size(nu);
    if ksize(1) == 1
     nu = nu';
    end
    cvec = nu' * Q;
    cvec = cvec/norm(cvec);
    d = cvec(M+1);
    cvec = cvec(1:M);
    % deflate away zeros at infinity
    count = 0;
    temp = 10^(-8);
    while abs(d) <= temp,
     d = cvec*bvec;
     cvec = cvec*Amat;
     count = count+1;
    end
    % now compute zeros; those which were at infinity are now at zero.
    zero = eig(Amat - bvec*cvec/d);
    if count > 0
     zero = sort(zero);
     zero(1:count) = inf * ones(count,1);
    end
    % end of lat2pz().
    %
    % Last modified: March 1997, Phil Regalia