function [nu,sth,cth,Tdl] = dir2lat(a,b)
    % [nu,sth,cth,Tdl] = dir2lat(a,b)
    % converts the direct form parameters to lattice form.
    % a and b are vectors which contain the direct form coefficients,
    % according to
    % A(z) = a(1) + a(2)*z + a(3)*z^2 + ... + a(M+1)*z^M
    % B(z) = b(1) + b(2)*z + b(3)*z^2 + ... + b(M+1)*z^M
    % as denominator and numerator, respectively.
    % If vectors a and b have different lengths,
    % the shorter of the two is appended with zeros, and the degree M is
    % taken as length(a)-1.
    % On output, nu contains the tap parameters, while sth and cth contain
    % the sines and cosines of the rotation angles.
    % Tdl is the triangular transformation matrix connecting
    % the direct form and lattice descriptions.
    %
    % [nu,sth] = dir2lat(a,b)
    % returns only the tap parameters in vector nu, and the sines
    % of the rotation angles in vector sth.

    [k1,k2] = size(a);
    if k1<k2,
     a = a';
    end
    [k1,k2] = size(b);
    if k1<k2,
     b = b';
    end
    k1 = length(a);
    k2 = length(b);
    if k1 > k2,
      b = [b; zeros(k1-k2,1)];
    elseif k1 < k2
     a = [a; zeros(k2-k1,1)];
    end
    a = a/a(1);
    b = b/a(1);

    M = length(a)-1;
    dL = a;
    scale = 1.0;
    Tdl = zeros(M+1,M+1);
    sth = zeros(M,1);
    cth = ones(M,1);
    for k=M:-1:1
     sth(k) = dL(k+1);
     drev = zeros(M+1,1);
     drev(1:k+1) = dL(k+1:-1:1);
     temp = 1.0 - sth(k)^2;
     dL = (dL - sth(k) * drev)/temp;
     Tdl(:,k+1) = scale * drev;
     cth(k) = sqrt(temp);
     scale = scale * cth(k);
    end
    Tdl(1,1) = scale;
    nu = Tdl\b;
    % end of dir2lat
    %
    % Last modified: March 1997, Phil Regalia