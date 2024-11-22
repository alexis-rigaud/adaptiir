function xout = lattice(xin,nu,sth,cth)
    % xout = lattice(xin,nu,sth,cth)
    % performs internal computations of lattice filter.
    % xin is a column vector whose first M components are the
    % state vector, and whose last component is the input sample.
    % nu is an M+1 vector containing the tap coefficients.
    % sth is an M vector containing the sines of the rotation angles,
    % while cth contains the corresponding cosines.
    % cth may be omitted from input list.
    % Filter order M is taken as length(sth).
    % xout contains the updated state vector in the first M components,
    % and the output sample in the M+1st component.

    M = length(sth);
    if nargin == 3,
     cth = sqrt(1.0 - sth.^2);
    end
    xout = xin;
    for k=[M:-1:1]
     xout(k:k+1) = [-sth(k) cth(k);cth(k) sth(k)] * xout(k:k+1);
    end
    ksize = size(nu);
    if ksize(2) == 1
     xout(M+1) = nu' * xout;
    elseif ksize(1) == 1
     xout(M+1) = nu * xout;
    else
     fprintf('size(nu) = %d   %d; nu should be a vector.\n',ksize(1),ksize(2));
     fprintf('nu = [0 ... 0 1] has been substituted\n');
    end
    %end of lattice()
    %
    % Last modified: March 1997, Phil Regalia