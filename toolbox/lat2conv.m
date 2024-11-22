function yout = lat2conv(uin,nu,sth,xinit)
    %
    % yout = lat2conv(uin,nu,sth,xinit)
    %
    % convolves the input sequence, contained in vector uin, with
    % the impulse response of a normalized lattice filter.
    % The output samples are returned to the vector yout, whose
    % length is that of uin.
    % nu and sth are the tap parameters and reflection coefficients,
    % respectively, of the filter.
    % xinit is the initial condition on the state vector; if omitted
    % from the input list, as in
    %
    %  yout = lat2conv(uin,nu,sth)
    %
    % the initial condition on the internal state is set to zero.
    % To obtain an output sequence longer than the input sequence,
    % simply zero-pad the input.
    %
    % Requires file "lattice.m".

    NN = length(nu);
    if nargin==3
     xstate = zeros(NN,1);
    else
     xstate = xinit;
    end
    yout = zeros(size(uin));
    limit = length(uin);
    cth = sqrt(1.0-sth.^2);
    for kk=1:limit
     xstate(NN) = uin(kk);
     xstate = lattice(xstate,nu,sth,cth);
     yout(kk) = xstate(NN);
    end
    % end of lat2conv()
    %
    % Last modified: March 1997, Phil Regalia