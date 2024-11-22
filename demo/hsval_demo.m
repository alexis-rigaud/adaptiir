% To run this demonstration, type
%
% hsval_demo
%
% from matlab.
% N.B.: This demonstration includes a comparison using some routines from
% the control toolbox of matlab.  If you do not have the control toolbox,
% this demo may not run.
fprintf('\nHankel singular values are key input-output invariants\n');
fprintf('in linear system theory, and are also known as second \n');
fprintf('order modes.  Among many applications in linear system\n');
fprintf('theory, the Hankel singular values give precisely the\n');
fprintf('distance, in Hankel norm, between a given function and\n');
fprintf('a reduced order approximant.\n\n');
fprintf('This demonstration file illustrates how to use HSVAL\n');
fprintf('and HSVAL2 to compute Hankel singular values, using\n');
fprintf('a normalized lattice description of a transfer function.\n\n');
fprintf('An all-pass function has all Hankel singular values equal\n');
fprintf('to one, and hence serves as a useful test example which\n');
fprintf('allows us to compare computed with theoretical values.\n\n');
fprintf('(Hit any key to continue)\n\n');
pause
fprintf('The filter parameters are set as:\n\n');
fprintf('sth = 2*rand(1,30) - 1;\n');
fprintf('nu = zeros(31,1); nu(31) = 1;\n\n');
sth = 2*rand(1,30) - 1;
nu = zeros(31,1); nu(31) = 1;
fprintf('Setting nu to the last unit vector ensures that our transfer\n');
fprintf('function is all-pass.  The vector sth contains the reflection\n');
fprintf('coefficients which set the poles, and is set to a random vector\n');
fprintf('having elements uniformly distributed between -1 and +1, so that\n');
fprintf('the resulting filter is stable; the results will change each time\n');
fprintf('you run this demonstration. ')
fprintf('The filter order is thirty = length(sth).\n\n');
fprintf('The computation of Hankel singular values is not necessarily\n');
fprintf('a well conditioned problem.  These quantities may be defined\n');
fprintf('as the the square roots of the eigenvalues of the product\n');
fprintf('of the controllability and observability gramian matrices.\n');
fprintf('These matrices, in turn, are the solutions to two Lyapunov\n');
fprintf('equations, and it is known that these Lyapunov equations can\n');
fprintf('be ill conditioned if the transfer function has poles\n');
fprintf('quite close to the unit circle.\n');
fprintf('To see the pole locations for this example,\n');
fprintf('strike any key.\n');
poles = lat2pz(nu,sth);
pause
figure(1)
clf
t = [0:0.01:2*pi];
plot(cos(t),sin(t),'r');
axis('equal')
hold on
plot(poles,'x');
title('Pole locations for this example');
xlabel('Real part')
ylabel('Imaginary part')
hold off
radius_max = max(abs(poles));
fprintf('\n Largest pole radius = %16.14f\n\n',radius_max);
format long;
format compact;
fprintf('\nWe can now compute the Hankel singular values with the command\n\n');
fprintf('sval = hsval(nu,sth);\n\n');
fprintf('The singular values should theoretically all equal one;\n');
fprintf('hit any key to see the results.\n');

sval1 = hsval(nu,sth);
bias1 = sum(sval1-1.0)/length(sval1);
mean1 = bias1+1.0;
stdev1 = sqrt(sum((sval1-mean1).^2)/length(sval1));
pause
figure(2)
plot(sval1)
title('Computed Hankel singular values using HSVAL')
xlabel('Index')
ylabel('Computed Hankel singular value')
sval = sval1'
fprintf('Empirical mean = %f\n', mean1);
fprintf('Empirical standard deviation = %e\n\n', stdev1);
fprintf('(Hit any key to continue)\n\n');
pause
fprintf('Since the computation of Hankel singular values is not always\n');
fprintf('a well conditioned problem, a second routine, called hsval2, is\n');
fprintf('also provided.  Whereas hsval uses hyperbolic rotations to\n');
fprintf('build up a key matrix, hsval2 uses exclusively plane rotations.\n');
fprintf('The command is:\n\n');
fprintf('sval = hsval2(nu,sth);\n\n');
fprintf('Hit any key to see the results\n');
sval = hsval2(nu,sth);
bias2 = sum(sval-1.0)/length(sval);
mean2 = bias2+1.0;
stdev2 = sqrt(sum((sval-mean2).^2)/length(sval));
pause
figure(2)
plot(sval)
title('Computed Hankel singular values using HSVAL2');
xlabel('Index')
ylabel('Computed Hankel singular value')
sval = sval'
fprintf('Empirical mean = %f\n',mean2);
fprintf('Empirical standard deviation = %e\n\n',stdev2);
fprintf('(Hit any key to continue)\n\n');
pause
fprintf('One can also compute Hankel singular values using any system\n')
fprintf('balancing routine, such as DBALREAL in the control toolbox\n');
fprintf('(if you have this).  To compare the results, I can generate\n');
fprintf('the same orthogonal internal description as that used for the\n');
fprintf('the normalized lattice filter in this example, and then pass it\n')
fprintf('to the DBALREAL routine.  The internal description so passed on\n');
fprintf('is already in balanced form, but the routine has no way of\n');
fprintf('of knowing this.  It will then attempt to solve the two Lyapunov\n')
fprintf('equations, which need not be well conditioned, in order to compute\n')
fprintf('a balancing transformation, from which the Hankel singular values\n')
fprintf('may be obtained as a byproduct.\n\n')
fprintf('Assuming you have the control toolbox installed, we can now run\n')
fprintf('the computations to see how the results compare.  If you do not\n')
fprintf('have the control toolbox, matlab will probably spit out some weird\n')
fprintf('error message.  Strike any key to continue, but be warned that\n')
fprintf('this problem might be too ill conditioned for that routine to\n')
fprintf('handle, depending on the seed used for the random number generator\n')
fprintf('which gave us the values in the vector sth.\n');

M = length(sth);
cth = sqrt(1.0 - sth.^2);
Q = eye(M+1);
for k = 1:M
 Q(:,k:k+1) = Q(:,k:k+1) * [-sth(k) cth(k);cth(k) sth(k)];
end
Amat = Q(1:M,1:M);
bvec = Q(1:M,M+1);
cvec = Q(M+1,1:M);
pause
fprintf('\n(DBALREAL is working now)\n')
[Ab,Bb,Cb,Msval] = dbalreal(Amat,bvec,cvec);
bias3 = sum(Msval-1.0)/length(Msval);
mean3 = bias3+1.0;
stdev3 = sqrt(sum((Msval-mean3).^2)/length(Msval));
figure(2)
sval = sort(Msval);
sval(1:M) = sval(M:-1:1)
plot(sval)
hold on
title('Computed Hankel singular values using DBALREAL')
xlabel('Index')
ylabel('Computed Hankel singular value')
plot(sval1,'-.')
if (sval(1) - 1) < (1 - sval(M))
 hght1 = 0.1*(sval(1) - sval(M)) + sval(M);
 hght2 = 0.05*(sval(1) - sval(M)) + sval(M);
 rnge = M*[0.1 0.2];
else
 hght1 = 0.9*(sval(1) - sval(M)) + sval(M);
 hght2 = 0.85*(sval(1) - sval(M)) + sval(M);
 rnge = M*[0.6 0.7];
end
plot(rnge,[hght1 hght1],'-')
text((rnge(2)+1),hght1,'DBALREAL')
plot(rnge,[hght2 hght2],'-.')
text((rnge(2)+1),hght2,'HSVAL')
hold off
fprintf('Empirical mean = %f\n',mean3)
fprintf('Empirical standard deviation = %e\n\n',stdev3)
fprintf('Since we have gotten this far, let us compare the results of the\n')
fprintf('three routines:\n\n');
fprintf('                      HSVAL             HSVAL2            DBALREAL\n')
fprintf('Bias in mean value:   %-e\t %-e\t %-e\n',bias1,bias2,bias3)
fprintf('Standard deviation:   %-e\t %-e\t %-e\n',stdev1,stdev2,stdev3)
%
% Clear registers
%
clear Ab bias2 mean1 stdev2 Amat bias3 mean2 stdev3
clear Bb bvec  mean3 sth Cb cth nu sval M cvec poles sval1
clear Msval hght1 radius_max t Q hght2 rnge bias1 k stdev1
%end of hsval_demo
%
% Last modified: March 1997, Phil Regalia