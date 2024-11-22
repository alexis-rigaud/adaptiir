% To run this demo, type
%
% smm_demo
%
% from matlab.
fprintf('\nThis demo illustrates an adaptive IIR filter in normalized\n');
fprintf('lattice form, using a Steiglitz-McBride adaptation algorithm.\n');
fprintf('\nYou may select between two choices for the unknown system\n')
fprintf('(henceforth, "plant"):\n\n')
fprintf('Choice 1: A fourth order system, "easy" to identify;\n\n')
fprintf('Choice 2: A fifth order system, much "harder" to identify;\n\n')
syschoice = input('Enter "1" or "2" to select your choice: ','s');
if isempty(syschoice)
 eehh = 1;
else
 eehh = sscanf(syschoice,'%d');
 while eehh > 2
   fprintf('Sorry, only two choices are available; try again.\n')
   syschoice = input('Enter "1" or "2" to select your choice: ','s');
   if isempty(syschoice)
    eehh = 1;
   else
    eehh = sscanf(syschoice,'%d');
   end
 end
end
fprintf('\n Choice %d has been selected.\n\n',eehh);
fprintf('Figure 1 plots the poles and zeros for the system to be\n');
fprintf('identified.\n');
%
% Hard:
if eehh==2
 nud = [0.56058; 0.02272; 0.75780; -0.24935; 0.22090; 0];
 sthd = sin([0.25909; 0.33660; -0.43754; -0.37564; 0.26247]);
else
% easy:
 nud = [  -0.48161405756148;
  -0.63968456340982;
   1.32674764590246;
  -2.04373417873129;
   0.98010000000000];
 sthd = [ -0.05593757847844;
   0.44823375834249;
  -0.12494092992543;
   0.6855840000000];
end
cthd = sqrt(1.0-sthd.^2);
[pole,zero] = lat2pz(nud,sthd);
fprintf('The largest pole radius is %f\n',max(abs(pole)));
sysord = length(sthd);
t = [0:0.01:2*pi];
figure(1)
plot(cos(t),sin(t),'g')
axis('equal');
hold on
plot(zero,'o')
plot(pole,'x')
title('Poles and Zeros of Plant')
xlabel('Real Part')
ylabel('Imaginary Part')
hold off
fprintf('\n(Hit any key to continue).\n\n')
pause
fprintf('We now must select what order filter we shall adapt.\n');
fprintf('Adapting a %dth order filter yields the sufficient order case,\n',sysord);
fprintf('because the system to be identified is also %dth order.\n',sysord);
fprintf('Adapting a lower order filter gives a reduced order case.\n');
fprintf('Type an integer indicating the order of the adaptive filter:\n')
fprintf('(If this is your first run, go ahead and type %d or <return>)\n',sysord);
ss = input('Enter the adaptive filter order: ','s');
if isempty(ss)
 M = sysord;
else
 M = sscanf(ss,'%d');
end
if eehh == 1
  if M >= sysord
    iter = 5000;
    mu = 0.2;
  elseif M==3
    iter = 4000;
    mu = 0.1;
  elseif M==2
    iter = 3000;
    mu = 0.07;
  else
    iter = 3000;
    mu = 0.01;
  end
else
  if M >= sysord
   iter=25000;
   mu = 0.4;
  elseif M==4
   iter = 2500;
   mu = 0.3;
  elseif M==3
   iter = 2500;
   mu = 0.3;
  elseif M==2
   iter = 2500;
   mu = 0.08;
  else
   iter = 4500;
   mu = 0.02;
  end
end
fprintf('\nAdaptive filter order = %d\n\n',M)
fprintf('\nThe experiment will use uniformly distributed white noise as the\n')
fprintf('input for simplicity, and the adaptive\n')
fprintf('filter computations will now run %d\n',iter)
fprintf('iterations.  This may take a minute or two, so why not\n')
fprintf('go get yourself a cup of coffee?\n')
fprintf('The results will be waiting for you when you come back.\n')
fprintf('(Much faster versions are available using mex,\n')
fprintf('but are less portable for the purposes of this demonstration).\n')

scale = sqrt(12); %normalize random number generator to unit variance.
nud = nud/norm(nud); %scale system to unit L_2 norm.
xstd = zeros(sysord+1,1); %state + input vector for lattice routine
%
% initialize coefficients and states of adaptive filter.
%
sth = zeros(M,1); % sines of rotation parameters.
cth = ones(M,1); % cosines of rotation parameters.
nu = zeros(M+1,1); % tap coefficients
nutemp = nu;
nutemp(M+1) = 1; %``dummy'' tap vector.
x1 = zeros(M+1,1); % state + input vector for adaptive filter.
xpost = zeros(M+1,1); % state + input vector for second lattice.
%
% registers for plotting data
%
errvec = zeros(1,iter);
rotvec = zeros(M,iter);
nuvec = zeros(M+1,iter);
%
% Run algorithm
%
for k=1:iter
 uin = scale*(rand - 0.5); %zero mean, unit variance white noise, uniform.
 xstd(sysord+1) = uin;
 xstd = lattice(xstd,nud,sthd,cthd); %larger order system.
%
% adaptive filter computations
%
 x1(M+1) = uin;
 x1 = lattice(x1,nutemp,sth,cth);
 error = xstd(sysord+1) - nu' * x1; %output error
 nu = nu + mu * error * x1/(1 + x1'*x1);
 scal2 = mu*error/(1 + xpost'*xpost);
 for kk=M:-1:1
   value = sth(kk) - scal2 * xpost(kk);
   if abs(value) >= 1.0
     fprintf('stability problem at iteration %d\n',k);
   else
     sth(kk) = value;
   end
   scal2 = scal2 * cth(kk);
 end
 cth = sqrt(1.0 - sth.^2);
 xpost(M+1) = xstd(sysord+1);
 xpost = lattice(xpost,nutemp,sth,cth);
 errvec(k) = error^2;
 nuvec(:,k) = nu;
 rotvec(:,k) = sth;
end
%
% now write results
%
fprintf('\fNow the computations are complete, so we can plot some results.\n')
fprintf('First comes the squared output error in Figure 2.\n\n')
figure(2)
plot([1:iter],10*log10(errvec))
xlabel('Iteration Number')
ylabel('Squared Error (dB)')
if M<sysord
 sval = hsval(nud,sthd); %Hankel singular values of larger order system.
 hold on
 plot([1 iter],(20*log10(sval(M+1))*[1 1]),'r')
 hold off
 title('Output Error Squared + A Priori Bound')
else
 title('Output Error Squared')
end
if M < sysord
 fprintf('Since you have chosen a reduced order case, I have also plotted the\n')
 fprintf('Hankel singular value of index %d (= %f) for the plant.\n',M+1,sval(M+1))
 fprintf('If the algorithm converges and misadjustment is adequately\n')
 fprintf('controlled, the mean square output error will be smaller than\n')
 fprintf('the square of this singular value.\n');
end
fprintf('\n(Hit any key to continue).\n\n')
pause
fprintf('Now come plots of the coefficient trajectories,\n')
fprintf('in Figures 3 and 4.\n')
figure(3)
plot([1:iter],rotvec)
title('Sines of Rotation angles')
xlabel('Iteration Number')
ylabel('Coefficient value')
if M >= sysord
 hold on
 pval = [1;1]*sthd';
 plot([1 iter],pval,'r-.')
 hold off
end
figure(4)
plot([1:iter],nuvec)
title('Tap Parameters')
xlabel('Iteration Number')
ylabel('Coefficient Value')
if M >= sysord
 hold on
 pval = [1;1]*nud';
 plot([1 iter],pval,'r-.')
 hold off
 fprintf('For this sufficient order case, the dash-dot lines show\n')
 fprintf('the true coefficient values, for comparison purposes.\n')
 if M > sysord
  fprintf('Since you have chosen an adaptive filter order greater\n')
  fprintf('than the plant order, the converged adaptive filter will\n')
  fprintf('theoretically have pole-zero cancellations, but we cannot\n')
  fprintf('predict where.  As a result, the coefficient values need not\n')
  fprintf('match exactly those of the plant.\n')
  nud = nuvec(:,iter);
  sthd = rotvec(:,k);
  [pole,zero] = lat2pz(nud,sthd);
  figure(1)
  plot(cos(t),sin(t),'g')
  axis('equal');
  hold on
  plot(zero,'o')
  plot(pole,'x')
  title('Poles and Zeros of Adapted Filter')
  xlabel('Real Part')
  ylabel('Imaginary Part')
  hold off
  fprintf('Figure 1 now plots the poles and zeros of the adapted filter,\n')
  fprintf('after convergence, to see the poles and zeros which cancel.\n')
 end
 fprintf('\nObserve that the coefficients converge in about %d iterations,\n',iter)
 fprintf('but that the squared error has not yet reached zero.\n')
 fprintf('The plant and converged filter now have the same parameters,\n')
 fprintf('but the converged filter has not yet attained stationarity\n')
 fprintf('in its internal signals. The output error will decay to zero\n')
 fprintf('in much the same way as a "transient" or initial condition on\n')
 fprintf('the internal state will decay. The "settling time" of the filter\n')
 fprintf('depends on the pole radii, which are close to one for this\n')
 fprintf('example.\n')
else
 fprintf('Observe the fluctuations in the coefficient values.\n')
 fprintf('This is typical of reduced order cases, because the\n')
 fprintf('output error variance cannot go to zero.\n')
 fprintf('The coefficient updates will then be perpetually driven\n')
 fprintf('by the mismodelling error.\n')
end
%
% Clear registers
%
clear errvec nuvec rotvec pval
clear M k pole sval x1 cth kk scal2 syschoice xpost
clear cthd mu scale sysord xstd eehh nu ss t zero
clear error nud sth uin iter nutemp sthd value
% end of smm_demo.m
%
% Last modified: March 1997, Phil Regalia