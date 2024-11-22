% To run this demonstration from matlab, type:
%
% notch_demo
%
fprintf('\nAn adaptive notch filter is useful for estimating the frequency\n');
fprintf('of a sinusoidal signal buried in broadband background noise.\n');
fprintf('Such a signal appears as\n\n');
fprintf('u(n) = ampl*cos(phi(n)) + b(n)\n\n');
fprintf('where "ampl" is the amplitude of the sinusoid,\n');
fprintf('"phi(n)" is the phase angle at time n, and "b(n)" is the\n');
fprintf('background noise.  In all simulations using this demo file,\n');
fprintf('the noise b(n) is white, with its variance equal to the\n');
fprintf('average sinusoid signal power, giving a signal-to-noise\n');
fprintf('ratio of 0dB.\n\n');
tt = pi*[0:0.0667:5];
tt = cos(tt) + 1.8171*(rand(size(tt)) - 0.5);
figure(1)
clf
plot(tt,'go')
hold on
plot([1;1]*[1:length(tt)],[0;1]*tt,'g')
xlabel('Sample number')
ylabel('Sample value')
title('Noisy sinusoid (SNR = 0dB)')
hold off;
fprintf('Figure 1 shows how such a signal can appear, using a sinusoid whose\n')
fprintf('period is 33 samples and adding uniformly distributed noise.\n\n');
fprintf('(Hit any key to continue)\n\n');
pause
fprintf('An adaptive notch filter adjusts its notch frequency in\n');
fprintf('an attempt to match the instantaneous signal frequency,\n');
fprintf('defined as the time derivative of the phase angle.\n');
fprintf('The sampled version of the phase angle is "phi(n)" above.\n');
fprintf('The instantaneous notch frequency is then used as\n');
fprintf('the instantaneous frequency estimate.\n\n');

fprintf('Four simulation possibilities are provided with this\n');
fprintf('demonstration:\n\n');

fprintf('Choice 1: Frequency hop experiment.  The signal frequency is\n');
fprintf('          abruptly changed every 1000 iterations.  The adaptive\n');
fprintf('          filter does not know when the frequency will shift,\n');
fprintf('          though, nor does it know the next frequency value.\n\n');
fprintf('Choice 2: A chirped sinusoid.  The instantaneous signal frequency\n');
fprintf('          changes linearly in time, and the adaptive filter must\n')
fprintf('          first locate it and then track it.\n\n');
fprintf('Choice 3: A "super chirp" sinusoid.  The instantaneous frequency\n');
fprintf('          changes quadratically in time, and the adaptive filter\n')
fprintf('          must first locate it and then track it.\n\n');
fprintf('Choice 4: A "mega-chirp" sinusoid.  The instantaneous frequency\n');
fprintf('          changes cubically in time, and the adaptive filter must\n')
fprintf('          first locate it and then track it.\n\n');

syschoice=input('Enter an integer (1, 2, 3 or 4) to select your choice: ','s');
if isempty(syschoice)
 switchh = 1;
else
 switchh = sscanf(syschoice,'%d');
 while switchh > 4
   fprintf('Sorry, only four choices are available; try again.\n')
   syschoice = input('Enter "1", "2", "3" or "4" to select your choice: ','s');
   if isempty(syschoice)
    switchh = 1;
   else
    switchh = sscanf(syschoice,'%d');
   end
 end
end
fprintf('\nChoice %d has been selected. ',switchh);
if switchh==1
 fprintf('(Frequency hop experiment)\n\n');
elseif switchh==2
 fprintf('(Chirped sinusoid experiment)\n\n');
elseif switchh==3
 fprintf('(Super chirped sinusoid experiment)\n\n');
elseif switchh==4
 fprintf('(Mega-chirped sinusoid experiment)\n\n');
else
 fprintf('(Default choice 4 is selected for you)\n\n');
end

iter = 4000;
insamples = zeros(1,iter);
%switchh=3;
if switchh==1
 freqs = 0.48*rand(1,4) + 0.01;
 fprintf('The frequency hop sequence is chosen randomly,\n')
 fprintf('and for this experiment the normalized frequencies are:\n\n')
 fprintf('%f\t%f\t%f\t%f\n\n',freqs);
 insamples(1:1000) = cos(2*pi*freqs(1)*[1:1000]);
 insamples(1001:2000) = cos(2*pi*freqs(2)*[1:1000]);
 insamples(2001:3000) = cos(2*pi*freqs(3)*[1:1000]);
 insamples(3001:4000) = cos(2*pi*freqs(4)*[1:1000]);
elseif switchh==2
 slope = 1/300 - 0.003*rand;
 insamples = cos(slope*([1:iter].^2));
 freqs = slope*[1:iter]/pi;
 freqs = 0.5*acos(cos(2*pi*freqs))/pi; %alias value to interval [0 0.5].
 fprintf('The instaneous signal frequency will sweep linearly from\n')
 fprintf('0 Hz to 0.5 Hz (normalized frequency) and back again.\n')
 fprintf('The sweep rate is set randomly, and for this experiment\n')
 fprintf('is %e Hz/sample.\n\n',freqs(2)-freqs(1));
elseif switchh==3
 slope = 4e-07*rand;
 insamples = cos(slope*(([1:iter]-2000).^3) + 0.02*pi*([1:iter]-2000));
 freqs = 1.5*slope*(([1:iter]-2000).^2)/pi + 0.01;
 freqs = 0.5*acos(cos(2*pi*freqs))/pi; %alias value to interval [0 0.5].
 slope = abs(freqs(4) - freqs(3) - freqs(3) + freqs(2));
 fprintf('The instantaneous signal frequency will sweep quadratically\n');
 fprintf('in time, corresponding to a constant acceleration of\n');
 fprintf('%e (Hz/sample)/sample; the acceleration value will\n',slope);
 fprintf('change each time you run this.\n\n');
else
 freqs = [1:iter] - 2000;
 slope = 3e-04*rand;
 mu = (slope + 1.2e-04)/(8e06);
 insamples = cos(pi*(0.5*freqs + slope*freqs.^2 - mu*freqs.^4));
 freqs = slope*([1:iter]-2000) - (mu+mu)*([1:iter]-2000).^3 + 0.25;
 freqs = 0.5*acos(cos(2*pi*freqs))/pi; %alias value to interval [0 0.5].
end
figure(1)
clf
if switchh==1
 xx = [0 1000 1000 2000 2000 3000 3000 4000];
 yy = [freqs(1) freqs(1) freqs(2) freqs(2) freqs(3) freqs(3) freqs(4) freqs(4)];
 plot(xx,yy,'g-.');
else
 plot([1:iter],freqs,'g-.')
end
axis([0 iter 0 0.5]);
xlabel('Iteration Number');
ylabel('Normalized Frequency (Hz)');
title('Instantaneous Frequency of Input Signal (SNR = 0dB)');
hold on;
fprintf('\nFigure 1 now shows the instantaneous frequency evolution\n')
fprintf('of the input signal.\n')
if switchh==4
 fprintf('The curve so plotted is a cubic function of the time index,\n')
 fprintf('and will change somewhat each time you run this simulation.\n');
end
fprintf('\n(Hit any key to continue)\n\n');
pause;
fprintf('We will compare two adaptive notch filter algorithms.\n');
fprintf('First comes a simplified lattice algorithm,\n');
fprintf('corresponding to Table 10.3, p. 579 of\n')
fprintf('"Adaptive IIR Filtering in Signal Processing and Control"\n');
fprintf('(Marcel Dekker, New York, 1995.)\n\n');
fprintf('The notch frequency is initialized to a midband value\n')
fprintf('(0.25 Hz in normalized frequency).\n')
fprintf('The adaptive filter computations are now running and,\n');
fprintf('after %d iterations, the frequency estimates from the\n',iter);
fprintf('adaptive notch filter will be superimposed on Figure 1.\n');
%
% First comes simplified lattice adaptive notch filter algorithm.
%
% Add noise to input samples:
%
stdev = 2*(0.75)^(1/3); %set variance to 0.5 = average signal power.
insamples = insamples + stdev*(rand(1,iter)-0.5); %noisy input signal
%
% initialize adaptive filter parameters
%
xst = zeros(2,1); % state of adaptive filter.
temp = zeros(2,1); % intermediate signals.
pihalf = 0.5*pi; % pi/2.
theta = 0; % initial value for notch frequency parameter.
sth = sin(theta);
cth = cos(theta);
bw = 0.20*pi; % bandwidth parameter for notch filter.
sth2 = sin(bw);
cth2 = cos(bw);
mu = 0.023; % adaptive filter step size.
%
% Run adaptive lattice notch filter:
%
freqold = 0.25; %initial notch frequency.
for kk=1:iter
 insig = insamples(kk);
 temp = [cth2 -sth2;sth2 cth2]*[insig; xst(2)];
 error = mu*(insig + temp(2)); %notch filter output times step size.
 theta = theta - error*xst(1); %coefficient update.
 freqnew = 0.5*acos(cos(theta+pihalf))/pi; %instantaneous freq. estimate.
 if kk>1
   plot([kk-1 kk], [freqold freqnew]);
 end
 freqold = freqnew;
 sth = sin(theta);
 cth = cos(theta);
 xst = [cth -sth;sth cth] * [temp(1);xst(1)];
end
title('Frequency Estimates of Lattice Adaptive Notch Filter');
hold off
%
% Now comes direct form adaptive notch filter algorithm
%
fprintf('\nWe will now feed the same noisy input sequence to a common\n')
fprintf('direct-form adaptive notch filter, having constrained poles\n')
fprintf('and zeros; the algorithm may be found in Table 10.1, p. 573 of\n')
fprintf('"Adaptive IIR Filtering in Signal Processing and Control".\n');
figure(2)
clf
if switchh==1
 plot(xx,yy,'g-.');
else
 plot([1:iter],freqs,'g-.')
end
axis([0 iter 0 0.5]);
xlabel('Iteration Number');
ylabel('Normalized Frequency (Hz)');
title('Instantaneous Frequency of Input Signal (SNR = 0dB)');
hold on;
fprintf('\nFigure 2 shows again the instantaneous frequency evolution\n')
fprintf('of the input signal.\n')
fprintf('\n(Hit any key to continue)\n');
pause
%
% Direct form adaptive filter computations
% Table 10.1, p. 573.
% r^2 -> sth2 : same pole radius
% Del_a -> cth2
% insig -> temporary delay register.
% a(n) -> theta : notch frequency parameter
% pihalf -> pi/4 (for alias function).
%
fprintf('\nThe direct form algorithm is now running and,\n');
fprintf('after %d iterations, the instantaneous fequency\n', iter);
fprintf('estimates will be superimposed on Figure 2.\n');
fprintf('For consistency, the same pole radius as for the lattice\n');
fprintf('algorithm is used, and the initial notch frequency is set\n');
fprintf('again to 0.25 Hz.\n')
%
% Initialize parameters
%
mu = 0.013; % adaptive filter step size.
sth2 = sqrt(sth2);
temp = zeros(2,1);
xst = zeros(2,1);
theta = 0; % set initial notch frequency to 0.25 Hz.
pihalf = 0.5*pihalf; % now equals pi/4.
freqold = 0.25; % initial notch frequency
for kk=1:iter
 insig  = xst(1) - sth2*theta*temp(1) - sth2*sth2*temp(2);
 cth2 = insig - sth2*temp(2);
 temp(2) = temp(1);
 temp(1) = insig; % delay line moved, insig is now free.
 insig = insamples(kk) - sth2*theta*xst(1) - sth2*sth2*xst(2);
 error = insig + theta*xst(1) + xst(2); % notch filter output.
 xst(2) = xst(1);
 xst(1) = insig; % delay line moved; insig is now free.
 insig = theta - mu*error*cth2;
 if abs(insig) <= 2.0 % check to see if |theta| <= 2 (for consistency).
  theta = insig;      % if so, no problem, do update.
 else
  theta = asin(sin(insig*pihalf))/pihalf; % if not, alias update into [-2,2].
 end
 freqnew = 0.5*acos(-0.5*theta)/pi; %instantaneous notch frequency.
 if kk>1
  plot([kk-1 kk], [freqold freqnew]);
 end
 freqold = freqnew;
end
title('Frequency Estimates of Direct-Form Adaptive Notch Filter');
hold off;
%
% final comments
%
fprintf('\nIn comparing Figures 1 and 2, you may observe that for frequencies\n')
fprintf('in the midband range (not too far from 0.25 Hz), the two algorithms\n')
fprintf('provide comparable tracking capabilities, and visually comparable\n')
fprintf('variances in the frequency estimates.\n\n');
fprintf('For frequencies near the band edges 0 Hz and 0.5 Hz, on the other hand,\n')
fprintf('the variance in the frequency estimate is noticeably more pronounced\n')
fprintf('with the direct-form algorithm, in Figure 2.\n\n')
fprintf('This is because the noise gain of the direct-form filter varies\n')
fprintf('with the notch frequency, both for the notch filter output signal\n');
fprintf('and the filtered regressor used to update the notch frequency\n')
fprintf('parameter.\n')
fprintf('The noise gain of the lattice filter, by contrast, is fairly\n')
fprintf('uniform versus frequency, such that the variance of the frequency\n')
fprintf('estimates do not change much over the frequency range.\n')
%
% clear registers
%
if switchh==1
 clear xx yy;
else
 clear slope;
end
clear freqnew iter sth theta bw freqold kk sth2 xst;
clear cth freqs mu switchh cth2 insamples pihalf syschoice;
clear error insig stdev temp tt;
% end of notch_demo.m
%
% Last modified: March 1997, Phil Regalia