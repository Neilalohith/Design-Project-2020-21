clear
clc
close all
%set sample time and signal frequency
t = (0:0.01:40); %microseond
f = 0.1; %MHz

%define variable parameters%

RF_signal_amplitude=0.5; % g=gain, t=time constant, c=DC offset
RF_signal_offset=0.5; % applied to signal after integrator,
g1= 3; % number indicates integrator param.
t1=0.4; % zh=quantizer on, z1=quantizer off,
zh=1.8; % yh=output high, y1=output low
zl=1.2;
yh=2;
yl=0;
%input RF signal
RF_in = RF_signal_amplitude*sin(2*pi*f*t) + RF_signal_offset;
%setup integrator and output signal arrangement
bit_out = ones(size(t));
int1_out = zeros(size(t));
%integrator transfer function
f1=1/t1;
n1=g1;
d1=[1 f1];
[A,B,C,D] = tf2ss(n1,d1);
syst_ss1=ss(A,B,C,D);

%simulate ADC system%

for i=1:length(t)-1
 %signal through integrator
 [z1] = lsim(syst_ss1,RF_in+bit_out,t);

 %save data point
 int1_out(i) = z1(i);

 %signal into quantizer
 if z1(i)>=zh
 bit_out(i+1)=yl; % NEW BOUNDS?
 elseif z1(i) <= zl
 bit_out(i+1)=yh;
 else
 bit_out(i+1)=bit_out(i);
 end
end
% create a 3rd -order LP filter
fs=1/(0.01);
f_cutoff=0.12;
f_norm=f_cutoff/(fs/2);
[b,a] = butter(5, f_norm, 'low');
%filter binary output
demod=filter(b,a,bit_out);

%PLOTTING%

subplot(5,1,1); plot(t,RF_in,'k')
title ('Input RF Signal')
xlim([10 40])
subplot(5,1,2); plot(t, RF_in+0.45*bit_out, 'k')
title('Error Signal')
xlim([10 40])
subplot(5,1,3); plot(t(1:end-1), int1_out(1:end-1), 'k')
title('Integrator Output')
xlim([10 40])
ylabel('Power (mW)')
subplot(5,1,4); plot(t,bit_out, 'k')
title('Binary Output of ADC')
xlim([10 40])
ylim([
-0.5 2.5])
subplot(5,1,5); plot(t, demod, 'k')
title('Demodulated Binary Output')

xlim([10 40])
xlabel('Time (\mus)')