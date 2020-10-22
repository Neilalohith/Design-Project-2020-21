clear
clc
close all
%set sample time and signal frequency
t = (0:0.1:40); %microsecond
f = 0.1; %MHz

%define variable parameters%

RF_signal_amplitude=0.1; % g=gain, t=time constant, c=DC offset
RF_signal_offset=0.15; % applied to signal after integrator,
g1=1; % controlled by SOA current
t1=1; % microsecond (given by loop length)
%input RF signal
RF_in = RF_signal_amplitude*square(2*pi*f*t) + RF_signal_offset;
%setup integrator and output signal arrangement
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
 [z1] = lsim(syst_ss1,RF_in,t);
 

 %save data point
 int1_out(i) = z1(i);

end
%%%%%%PLOTTING%

subplot(2,1,1); plot(t,RF_in,'k')
title ('Input RF Signal','Fontsize',12)
xlabel('Time(\mus)','Fontsize',12)
ylabel('Power (mW)','Fontsize',12)
xlim([0 40])
ylim([-0.09 0.351])
subplot(2,1,2); plot(t(1:end-1), int1_out(1:end-1), 'k')
title('Integrator Output','Fontsize',12)
xlabel('Time(\mus)','Fontsize',12)
ylabel('Power (mW)','Fontsize',12)
xlim([0 40])