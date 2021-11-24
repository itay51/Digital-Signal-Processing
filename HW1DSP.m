v
%% Q1 :IIR Digital filters
%%define constants
T=1;
w_s=0.35*pi;
w_p=0.25*pi;
delta_s=0.175;
delta_p=0.1;
%impulse invariance analog filter
[b1,a1]=butter(8,0.2735*pi,'s');
[H1,w1]=freqs(b1,a1,1024);
%bilinear transformation analog filter
[b2,a2]=butter(7,0.291*pi,'s');
[H2,w2]=freqs(b2,a2,1024);
figure;
plot(w1,abs(H1),w2,abs(H2));
xlabel('\Omega [rad/sec]');
ylabel('Amplitude');
title('Analog Low Pass Filter');
legend('Impulce Invariance','Bilinear Transformation');
xlim([0 3]);
%impulse invariance digital filter
[b_z1,a_z1]=impinvar(b1,a1,1/T);
%bilinear transformation digital filter
[b_z2,a_z2]=bilinear(b2,a2,1/T);
x=0:pi/1024:pi-pi/1024;
[H_z1,w_z1]=freqz(b_z1,a_z1,1024);
[H_z2,w_z2]=freqz(b_z2,a_z2,1024);
figure;
plot(w_z1,abs(H_z1),w_z2,abs(H_z2));
xlabel('\omega [rad]');
ylabel('Amplitude');
title('Digital Low Pass Filter');
legend('Impulce Invariance','Bilinear Transformation');
%% Q2: FIR Digital filters using windowing
x=0:pi/1024:pi-pi/1024;
M=22;
w_c=0.3;
W1=rectwin(M);
b_rec=fir1(M-1,w_c,W1);
[H_rec,w]=freqz(b_rec,1,1024);
%M/2 Delay
H_rec=H_rec.*exp(sqrt(-1)*M/2.*w);
W2=hamming(M);
b_ham=fir1(M-1,w_c,W2);
[H_ham,w]=freqz(b_ham,1,1024);
H_ham=H_ham.*exp(sqrt(-1)*M/2.*w);
W3=hann(M);
b_han=fir1(M-1,w_c,W3);
[H_han,w]=freqz(b_han,1,1024);
H_han=H_han.*exp(sqrt(-1)*M/2.*w);
W4=blackman(M);
b_black=fir1(M-1,w_c,W4);
[H_black,w]=freqz(b_black,1,1024);
H_black=H_black.*exp(sqrt(-1)*M/2.*w);
W5=bartlett(M);
b_bar=fir1(M-1,w_c,W5);
[H_bar,w]=freqz(b_bar,1,1024);
H_bar=H_bar.*exp(sqrt(-1)*M/2.*w);
H_ideal=cat(1,ones(306,1),zeros(1024-306,1));
figure;
plot(x,abs(H_rec),x,abs(H_ham),x,abs(H_han),x,abs(H_black),x,abs(H_bar),x,H_ideal);
xlabel('\Omega [rad/sec]');
ylabel('Amplitude');
title('Digital FIR filters');
legend('rectangle windowing','hamming windowing','hann windowing','blackman windowing','barlett windowing','ideal filter');
figure;
plot(x,phase(H_rec),x,phase(H_ham),x,phase(H_han),x,phase(H_black),x,phase(H_bar),x,phase(H_ideal));
xlabel('\Omega [rad/sec]');
ylabel('Phase');
title('Phase of window filters');
legend('rectangle windowing','hamming windowing','hann windowing','blackman windowing','barlett windowing','ideal filter');
%find peak side lobe amplitude & aprox width odf main lobe
%wvtool(W1);
%wvtool(W2);
%wvtool(W3);
%wvtool(W4);
%wvtool(W5);

%find error
error=zeros(1,5);
error_rec=(abs(H_ideal-H_rec)).^2;
error_ham=(abs(H_ideal-H_ham)).^2;
error_han=(abs(H_ideal-H_han)).^2;
error_black=(abs(H_ideal-H_black)).^2;
error_bar=(abs(H_ideal-H_bar));
error(1)=trapz(error_rec);
error(2)=trapz(error_ham);
error(3)=trapz(error_han);
error(4)=trapz(error_black);
error(5)=trapz(error_bar);

%% Q3: Kaiser
w_k_s=0.46*pi;
w_k_p=0.3*pi;
delta_k_s=0.009;
delta_k_p=0.03;
A_k=-20*(log(delta_k_s)/log(10));
beta_k=0.5842*(A_k-21)^0.4+0.7886*(A_k-21);
M_k=29;
alpha_k=M_k/2;
w_k_c=(w_k_s+w_k_p)/2;
w_k=kaiser(M_k,beta_k);
figure
n=1:29;
plot(n,w_k)
title('Kaiser Window');
xlabel('n');
ylabel('Amplitude');
b_k=fir1(M_k-1,w_k_c/pi,w_k);
figure;
plot(n,b_k);
title('Parameters of Kaiser Window');
xlabel('n');
ylabel('Amplitude');
[H_k,w]=freqz(b_k,1,1024);
H_k=H_k.*exp(sqrt(-1)*M/2.*w);
figure;
plot(w,log(abs(H_k)));
title('Kaiser Window in dB');
xlabel('\omega [rad/sec]');
ylabel('amp dB');
figure;
plot(w,phase(H_k));
title('Kaiser Window-Phase');
xlabel('\omega [rad/sec]');
ylabel('phase');
E_k=zeros(1,1024);
for i=1:1024
if(w(i)<w_k_p)
    E_k(i)=1-abs(H_k(i));
end
 if (w(i)>w_k_s)
        E_k(i)=abs(H_k(i));
 end
end
    figure
    plot(w,E_k)
    title('Kaiser Window-Error');
xlabel('\omega [rad/sec]');
ylabel('Error');



