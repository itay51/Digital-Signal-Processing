%Q1
clear all
clc
b=1;
%Subsection 1
N=256;
a0=[1, -.5, .5 ];
a1=[1,0,-.5];
a=conv(a0,a1);
[h,w]=freqz(b,a,N);
spec=h.*conj(h);
p=4;
emse=0;
for l=1:1000
    u=randn(1,N);
    x=filter(b,a,u);
    [ar_spec,e]=aryule(x,p);
    estimated_h=freqz(1,ar_spec,N);
    estimated_s_spec=e.*estimated_h.*conj(estimated_h);
    emse=emse+(estimated_s_spec-spec).^2;
end
emse=emse./1000;
figure;
plot(w,spec,w,real(estimated_s_spec));
xlim ([0 pi])
title 'Yule Walker - p=4'
xlabel('\omega')
ylabel('Spectrum')
legend('Real','Yule walker')
%subsection 2
figure;
plot(w,10.*log(abs(emse))./log(10))
xlim ([0 pi])
title 'EMSE(\omega)- p=4'
xlabel('\omega')
ylabel('10 log(EMSE)[dB]')
%subsection 2
figure;
hold on;
plot(w,10.*log(abs(emse))./log(10));
for p=[8,16,64]
    emse=0;
    for l=1:1000
       u=randn(1,N);
       x=filter(b,a,u);
        [ar_spec,e]=aryule(x,p);
        estimated_h=freqz(1,ar_spec,N);
        estimated_s_spec=e.*estimated_h.*conj(estimated_h);
        emse=emse+(estimated_s_spec-spec).^2;
    end
    emse=emse./1000; 
    plot(w,10.*log(abs(emse))./log(10))    
end
xlim ([0 pi])
title 'EMSE(\omega)'
xlabel('\omega')
ylabel('10 log(EMSE)[dB]')
legend ('p=4','p=8','p=16','p=64')
hold off
%Q2
clear all
clc
a=[1,-1.3435,0.9025];
b=[1,1.3435,0.9025];
%Subsection 1+2
for N=[128,512,2048,2048]
    [h,w]=freqz(b,a,N);
    spec=h.*conj(h);
    u=randn(1,N);
    x=filter(b,a,u);
    r=zeros(1,N);
    s_per=zeros(N,1);
    s_check=zeros(N,1);
    for k=0:N-1
       r(k+1)=cat(2,zeros(1,k),x)*(cat(2,x,zeros(1,k)))';
       r(k+1)=r(k+1)/N;
       s_per=s_per+((x(k+1).*exp(-1i*k.*w))./sqrt(N));
    end
    s_per=s_per.*conj(s_per);
    figure
    plot(w,10.*log(abs(s_per))./log(10));
    xlim ([0 pi])
    title(['Periodogram Estimation, N=',num2str(N)])
    xlabel ('\omega [rad]')
    ylabel ('Spectrum')
end
plot(w,10.*log(abs(spec))./log(10));
xlim ([0 pi])
title 'Real Spectrum'
xlabel ('\omega [rad]')
ylabel 'Spectrum'
%Subsection 3- bartlett
N=2048;
for K=[32,128,256]
    L=2*K+1;
    w_b=bartlett(L);
    S_bartlett=zeros(N,1);
    for k=0:K
       if (k>=1)
         S_bartlett=S_bartlett+w_b(k+K+1).*r(k+1).*(exp(-k.*1i.*w)+exp(k.*1i.*w));
       else
         S_bartlett=r(k+1);
       end  
    end
    figure
    plot(w,10.*log(abs(S_bartlett))./log(10));
    xlim ([0 pi])
    title(['Bartlett Estimation, K=',num2str(K)])
    xlabel ('\omega [rad]')
    ylabel 'Spectrum'
end
figure
plot(w,10.*log(abs(spec))./log(10));
xlim ([0 pi])
title 'Real Spectrum'
xlabel ('\omega [rad]')
ylabel ('Spectrum')





