clc; 
clear all;
close all;

% Nakagami-m Shape Parameter
md    = 2.3;

% Rician k Parameter
k     = 2.4;

% Square of RMS power of Rician Parameter (x)
x_rms = 3.4;

% AWGN variables
alpha = 0.5;
beta  = 0.5;
p1    = 0.5;
p2    = 0.5;

% Considering SNR and SNI values in dB
sni_dB = 10;
snr_dB = 10;

% Conversion of dB to normal
sni = 10^(sni_dB/10);
snr = 10^(snr_dB/10);

% SNR and SNI in terms of Ng and Ni respectively
Ng = 2*(snr^2);
Ni = 2*(sni^2);

% Additive Noise Parameter
r = linspace(0,100,10);

% PDF of Additive Noise
final = zeros(length(r),1);

for r_len = 1:length(r)
    
    % 1st part of equation
    term1 = alpha*beta*p1*p2/sqrt(pi*(Ng+Ni));
    term2 = ((md/(md+k))^md) * ((1+k)/x_rms) * gamma(1.5);
    term3 = ((r(r_len).^2/(Ng+Ni))+((1+k)/x_rms)).^-(3/2);
    term4 = hypergeom([md,1.5],1,(k/(md+k)) * ((1+k)/x_rms) / ((r(r_len).^2/(Ng+Ni))+((1+k)/x_rms)));
    
    head  = term1*term2*term3*term4;
    
    % 2nd part of equation
    term5 = (1-alpha*beta*p1*p2)/sqrt(pi*(Ng));
    term6 = ((md/(md+k))^md) * ((1+k)/x_rms) * gamma(1.5);
    term7 = ((r(r_len).^2/(Ng))+((1+k)/x_rms)).^-(3/2);
    term8 = hypergeom([md,1.5],1,(k/(md+k)) * ((1+k)/x_rms) / ((r(r_len).^2/(Ng))+((1+k)/x_rms)));
    
    tail  = term5*term6*term7*term8;
    
    % Adding 2 parts of equation
    final(r_len,:) = head + tail;
    
end

% Calculating CDF of Additive Noise
P_pdf = final./sum(final);
P_cdf = cumsum(P_pdf);

% Plotting CDF of Additive Noise
figure(1), clf
plot(r,P_cdf,'k','linew',2)
xlabel('Additive Noise','FontSize',16)
title('CDF OF ADDITIVE NOISE','FontSize',16)




