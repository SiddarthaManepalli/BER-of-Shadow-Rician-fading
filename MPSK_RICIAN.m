function pe = MPSK_RICIAN(M)
% ========================================================================================================================
% BIT ERROR PROBABILITY OF M-PSK OVER DOUBLE SHADOWED RICIAN FADING CHANNEL SUBJECTED TO DOUBLE GATED NOISE
% The Bit Error Probability of M-PSK scheme is obtained by a formula where
% Error Function erfc can be obtained by 2-(2*CDF(Additive Noise)) then by
% Multiplying it with weights W(i,k,M) followed by summation with i and k
% Hence CDF of Additive noise is obtained and called into the main function 


% Probability of Error is validated with respect to
snr_db = linspace(0,30,30);

% BER of M-PSK scheme

final_count  = 0;
for K        = 1:log2(M/2)
    
    count    = 0;
    i        = 0:1:(((1-2.^(-K))*(M/2))-1);
   for i_len = 1:length(i)
      
      % Weights w(i,k,M) is defined as
      W     = ((-1).^(floor((i(i_len).*(2.^(K-1)))/(M/2)))).*((2.^(K-1))-(floor(((i(i_len).*(2.^(K-1)))/(M/2))+0.5)));
      % a(i,M) is defined as
      a     = (4.*(((2.*i(i_len))+1).^2).*log2(M).*((sin(((2.*i(i_len))+1)*(pi/M))).^2));
      % Error Function for M-PSK scheme is defined as
      erfc  =  2-(2.*CDF(a));
      % Validating inner Summation
      count =  count + (W.*erfc);
   end % for i_len
   
    % Validating Outer Summation 
    final_count = final_count + count;
end % for K

% PROBABILITY OF ERROR OVER M-PSK SCHEME
pe = abs((1/((M/2).*log2(sqrt(M/2)))).*final_count);
  
% ----------------- Plotting BER vs SNR ------------------------------------
figure(2), clf
semilogy(snr_db,pe(1,:),'->b',...
         snr_db,pe(2,:),'-->r',...
         snr_db,pe(3,:),'-.>g',...
         snr_db,pe(4,:),':>k',...
         'LineWidth',2);

xlim([0 25])
set(gcf,'color','white')
xlabel('Signal to Noise Ratio (dB)','FontSize',16)
ylabel('Bit Error Probability (Pe)','FontSize',16)
title('BEP of M-PSK over Double Shadowed Rician Fading Channel subjected to Double Gated Noise','FontSize',16)
legend({'SNI(dB)=10, md=2.3, k=2.4',...
        'SNI(dB)=15, md=2.3, k=2.4',...
        'SNI(dB)=20, md=2.3, k=2.4',...
        'SNI(dB)=25, md=2.3, k=2.4'},...
        'FontSize',16,'Location','southwest');
    
legend('boxoff')
% ---------------------- End of Plotting -----------------------------------

% =========================================================================================================================


function final = CDF(A)
% =========================================================================================================================
% CUMMULATIVE DISTRIBUTION FUNCTION OF ADDITIVE NOISE F(r)
% The CDF of Additive noise is obtained by normalizing and finding
% Cummulative summation of PDF of Additive Noise 
% PDF of Additive Noise is obtained by integrating the product of 
% The PDF's of Double Shadowed Rician fading model and G2AWGN noise model 


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
sni_dB = [10 15 20 25];
snr_dB = linspace(0,30,30);

% Conversion of dB to normal
sni = 10.^(sni_dB/10);
snr = 10.^(snr_dB/10);

% PDF of Additive Noise
 pre_final = zeros(length(snr),1);
 final     = zeros(length(sni),length(snr));

for sni_len = 1:length(sni)
    
    for snr_len = 1:length(snr)
        
        % SNR and SNI in terms of Ng and Ni respectively
        Ng = 2*(snr(snr_len)^2);
        Ni = 2*(sni(sni_len)^2);
        
        % Substituting a(i,M)*Eb in r
        % Eb = Ng*SNR or Eb = Ni*SNI
        r = A.*Ng.*snr(snr_len);
    
        % 1st part of equation
        term1 = alpha*beta*p1*p2/sqrt(pi*(Ng+Ni));
        term2 = ((md/(md+k))^md) * ((1+k)/x_rms) * gamma(1.5);
        term3 = ((r.^2/(Ng+Ni))+((1+k)/x_rms)).^-(3/2);
        term4 = hypergeom([md,1.5],1,(k/(md+k)) * ((1+k)/x_rms) / ((r.^2/(Ng+Ni))+((1+k)/x_rms)));
    
        head  = term1*term2*term3*term4;
    
        % 2nd part of equation
        term5 = (1-alpha*beta*p1*p2)/sqrt(pi*(Ng));
        term6 = ((md/(md+k))^md) * ((1+k)/x_rms) * gamma(1.5);
        term7 = ((r.^2/(Ng))+((1+k)/x_rms)).^-(3/2);
        term8 = hypergeom([md,1.5],1,(k/(md+k)) * ((1+k)/x_rms) / ((r.^2/(Ng))+((1+k)/x_rms)));
    
        tail  = term5*term6*term7*term8;
    
    
        % Adding 2 parts of equation
        pre_final(snr_len,:) = head + tail;
    
        % Calculating CDF of Additive Noise
        P_pdf = pre_final./sum(pre_final);
        P_cdf = cumsum(P_pdf);
    
    end % for snr_len
    final(sni_len,:) = P_cdf;
    
end     % for sni_len

% ---------------- Plotting CDF of Additive Noise --------------------------
figure(1), clf
semilogy(snr_dB,final(1,:),'->b',...
         snr_dB,final(2,:),'-->r',...
         snr_dB,final(3,:),'-.>g',...
         snr_dB,final(4,:),':>k',...
         'LineWidth',2);

xlim([0 5])
set(gcf,'color','white')
xlabel('Additive Noise','FontSize',16)
title('CDF OF ADDITIVE NOISE','FontSize',16)
legend({'SNI(dB)=10, md=2.3, k=2.4',...
        'SNI(dB)=15, md=2.3, k=2.4',...
        'SNI(dB)=20, md=2.3, k=2.4',...
        'SNI(dB)=25, md=2.3, k=2.4'},...
        'FontSize',16,'Location','southwest');
    
legend('boxoff')
% ------------------ End of Plotting ---------------------------------------
end     % function sub

% =========================================================================================================================


end     % function main
