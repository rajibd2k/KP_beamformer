%% Evaluate the BROADBAND Metrics of Beamformer

function Perf = Metrics(h, d, phiX1, phiV1, phiV)

phiX1 = abs( phiX1 ) ;
phiV1 = abs( phiV1 ) ;

num_f = size(phiX1,1) ;
f = [0 : (num_f - 1)]' / num_f ;
df = f(2) - f(1) ;

phiX1_fd = phiX1 .* abs(sum(conj(h).*d,1)').^2 ;
phiE_fd = phiX1 .* (abs(sum(conj(h).*d,1)'- 1).^2) ;

phiV1_rn = zeros( num_f , 1 ) ;
for idxf=1:num_f
    phiV1_rn(idxf) = abs( h(:,idxf)'*phiV(:,:,idxf)*h(:,idxf) ) ;
end

% input SINR
iSINR = nansum( phiX1 ) / nansum( phiV1 )  ;
iSINR_dB = 10*log10(iSINR) ;
Perf.iSINR_dB = iSINR_dB ;

% output SINR
oSINR = nansum( phiX1_fd ) / nansum( phiV1_rn ) ;
oSINR_dB = 10*log10(oSINR) ;
Perf.oSINR_dB = oSINR_dB ;

% Gain
Perf.Gain_dB = oSINR_dB - iSINR_dB ;

% Noise Rejection factor
xi_n = nansum( phiV1 ) / nansum( phiV1_rn ) ;
xi_n_dB = 10*log10(xi_n) ;
Perf.xi_n_dB = xi_n_dB ;

% Desired-Signal Reduction factor
xi_d = nansum( phiX1 ) / nansum( phiX1_fd ) ;
xi_d_dB = 10*log10(xi_d) ;
Perf.xi_d_dB = xi_d_dB ;

% Desired-Signal Distortion Index
nu = nansum( phiE_fd ) / nansum( phiX1 ) ;
nu_dB = 10*log10(nu) ;
Perf.nu_dB = nu_dB ;

% Mean-Squared Error
J = nansum(  phiE_fd + phiV1_rn  ) * df  ;
J_dB = 10*log10(J) ;
Perf.J_dB = J_dB ;

end