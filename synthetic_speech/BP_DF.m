%% Calculate BeamPattern (BP) and Directivity factor (DF)
%% f : digital frequency [0, 0.5] 
%% h_f : filter at frequency 'f' 
%% theta_d : angle (in degrees) of the source signal impinging on the ULA
%% BP_f_theta_d (dB): BeamPattern value at frequency 'f' and angle 'theta_d'
%% BP_f : BeamPattern values at frequency 'f' for all angles [0, 180]
%% DF_f (dB): Directivity Factor value at frequency 'f' 

function [BP_f_theta_d, BP_f, DF_f] = BP_DF(f, h_f , theta_d)

M = length(h_f) ;

theta_d = theta_d/180*pi ; 

% Beampattern
phi=(0:1:180)';
phi_rad=phi/180*pi;
[phi_mat,m_mat]=meshgrid(phi_rad,0:M-1);
D_phi=exp(-1i*2*pi*f*m_mat.*cos(phi_mat));
BP_f=D_phi'*h_f;
BP_f=abs(BP_f);
BP_f = BP_f / max(BP_f) ;
BP_f=20*log10(BP_f);

m_mat=[0:M-1]' ;
d=exp(-1i*2*pi*f*m_mat.*cos(theta_d));

% Beampattern at (f ,theta_d)
BP_f_theta_d=d'*h_f;
BP_f_theta_d=abs(BP_f_theta_d);
BP_f_theta_d=20*log10(BP_f_theta_d);

% Directivity Factor
[j_mat,i_mat]=meshgrid(1:M,1:M);
Gamma0=sinc( 2 * f * (j_mat-i_mat) );
tmp_num = ( abs( h_f'*d ) )^2 ;
tmp_den = abs( h_f'*Gamma0*h_f ) ;
DF_f = 10*log10( tmp_num / tmp_den ) ;

end
