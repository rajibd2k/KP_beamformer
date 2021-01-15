%% Calculate PowerPattern (PP) for polar plot
%% f : digital frequency [0, 0.5] 
%% h_f : filter at frequency 'f' 
%% theta_d : angle (in degrees) of the source signal impinging on the ULA
%% delta = 10^(-2) ; % 1 cm (inter-sensor-distance)
%% c = 340 ; % m/s
%% Ts = 1/8000 ; % s
%% PP_f : Square of magnitude of BeamPattern values, at frequency 'f', for all angles [0, 360]

function [phi_rad, PP_f, PP_f_norm] = PP(f, h_f , delta, c, Ts)

M = length(h_f) ;

% Powerpattern
phi=(0:1:180)'; % phi=-90:1:90 ;
phi_rad=phi/180*pi;
[phi_mat,m_mat]=meshgrid(phi_rad,0:M-1);
D_phi=exp(-1i*2*pi*f*(delta/c/Ts)*m_mat.*cos(phi_mat));
PP_f=D_phi'*h_f;
PP_f=abs(PP_f);
PP_f_norm = PP_f / max(PP_f) ;

end
