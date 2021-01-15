%% Calculate PowerPattern (PP) for polar plot
%% f : digital frequencies within [0, 1) 
%% h : filters at different frequencies 
%% delta = 10^(-2) ; % 1 cm (inter-sensor-distance)
%% c = 340 ; % m/s
%% Ts = 1/8000 ; % s
%% PP_norm : Magnitude of BeamPattern values, averaged over all frequencies, for all angles [0, 180]

function [phi_rad, PP_norm, PP_f_norm] = PP(h , delta, c, Ts)

M = size(h,1) ;
num_f = size(h,2) ;
f = [0:(num_f-1)]' / num_f ;

% Powerpattern
phi = [0:1:180]'; % phi=-90:1:90 ;
phi_rad = phi/180*pi;
[phi_mat,m_mat] = meshgrid(phi_rad,[0:M-1]);

PP_f_norm = zeros( length(phi_rad) , num_f) ;

for idx_f = 1 : num_f

    D_phi = exp(-1i*2*pi* f(idx_f) *(delta/c/Ts)*m_mat.*cos(phi_mat));
    PP_f = D_phi' * h(:,idx_f);
    PP_f = abs(PP_f);
    PP_f = PP_f / max(PP_f) ;
    
    PP_f_norm(:,idx_f) = PP_f ; 

end

PP_f_norm( find(isnan(PP_f_norm))  ) = 0 ;
PP_f_norm( find(isinf(PP_f_norm))  ) = 1 ;

PP_norm = nanmean(PP_f_norm,2) ;
PP_norm = PP_norm / max(PP_norm) ;

end
