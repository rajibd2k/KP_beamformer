%% MVDR (Minimum Variance Distortionless Response) Filter Beamforming 
%% est_d : Steering vector of the SOI at each frequency
%% est_phiY : estimated covariance matrix of the data, at each frequency

function [h] = MVDR_error( est_phiY, est_d )

num_f = size( est_phiY, 3 ) ;
f = [0 : (num_f - 1)]' / num_f ;

M = size( est_phiY , 1 ) ;

h=zeros(M,length(f));

for idxf=1:length(f)
%     h(:,idxf) = inv( est_phiY(:,:,idxf) )*est_d(:,idxf)/(est_d(:,idxf)'*inv( est_phiY(:,:,idxf) )*est_d(:,idxf));
    try
        h(:,idxf) = pinv( est_phiY(:,:,idxf) )*est_d(:,idxf)/(est_d(:,idxf)'*pinv( est_phiY(:,:,idxf) )*est_d(:,idxf));
    catch
        h(:,idxf) = inv( est_phiY(:,:,idxf) )*est_d(:,idxf)/(est_d(:,idxf)'*inv( est_phiY(:,:,idxf) )*est_d(:,idxf));
    end
end

end
