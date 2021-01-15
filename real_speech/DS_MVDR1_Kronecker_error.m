%% Combined DS (Delay Sum) + MVDR (Minimum Variance Distortionless Response) Kronecker Filter Beamforming
%% est_phiY : estimated covariance matrix of the data, at each frequency
%% d1, d2 : Steering vectors of the SOI at each frequency, using the virtual ULAs
%% n_iter : number of iterations

function [h] = DS_MVDR1_Kronecker_error( est_phiY, d_1, d_2, n_iter )

num_f = size( est_phiY, 3 ) ;
f = [0 : (num_f - 1)]' / num_f ;

M1 = size( d_1 , 1 ) ; 
I_M1 = eye(M1) ;
h_1 = zeros(M1,length(f)) ;

M2 = size( d_2 , 1 ) ; 
I_M2 = eye(M2) ;
h_2 = zeros(M2,length(f)) ;

phiY_2 = est_phiY([1:M2],[1:M2],:) ;
for idxf=1:length(f)
%     h_2(:,idxf) = inv( phiY_2(:,:,idxf) )*d_2(:,idxf)/(d_2(:,idxf)'*inv( phiY_2(:,:,idxf) )*d_2(:,idxf));
    try
        h_2(:,idxf) = pinv( phiY_2(:,:,idxf) )*d_2(:,idxf)/(d_2(:,idxf)'*pinv( phiY_2(:,:,idxf) )*d_2(:,idxf));
    catch
        h_2(:,idxf) = inv( phiY_2(:,:,idxf) )*d_2(:,idxf)/(d_2(:,idxf)'*inv( phiY_2(:,:,idxf) )*d_2(:,idxf));
    end
end

%phiY_2 = zeros(M1,M1,length(f)) ;
phiY_1 = zeros(M2,M2,length(f)) ;

if n_iter == 0
    for idxf=1:length(f)
        h_1(:,idxf)= d_1(:,idxf) / M1 ;
        h(:,idxf) = kron( h_1(:,idxf),h_2(:,idxf) ) ;
    end
    return ;
end

for idx_iter = 1 : n_iter
    for idxf=1:length(f)

        h_1(:,idxf)= d_1(:,idxf) / M1 ;
        
        phiY_1(:,:,idxf) = kron( h_1(:,idxf),I_M2 )' * est_phiY(:,:,idxf) * kron( h_1(:,idxf),I_M2 );
%         h_2(:,idxf) = inv( phiY_1(:,:,idxf) )*d_2(:,idxf)/(d_2(:,idxf)'*inv( phiY_1(:,:,idxf) )*d_2(:,idxf));
        try
            h_2(:,idxf) = pinv( phiY_1(:,:,idxf) )*d_2(:,idxf)/(d_2(:,idxf)'*pinv( phiY_1(:,:,idxf) )*d_2(:,idxf));
        catch
            h_2(:,idxf) = inv( phiY_1(:,:,idxf) )*d_2(:,idxf)/(d_2(:,idxf)'*inv( phiY_1(:,:,idxf) )*d_2(:,idxf));
        end
                
        h(:,idxf) = kron( h_1(:,idxf),h_2(:,idxf) ) ;
    end
end

end
