function [SE_MR,A] = functionComputeMonteCarloSE_UL(Hhat,H,tau_c,tau_p,nbrOfRealizations,N,K,M,p)
%Compute uplink SE for Cell-free mMIMO for the four different receiver
%cooperation levels, using either MR or MMSE/L-MMSE combining
%
%This function was developed as a part of the paper:
%
%Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
%Competitive With MMSE Processing and Centralized Implementation,"
%IEEE Transactions on Wireless Communications, To appear.
%
%Download article: https://arxiv.org/abs/1903.10611
%
%This is version 1.0 (Last edited: 2019-03-19)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.
%
%INPUT:
%Hhat              = Matrix with dimension N*L x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel from
%                    all BSs to UE k at channel realization n.
%H                 = Matrix with dimension N*L x nbrOfRealizations x K
%                    where (:,n,k) is the true collective channel from all
%                    BSs to UE k at channel realization n.
%R                 = Matrix with dimension N x N x L x K where
%                    (:,:,l,k) is the spatial correlation matrix between BS
%                    l and UE k in setup n, normalized by the noise power
%B                 = Matrix with dimension N x N x L x K where
%                    (:,:,l,k) is the spatial correlation matrix of the
%                    estimate between BS l and UE k in setup n, normalized
%                    by the noise power
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences and number of UEs per cell
%nbrOfRealizations = Number of channel realizations
%N                 = Number of antennas per AP
%K                 = Total number of UEs
%L                 = Number of APs
%p                 = Matrix K x 1 where element k is the uplink transmit
%                    power of UE k (If it is a scalar, the same value is
%                    used for all users)
%p1                = (Optional) Same as p but only for Level 1
%
%OUTPUT:
%SE_MR     = K x 4 matrix where the (k,n):th element is the uplink SE of 
%            UE k achieved with MR combining at cooperation level n
%SE_MMMSE  = Same as SE_MR but with MMSE or L-MMSE combining
%sumSE_SIC = Scalar with the sum SE achieved using MMSE-SIC combining



%If only one transmit power is provided, use the same for all the UEs
if length(p) == 1
   p = p*ones(K,1);
end

%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);
%Prepare to store simulation results
SE_MR = zeros(K,1);
%Diagonal matrix with transmit powers and its square root
Dp12 = diag(sqrt(p));


signal_MR = zeros(M,K);
scaling_MR = zeros(M,K);
Gp_MR = zeros(M,M,K);
A = zeros(M,M,K);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    
    
    %-----------------Levels 1-3
    gp_MR = zeros(M,K,K);
    
    
    %Go through all APs
    for m = 1:M
        
        %Extract channel realizations from all UEs to AP m
        Hallj = reshape(H(1+(m-1)*N:m*N,n,:),[N K]);
        
        %Extract channel estimate realizations from all UEs to AP m
        Hhatallj = reshape(Hhat(1+(m-1)*N:m*N,n,:),[N K]);

        %Compute MR combining
        V_MR = Hhatallj;
 
       for k=1:K
            v = V_MR(:,k); %Extract combining vector
            signal_MR(m,k) = signal_MR(m,k) + (v'*Hallj(:,k))/nbrOfRealizations; % (v_mk)'h_mk
            gp_MR(m,:,k) = gp_MR(m,:,k) + (v'*Hallj)*Dp12;                       % 
            scaling_MR(m,k) = scaling_MR(m,k) + norm(v).^2/nbrOfRealizations;    % ||v_mk||^2
       end
    end
for k = 1:K
        
        Gp_MR(:,:,k) = Gp_MR(:,:,k) + gp_MR(:,:,k)*gp_MR(:,:,k)'/nbrOfRealizations;
        
    end
    
end
    
%Compute the SE
for k = 1:K
    a=ones(M,1);
    %With MR combining
    b = signal_MR(:,k);
    A(:,:,k) = Gp_MR(:,:,k) + diag(scaling_MR(:,k))- p(k)*(b*b');
%     SE_MR(k) =  prelogFactor*real(log2(1+p(k)*(b'*b)/(trace(A(:,:,k)))));
        SE_MR(k) = prelogFactor*real(log2(1+p(k)*(b'*b)/(a'*A(:,:,k)*a))); 
end