function [coherentx,nonCoherentx] = functionMMSE_interferenceLevels( R_AP,HMean_Withoutphase,A_singleLayer,M,K,N,tau_p,p,Pset)
  
%This function is used to Check the levels of coherent and non-coherent interference levels (used
%for pilot allocation)
%And each AP is equipped with N antennas.
%This is version 1.0 (Last edited: 2020-04-10)

%INPUT:
%R_AP                 = Matrix with dimension N x N x M x K  where(:,:,m,k) is
%                       the spatial correlation matrix between AP m and UE k 
%                       in setup n, normalized by the noise power
%HMean                = Matrix with dimension MN x K ,where (mn,k) is the
%                       channel mean between the n^th antenna of AP m and UE k, normalized by
%                       noise power and with random phase shifts
%HMeanWithoutPhase    = Matrix with dimension MN x K ,where (mn,k) is the
%                       channel mean between the n^th antenna of AP m and UE k, normalized by
%                       noise power and without random phase shifts
%A                    = large-scale fading coefficients. In this setup,LSFD
%                       coefficients are set to 1 for single layer decoding                 
%M                    = Number of APs
%K                    = Number of UEs 
%p                    = Column vector with dimension K x 1, where p(k) is
%                       the uplink transmit power of UE k
%tau_p                = Pilot length 
%Pset                 = Pilot set
%N                    = Number of antennas per AP
%
%OUTPUT:
%
%coherentx            = Coherent interference of UE k
%                     
%nonCoherentx         = Non-Coherent interference of UE k

% M = CellFreeParameter.M;
% K = CellFreeParameter.K;
% N = CellFreeParameter.N;
% tau_p = CellFreeParameter.tau_p;

%Prepare to store the results
Ksi = zeros(M,M,K,K); %Ksi
X_p1 = zeros(M,M,K,K); %X(1)
% X_p2 = zeros(M,M,K,K); %X(2)
nonCoherent = zeros(K,K);
coherent = zeros(K,K);



%Generate matrix used in this setup

Lk = zeros(N,N,M,K);
Phi = zeros(N,N,M,K);
Omega = zeros(N,N,M,K);

for m = 1:M
    for k = 1:K
       
         Lk(:,:,m,k) = HMean_Withoutphase((m-1)*N+1:m*N,k)*(HMean_Withoutphase((m-1)*N+1:m*N,k))';
         
    end
end

%Go through all APs          
for m = 1:M
    
    %Go through all UEs
    for k = 1:K
        
        %Compute the UEs indexes that use the same pilot as UE k
        inds = Pset(:,k);
        PsiInv = zeros(N,N);
        
        %Go through all UEs that use the same pilot as UE k 
        for z = 1:length(inds)   
            
            PsiInv = PsiInv + p(inds(z))*tau_p*R_AP(:,:,m,inds(z));
        
        end
            PsiInv = PsiInv + eye(N);
            
            for z = 1:length(inds)
                
                Phi(:,:,m,inds(z)) = PsiInv;
            
            end
            
            Omega(:,:,m,k) = R_AP(:,:,m,k)/PsiInv*R_AP(:,:,m,k);
            
    end
end



% Go through all APs
for m = 1:M
    
    %Go through all UEs
    for k = 1:K
        
        for l=1:K  %Non-coherent interference (i=k')
            
            Ksi(m,m,k,l) = p(k)*tau_p*trace(R_AP(:,:,m,l)*Omega(:,:,m,k))+...
                trace(Lk(:,:,m,k)*R_AP(:,:,m,l)) + p(k)*tau_p*trace(Lk(:,:,m,l)*Omega(:,:,m,k))+...
                trace(Lk(:,:,m,k))*trace(Lk(:,:,m,k));
            
           
           if any(l==Pset(:,k)) %Coherent interference (If there is pilot contamination)
           
           X_p1(m,m,k,l) = sqrt(p(k)*p(l))*tau_p*trace(R_AP(:,:,m,l)/Phi(:,:,m,k)*R_AP(:,:,m,k));
           

           end
           
        end
        
    end
    
end

%Go through all UEs
for k = 1:K
    
    for l=1:K  %Non-coherent interference (i=k')
    
        nonCoherent(k,l) = p(l)*trace(A_singleLayer(:,:,k)'*Ksi(:,:,k,l)*A_singleLayer(:,:,k));
        
        if any(l==Pset(:,k)) 
            
            coherent(k,l)=  p(l)*abs(trace(A_singleLayer(:,:,k)*X_p1(:,:,k,l)))^2; 
            
        end
        
    end
    
end

coherentx = sum(coherent,2);
nonCoherentx = sum(nonCoherent,2);
              
end
    
