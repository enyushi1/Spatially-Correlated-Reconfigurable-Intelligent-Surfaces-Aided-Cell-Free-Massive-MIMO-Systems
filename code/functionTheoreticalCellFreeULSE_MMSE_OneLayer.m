function [SE_CC] = functionTheoreticalCellFreeULSE_MMSE_OneLayer(HMean_Withoutphase,R_AP,M,K,N,tau_p,tau_c,p,Pset)     

%---This function is used to computes the theoretical uplink SE for
%phase-aware MMSE estimator (OneLayer).
%And each AP is equipped with N antennas.
%This is version 1.0 (Last edited: 2020-04-17)

% M = CellFreeParameter.M;
% K = CellFreeParameter.K;
% N = CellFreeParameter.N;
% tau_p = CellFreeParameter.tau_p;
% tau_c = CellFreeParameter.tau_c;

%INPUT:
%R_AP                 = Matrix with dimension N x N x M x K  where(:,:,m,k) is
%                       the spatial correlation matrix between AP m and UE k 
%                       in setup n, normalized by the noise power 
%                       normalized by the noise power
%Lk                   = Matrix with dimension MN x MN x K 
%  
%Phi                  = Matrix with dimension MN x MN x K
%                    
%Omega                = Matrix with dimension MN x MN x K
%
%A                    = Diagonal matrix with dimension M x M x K where (:,:,k)
%                       is the LSFD coefficients of UE k (when MMSE
%                       estimator is used.)
%M                   = Number of APs
%K                   = Number of UEs 
%p                   = 1xK vector, uplink power at each UE
%tau_p               = Pilot length
%tau_c               = Length of the coherence block
%Pset                = Pilot allocation set
%
%
%OUTPUT:
%
%SE_CC              = Vector with dimension K x 1 where (k) is the SE of UE k

SE_CC=zeros(K,1);
% for m = 1:M
%     for k = 1:K
%         
% 
%          Lk(:,:,m,k) = HMean_Withoutphase((m-1)*N+1:m*N,k)*(HMean_Withoutphase((m-1)*N+1:m*N,k))';
%          Dk(:,:,m,k) = diag(diag(R_AP(:,:,m,k)));
%          
%     end
% end
%  Rp = R_AP + Lk;

Phi = zeros(N,N,M,K);
Omega = zeros(N,N,M,K);
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

% PHI_MMSE=zeros(M,M,K);
% Psi = zeros(M,M,K);
% G=zeros(M,M,K);
% Sigma=zeros(M,M,K);
% for k=1:K
%     for m=1:M
%        PHI_MMSE(m,m,k)= trace(Phi(:,:,m,k));
%        Psi(m,m,k)=trace(Omega(:,:,m,k));
%        G(m,m,k)=HMean_Withoutphase((m-1)*N+1:m*N,k)'*(HMean_Withoutphase((m-1)*N+1:m*N,k));
%        Sigma(m,m,k)=trace(R_AP(:,:,m,k));
%     end
% end
G=zeros(N,N,M,K);
for m = 1:M
    for k = 1:K
         G(:,:,m,k) = HMean_Withoutphase((m-1)*N+1:m*N,k)*(HMean_Withoutphase((m-1)*N+1:m*N,k))';
    end
end

prelogFactor = (tau_c-tau_p)/(tau_c);

CCterm1 = zeros(K,1); %Store E{v^H_k h_k}= E{|v_k|^2}
CCterm3 = zeros(K,1); %Store E{v^H_k h_k}= E{|v_k|^2}
CCterm2_P1 = zeros(K,K);
Zk = zeros(M,M,K);
Ksi = zeros(M,M,K,K); %Ksi
X_p1 = zeros(M,M,K,K); %X(1);
X_p2 = zeros(M,M,K,K);
LK = zeros(M,M,K);
for k = 1:K
    for m = 1:M
       Zk(m,m,k) = trace(p(k)*tau_p*Omega(:,:,m,k)+G(:,:,m,k)); 
       for l=1:K
          Ksi(m,m,k,l) =p(k)*tau_p*trace(R_AP(:,:,m,k)*Omega(:,:,m,k))...
          +p(k)*tau_p*(HMean_Withoutphase((m-1)*N+1:m*N,l))'*Omega(:,:,m,k)*HMean_Withoutphase((m-1)*N+1:m*N,l)...
          +(HMean_Withoutphase((m-1)*N+1:m*N,k))'*R_AP(:,:,m,l)*HMean_Withoutphase((m-1)*N+1:m*N,k)...
          +abs((HMean_Withoutphase((m-1)*N+1:m*N,k))'*HMean_Withoutphase((m-1)*N+1:m*N,l))^2;
      
       if any(l==Pset(:,k)) && l~=k
           
         X_p1(m,m,k,l) = sqrt(p(k)*p(l))*tau_p*trace(R_AP(:,:,m,l)*R_AP(:,:,m,k)/Phi(:,:,m,k));
         X_p2(m,m,k,l) =2*sqrt(p(k)*p(l))*tau_p*real(trace(R_AP(:,:,m,l)*R_AP(:,:,m,k)/Phi(:,:,m,k))*(HMean_Withoutphase((m-1)*N+1:m*N,k))'*(HMean_Withoutphase((m-1)*N+1:m*N,l)));
       end
       end
       
      LK(m,m,k) = trace(G(:,:,m,k)); 
    end
end


for k = 1:K
    
    CCterm1(k) = trace(Zk(:,:,k));
    CCterm3(k) = trace(Zk(:,:,k));
     for l=1:K  %Non-coherent interference (i=k')
    
        CCterm2_P1(k,l) = p(l)*trace(Ksi(:,:,k,l));
        
        if any(l==Pset(:,k)) && l~=k
            
            CCterm2_P1(k,l)=  CCterm2_P1(k,l)+p(l)*(abs(trace(X_p1(:,:,k,l)))^2+abs(trace(X_p2(:,:,k,l)))); 
            
        end
         if l == k
            
            CCterm2_P1(k,l)= CCterm2_P1(k,l) -  p(k)*trace(LK(:,:,k)*LK(:,:,k));
        
        end
        
    end
    
end
CCterm2 = sum(CCterm2_P1,2);
for k=1:K
    
    SE_CC(k) = prelogFactor*real(log2(1+(p(k)*abs(CCterm1(k))^2)/((CCterm2(k))+CCterm3(k))));
    
end
% for k=1:K
%     CCterm1(k)=p(k)*abs(p(k)*tau_p*trace(Psi(:,:,k))+trace(G(:,:,k)))^2;
%     CCterm3(k)=p(k)*tau_p*trace(Psi(:,:,k))+trace(G(:,:,k));
% end
% for k=1:K
%     for l=1:K
%         CCterm2_p1(k,l)=p(l)*(p(k)*tau_p*trace(Psi(:,:,k)*Sigma(:,:,k))+p(k)*tau_p*trace(Psi(:,:,k)*G(:,:,l))...
%             +trace(Sigma(:,:,l)*G(:,:,k))+trace(G(:,:,k)*G(:,:,l)));
%     if  any(l==Pset(:,k)) && l~=k 
%         CCterm2_p1(k,l)=CCterm2_p1(k,l)+p(l)*p(k)*tau_p*tau_p*trace(Sigma(:,:,k))*trace(Sigma(:,:,l))...
%             +2*sqrt(p(l)*p(k))*tau_p*real(trace(Sigma(:,:,k)*Sigma(:,:,l)/PHI_MMSE(:,:,k))*HMean_Withoutphase(:,k)'*HMean_Withoutphase(:,l));
%     end
%         if l==k
%             CCterm2_p1(k,l)= CCterm2_p1(k,l) - p(k)*trace(G(:,:,k)*G(:,:,k));
%         end
%     end
% end
%     
%   CCterm2=sum(CCterm2_p1,2);  
% for k=1:K
%     SE_CC(k)= prelogFactor*log2(1+ CCterm1(k)/(abs(CCterm2(k)) +  CCterm3(k) )  )  ;
% end


end    

    
    
    
    








