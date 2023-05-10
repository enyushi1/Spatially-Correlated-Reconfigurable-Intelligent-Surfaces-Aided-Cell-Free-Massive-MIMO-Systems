function [Hhat_MMSE] = functionChannelEstimates_MMSE(R_AP,RHat_AP,HMean,H,nbrOfRealizations,M,K,N,tau_p,p,Pset)

%MMSE channel estimator for Cell-Free setup. The estimation is locally
%performed at the APs.
%And each AP is equipped with N antennas.
%This is version 1.0 (Last edited: 2020-04-17)

%INPUT:
%R_AP                 = Matrix with dimension N x N x M x K x nbrOfRealzations where (:,:,m,k,i) is
%                       the spatial correlation matrix between AP l and UE k 
%                       in i^th realization, normalized by the noise power
%HMean                = Matrix with dimension MN  x K x nbrOfRealzations
%                       where (mn,i,k) is the i^th realization of the channel mean
%                       between the n^th antenna of AP m and UE k (with phase shift)                  
%H                    = Matrix with dimension MN x K x nbrOfRealzations 
%                       where (mn,i,k) is the i^th channel realization
%                       between the n^th antenna of AP m and UE k
%nbrOfRealizations    = Number of realizations
%M                    = Number of APs
%K                    = Number of UEs 
%N                    = Number of AP antennas
%p                    = 1xK vector, uplink power at each UE
%tau_p                = Pilot length
%Pset                 = Pilot allocation set
%
%OUTPUT:
%Hhat_MMSE            = Matrix with dimension MN x nbrOfRealzations x K
%                       where (mn,k,i) is the i^th  realization of phase-aware
%                       MMSe channel estimate between the n^th antenna of AP m and UE k
%                     

% M = CellFreeParameter.M;
% K = CellFreeParameter.K;
% N = CellFreeParameter.N;
% tau_p = CellFreeParameter.tau_p;



%Prepare to store MMSE channel estimates
Hhat_MMSE = zeros(M*N,nbrOfRealizations,K);

%Store identity matrix of size N x N
eyeN = eye(N);

%Generate realizations of normalized noise 
Np = sqrt(0.5)*(randn(N,nbrOfRealizations,M,K) + 1i*randn(N,nbrOfRealizations,M,K));


for m = 1:M
    for k = 1:K
        
        yp = zeros(N,nbrOfRealizations);
        yMean = zeros(N,nbrOfRealizations);
        PsiInv = zeros(N,N);
        inds = Pset(:,k); 
        
        for z = 1:length(inds)
            
            yp = yp + sqrt(p(inds(z)))*tau_p*H((m-1)*N+1:m*N,:,inds(z));
            yMean = yMean + sqrt(p(inds(z)))*tau_p*HMean((m-1)*N+1:m*N,:,inds(z));
            PsiInv = PsiInv + p(inds(z))*tau_p*RHat_AP(:,:,m,inds(z));
            
        end
        yp = yp + sqrt(tau_p)*Np(:,:,m,k);
        PsiInv = PsiInv + eyeN;
        
      
        for z = 1:length(inds)
            
            RPsi = RHat_AP(:,:,m,inds(z))/PsiInv;
            Hhat_MMSE((m-1)*N+1:m*N,:,inds(z)) = HMean((m-1)*N+1:m*N,:,inds(z)) + sqrt(p(inds(z)))*RPsi*(yp-yMean);
            
        end
    end
end
