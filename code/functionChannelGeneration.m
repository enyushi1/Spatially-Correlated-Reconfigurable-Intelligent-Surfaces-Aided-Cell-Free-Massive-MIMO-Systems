function [H,HMean] = functionChannelGeneration( R_AP,HMean_Withoutphase,M,K,N,nbrOfRealizations)

%---This function is used to generate the channel realizations for given covariance and
%mean matrices. The outputs are channel realizations and channel means with
%random phase shifts at a coherence block.
%And each AP is equipped with N antennas.
%This is version 1.0 (Last edited: 2020-04-14)

%INPUT:
%channelGain_LoS       = Matrix with dimension M x K  where(m,k) is
%                       LoS channel gain between AP m and UE k 
%                       in setup n, normalized by the noise power
%channelGain_NLoS      = Matrix with dimension M x K ,where (m,k) is the
%                        NLos channel mean between the AP m and UE k, normalized by
%                        noise power
%M                     = Number of APs
%nbrOfRealizations     = Number of realizations
%K                     = Number of UEs 
%N                     = Number of antennas per AP
%
%OUTPUT:
%
%R_AP                  = Matrix with dimension N x N x M x K x nbrOfRealizations where (:,:,m,i,k) is 
%                        the spatial correlation matrix between AP m and UE k in i^th channel realization, normalized by noise
%                        power
%H                     = Matrix with dimension MN x K x nbrOfRealzations
%                        where (mn,k,i) is the i^th channel realization
%                        between the n^th antenna of AP m and UE k
%                     
%HMean                 = Matrix with dimension MN x K x nbrOfRealzations 
%                        where (mn,k,i) is the i^th realization of the channel mean
%                        between the n^th antenna of AP m and UE k

% M = CellFreeParameter.M;
% K = CellFreeParameter.K;
% N = CellFreeParameter.N;



%Prepare to store the results   
R = zeros(M*N,M*N,K);
W = (randn(M*N,nbrOfRealizations,K)+1i*randn(M*N,nbrOfRealizations,K));
H = zeros(M*N,nbrOfRealizations,K);

%Same channelGain_LoS and channelGain_NLoS for all realizations (at each setup) but phases shift of LoS are different at each
%coherence block
HMean = zeros(M*N,nbrOfRealizations,K); 
HMeanx = reshape(repmat(HMean_Withoutphase,nbrOfRealizations,1),M*N,nbrOfRealizations,K); 

%Create phase shifts and store them in a matrix
%uniformly distributed random variables
% angles = -pi + 2*pi*rand(M*N,nbrOfRealizations,K);
% phaseMatrix = exp(1i*angles);

%--Phase shift same for all antennas
angles = -pi + 2*pi*rand(M,nbrOfRealizations,K);
phaseMatrix = exp(1i*angles*0);
v_kron = ones(N,1);
PhaseMatrix = zeros(M*N,nbrOfRealizations,K);

for k = 1:K
    PhaseMatrix(:,:,k) = kron(phaseMatrix(:,:,k),v_kron);
end
    
    


for m = 1:M
    for k = 1:K
        
        R((m-1)*N+1:m*N,(m-1)*N+1:m*N,k) = R_AP(:,:,m,k);
        
    end
end


%Go through all UEs and apply the channel gains to the spatial
%correlation and mean matrices and introduce the phase shifts 
for k = 1:K
    
%     HMean(:,:,k)= phaseMatrix(:,:,k).*HMeanx(:,:,k);
    HMean(:,:,k)= HMeanx(:,:,k);
    Rsqrt = sqrtm(R(:,:,k));
    H(:,:,k) = sqrt(0.5)*Rsqrt*W(:,:,k) + HMean(:,:,k);
       
end








