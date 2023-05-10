%This Matlab script can be used to generate Figure 3 in the article:
%
%? Özdogan, E. Björnson and J. Zhang, "Performance of Cell-Free Massive MIMO with Rician Fading 
%and Phase Shifts," in IEEE Transactions on Wireless Communications.
%doi: 10.1109/TWC.2019.2935434
%
%Download article: https://arxiv.org/abs/1903.07335
%
%This is version 1.0 (Last edited: 2019-11-01)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.


%Empty workspace and close figures
clear
close all
clc

%Define the range of number of Access Points (APs)
M=20:20:100;
NN= [4 9 16];
%Number of UEs
K=10; 
delement=0.5;
%Pilot length
tau_p=10;
%Select length of coherence block
tau_c=200;

%The cell size 1km x 1km
cellRange=1000;

%Compute the noise power
noiseFigure = 7;
B=20e6; %Communication bandwidth
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
sigma2=db2pow(noiseVariancedBm-30); %(in Watt)

%Uplink transmit power per UE (W)
p=0.2; %200 mW
%Create the power vector for all UEs (The uplink power is the same
%(p)at each UE)
pv=p*ones(1,K);



%Select the number of setups with random AP/UE  locations
nbrOfSetups=200;
%Select the number of channel realizations per setup
nbrOfRealizations=5e2;

%Single-layer decoding
%Prepare to save simulation results for Theoretical and Monte-Carlo
%MR for Monte-Carlo and CC for theoretical
sumSE_MR_MMSE= zeros(length(M),nbrOfSetups,3);
sumSE_CC_MMSE= zeros(length(M),nbrOfSetups,3);
%sumSE_MR_LMMSE= zeros(length(M),nbrOfSetups);
%sumSE_CC_LMMSE= zeros(length(M),nbrOfSetups);
%sumSE_MR_LS= zeros(length(M),nbrOfSetups);
%sumSE_CC_LS= zeros(length(M),nbrOfSetups);


%Two-layer decoding
%Prepare to save simulation results for Theoretical and Monte-Carlo
% sumSE_MR_MMSE_LSFD= zeros(length(M),nbrOfSetups);
% sumSE_CC_MMSE_LSFD= zeros(length(M),nbrOfSetups);
% sumSE_MR_LMMSE_LSFD= zeros(length(M),nbrOfSetups);
% sumSE_CC_LMMSE_LSFD= zeros(length(M),nbrOfSetups);
% sumSE_MR_LS_LSFD= zeros(length(M),nbrOfSetups);
% sumSE_CC_LS_LSFD= zeros(length(M),nbrOfSetups);

%%
for ii=1:3
    N=NN(ii);
%Go through the number of APs
for m=1:length(M)
    
    %Deploy APs randomly
    APpositions=cellRange*(rand(M(m),1) + 1i*rand(M(m),1));
    
    %For single layer decoding set all large-scale fading coefficients to 1
    A_singleLayer=reshape(repmat(eye(M(m)),1,K),M(m),M(m),K);
    
    %Go through all setups
    for n=1:nbrOfSetups

       %Deploy UEs and generate the covariance and mean matrices
       [R_AP,RHat_AP,HMean_Withoutphase,~] = RandomAP_generateSetup_Rician_Multi_Antenna(M(m),K,N,delement,1,1);
       
       
       %Create channel generations for each UE-AP pair
     [H,HMean] = functionChannelGeneration( R_AP,HMean_Withoutphase,M(m),K,N,nbrOfRealizations);
%        [H,HMean,R_xg] = functionChannelGeneration( R,HMean_Withoutphase,M(m),nbrOfRealizations,K,L );
       
       %Pilot allocation
        [Pset] = functionPilotAllocation( R_AP,HMean_Withoutphase,A_singleLayer,M(m),K,N,tau_p,pv);
%        [Pset] = functionPilotAllocation( R,HMean_Withoutphase,A_singleLayer,K,M(m),L,pv,tau_p);
%         [Rp,Dk,Lk,Phi,Phip,Phi_EW,Omega,Omegap,Sigma,C_MMSE,C_EW_MMSE,C_LMMSE,C_LS] = functionMatrixGeneration(R_AP,HMean_Withoutphase,M(m),K,N,tau_p,pv,Pset);
       %Second Layer Decoding
       %Generate Large-scale fading coefficients for phase-aware MMSE,
       %Linear MMSE and LS estimators
%        [ A_MMSE,A_LMMSE,A_LS] = functionLSFD( R,HMeanWithoutPhase,M(m),K,pv,tau_p,Pset);
       
       %Channel estimation 
       %Channel estimation with phase-aware MMSE estimator
        [Hhat_MMSE] = functionChannelEstimates_MMSE(R_AP,RHat_AP,HMean,H,nbrOfRealizations,M(m),K,N,tau_p,pv,Pset);
%        [Hhat_MMSE] =functionCellFreeMMSE(R_xg,HMean,H,nbrOfRealizations,M(m),K,L,pv,tau_p,Pset);
       %Channel estimation with LMMSE estimator
%        [Hhat_LMMSE] =functionCellFreeLMMSE(R,HMeanWithoutPhase,H,nbrOfRealizations,M(m),K,pv,tau_p,Pset);
       %Channel estimation with LS estimator
%        [Hhat_LS] = functionCellFreeLS( H,nbrOfRealizations,M(m),K,pv,tau_p,Pset);
%        [Rp,Dk,Lk,Phi,Phip,Phi_EW,Omega,Omegap,Sigma,C_MMSE,C_EW_MMSE,C_LMMSE,C_LS] = functionMatrixGeneration(R_AP,HMean_Withoutphase,M,K,N,tau_p,p,Pset);

       
       
       
       %Second layer decoding Spectral Efficiency (SE) computation
       %SE with MMSE estimator
%        [SE_MR_MMSE_LSFD ]= functionMonteCarloSE_UL(Hhat_MMSE,H,A_MMSE,tau_c,tau_p,nbrOfRealizations,M(m),K,pv);
%        [SE_CC_MMSE_LSFD] = functionTheoreticalCellFreeULSE_MMSE( R,HMeanWithoutPhase,A_MMSE,M(m),K,pv,tau_p,tau_c,Pset);
       
       %SE with LMMSE estimator
%        [SE_MR_LMMSE_LSFD ]= functionMonteCarloSE_UL(Hhat_LMMSE,H,A_LMMSE,tau_c,tau_p,nbrOfRealizations,M(m),K,pv);
%        [SE_CC_LMMSE_LSFD] = functionTheoreticalCellFreeULSE_LMMSE( R,HMeanWithoutPhase,A_LMMSE,M(m),K,pv,tau_p,tau_c,Pset);
      
       %SE with LS estimator
%        [SE_MR_LS_LSFD ]= functionMonteCarloSE_UL(Hhat_LS,H,A_LS,tau_c,tau_p,nbrOfRealizations,M(m),K,pv);
%        [SE_CC_LS_LSFD] = functionTheoreticalCellFreeULSE_LS( R,HMeanWithoutPhase,A_LS,M(m),K,pv,tau_p,tau_c,Pset);
       
      
    
       %One-layer decoding
       %SE with MMSE estimator LSFD coefficients are set to 1
       [SE_MR_MMSE] = functionComputeMonteCarloSE_UL(Hhat_MMSE,H,tau_c,tau_p,nbrOfRealizations,N,K,M(m),pv);
%        [SE_MR_MMSE ]= functionMonteCarloSE_UL(Hhat_MMSE,H,A_singleLayer,tau_c,tau_p,nbrOfRealizations,M(m),K,L,pv);
%   [SE_CC_MMSE] = functionTheoreticalCellFreeULSE_MMSE(HMean_Withoutphase,R_AP,Lk,Omega,Phi,M(m),K,N,tau_p,tau_c,pv,Pset); 
%        [SE_CC_MMSE] = functionTheoreticalCellFreeULSE_MMSE( R_xg,HMean_Withoutphase,A_singleLayer,M(m),K,L,pv,tau_p,tau_c,Pset);
         [SE_CC_MMSE] = functionTheoreticalCellFreeULSE_MMSE_OneLayer(HMean_Withoutphase,R_AP,M(m),K,N,tau_p,tau_c,pv,Pset); 
       %SE with LMMSE estimator LSFD coefficients are set to 1
%        [SE_MR_LMMSE]= functionMonteCarloSE_UL(Hhat_LMMSE,H,A_singleLayer,tau_c,tau_p,nbrOfRealizations,M(m),K,pv);
%        [SE_CC_LMMSE] = functionTheoreticalCellFreeULSE_LMMSE( R,HMeanWithoutPhase,A_singleLayer,M(m),K,pv,tau_p,tau_c,Pset);
       
       %SE with LS estimator LSFD coefficients are set to 1
%        [SE_MR_LS]= functionMonteCarloSE_UL(Hhat_LS,H,A_singleLayer,tau_c,tau_p,nbrOfRealizations,M(m),K,pv);
%        [SE_CC_LS] = functionTheoreticalCellFreeULSE_LS( R,HMeanWithoutPhase,A_singleLayer,M(m),K,pv,tau_p,tau_c,Pset);

       
       %Average of all users 
       %One layer
       sumSE_MR_MMSE(m,n,ii) = mean(SE_MR_MMSE,1);
       sumSE_CC_MMSE(m,n,ii) = mean(SE_CC_MMSE,1);
%        sumSE_MR_LMMSE(m,n) = mean(SE_MR_LMMSE,1);
%        sumSE_CC_LMMSE(m,n) = mean(SE_CC_LMMSE,1);
%        sumSE_MR_LS(m,n) = mean(SE_MR_LS,1);
%        sumSE_CC_LS(m,n) = mean(SE_CC_LS,1);
       %Second layer
%        sumSE_MR_MMSE_LSFD(m,n) = mean(SE_MR_MMSE_LSFD,1);
%        sumSE_CC_MMSE_LSFD(m,n) = mean(SE_CC_MMSE_LSFD,1);
%        sumSE_MR_LMMSE_LSFD(m,n) = mean(SE_MR_LMMSE_LSFD,1);
%        sumSE_CC_LMMSE_LSFD(m,n) = mean(SE_CC_LMMSE_LSFD,1);
%        sumSE_MR_LS_LSFD(m,n) = mean(SE_MR_LS_LSFD,1);
%        sumSE_CC_LS_LSFD(m,n) = mean(SE_CC_LS_LSFD,1);
%        
       
       
    %Output simulation progress
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    end
   clear R HMean Hhat_MMSE Hhat_LMMSE Hhat_LS
   %Output simulation progress
   disp([num2str(M(m)) ' APs out of ' num2str(M(end))]);
end
end
% Plot simulation results
%%
figure;
m1=plot(M,mean(sumSE_CC_MMSE(:,:,1),2),'ko');
hold on
m2=plot(M,mean(sumSE_CC_MMSE(:,:,2),2),'blueo');
hold on
m3=plot(M,mean(sumSE_CC_MMSE(:,:,3),2),'redo');
hold on
c1=plot(M,mean(sumSE_CC_MMSE(:,:,1),2),'k:','linewidth',2);
hold on
c2=plot(M,mean(sumSE_CC_MMSE(:,:,2),2),'blue-.','linewidth',2);
hold on
c3=plot(M,mean(sumSE_CC_MMSE(:,:,3),2),'red-','linewidth',2);
hold on
% c2=plot(M,mean(sumSE_CC_LMMSE,2),'k-x');
% hold on
% plot(M,mean(sumSE_MR_LS,2),'ks');
% hold on
% c3=plot(M,mean(sumSE_CC_LS,2),'k-o');
% hold on
% plot(M,mean(sumSE_MR_MMSE_LSFD,2),'rs');
% hold on
% plot(M,mean(sumSE_CC_MMSE_LSFD,2),'r-+');
% hold on
% plot(M,mean(sumSE_MR_LMMSE_LSFD,2),'rs');
% hold on
% plot(M,mean(sumSE_CC_LMMSE_LSFD,2),'r-x');
% hold on
% plot(M,mean(sumSE_MR_LS_LSFD,2),'rs');
% hold on
% plot(M,mean(sumSE_CC_LS_LSFD,2),'r-o');
xlabel('Number of APs, M');
ylabel('Average UL SE [bit/s/Hz]');
legend([c1 c2 c3 m1 m2 ],'L=4','L=9','L=16','Monte-Carlo','Location','Northwest')
 grid on;%set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
axis([20 100 1.5 5.5])


% save('SE.mat','sumSE_CC_MMSE');