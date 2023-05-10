function [ R_AU,R_AR,R_RU,R,HMean] = functionRISCellFreeSetup( M,K,NN,cellRange,APpositions,sigma2 )


%The minimum allowed distance (Access Point to User Equipment)
dmin=10; 
%The maximum distance of LoS occurence
%maxDistanceLoS=300; 
%The height of Access Point
APheigth=12.5;  
%The height of User Equipment
UEheigth=1.5;  
RISheigth=15;
%The height difference between AP and UE 
verticalDistance_AU=APheigth-UEheigth;
% verticalDistance_AR=APheigth-RISheigth;
% verticalDistance_RU=RISheigth-UEheigth;

%Dropping all UEs while minimum distance requirement is satisfied.
droppedUEs=0;
%Preparing to save the distances and UE positions
distance_AU=zeros(M,K);
UEpositions=zeros(K,1);

%Compute alternative AP locations by using wrap around
wrapHorizontal = repmat([-cellRange 0 cellRange],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[M 1]);
%Dropping all UEs
while droppedUEs <K
    
    UEposition=cellRange*(rand(1,1) + 1i*rand(1,1));
    horizontalDistance_AU=sqrt(real(APpositions-UEposition).^2+imag(APpositions-UEposition).^2);
    distanceAU=sqrt(verticalDistance_AU.^2 + horizontalDistance_AU.^2);
    
    if isempty(distanceAU(distanceAU<dmin))
        droppedUEs=droppedUEs+1;
        distance_AU(:,droppedUEs)=distanceAU;
        
        %Store UE positions
        UEpositions(droppedUEs)=UEposition;
    end
    
end


%Calculating the distances between all AP and UE pairs
for k = 1:K
    distance_AU(:,k) = min(abs(APpositionsWrapped - repmat(UEpositions(k),size(APpositionsWrapped))),[],2);
end

RISpositions=cellRange*(0.5 + 1i*0.5);
distance_AR=sqrt(real(APpositions-RISpositions).^2+imag(APpositions-RISpositions).^2+(APheigth-RISheigth).^2);
distance_RU=sqrt(real(RISpositions-UEpositions).^2+imag(RISpositions-UEpositions).^2+(RISheigth-UEheigth).^2);

%三段式求大尺度beta
d0 = 10; %meter
d1 = 50; %meter
%Constant term in the model from (53) in [15]
L = 140.7151;
%Compute the pathloss using the three-slope model in (52) in [15]
PL_AU = zeros(M,K);
PL_AR = zeros(M,1);
PL_RU = zeros(K,1);
for m=1:M
    for k=1:K
    if distance_AU(m,k)<=d0
        PL_AU(m,k) = -L -15*log10(d1/1000) -20*log10(d0/1000);
    elseif distance_AU(m,k)<=d1
        PL_AU(m,k) = -L -15*log10(d1/1000) -20*log10(distance_AU(m,k)/1000);
    else
        PL_AU(m,k) = -L -35*log10(distance_AU(m,k)/1000);
    end
    end
end
for m=1:M
    if distance_AR(m,1)<=d0
        PL_AR(m,1) = -L -15*log10(d1/1000) -20*log10(d0/1000);
    elseif distance_AR(m,1)<=d1
        PL_AR(m,1) = -L -15*log10(d1/1000) -20*log10(distance_AR(m,1)/1000);
    else
        PL_AR(m,1) = -L -35*log10(distance_AR(m,1)/1000);
    end
end
for k=1:K
    if distance_RU(k,1)<=d0
        PL_RU(k,1) = -L -15*log10(d1/1000) -20*log10(d0/1000);
    elseif distance_RU(k,1)<=d1
        PL_RU(k,1) = -L -15*log10(d1/1000) -20*log10(distance_RU(k,1)/1000);
    else
        PL_RU(k,1) = -L -35*log10(distance_RU(k,1)/1000);
    end
end

sigma_sf=8;
shadowFading_AU=10.^(sigma_sf*randn(M,K)/10); 
shadowFading_AR=10.^(sigma_sf*randn(M,1)/10); 
shadowFading_RU=10.^(sigma_sf*randn(K,1)/10); 
channelGain_AU=(db2pow(PL_AU).*shadowFading_AU)/sigma2;
channelGain_AR=(db2pow(PL_AR).*shadowFading_AR)/sigma2;
channelGain_RU=(db2pow(PL_RU).*shadowFading_RU)/sigma2;

R_AU=zeros(M,M,K);
for k=1:K
R_AU(:,:,k)=diag(channelGain_AU(:,k));
end

R_AR=diag(channelGain_AR(:,1));

R_RU=diag(channelGain_RU(:,1));

R=zeros(M,M,K);
r=zeros(M,K);
for m=1:M
    for k=1:K
        r(m,k)=R_AU(m,m,k)+NN*R_AR(m,m)*R_RU(k,k);
    end
end
for k=1:K
    R(:,:,k)=diag(r(:,k));

end
    HMean=0*ones(M,K);




 end