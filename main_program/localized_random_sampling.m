%% Source code for
%% [1] N.J.Myers, and R. W. Heath Jr.,"Localized random sampling for robust compressive beam alignment",
%% to appear in proc. of the IEEE ICASSP 2019
%% If you use this code, please cite our ICASSP paper [1] (or)
%% [2] N.J.Myers, A. Mezghani and R. W. Heath Jr., "Swift-Link: A compressive beam alignment algorithm for practical mmWave radios" IEEE Lrxans. on Signal Process., vol. 67, no. 4, pp. 1104?1119, 2019.
%% In this paper, we show that by a clever design of compressed sensing (CS) matrix, standard CS for mmWave beam alignment can be made robust to carrier frequency offset
clear all;
close all;
clc;
addpath('../channel_data_20MHz_NB_38GHz_64_64/')
tic
rng(24);
N=64;                             % No. of antennas at the TX= No. of antennas at the RX
UN=dftmtx(N)/sqrt(N);             % Standard DFT matrix - for sparsifying dictionary at the TX / RX

Meas=126;                         % Use 2*Meas so that Meas are acquired along p-trajectory
qb=3;                             % Resolution of phase shifters at the TX and the RX
Zrx=quantz(exp(1i*pi*11*([0:1:N-1].^2)/N),qb).'/sqrt(N);    % Zadoff-Chu sequence for the RX. Root 11 is chosen such that Zrx is a perfect sequence even after quantization.
Ztx=quantz(exp(1i*pi*11*([0:1:N-1].^2)/N),qb).'/sqrt(N);    % Zadoff-Chu sequence for the RX. Root 11 is chosen such that Ztx is a perfect sequence even after quantization.
% quantz does qbit quantization

SNRvec_new=[];                                             %  New procedure- CS using the training defined in Algorithm 1 of [1]
SNRvec_old=[];                                             %  Old procedure- CS using standard random IID phase shift-based training
SNRvec_old_broad=[];                                       %  Beam broadening after CS with random IID phase shift-based training
SNRvec_max=[];
resultid=2;                                                %  Use 1 for experiment in Fig. 1, and 2 for experiment in Fig. 2.

if(resultid==1)
    cfolist=[20.5579];                                     % Maximally off-grid CFO within the 40ppm limit
% For a faster simulation, I set Nruns=5; The results in [1] use Nruns=20
    Nreal=100;                                             % Number of channel realizations
    Nruns=5;                                             % Number of CS matrix realizations
else
    cfolist=[0,5,15,25,40];    % Results in [1] are for cfolist=[-40:5:40], Nreal=100 and Nruns=20;
    Nreal=30;
    Nruns=5;
end

SNRgain_max=0;
for cfoppm=cfolist
    Ind_err_old=[];                                             % Angle error induced by CFO with the old procedure
    Ind_err_new=[];                                             % Angle error induced by CFO with the new procedure
    SNRnew=[];
    SNRold=[];
    SNRold_broad=[];
    BW=100;                                                % Bandwidth in MHz
    f_e=38*cfoppm*1e-3;                                    % CFO in MHz
    we=2*pi*f_e/BW;                                        % CFO in radians
    SNR=0;                                                 % SNR- set it to 0dB in [1]
    
    OMPiter=1;                                             % Number of OMP iterations- just 1 here as we test the best direction given channel measurements
    cfoerr=exp(1i*we*[1:1:Meas]).';                        % phase error induced by CFO in samples
    
    for realz=1:1:Nreal
        Hv=dlmread(strcat('it_',num2str(realz),'time_dom_channel.txt')); % Channel derived from NYU simulator
        Hmat=reshape(Hv,[N,N]);
        for runs= 1:1:Nruns
            noise=(randn(Meas,1)+1i*randn(Meas,1))/sqrt(2);              % noise for AWGN in channel measurements, scaled according to SNR
            %% Training Generation with the new design, i.e.,   Choose Meas number of coordinates  according to Algorithm 1
            Lrx=[];
            Ltx=[];
            indx=[N-1-(Meas/2):1:N+(Meas/2)-2];                           % In [1], we consider a specific case of Meas=126.
            for ms=1:1:Meas
                if(indx(ms)< N-1)
                    if(indx(ms)==0)
                        rs=[0];
                        cs=[0];
                    else
                        rs=binornd(indx(ms),0.5,[1,1]);  %indx(ms)+1 points  0--->indx(ms)
                        cs=indx(ms)-rs;
                    end
                else
                    if(indx(ms)==2*N-2)
                        rs=[N-1];
                        cs=[N-1];
                    else
                        rs=indx(ms)-N+1+binornd(2*(N-1)-indx(ms),0.5,[1,1]);
                        cs=indx(ms)-rs;
                    end
                end
                Lrx=[Lrx,rs];                   % Notice that Lrx(ms)+Ltx(ms)=ms according to Algorithm 1 in [1]
                Ltx=[Ltx,cs];
            end
            Wrp=zeros(N,Meas);
            Ftp=zeros(N,Meas);
            for i=1:1:Meas
                Wrp(:,i)=circshift(Zrx,[Lrx(1,i),0]);
                Ftp(:,i)=circshift(Ztx,[Ltx(1,i),0]);
            end
            % Training design according to the proposed technique complete-Wrp,Ftp
            % The m^th column of Wrp represents the beam training vector used at the RX in the m^th measurement interval
            % The m^th column of Ftp represents the beam training vector used at the TX in the m^th measurement interval
            %% Random IID phase shift-based training
            
            Wrr= (pi/2^qb)+ (randi(2^qb,[N,Meas])*2*pi/(2^qb));         %random IID phase shifts
            Ftr= (pi/2^qb)+ (randi(2^qb,[N,Meas])*2*pi/(2^qb));
            Wrr=exp(1i*Wrr)/sqrt(N); Ftr=exp(1i*Ftr)/sqrt(N);
            
            %% Perfect CSI setting, DFT-based codebook
            if(runs==1 & cfoppm==cfolist(1))                                                % Calculate once per channel
            Xbest=ifft2(Hmat);
            [val,indbest]=max(abs(Xbest(:)));                          %Pick max location corresponding to DFT-based beamspace
            [rx_best,tx_best]=ind2sub(size(Xbest),indbest);
            arx_best=quantz(conj(UN(:,rx_best)),qb);                   % find beamforming vectors corresponding to the best beamspace location
            atx_best=quantz(conj(UN(:,tx_best)),qb);
            SNRgain_best=10*log10(abs(arx_best.'*Hmat*atx_best)^2);    % SNR obtained after DFT codebook-based beam alignment under perfect CSI
            SNRgain_max=SNRgain_max+SNRgain_best;
            end
            %% CS matrix construction with New and Old procedures and channel measurement definition
            Yp=zeros(Meas,1); 
            Ysimp=zeros(Meas,1);
            for ms=1:1:Meas
                PG=kron(Ftp(:,ms).',Wrp(:,ms).')*Hv;                    % Pure sample acquired with new design- not corrupted by CFO / AWGN
                Yp(ms)= sum(PG)*cfoerr(ms);                             % CFO error introduced
                Ap(ms,:)=kron(Ftp(:,ms).'*UN,Wrp(:,ms).'*UN);           % CS matrix corresponding to the new design
            end
            for ms=1:1:Meas
                Gain=kron(Ftr(:,ms).',Wrr(:,ms).')*Hv;                  % The same stuff as above- but with random IID phase shift-based training
                Ysimp(ms)= sum(Gain)*cfoerr(ms);
                Aold(ms,:)=kron(Ftr(:,ms).'*UN,Wrr(:,ms).'*UN);
            end
            
            sig=10^(-SNR/20);
            Ynew=Yp+(sig*noise([1:1:Meas]));                          % Can be verified here that norm(Yp(k)) is 1 in expectation, so SNR defn. is consistent
            % Now, Ynew is the channel measurement with proposed beam training, perturbed by CFO and AWGN
            Yold=Ysimp+(sig*noise([1:1:Meas]));
            % Yold is the channel measurement with random phase shift-based beam training, perturbed by CFO and AWGN
         
            %% Run a single step of matching pursuit algorithm to find the best beamspace location that best explains the measurements
            % Proposed training
            Xe1=OMPf(Ynew,Ap,sig,OMPiter);
            Xin=reshape(Xe1,[N,N]);
            [val,indnew]=max(abs(Xin(:)));
            [rx_new,tx_new]=ind2sub(size(Xin),indnew); 
            Ind_err_new=[Ind_err_new,abs(rx_new-rx_best)];              % Find beam misalignment induced by CFO at the receiver, for CS with the proposed training              
            % IID phase shift-based training
            Xold=OMPf(Yold(1:Meas),Aold(1:Meas,:),sig,OMPiter);
            Xold=reshape(Xold,[N,N]);
            [val,indsimp]=max(abs(Xold(:)));
            [rx_old,tx_old]=ind2sub(size(Xold),indsimp);
            Ind_err_old=[Ind_err_old,abs(rx_old-rx_best)];              % Find beam misalignment induced by CFO at the receiver, for CS with random IID phase shift-based training 
            
            %% Beam alignment with the mismatched beamformer-old procedure
            arx_old=quantz(conj(UN(:,rx_old)),qb);
            atx_old=quantz(conj(UN(:,tx_old)),qb);
            SNRgain_simp=10*log10(abs(arx_old.'*Hmat*atx_old)^2);
            
            %% Beam alignment with a broadened beamformer-old procedure
            energ_rx=spreadenerg(N,rx_old);                                 % spreads out the beam by 1 unit
            arx_old=quantz(designawv(N,energ_rx,qb, 100),3).';              % broadened beamforming vector RX
            energ_tx=spreadenerg(N,tx_old);                                 % spreads it out by 1 unit
            atx_old=quantz(designawv(N,energ_tx,qb, 100),3).';              % broadened beamforming vector RX
            SNRgain_simp_broad=10*log10(abs(arx_old.'*Hmat*atx_old)^2);

             %% Beam alignment with a broadened beamformer-new procedure
            energ_rx=spreadenerg(N,rx_new);                                 %  spreads out the beam by 1 unit
            arx_new=quantz(designawv(N,energ_rx,qb, 100),3).';              % broadened beamforming vector RX
            energ_tx=spreadenerg(N,tx_new);                                 % spreads out the beam by 1 unit
            atx_new=quantz(designawv(N,energ_tx,qb, 100),3).';              % broadened beamforming vector TX
            SNRgain_new=10*log10(abs(arx_new.'*Hmat*atx_new)^2);
            SNRnew=[SNRnew,SNRgain_new];
            SNRold=[SNRold,SNRgain_simp];
            SNRold_broad=[SNRold_broad,SNRgain_simp_broad];
        end
    end
    SNRvec_new=[SNRvec_new,SNRnew.'];
    SNRvec_old=[SNRvec_old,SNRold.'];
    SNRvec_old_broad=[SNRvec_old_broad,SNRold_broad.'];
end
SNRgain_max=SNRgain_max/Nreal;

%% Plotting results
if(resultid==1)  % plot CDF of beam misalignment induced by CFO
    Ind_err_old=Ind_err_old*2*pi/N;             % To make Angle errors in [-pi,pi]
    Ind_err_old(Ind_err_old>pi)=2*pi-Ind_err_old(Ind_err_old>pi);
    Ind_err_new=Ind_err_new*2*pi/N;
    Ind_err_new(Ind_err_new>pi)=2*pi-Ind_err_new(Ind_err_new>pi);
    Ind_err_old=abs(Ind_err_old);            % Absolute of angle error in [0,pi]
    Ind_err_new=abs(Ind_err_new);
    cdfy=[1:1:length(Ind_err_new)]/length(Ind_err_new);
    figure(1)
    plot(sort(Ind_err_new),cdfy,'b--','LineWidth',2)
    hold on;
    plot(sort(Ind_err_old),cdfy,'k-','LineWidth',2)
    h_legend=legend('Proposed training', 'IID random training');
    set(h_legend,'Position',[0.602687971932548 0.115238096600487 0.297312028067453 0.0907142843518938],'Interpreter','Latex','FontSize',14);
    xlab=xlabel('Angle (rad) ','Interpreter','Latex')
    xlim([0,pi])
    set(xlab,'FontSize',14);
    ylab=ylabel('Beam misalignment CDF','Interpreter','Latex');
    set(ylab,'FontSize',14);
    set(gca,'fontsize',14);
end

if(resultid==2)
   PSNRvec_old=mean(max(SNRvec_old,0));
   PSNRvec_new=mean(max(SNRvec_new,0));
   PSNRvec_obd=mean(max(SNRvec_old_broad,0));
   AL=6.8*ones(size(PSNRvec_obd));                  % Evaluated for Agile Link(non-coherent algorithm) at SNR=0dB. 
   MaxSNR=SNRgain_max*ones(size(PSNRvec_obd));       
   plot(cfolist,MaxSNR,'r-','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor','r')
   hold on;
   plot(cfolist,PSNRvec_new,'bd--','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','r')
   plot(cfolist,PSNRvec_old,'k^-','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','g')
   plot(cfolist,PSNRvec_obd,'ko-','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','y')
   plot(cfolist,AL,'mv-','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','c')
   set(gca,'Fontsize',14)
   xlab=xlabel('CFO in ppm ','Interpreter','Latex')
   set(xlab,'FontSize',14);
   ylab=ylabel('SNR (dB)','Interpreter','Latex');
   set(ylab,'FontSize',14);
   h_legend=legend('Perfect CSI','Proposed training, broaden', 'Random training','Random training, broaden','Agile-Link')
   set(h_legend,'Position',[0.528372682843889 0.283333335036324 0.371826062883649 0.214285710879735],'Interpreter','Latex','FontSize',14);
end
toc