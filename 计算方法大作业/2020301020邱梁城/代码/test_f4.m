close all; clc; clear all;
linewd = 0.8;
hcfontsize = 15;

MarkerSize=5;
A0=0.5;
f0=0.1;
L=30;
M_values=[1 2 5 10]; % Number of sensors
[t1,t2]=meshgrid(1:L,1:L);
Rx1=0.5*A0^2*cos(2*pi*f0*(t1-t2));

iSNR_db_values=-5:15;  
Gain_dB_values=zeros(length(iSNR_db_values),length(M_values));
Gain_dB_values_my=zeros(length(iSNR_db_values),length(M_values));
MSE_values=zeros(length(iSNR_db_values),length(M_values));
MSE_values_my=zeros(length(iSNR_db_values),length(M_values));
xi_n_dB_values=zeros(length(iSNR_db_values),length(M_values));
xi_n_dB_values_my=zeros(length(iSNR_db_values),length(M_values));
xi_d_dB_values=zeros(length(iSNR_db_values),length(M_values));
xi_d_dB_values_my=zeros(length(iSNR_db_values),length(M_values));
for idxM=1:length(M_values)
    M=M_values(idxM);
    
    Rvw=0.01*eye(L*M);
    Rvi=eye(L*M);
    for m1=0:M-1
        for m2=0:M-1
            if L>abs(m2-m1)
                Rvi((m1*L+1):(m1*L+L),(m2*L+1):(m2*L+L))=diag(ones(L-abs(m2-m1),1),m1-m2);
            end
        end
    end
    Rv0=Rvi+Rvw;
    Rx=kron(ones(M),Rx1);
    for idxS=1:length(iSNR_db_values)
        W=zeros(M*L,L);
        iSNR_dB=iSNR_db_values(idxS);
        iSNR=10^(iSNR_dB/10);%求的isnr
        
        sigma_i2=A0^2/2/(1.01*iSNR);
        Rv=sigma_i2*Rv0;%求RV
        Ry=Rx+Rv;
        Q=[eye(L), zeros(L, (M-1)*L)];
        B=(Ry-Rv)*Q';
        tic
        for i=1:L
            b=B(:,i);
            x0 = b;
            x = x0;
            r = b - Ry*x0;
            rho0 = conj(r')*r;
            rho00 = rho0;
            p=r;
            while(rho0>1e-9)
                w= Ry*p;
                alpha = rho0/(conj(p')*w);
                x = x+alpha*p;
                r = r- alpha*w;
                rho1 = conj(r')*r;
                beta = rho1/rho0;
                p = r + beta*p;
                rho0 = rho1;   
            end
            W(:,i)=x;
        end
        elapsed_time = toc;
        %生成系数矩阵
        W=W';
        Rxfd=W*Rx*W';
        Rvrn=W*Rv*W';
        tic
        [T,Lambda]=jeig(Rx,Rv);
        [diagL,idx1]=sort(diag(Lambda),'descend');
        Lambda=diag(diagL);
        T=T(:,idx1);
        
        I_i=[eye(L) zeros(L,(M-1)*L)];
        A=I_i*Rx*T/(Lambda+eye(L*M));%求得滤波器 A T的值
        esssss=toc;
        time(idxS)=esssss-elapsed_time;
        oSNR_MY=trace(Rxfd)/trace(Rvrn);
        oSNR_dB_MY=10*log10(oSNR_MY);
        Gain_dB_values_my(idxS,idxM)=oSNR_dB_MY-iSNR_dB;
        oSNR=trace(A*Lambda*A')/trace(A*A');
        oSNR_dB=10*log10(oSNR);
        Gain_dB_values(idxS,idxM)=oSNR_dB-iSNR_dB;
        MSE_values_my(idxS,idxM)=10*log10((trace(Rx1)+trace(W*Ry*W')-2*trace(W*((Ry-Rv)*Q'))));
        MSE_values(idxS,idxM)=10*log10( trace(Rx1-2*A*T'*Rx*I_i'+A*(Lambda+eye(L*M))*A'));
        xi_n_my=trace(Rv(1:L,1:L))/trace(Rvrn);
        xi_n_dB_values_my(idxS,idxM)=10*log10(xi_n_my);
        xi_n=trace(Rv(1:L,1:L))/trace(A*A');
        xi_n_dB_values(idxS,idxM)=10*log10(xi_n);
        xi_d_my=trace(Rx1)/trace(Rxfd);
        xi_d_dB_values_my(idxS,idxM)=10*log10(xi_d_my);
        xi_d=trace(Rx1)/trace(A*Lambda*A');
        xi_d_dB_values(idxS,idxM)=10*log10(xi_d);
    end
end
% figure
% plot(iSNR_db_values(1:2:end),Gain_dB_values_my(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
% hold on
% plot(iSNR_db_values(1:2:end),Gain_dB_values_my(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
% plot(iSNR_db_values(1:2:end),Gain_dB_values_my(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
% plot(iSNR_db_values(1:2:end),Gain_dB_values_my(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
% hold off
% set(gca, 'Color', [1, 1, 1]); 
% set(gca, 'FontName', 'Times New Roman');
% set(gca, 'FontSize', hcfontsize);
% set(gca, 'LineWidth', linewd); 
% box on; grid on;
% 
% figure
% plot(iSNR_db_values(1:2:end),MSE_values_my(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
% hold on
% plot(iSNR_db_values(1:2:end),MSE_values_my(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
% plot(iSNR_db_values(1:2:end),MSE_values_my(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
% plot(iSNR_db_values(1:2:end),MSE_values_my(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
% hold off
% set(gca, 'Color', [1, 1, 1]); 
% set(gca, 'FontName', 'Times New Roman');
% set(gca, 'FontSize', hcfontsize);
% set(gca, 'LineWidth', linewd); 
% box on; grid on;

figure
plot(iSNR_db_values(1:2:end),xi_n_dB_values_my(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),xi_n_dB_values_my(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_n_dB_values_my(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_n_dB_values_my(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
ylabel('\xi _n(H)')
xlabel('iSNR(dB)')
legend('N=1','N=2','N=5','N=10')
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_db_values(1:2:end),xi_d_dB_values_my(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),xi_d_dB_values_my(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_d_dB_values_my(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_d_dB_values_my(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
ylabel('\xi _d(H)')
xlabel('iSNR(dB)')
legend('N=1','N=2','N=5','N=10')
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_db_values(1:2:end),xi_n_dB_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),xi_n_dB_values(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_n_dB_values(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_n_dB_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;

figure
plot(iSNR_db_values(1:2:end),xi_d_dB_values(1:2:end,1),'-bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(iSNR_db_values(1:2:end),xi_d_dB_values(1:2:end,2),'--g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_d_dB_values(1:2:end,3),':rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(iSNR_db_values(1:2:end),xi_d_dB_values(1:2:end,4),'-.c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
set(gca, 'Color', [1, 1, 1]); 
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
box on; grid on;


