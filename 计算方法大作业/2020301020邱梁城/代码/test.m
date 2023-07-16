close all; clc; clear all;

linewd = 0.6;
hcfontsize = 10;
MarkerSize=2;

A0=0.5;
f0=0.1;
sigma_i2=0.5; %高斯白噪声
sigma_w2=sigma_i2/100;  %传感器里的噪声

L_values=1:60;  % length of the filter     
M_values=[1 2 5 10]; % Number of sensors   相当于N的个数
Gain_dB_values=zeros(length(L_values),length(M_values));
Gain_dB_values_my=zeros(length(L_values),length(M_values));
MSE_values=zeros(length(L_values),length(M_values));
MSE_values_my=zeros(length(L_values),length(M_values));
time=zeros(4,2)
for idxM=1:length(M_values)
    M=M_values(idxM);
    for idx=1:length(L_values)
        L=L_values(idx); %变换滤波器的阶数  从1到60
        W=zeros(M*L,L);
        Rvw=sigma_w2*eye(L*M);
        Rvi=eye(L*M);
        for m1=0:M-1
            for m2=0:M-1
                if L>abs(m2-m1)
                    Rvi((m1*L+1):(m1*L+L),(m2*L+1):(m2*L+L))=diag(ones(L-abs(m2-m1),1),m1-m2);
                end
            end
        end
        Rvi=sigma_i2*Rvi;
        Rv=Rvi+Rvw; %噪声的自相关矩阵 直接拿来用？
        [t1,t2]=meshgrid(1:L,1:L);
        %Rx(i,j) = 0.5*A0^2*cos(2*pi*f0*(t1(i,j)-t2(i,j)))
        Rx1=0.5*A0^2*cos(2*pi*f0*(t1-t2));
        Rx=kron(ones(M),Rx1);
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
            while(rho0>1e-6)
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
        time(idxM,1)=elapsed_time;
        %生成RX的自相关矩阵
        W=W';
%     figure(1)
%     plot(L_eig);
%     hold on
%    plot(L_eig_Rx([1:60]));
        tic
        [T,Lambda]=jeig(Rx,Rv);
        [diagL,idx1]=sort(diag(Lambda),'descend');
        Lambda=diag(diagL);
        T=T(:,idx1);
        Rxfd=W*Rx*W';
        Rvrn=W*Rv*W';
        iSNR_my=trace(Rx1)/trace(Rv(1:L,1:L));
        iSNR_dB_my=10*log10(iSNR_my);
        oSNR_MY=trace(Rxfd)/trace(Rvrn);
        oSNR_dB_MY=10*log10(oSNR_MY);
        Gain_dB_values_my(idx,idxM)=oSNR_dB_MY-iSNR_dB_my;
        I_i=[eye(L) zeros(L,(M-1)*L)];
        A=I_i*Rx*T/(Lambda+eye(L*M));
        real_filter=A*T';
        elapsed_time1=toc;
        time(idxM,2)=elapsed_time1;
        iSNR=trace(Rx1)/trace(Rv(1:L,1:L));
        iSNR_dB=10*log10(iSNR);
        oSNR=trace(A*Lambda*A')/trace(A*A');
        oSNR_dB=10*log10(oSNR);
        Gain_dB_values(idx,idxM)=oSNR_dB-iSNR_dB;
        MSE_values(idx,idxM)=10*log10(1/L*trace(Rx1-2*A*T'*Rx*I_i'+A*(Lambda+eye(L*M))*A'));
        MSE_values_my(idx,idxM)=10*log10(1/L*(trace(Rx1)+trace(W*Ry*W')-2*trace(W*((Ry-Rv)*Q'))));
    end
end
L_eig=eig(Rx1);
L_eig_Rx=eig(Rx);
figure
plot(L_eig,'r') 
figure
plot(L_eig_Rx([1:60]),'y');
figure(1)
subplot(121)
plot(L_values,Gain_dB_values(:,1),'bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(L_values,Gain_dB_values(:,2),'g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(L_values,Gain_dB_values(:,3),'rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(L_values,Gain_dB_values(:,4),'c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
title('Tradition GAIN')
set(gca, 'Color', [1, 1, 1]); 
ylabel('G(H_W)(dB)')
xlabel('filter:L')
legend('N=1','N=2','N=5','N=10')
set(gca, 'FontName','Microsoft YaHei');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 max(L_values)]);
set(gca,'XTick',10:10:max(L_values)); 
box on; grid on;

figure(2)
subplot(121)
plot(L_values,MSE_values(:,1),'bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(L_values,MSE_values(:,2),'g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(L_values,MSE_values(:,3),'rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(L_values,MSE_values(:,4),'c^','linewidth',linewd, 'MarkerSize',MarkerSize);
hold off
title('Tradition MSE')
set(gca, 'Color', [1, 1, 1]); 
ylabel('J(H_W)/L(dB)')
xlabel('filter:L')
legend('N=1','N=2','N=5','N=10')
% set(gca, 'FontName', 'Microsoft YaHei');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 max(L_values)]);
set(gca,'XTick',10:10:max(L_values)); 
box on; grid on;

%CG
%figure
L_eig=eig(Rx1);
% L_eig_Rx=eig(Rx);
% figure
% plot(L_eig,'r') 
% figure
% plot(L_eig_Rx([1:60]),'y');
figure(1)
subplot(122)
plot(L_values,Gain_dB_values_my(:,1),'bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(L_values,Gain_dB_values_my(:,2),'g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(L_values,Gain_dB_values_my(:,3),'rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(L_values,Gain_dB_values_my(:,4),'c^','linewidth',linewd, 'MarkerSize',MarkerSize);
title('CG GAIN')
ylabel('G(H_W)(dB)')
xlabel('filter:L')
hold off
legend('N=1','N=2','N=5','N=10')
set(gca, 'Color', [1, 1, 1]); 
% set(gca, 'FontName','Microsoft YaHei');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 max(L_values)]);
set(gca,'XTick',10:10:max(L_values)); 
box on; grid on;


%CG
figure(2)
subplot(122)
plot(L_values,MSE_values_my(:,1),'bo','linewidth',linewd, 'MarkerSize',MarkerSize);
hold on
plot(L_values,MSE_values_my(:,2),'g*','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(L_values,MSE_values_my(:,3),'rs','linewidth',linewd, 'MarkerSize',MarkerSize);
plot(L_values,MSE_values_my(:,4),'c^','linewidth',linewd, 'MarkerSize',MarkerSize);
title('CG MSE')
legend('N=1','N=2','N=5','N=10')
ylabel('J(H_W)/L(dB)')
xlabel('filter:L')
hold off
set(gca, 'Color', [1, 1, 1]); 
% set(gca, 'FontName','Microsoft YaHei');
set(gca, 'FontSize', hcfontsize);
set(gca, 'LineWidth', linewd); 
set(gca,'XLim',[1 max(L_values)]);
set(gca,'XTick',10:10:max(L_values)); 
box on; grid on;