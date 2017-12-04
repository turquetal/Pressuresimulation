clear all
close all
clc 
addpath ~/Dropbox/Super_Mega_Script_Folder/
addpath C:\Users\usr1\Dropbox\Super_Mega_Script_Folder
addpath D:\Dropbox\Dropbox\Super_Mega_Script_Folder
addpath C:\Users\User1\Dropbox\Super_Mega_Script_Folder
% svloc= 'C:\Users\User1\Dropbox\Work\Oct 16\JGR Code\Calisan\'
FNAMEFMT = '2exp_C001H001S0001000%03d.bmp';
imglimitx=23:1002;
imglimity=186:560;
num_images=1;
% im=imread(sprintf(FNAMEFMT, 1));
% IMG=im2bw(im(186:560,23:1002));
rati=0.1;      %resize factos
pixels=870*rati;

row=floor(linspace(23,1002,6));
thr=[0.245 0.36 0.37 0.47 0.55];
pix=0.0077;


% for  i=1:350
% im=imread(sprintf(FNAMEFMT, (i)));
% imr(:,:,i)=im(186:560,23:1002);
% 
% 
%  for j=1:length(row)-1
% img(:,row(j)-22:row(j+1)-22)=im2bw(im(186:560,row(j):row(j+1)),thr(j));
% 
%  end
% 
% 
% B=1-imresize(img,0.1);
% 
% B1=zeros(size(B,1)+2,size(B,2)+2);
% bound=zeros(size(B,1)+2,size(B,2)+2);
% bound(2:size(B,1)+1,2:size(B,2)+1)=1;
% B1(2:size(B,1)+1,2:size(B,2)+1)=B;
% 
% 
% fin_img(:,:,i)=B1;f
% end

load imdata.mat

B1=fin_img(:,:,1);

B=zeros(size(B1,1)-2,size(B1,2)-2);
bound=zeros(size(B,1)+2,size(B,2)+2);
bound(1:size(B,1)+2,2:size(B,2)+2)=1;


    



% % % % % % % % % % % % % % % % % % % % %            Data for lamb waves
addpath C:\Users\usr1\Dropbox\Super_Mega_Script_Folder
% Xs=(-0.2:0.005:0.2);
% Ys=(-0.05:0.005:0.05);
% rhob = 7800; %kg/m^3 densite 2500 verre / 7800 acier / 1140 polyamide / 2600 granite
% Eb = 203e9; %Pa module d'Young 70e9 Verre / 203e9 acier / 4e9 polyamide / 60e9 granite
% nub = 0.3; % coef Poisson 0.2 verre / 0.3 acier / 0.4 polyamide / 0.27 granite
% 
% rhop = 2500; %kg/m^3 2500 verre / 2800 marbre / 1180 gplexi
% Ep = 70e9; %Pa 70e9 verre / 26e9 marbre / 2.4e9 gplexi
% nup = 0.2; % 0.2 verre / 0.3 marbre / 0.375 gplexi
% h = 0.01; %m % plate thickness
% B = h^3*Ep/(12*(1-nup^2)); % bending stiffness m^3 Pa
% E = 1/((1 - nub^2)/Eb + (1 - nup^2)/Ep); % (Pa) module equivalent
% G=Ep/(2*(1+nup));  %shear modulus
% Vray=(0.87+1.12*nup)/(1+nup).*sqrt(G/rhop);
% omega = 2*pi*f0;
% k = ((rhop*h/B)^(1/4).*sqrt(omega));% wave number
% Vph=omega./k;
% c=zeros(length(Vph));
% c=(diff(omega)./diff(k));
% vind=find(c>Vray,1); %cut off frequency for non-dispersive phase
% c(vind:end)=Vray;
% k(vind:end)=k(vind)/omega(vind)*omega(vind:end);
load dispexp
w1=omega;
k1=k;
load Parameters
k=interp1(w1,k1,omega,'linear','extrap');

xrec=19;
yrec=19;
% % % % % 
dt=1e-6;
winsize=8e-3;
winval=7999;

timewin=600;
tslot=0:dt:7999*dt;
% GG=0;
GG=zeros(timewin,4097);
GG1=zeros(timewin,4097);
GGp=zeros(timewin,4097);
GGs=zeros(timewin,4097);
% coco=1;
twin=0:8e-3:8e-3*599;
% SS=0;
Psave=zeros(40,100,1);
% comparison of different time windows in long time window
% flist=[1 10 50 100 130 170];
flist=1:599;
for z=1:length(flist)
P2=zeros(40,100,1);
P=zeros(40,100,1);
for s=1:2
    
 svname=['neweq_efsig_U_test',num2str(flist(z)+s-1),'.mat'];
 load(svname)
sigxx=repmat(sigxx,[1,1,size(U(:,:,2:end-1),3)]);
 P2(:,:,end+1:end+size(U(:,:,2:end-1),3))=sigxx+U(:,:,2:end-1);
 P(:,:,end+1:end+size(U(:,:,2:end-1),3))=U(:,:,2:end-1);
 Bcheck=fin_img(:,:,flist(z));
end
        t1=0:1e-5:1e-5*(length(P)-1);
        time=0:1e-6:max(t1);

 Psol=P2-P;

 clear U efsig
 z
 
for xi=2:98
    for yi=2:39
    
aa2=squeeze(P2(yi,xi,:)-Psave(yi,xi,1));
        aa2=interp1(t1,aa2,time);
  
if aa2<2e3  
aa2=squeeze(P2(yi,xi,:)-P2(yi,xi,1));
 aa2=interp1(t1,aa2,time);
end
aa=squeeze(P(yi,xi,:)-P(yi,xi,1));
asol=squeeze(Psol(yi,xi,:)-Psol(yi,xi,1));
 asol=interp1(t1,asol,time);
  aa=interp1(t1,aa,time);
if sum(abs(asol))>0
[w0,f0]=fftrl(fliplr(asol),time,0,2^nextpow2(length(tslot)));
r=sqrt((xi-xrec)^2+(yi-yrec)^2);
r=r./87.*0.66;
GGsol=-1i.*w0./(8*B.*k.^2).*sqrt(2./(r.*k*pi)).*exp(-1i*(r.*k-pi/4)).*omega.^2;
GGsol(1)=0;
GGsol(isinf(GGsol))=0;
GGs(z,:)=GGs(z,:)+GGsol;
clear aasol GGsol
end    
if sum(abs(aa2))>0
     if Bcheck(yi,xi)
% plot(aa2)
% pause
if sum(abs(aa))>0
[w0,f0]=fftrl(fliplr(aa),time,0,2^nextpow2(length(tslot)));
r=sqrt((xi-xrec)^2+(yi-yrec)^2);
r=r./87.*0.66;
GG2=-1i.*w0./(8*B.*k.^2).*sqrt(2./(r.*k*pi)).*exp(-1i*(r.*k-pi/4)).*omega.^2;
GG2(1)=0;
GG2(isinf(GG2))=0;
GG(z,:)=GG(z,:)+GG2;
end
     else
if sum(abs(aa))>0
[w0,f0]=fftrl(fliplr(aa),time,0,2^nextpow2(length(tslot)));
r=sqrt((xi-xrec)^2+(yi-yrec)^2);
r=r./87.*0.66;
GG2p=-1i.*w0./(8*B.*k.^2).*sqrt(2./(r.*k*pi)).*exp(-1i*(r.*k-pi/4)).*omega.^2;
GG2p(1)=0;

GG2p(isinf(GG2p))=0;
GGp(z,:)=GGp(z,:)+GG2p;
end
     end

[w0,f0]=fftrl(aa2,time,0,2^nextpow2(length(tslot)));
r=sqrt((xi-xrec)^2+(yi-yrec)^2);
r=r./87.*0.66;
GG3=-1i.*w0./(8*B.*k.^2).*sqrt(2./(r.*k*pi)).*exp(-1i*(r.*k-pi/4)).*omega.^2;
GG3(1)=0;
GG3(isinf(GG3))=0;
GG1(z,:)=GG1(z,:)+GG3;
% [acc,t0]=ifftrl(GG2,f0);
% acc(isnan(acc))=0;
% SS(:,coco)=SS(:,coco)+acc;
% GG2(z,:)=GG2(z,:)+abs(-1i.*w0'./(8*B.*k.^2).*sqrt(2./(r.*k*pi)).*exp(-1i*(r.*k-pi/4)).*omega.^2).^2;
% coco=coco+1;
clear aa2 GG2 GG3 aa
end



     
end
    end
Psave=P2(:,:,end);
end


NN=floor(size(GG,2)/2+1);
GGup=GG(:,1:NN)+fliplr(conj(GG(:,NN:end)));
GGpup=GGp(:,1:NN)+fliplr(conj(GGp(:,NN:end)));
GGup1=GG1(:,1:NN)+fliplr(conj(GG1(:,NN:end)));
GGups=GGs(:,1:NN)+fliplr(conj(GGs(:,NN:end)));
powabs=abs(GGup(:,1:NN)).^2;
powabsp=abs(GGpup(:,1:NN)).^2;
powabs1=abs(GGup1(:,1:NN)).^2;
powabss=abs(GGups(:,1:NN)).^2;

powabsRev=abs(GGup(:,1:NN)).^2;
powabspRev=abs(GGpup(:,1:NN)).^2;
powabs1Rev=abs(GGup1(:,1:NN)).^2;
powabssRev=abs(GGups(:,1:NN)).^2;
f_up=linspace(f0(1),f0(end),size(powabs,2));
t=0:8e-3:(timewin-1)*8e-3;
figure
imagesc(t,f_up,log10(powabs'))
colorbar
set(gca,'YDir','normal')
set(gca,'fontsize',18)
title('Total Stress')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
% caxis([10 22])
caxis([13 19])
figure
imagesc(t,f_up,log10(powabss'))
set(gca,'YDir','normal')
set(gca,'fontsize',18)
title('Solid Stress')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
% caxis([10 22])
caxis([15 17])
figure
imagesc(t,f_up,log10(powabsp'))
set(gca,'fontsize',18)
title('Air in Pores')
set(gca,'YDir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
% caxis([10 22])
caxis([13 19])
figure
imagesc(t,f0,log10(powabs1'))
set(gca,'fontsize',18)
title('Air in Empty Space')
set(gca,'YDir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
% caxis([10 22])
caxis([13 19])







load expnocorr.mat
figure
imagesc(t,f0,log10(pow))
colorbar
set(gca,'YDir','normal')
set(gca,'fontsize',18)
title('Experimental Result')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
caxis([-5 1])

figure
imagesc(t,f0(1:411),log10(powabs1(:,1:411)'))
title('Air in Empty Space')
set(gca,'YDir','normal')
caxis([10 22])

figure
imagesc(t,f0(1:411),log10(powabsp(:,1:411)'))
title('Air in Pores')
set(gca,'YDir','normal')
caxis([10 22])

figure
imagesc(t,f0(1:411),log10(powabs(:,1:411)'))
title('Total Stress')
set(gca,'YDir','normal')
caxis([10 22])

figure
imagesc(t,f0(1:411),log10(powabss(:,1:411)'))
title('Solid Stress')
set(gca,'YDir','normal')
caxis([10 22])

save('Resultsall.mat')

