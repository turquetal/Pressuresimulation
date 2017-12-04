clear all
close all
clc 
addpath C:\Users\usr1\Dropbox\Super_Mega_Script_Folder
addpath D:\Dropbox\Dropbox\Super_Mega_Script_Folder
addpath C:\Users\User1\Dropbox\Super_Mega_Script_Folder
addpath C:\Users\Dhareion\Dropbox\Super_Mega_Script_Folder
addpath /Users/turquetal/Dropbox/Super_Mega_Script_Folder
% load '2bar_acceldata_161213.mat'
addpath ('/Users/turquetal/Dropbox/ExperimentalPics/2exp_C001H001S0001/pics');
FNAMEFMT = '2exp_C001H001S0001000%03d.bmp';
imglimitx=23:1002;
imglimity=186:560;
num_images=1;
% im=imread(sprintf(FNAMEFMT, 1));
% IMG=im2bw(im(186:560,23:1002));
fig=0; %1 if you want to plot figure
pts=[772 182;769 175; 775 167];% radiusofcurv of a finger experimentally
Rc = radiusCurvature(pts(:,1), pts(:,2));


pts4=[861 366;847 186; 855 8];% radiusofcurv of a finger experimentally
Rc4 = radiusCurvature(pts4(:,1), pts4(:,2));


pts20=[536 216;534 205; 539 193];% radiusofcurv of a finger experimentally
Rc20 = radiusCurvature(pts20(:,1), pts20(:,2));

pts198=[344 43;346 41; 349 44];% radiusofcurv of a finger experimentally
Rc198 = radiusCurvature(pts198(:,1), pts198(:,2));
L=Rc*(0.48/0.6)/(1-0.48/0.6)*0.1; % b.sandnes nature comm with a resize factor 0.1
L=ceil(L);
range=L;

rati=0.1;      %resize factor
pixels=870*rati;
L=L/pixels*0.66; %pixel to m correction;
row=floor(linspace(23,1002,6));
thr=[0.245 0.36 0.37 0.47 0.55];

% sound wave in porous medium (see ref winkler95porous)
Vair=343;
Vpoly1=2350; %longitudinal wave
Vpoly2=1120; %shear wave
Vpor=(0.4/Vair+0.6/Vpoly1)^-1;
load imdata.mat

siggen=0; % 0=no signal generation
% for  i=1:600
% im=imread(sprintf(FNAMEFMT, (i)));
% imr(:,:,i)=im(186:560,23:1002);
% 
% 
%  for j=1:length(row)-1
% img(:,row(j)-22:row(j+1)-22)=im2bw(im(186:560,row(j):row(j+1)),thr(j));
% 
%  end
% imbin(:,:,i)=img;
% 
% B=1-imresize(img,rati);
% 
% B1=zeros(size(B,1)+2,size(B,2)+2);
% bound=zeros(size(B,1)+2,size(B,2)+2);
% bound(1:size(B,1)+2,2:size(B,2)+2)=1;
% B1(2:size(B,1)+1,2:size(B,2)+1)=B;
% 
% 
% fin_img(:,:,i)=B1;
% end

% % % Diffusion Boundaries set to 2 in the Grid
B1=fin_img(:,:,1);
% B1= imresize(B1,2);
% B1=im2bw(B1);

B=zeros(size(B1,1)-2,size(B1,2)-2);
bound=ones(size(B,1)+2,size(B,2)+2);
% bound(1:size(B,1)+2,2:size(B,2)+2)=1;
% Grid Planning
% define pixel size
% in meters
c=343;
X=linspace(-0.4,0.35,size(B,2));
pix=abs(X(2)-X(1));
Y=linspace(-0.15,0.15,size(B,1));
piy=abs(Y(2)-Y(1));
% diffusion coefficient
fi_void=1;
fi_solid=0.55;
phi=0.52; %initial porosity
K_void=(0.001)^2/12;
K_loose=(80e-6)^2/180*0.52^3/0.48^2;
K_pore=(80e-6)^2/180*0.4^3/0.6^2;
P0=2e5;
mu=1.810e-5;
D1=1/mu/fi_void*K_void*P0;
D2=1/mu/fi_solid*K_pore*P0;
rho=1005;
g=10;
rho_air=1.225;
comp_air=1/(c^2*rho_air);
K=1; %Janssen's coefficient
dt=1e-5;
t=0:dt:8e-3;

r=(t(2)-t(1))/pix;
r2=(t(2)-t(1))/pix^2;





%%%%
U=zeros(size(B1,1),size(B1,2),2);
A=U;
U(2:end-1,end-1,2)=P0;
Pres=U;
cntr=2;
saver = 'airvib%02d.mat';
con=1;
soro=1;
cono=1;
% P1=ones(length(t),1);
% P1=awgn(P1,40).*P0;

% h = waitbar(0,'Hacking the system...');


% new = VideoWriter('airvib2','mpeg-4'); 
% open(new);
Pop=0;
for i=1:600

% waitbar(i / 600)
    clear B1
    B1=fin_img(:,:,i);
%     B1=interp2(fin_img(:,:,i),1);
% B1= imresize(B1,2);
%  B1=im2bw(B1);
    
    [ distmap ] = skindepth(1-B1);
Lmap=distmap./pixels.*0.66; %pixel to m correction
    if range
    [row_down,col_down]=find(distmap<=range);
    AA=zeros(size(distmap,1),size(distmap,2));
    rhomap=zeros(size(distmap,1),size(distmap,2));
    for rowi=1:length(row_down)
    AA(row_down(rowi),col_down(rowi))=distmap(row_down(rowi),col_down(rowi));
    end
    [row_2,col_2]=find(AA>0);
    [row_up,col_up]=find(distmap>range);
    for rowi=1:length(row_2)
    rhomap(row_2(rowi),col_2(rowi))=0.6;
    end
   for rowi=1:length(row_up)
    rhomap(row_up(rowi),col_up(rowi))=0.48;
   end
phimap=1-rhomap;
rhomap=rhomap.*1005/0.6;
rhorho(:,:,i)=rhomap;
Lnew=L-Lmap;
for xxi=1:size(B1,2)
   for yyi=1:size(B1,1)
       LL(yyi,xxi)=(xxi-1)./pixels.*0.66;
   end
end
    
    [row_up,col_up]=find(Lnew<0);
    for rowi=1:length(row_up)
    Lnew(row_up(rowi),col_up(rowi))=0;
    end
    [row_up1,col_up1]=find(Lmap<=L & Lmap>0);
    [row_up2,col_up2]=find(Lmap==0);
%     Lmap2=L-Lmap;
%     for rowi=1:length(row_up1)
%     Lmap2(row_up1(rowi),col_up1(rowi))=L;
%     end
%     for rowi=1:length(row_up2)
%     Lmap2(row_up2(rowi),col_up2(rowi))=0;
%     end
%     distmap=skindepth2(Lmap2,L);
%     Lmap2=-(distmap./87.*0.66); 
%         for rowi=1:length(row_up1)
%     Lmap2(row_up1(rowi),col_up1(rowi))=L;
%     end
%     for rowi=1:length(row_up2)
%     Lmap2(row_up2(rowi),col_up2(rowi))=0;
%     end
%     
    if i==1
        V(:,:,i)=zeros(size(B1,1),size(B1,2),1);
        V1(:,:,i)=zeros(size(B1,1),size(B1,2),1);
           delV(:,:,i)=gradient(V(:,:,i))./pix.*2;
        delV(delV<0)=0;
    else
        if i<600
%         Va=im2bw(interp2(fin_img(:,:,i+1),1));
%         Vb=im2bw(interp2(fin_img(:,:,i),1));
        Va=fin_img(:,:,i+1);
        Vb=fin_img(:,:,i);
        Vk=(Va-Vb);
            
            Vk(Vk<0)=0;
            [oy,ox]=find(Vk==1);
            
            for ceku=1:length(ox)
            [invx]=find(Vk(oy(ceku),ox(ceku):end)==0,1);
            if isempty(invx)
            else
            Vk(oy(ceku),ox(ceku))=invx-1;
            end
            end
           V(:,:,i)=Vk./0.008/pixels*0.66;
                delV(:,:,i)=gradient(V(:,:,i))./pix.*2;
        delV(delV<0)=0;
        else
        V(:,:,i)=zeros(size(B1,1),size(B1,2),1);
        delV(:,:,i)=gradient(V(:,:,i))./pix.*2;
        delV(delV<0)=0;
        end
    end
    
    con=con+1;
%     clear rhomap
    end
    sigxx=K.*rhomap.*g.*0.0015./(2*K).*((K*.30+1).*exp(2*0.30*K.*(Lnew)./0.0015)-1); %from sandnes nature

    

cntr=2;
BB=1-B1;
Vtrig=0;
for ti=2:length(t)-1
        delP(:,:,ti)=zeros(size(B1,1),size(B1,2),1);
        delP2(:,:,ti)=zeros(size(B1,1),size(B1,2),1);
         delP(2:end-1,2:end-1,ti)=gradient(U(2:end-1,2:end-1,ti).*bound(2:end-1,2:end-1))*2/pix;
         delP2(2:end-1,2:end-1,ti)=gradient(gradient((U(2:end-1,2:end-1,ti).*bound(2:end-1,2:end-1))))*4/pix^2;
    
        

for xi=2:size(B1,2)-1
for yi=2:size(B1,1)-1
       
       if B1(yi,xi)
                  
           KK=K_void;
           c=Vair;
           
% delP(yi,xi,ti)=(bound(yi,xi+1)*U(yi,xi+1,ti)...
%              +bound(yi+1,xi)*U(yi+1,xi,ti)...
%             -2*bound(yi,xi)*U(yi,xi,ti))/pix;
% % dP1(yi,xi)=delP;        
%  delP2(yi,xi,ti)=(bound(yi,xi+1)*U(yi,xi+1,ti)...
%             +bound(yi,xi-1)*U(yi,xi-1,ti)...
%             +bound(yi+1,xi)*U(yi+1,xi,ti)...
%             +bound(yi-1,xi)*U(yi-1,xi,ti)...
%             -4*bound(yi,xi)*U(yi,xi,ti))/pix^2;
        
Phat=U(yi,xi,ti)-P0+1/comp_air;
rho_f=Phat*rho_air*comp_air;
U(yi,xi,ti+1)=dt^2*Phat/rho_f*delP2(yi,xi,ti)+2*U(yi,xi,ti)-U(yi,xi,ti-1);
%         D=Phat*KK/mu/phimap(yi,xi);
%         omom=rho_f*KK/mu/phimap(yi,xi);
%         TTT=(-V(yi,xi)*delP+D*(delP2)-Phat/(phimap(yi,xi))*delV);
%         U(yi,xi,ti+1)=1/(dt+omom)*(TTT*dt^2+omom*(2*U(yi,xi,ti)-U(yi,xi,ti-1))+U(yi,xi,ti)*dt);
% %         

       if isnan(U(yi,xi,ti+1))
    pause
    end
       if isinf(U(yi,xi,ti+1))
    pause
       end
       
        else
                    
           
           c=Vpor;
           
           if rhorho(yi,xi,i)==1005
           KK=K_pore;
           else
              KK=K_loose; 
           end
% delVt=BB(yi,xi)*rhomap(yi,xi)*(1-phi)/phi*(dv1+dvx1-dva-dvxa)/0.008;

% Putting the new equation in 11 Oct 2016
% delP(yi,xi,ti)=(bound(yi,xi+1)*U(yi,xi+1,ti)...
%              +bound(yi+1,xi)*U(yi+1,xi,ti)...
%             -2*bound(yi,xi)*U(yi,xi,ti))/pix;
% % dP1(yi,xi)=delP;        
%  delP2(yi,xi,ti)=(bound(yi,xi+1)*U(yi,xi+1,ti)...
%             +bound(yi,xi-1)*U(yi,xi-1,ti)...
%             +bound(yi+1,xi)*U(yi+1,xi,ti)...
%             +bound(yi-1,xi)*U(yi-1,xi,ti)...
%             -4*bound(yi,xi)*U(yi,xi,ti))/pix^2;
% dP2(yi,xi)=delP; 
% Phat=U(yi,xi,ti)-P0+1/comp_air;
% find nearest pressure in gap
Phat=U(yi,xi+1,ti-1);
% cheko=U(:,:,ti).*B1;
% icx=find(abs(cheko(yi,xi:end))>0,1);
% 
% Phat=cheko(yi,xi+icx-1);
rho_f=Phat*rho_air*comp_air;
if Vtrig==0
%     V(yi,xi,i)=0;
    delVtr=delV(yi,xi,i);
    Vtr=V(yi,xi);
    Vtrig=1;
else
    Vtr=0;
    delVtr=0;
end
        D=Phat*KK/mu/phimap(yi,xi);       
%         TTT=(-V(yi,xi)*delP+D*(delP2)-Phat/phimap(yi,xi)*delV);
%         U(yi,xi,ti+1)=(TTT*dt+U(yi,xi,ti));
        omom=rho_f^2/rho_air*KK/mu/phimap(yi,xi);

        
        TTT=(-Vtr*delP(yi,xi,ti)+D*(delP2(yi,xi,ti))-Phat/(phimap(yi,xi))*delVtr);
%         contV(yi,xi,i)=D*(delP2);
%         contD(yi,xi,i)=-V(yi,xi)*delP-Phat/(phimap(yi,xi))*delV;
        U(yi,xi,ti+1)=1/(dt+omom)*(TTT*dt^2+omom*(2*U(yi,xi,ti)-U(yi,xi,ti-1))+U(yi,xi,ti)*dt);
%         

       if isnan(U(yi,xi,ti+1))
    pause
    end
       if isinf(U(yi,xi,ti+1))
    pause
    end       
           
        
        
        end

%         U(yi,xi,ti+1)=u;
%         clear u
    end
end
% imagesc(U(:,:,ti+1))
% pause
% close all
U(1,:,ti+1)=U(2,:,ti+1);
U(end,:,ti+1)=U(end-1,:,ti+1);
U(:,end,ti+1)=U(:,end-1,ti+1);
% cntr=cntr+1;
% P2=ones(38,1).*P0;
U(2:end-1,end-1,end)=P0;
U(2:end-1,1:2,end)=0;
Pres(:,:,ti)=U(:,:,end);
    for iii=1:size(Pres,3)
    efsig(:,:,iii)=sigxx+Pres(:,:,iii);
    end
end
% close all
% figure
% imagesc(U(:,:,200))
% pause(0.1)
% close all


svname=['neweq_efsig_U_test',num2str(i)];
save(svname,'sigxx','U')
inU(:,:,1:2)=U(:,:,end-1:end);
clear U
U=inU;

clear inU

if siggen
% pause
xrec=19;
yrec=19;
% tic

signal=zeros(8193,1);
coco=1;
for yi=2:size(B1,1)-1
    for xi=2:size(B1,2)-1
  
        source=squeeze(efsig(yi,xi,:));
        source(1)=source(2);
        time=0:1e-5:(length(source)-1)*1e-5;
        radd=sqrt((xi-xrec)^2+(yi-yrec)^2);
        radd=radd./pixels.*0.66;
        inin=find(source<0);
        source(inin)=0;
         source=source.*pix^2;
        if sum(source)
       
        [acc,t0]=lambgeneratornoref(source,time,radd);
        if isnan(sum(acc))
        else
           signal(1:length(acc),coco)=acc;
           
        coco=coco+1;
      
        end
           clear source time
        else
           clear source time
        end
        

    end
end

if sum(sum(signal))
  Signal=sum(signal,2);
         svname=['completesig_',num2str(i)];
    save(svname,'Signal','t0')
end
end
% plot(Signal)
% pause
cntr=2;
% toc
% contourf(log(efsig),log(c));
% caxis(log([c(1) c(length(c))]));
% colorbar('FontSize',11,'YTick',log(c),'YTickLabel',c);
if fig
figure
set(gcf,'units','normalized','outerposition',[0 0 1 1])

subplot(2,2,1)
set(gca,'FontSize',30)


imshow(imr(:,:,i))
%  caxis([-15e4 15e4])
% colormap('gray')
freezeColors
% cbfreeze
title(['Experimental Pictures t=',num2str(i*0.008),'s'])
xlabel('Length(m)')
ylabel('Length(m)')

% subplot(3,2,4)

subplot(2,2,2)
set(gca,'FontSize',30)

imagesc(sigxx)
caxis([-15e4 15e4])
% caxis([-7e4 15e4])
% contourf(log(sigxx),log(c));
% caxis(log([c(1) c(length(c))]));
colormap('jet')
title('Total Stress of The Medium');
% title(num2str(i+1))
% colorbar
xlabel('Length(m)')
ylabel('Length(m)')
% caxis([2e-3 3e-3])

freezeColors
cbfreeze(colorbar)

subplot(2,2,3)
set(gca,'FontSize',30)
imagesc(efsig(:,:,i))
caxis([-15e4 15e4])
% contourf(log(efsig),log(c));
% caxis(log([c(1) c(length(c))]));

title('Effective Stress of The Medium');
% title(num2str(i+1))
% colorbar
xlabel('Length(m)')
ylabel('Length(m)')
% caxis([2e-3 3e-3])
% colormap('pink')
freezeColors
cbfreeze(colorbar)

subplot(2,2,4)
set(gca,'FontSize',30)


% contourf(log(Pres),log(c));
% caxis(log([c(1) c(length(c))]));
imagesc(Pres(:,:,end))
caxis([-15e4 15e4])
title('Pore Fluid Pressure');
% title(num2str(i+1))
% colorbar
xlabel('Length(m)')
ylabel('Length(m)')
% caxis([2e-3 3e-3])
% colormap('pink')
freezeColors
cbfreeze(colorbar)

sname=['movie',num2str(i+1),'png'];
% export_fig(gcf,sname)
saveas(1,num2str(i+1),'png')
close all
end
end
% close(new)
% close(h) 
% save('efsig.mat','efsig')
% 
% c=logspace(1,7);
% contourf(log(efsig),log(c));
% caxis(log([c(1) c(length(c))]));
% colorbar('FontSize',11,'YTick',log(c),'YTickLabel',c);

