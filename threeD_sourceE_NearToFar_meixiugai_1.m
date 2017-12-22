clear;clc;
%**********************************************************
pi=3.1415926;
muz=4*pi*1e-7;
eps0=8.851*1e-12;
c=1/sqrt(muz*eps0);
f=1e8;%%时谐激励源的频率
k=2*pi*f*sqrt(muz*eps0);%%波数
Z=sqrt(muz/eps0);%%阻抗
wave_length=c/f;%%波长
ds=wave_length/80.0;
dt=ds/c/2.0;
TT=1.0/f/dt;%%一个周期的时间步数
%**********************************************************
NX=80;NY=80;NZ=80;AA=2;nT=5;nmax=round(TT*nT);      %外推边界距离吸收边界的距离为AA,时间循环为nT个周期，一般3-5个周期达到稳定
%**********************************************************     
Hx=zeros(NX+1,NY,NZ);Hy=zeros(NX,NY+1,NZ);Hz=zeros(NX,NY,NZ+1);
Ex=zeros(NX,NY+1,NZ+1);Ey=zeros(NX+1,NY,NZ+1);Ez=zeros(NX+1,NY+1,NZ);
Ez_near_analytical(nmax)=0;
Ez_near_save(nmax)=0;
%**********************************************************
hx_1=zeros(NX+1,NY,NZ);hy_1=zeros(NX,NY+1,NZ);hz_1=zeros(NX,NY,NZ+1);  %三个时刻的场（边界问题）
hx_2=zeros(NX+1,NY,NZ);hy_2=zeros(NX,NY+1,NZ);hz_2=zeros(NX,NY,NZ+1);
hx_3=zeros(NX+1,NY,NZ);hy_3=zeros(NX,NY+1,NZ);hz_3=zeros(NX,NY,NZ+1);
ex_1=zeros(NX,NY+1,NZ+1);ey_1=zeros(NX+1,NY,NZ+1);ez_1=zeros(NX+1,NY+1,NZ);
ex_2=zeros(NX,NY+1,NZ+1);ey_2=zeros(NX+1,NY,NZ+1);ez_2=zeros(NX+1,NY+1,NZ);
ex_3=zeros(NX,NY+1,NZ+1);ey_3=zeros(NX+1,NY,NZ+1);ez_3=zeros(NX+1,NY+1,NZ);
%**********************************************************
midx=fix(NX/2+1);midy=fix(NY/2+1);midz=fix(NZ/2+1);%%中央位置
i0=1+AA;j0=1+AA;k0=1+AA;ia=NX+1-AA;jb=NY+1-AA;kc=NZ+1-AA;  %各个方向外推（输出）边界的位置
%%%%_x1,_x2分别代表x向的两个输出面；_y1,_y2分别代表y向的两个输出面；_z1,_z2分别代表z向的两个输出面
FMz_x1=0;FMy_x1=0;FMz_y1=0;FMx_y1=0;FMy_z1=0;FMx_z1=0;%%各个输出面上磁流的各个分量（沿面的切向）的复数形式
FMz_x2=0;FMy_x2=0;FMz_y2=0;FMx_y2=0;FMy_z2=0;FMx_z2=0;
FJz_x1=0;FJy_x1=0;FJz_y1=0;FJx_y1=0;FJy_z1=0;FJx_z1=0;%%各个输出面上电流的各个分量（沿面的切向）的复数形式
FJz_x2=0;FJy_x2=0;FJz_y2=0;FJx_y2=0;FJy_z2=0;FJx_z2=0;
%**********************************************************
r0=1/3*ds;c1=dt/muz/(ds*ds-pi*r0*r0/4);c2=ds/2*log(ds/r0);ll=9;
is=midx;js=midy;ks=midz;%%时谐激励源的位置
%**********************************************************
tic

for n=1:nmax
      %-------------------电磁场循环---------------------------
    hx_3=hx_2;hx_2=hx_1;hx_1=Hx;hy_3=hy_2;hy_2=hy_1;hy_1=Hy;hz_3=hz_2;hz_2=hz_1;hz_1=Hz;
    Hx=Hx+dt/muz/ds*(Ey(:,:,2:NZ+1)-Ey(:,:,1:NZ)+Ez(:,1:NY,:)-Ez(:,2:NY+1,:));
    Hy=Hy+dt/muz/ds*(Ez(2:NX+1,:,:)-Ez(1:NX,:,:)+Ex(:,:,1:NZ)-Ex(:,:,2:NZ+1));
    Hz=Hz+dt/muz/ds*(Ex(:,2:NY+1,:)-Ex(:,1:NY,:)+Ey(1:NX,:,:)-Ey(2:NX+1,:,:));
      
    ex_3=ex_2;ex_2=ex_1;ex_1=Ex;ey_3=ey_2;ey_2=ey_1;ey_1=Ey;ez_3=ez_2;ez_2=ez_1;ez_1=Ez;
    Ex(:,2:NY,2:NZ)=Ex(:,2:NY,2:NZ)+dt/eps0/ds*(Hz(:,2:NY,2:NZ)-Hz(:,1:NY-1,2:NZ)+Hy(:,2:NY,1:NZ-1)-Hy(:,2:NY,2:NZ));
    Ey(2:NX,:,2:NZ)=Ey(2:NX,:,2:NZ)+dt/eps0/ds*(Hx(2:NX,:,2:NZ)-Hx(2:NX,:,1:NZ-1)+Hz(1:NX-1,:,2:NZ)-Hz(2:NX,:,2:NZ));
    Ez(2:NX,2:NY,:)=Ez(2:NX,2:NY,:)+dt/eps0/ds*(Hy(2:NX,2:NY,:)-Hy(1:NX-1,2:NY,:)+Hx(2:NX,1:NY-1,:)-Hx(2:NX,2:NY,:));
    
    r_near=38*ds;
    Ez(is,js,ks)=Ez(is,js,ks)+sin(2*pi*f*n*dt);  %加的是正弦信号源
    Ez_near_analytical(n)=muz/(r_near*4*pi)*(-2*c/r_near*eps0*ds^3/dt*sin(2*pi*f*(n*dt-r_near/c))+2*c^2/r_near^2*eps0*ds^3/dt/2/pi/f*cos(2*pi*f*(n*dt-r_near/c)));
    Ez_near_save(n)=Ez(is,js,ks+38);
    %************吸收边界条件*************
    Ey(1,:,:)=3*ey_1(1+1,:,:)-3*ey_2(1+2,:,:)+ey_3(1+3,:,:);
    Ey(NX+1,:,:)=3*ey_1(NX+1-1,:,:)-3*ey_2(NX+1-2,:,:)+ey_3(NX+1-3,:,:);
    Ez(1,:,:)=3*ez_1(1+1,:,:)-3*ez_2(1+2,:,:)+ez_3(1+3,:,:);
    Ez(NX+1,:,:)=3*ez_1(NX+1-1,:,:)-3*ez_2(NX+1-2,:,:)+ez_3(NX+1-3,:,:);
    Ex(:,1,:)=3*ex_1(:,1+1,:)-3*ex_2(:,1+2,:)+ex_3(:,1+3,:);
    Ex(:,NY+1,:)=3*ex_1(:,NY+1-1,:)-3*ex_2(:,NY+1-2,:)+ex_3(:,NY+1-3,:);
    Ez(:,1,:)=3*ez_1(:,1+1,:)-3*ez_2(:,1+2,:)+ez_3(:,1+3,:);
    Ez(:,NY+1,:)=3*ez_1(:,NY+1-1,:)-3*ez_2(:,NY+1-2,:)+ez_3(:,NY+1-3,:);
    Ex(:,:,1)=3*ex_1(:,:,1+1)-3*ex_2(:,:,1+2)+ex_3(:,:,1+3);
    Ex(:,:,NZ+1)=3*ex_1(:,:,NZ+1-1)-3*ex_2(:,:,NZ+1-2)+ex_3(:,:,NZ+1-3);
    Ey(:,:,1)=3*ey_1(:,:,1+1)-3*ey_2(:,:,1+2)+ey_3(:,:,1+3);
    Ey(:,:,NZ+1)=3*ey_1(:,:,NZ+1-1)-3*ey_2(:,:,NZ+1-2)+ey_3(:,:,NZ+1-3);
    
%***********存储最后一个周期输出面上的切向电磁场分量**************
%%%%保证各个面上的切向电、磁场量均采在一个点上，如在x1和x2两个面（位于网格边）上，y和z向均采在网格中央，即(i,j+1/2,k+1/2)
    if n>TT*(nT-1)        
        ey_x1=(Ey(i0,j0:jb-1,k0:kc-1)+Ey(i0,j0:jb-1,k0+1:kc))/2;%%%(i,j+1/2,k+1/2)
        ez_x1=(Ez(i0,j0:jb-1,k0:kc-1)+Ez(i0,j0+1:jb,k0:kc-1))/2;%%%(i,j+1/2,k+1/2)
        ex_y1=(Ex(i0:ia-1,j0,k0:kc-1)+Ex(i0:ia-1,j0,k0+1:kc))/2;%%%(i+1/2,j,k+1/2)
        ez_y1=(Ez(i0:ia-1,j0,k0:kc-1)+Ez(i0+1:ia,j0,k0:kc-1))/2;%%%(i+1/2,j,k+1/2)
        ex_z1=(Ex(i0:ia-1,j0:jb-1,k0)+Ex(i0:ia-1,j0+1:jb,k0))/2;%%%(i+1/2,j+1/2,k)
        ey_z1=(Ey(i0:ia-1,j0:jb-1,k0)+Ey(i0+1:ia,j0:jb-1,k0))/2;%%%(i+1/2,j+1/2,k)
        
        ey_x2=(Ey(ia,j0:jb-1,k0:kc-1)+Ey(ia,j0:jb-1,k0+1:kc))/2;
        ez_x2=(Ez(ia,j0:jb-1,k0:kc-1)+Ez(ia,j0+1:jb,k0:kc-1))/2;
        ex_y2=(Ex(i0:ia-1,jb,k0:kc-1)+Ex(i0:ia-1,jb,k0+1:kc))/2;
        ez_y2=(Ez(i0:ia-1,jb,k0:kc-1)+Ez(i0+1:ia,jb,k0:kc-1))/2;
        ex_z2=(Ex(i0:ia-1,j0:jb-1,kc)+Ex(i0:ia-1,j0+1:jb,kc))/2;
        ey_z2=(Ey(i0:ia-1,j0:jb-1,kc)+Ey(i0+1:ia,j0:jb-1,kc))/2;
        
        hy_x1=(Hy(i0,j0:jb-1,k0:kc-1)+Hy(i0-1,j0:jb-1,k0:kc-1)+Hy(i0,j0+1:jb,k0:kc-1)+Hy(i0-1,j0+1:jb,k0:kc-1))/4;%%%(i,j+1/2,k+1/2)
        hz_x1=(Hz(i0,j0:jb-1,k0:kc-1)+Hz(i0-1,j0:jb-1,k0:kc-1)+Hz(i0,j0:jb-1,k0+1:kc)+Hz(i0-1,j0:jb-1,k0+1:kc))/4;%%%(i,j+1/2,k+1/2)
        hx_y1=(Hx(i0:ia-1,j0,k0:kc-1)+Hx(i0:ia-1,j0-1,k0:kc-1)+Hx(i0+1:ia,j0,k0:kc-1)+Hx(i0+1:ia,j0-1,k0:kc-1))/4;%%%(i+1/2,j,k+1/2)
        hz_y1=(Hz(i0:ia-1,j0,k0:kc-1)+Hz(i0:ia-1,j0-1,k0:kc-1)+Hz(i0:ia-1,j0,k0+1:kc)+Hz(i0:ia-1,j0-1,k0+1:kc))/4;%%%(i+1/2,j,k+1/2)
        hx_z1=(Hx(i0:ia-1,j0:jb-1,k0)+Hx(i0:ia-1,j0:jb-1,k0-1)+Hx(i0+1:ia,j0:jb-1,k0)+Hx(i0+1:ia,j0:jb-1,k0-1))/4;%%%(i+1/2,j+1/2,k)
        hy_z1=(Hy(i0:ia-1,j0:jb-1,k0)+Hy(i0:ia-1,j0:jb-1,k0-1)+Hy(i0:ia-1,j0+1:jb,k0)+Hy(i0:ia-1,j0+1:jb,k0-1))/4;%%%(i+1/2,j+1/2,k)
        
        hy_x2=(Hy(ia,j0:jb-1,k0:kc-1)+Hy(ia-1,j0:jb-1,k0:kc-1)+Hy(ia,j0+1:jb,k0:kc-1)+Hy(ia-1,j0+1:jb,k0:kc-1))/4;
        hz_x2=(Hz(ia,j0:jb-1,k0:kc-1)+Hz(ia-1,j0:jb-1,k0:kc-1)+Hz(ia,j0:jb-1,k0+1:kc)+Hz(ia-1,j0:jb-1,k0+1:kc))/4;
        hx_y2=(Hx(i0:ia-1,jb,k0:kc-1)+Hx(i0:ia-1,jb-1,k0:kc-1)+Hx(i0+1:ia,jb,k0:kc-1)+Hx(i0+1:ia,jb-1,k0:kc-1))/4;
        hz_y2=(Hz(i0:ia-1,jb,k0:kc-1)+Hz(i0:ia-1,jb-1,k0:kc-1)+Hz(i0:ia-1,jb,k0+1:kc)+Hz(i0:ia-1,jb-1,k0+1:kc))/4;
        hx_z2=(Hx(i0:ia-1,j0:jb-1,kc)+Hx(i0:ia-1,j0:jb-1,kc-1)+Hx(i0+1:ia,j0:jb-1,kc)+Hx(i0+1:ia,j0:jb-1,kc-1))/4;
        hy_z2=(Hy(i0:ia-1,j0:jb-1,kc)+Hy(i0:ia-1,j0:jb-1,kc-1)+Hy(i0:ia-1,j0+1:jb,kc)+Hy(i0:ia-1,j0+1:jb,kc-1))/4;
        
        ey_x1=squeeze(ey_x1);ey_x2=squeeze(ey_x2);
        ez_x1=squeeze(ez_x1);ez_x2=squeeze(ez_x2);
        ex_y1=squeeze(ex_y1);ex_y2=squeeze(ex_y2);
        ez_y1=squeeze(ez_y1);ez_y2=squeeze(ez_y2);
        ex_z1=squeeze(ex_z1);ex_z2=squeeze(ex_z2);
        ey_z1=squeeze(ey_z1);ey_z2=squeeze(ey_z2);
        
        hy_x1=squeeze(hy_x1);hy_x2=squeeze(hy_x2);
        hz_x1=squeeze(hz_x1);hz_x2=squeeze(hz_x2);
        hx_y1=squeeze(hx_y1);hx_y2=squeeze(hx_y2);
        hz_y1=squeeze(hz_y1);hz_y2=squeeze(hz_y2);
        hx_z1=squeeze(hx_z1);hx_z2=squeeze(hx_z2);
        hy_z1=squeeze(hy_z1);hy_z2=squeeze(hy_z2);
        %************************* 各个输出面上的磁流密度（表7-1）瞬时值      
        Mz_x1=ey_x1;Mz_x2=-ey_x2;
        My_x1=-ez_x1;My_x2=ez_x2;
        Mz_y1=-ex_y1;Mz_y2=ex_y2;
        Mx_y1=ez_y1;Mx_y2=-ez_y2;
        My_z1=ex_z1;My_z2=-ex_z2;
        Mx_z1=-ey_z1;Mx_z2=ey_z2;
        %************************* 各个输出面上的电流密度（表7-1）瞬时值    
        Jz_x1=-hy_x1;Jz_x2=hy_x2;
        Jy_x1=hz_x1;Jy_x2=-hz_x2;
        Jz_y1=hx_y1;Jz_y2=-hx_y2;
        Jx_y1=-hz_y1;Jx_y2=hz_y2;
        Jy_z1=-hx_z1;Jy_z2=hx_z2;
        Jx_z1=hy_z1;Jx_z2=-hy_z2;
        %*************************%傅立叶变换，得到各个输出面上的电流/磁流密度的幅值和相位 
        phex=exp(-1i*2*pi*(n-TT*(nT-1))*f*dt);
        FMz_x1=FMz_x1+Mz_x1*phex;FMz_x2=FMz_x2+Mz_x2*phex;             
        FMy_x1=FMy_x1+My_x1*phex;FMy_x2=FMy_x2+My_x2*phex;
        FMz_y1=FMz_y1+Mz_y1*phex;FMz_y2=FMz_y2+Mz_y2*phex;
        FMx_y1=FMx_y1+Mx_y1*phex;FMx_y2=FMx_y2+Mx_y2*phex;
        FMy_z1=FMy_z1+My_z1*phex;FMy_z2=FMy_z2+My_z2*phex;
        FMx_z1=FMx_z1+Mx_z1*phex;FMx_z2=FMx_z2+Mx_z2*phex;
        
        phex=exp(-1i*2*pi*(n-1/2-TT*(nT-1))*f*dt);
        FJz_x1=FJz_x1+Jz_x1*phex;FJz_x2=FJz_x2+Jz_x2*phex;
        FJy_x1=FJy_x1+Jy_x1*phex;FJy_x2=FJy_x2+Jy_x2*phex;
        FJz_y1=FJz_y1+Jz_y1*phex;FJz_y2=FJz_y2+Jz_y2*phex;
        FJx_y1=FJx_y1+Jx_y1*phex;FJx_y2=FJx_y2+Jx_y2*phex;
        FJy_z1=FJy_z1+Jy_z1*phex;FJy_z2=FJy_z2+Jy_z2*phex;
        FJx_z1=FJx_z1+Jx_z1*phex;FJx_z2=FJx_z2+Jx_z2*phex;
    end
    
end
iii=1:nmax;
plot(iii,Ez_near_analytical,iii,Ez_near_save,'-.');
%**********************************************************%积分变量设置
%%%各个面均采在网格边上
x1_x=repmat(i0*ds,[jb-j0 kc-k0]);x2_x=repmat(ia*ds,[jb-j0 kc-k0]);
y1_y=repmat(j0*ds,[ia-i0 kc-k0]);y2_y=repmat(jb*ds,[ia-i0 kc-k0]);
z1_z=repmat(k0*ds,[ia-i0 jb-j0]);z2_z=repmat(kc*ds,[ia-i0 jb-j0]);

%%%各个面上的其它两向均采在网格中央，即x的两个面上的节点数为（YY-SY）*（ZZ-SZ）
nx=((i0+1/2):(ia-1/2))*ds;     
ny=((j0+1/2):(jb-1/2))*ds;
nz=((k0+1/2):(kc-1/2))*ds;

x1_y=repmat(ny',[1 kc-k0]);x2_y=repmat(ny',[1 kc-k0]);
y1_x=repmat(nx',[1 kc-k0]);y2_x=repmat(nx',[1 kc-k0]);
z1_x=repmat(nx',[1 jb-j0]);z2_x=repmat(nx',[1 jb-j0]);

x1_z=repmat(nz,[jb-j0 1]);x2_z=repmat(nz,[jb-j0 1]);
y1_z=repmat(nz,[ia-i0 1]);y2_z=repmat(nz,[ia-i0 1]);
z1_y=repmat(ny,[ia-i0 1]);z2_y=repmat(ny,[ia-i0 1]);
%************************
r=500;
step_angel=5;
elevation=0:step_angel:180;      
elevation=elevation*pi/180;
azimuth=0:step_angel:360;
azimuth=azimuth*pi/180;
Esita_far=zeros(length(elevation),length(azimuth));
Esita_far_analytical=zeros(length(elevation),length(azimuth));
x=zeros(length(elevation),length(azimuth));
y=zeros(length(elevation),length(azimuth));
z=zeros(length(elevation),length(azimuth));

for i_an=1:length(elevation)
    theta=elevation(i_an);
    for j_an=1:length(azimuth)
        fai=azimuth(j_an);    
        index1=1i*k*sin(theta)*cos(fai);
        index2=1i*k*sin(theta)*sin(fai);
        index3=1i*k*cos(theta);
       
       fz_x1= sum(sum(FJz_x1.*exp(index1*x1_x+index2*x1_y+index3*x1_z)*ds*ds));        %x向第一个面
       fy_x1= sum(sum(FJy_x1.*exp(index1*x1_x+index2*x1_y+index3*x1_z)*ds*ds)); 
       fmz_x1=sum(sum(FMz_x1.*exp(index1*x1_x+index2*x1_y+index3*x1_z)*ds*ds));
       fmy_x1=sum(sum(FMy_x1.*exp(index1*x1_x+index2*x1_y+index3*x1_z)*ds*ds));
       
       fz_x2= sum(sum(FJz_x2.*exp(index1*x2_x+index2*x2_y+index3*x2_z)*ds*ds));       %x向第二个面
       fy_x2= sum(sum(FJy_x2.*exp(index1*x2_x+index2*x2_y+index3*x2_z)*ds*ds)); 
       fmz_x2=sum(sum(FMz_x2.*exp(index1*x2_x+index2*x2_y+index3*x2_z)*ds*ds));
       fmy_x2=sum(sum(FMy_x2.*exp(index1*x2_x+index2*x2_y+index3*x2_z)*ds*ds));
       
       fx_y1= sum(sum(FJx_y1.*exp(index1*y1_x+index2*y1_y+index3*y1_z)*ds*ds));        %y
       fz_y1= sum(sum(FJz_y1.*exp(index1*y1_x+index2*y1_y+index3*y1_z)*ds*ds));
       fmx_y1=sum(sum(FMx_y1.*exp(index1*y1_x+index2*y1_y+index3*y1_z)*ds*ds));
       fmz_y1=sum(sum(FMz_y1.*exp(index1*y1_x+index2*y1_y+index3*y1_z)*ds*ds));
       
       fx_y2= sum(sum(FJx_y2.*exp(index1*y2_x+index2*y2_y+index3*y2_z)*ds*ds));  
       fz_y2= sum(sum(FJz_y2.*exp(index1*y2_x+index2*y2_y+index3*y2_z)*ds*ds));
       fmx_y2=sum(sum(FMx_y2.*exp(index1*y2_x+index2*y2_y+index3*y2_z)*ds*ds));
       fmz_y2=sum(sum(FMz_y2.*exp(index1*y2_x+index2*y2_y+index3*y2_z)*ds*ds));
       
       fx_z1= sum(sum(FJx_z1.*exp(index1*z1_x+index2*z1_y+index3*z1_z)*ds*ds));         %z
       fy_z1= sum(sum(FJy_z1.*exp(index1*z1_x+index2*z1_y+index3*z1_z)*ds*ds));
       fmx_z1=sum(sum(FMx_z1.*exp(index1*z1_x+index2*z1_y+index3*z1_z)*ds*ds));
       fmy_z1=sum(sum(FMy_z1.*exp(index1*z1_x+index2*z1_y+index3*z1_z)*ds*ds));
       
       fx_z2= sum(sum(FJx_z2.*exp(index1*z2_x+index2*z2_y+index3*z2_z)*ds*ds));
       fy_z2= sum(sum(FJy_z2.*exp(index1*z2_x+index2*z2_y+index3*z2_z)*ds*ds));
       fmx_z2=sum(sum(FMx_z2.*exp(index1*z2_x+index2*z2_y+index3*z2_z)*ds*ds));
       fmy_z2=sum(sum(FMy_z2.*exp(index1*z2_x+index2*z2_y+index3*z2_z)*ds*ds));
       
       fx=fx_y1+fx_y2+fx_z1+fx_z2;                                     
       fmx=fmx_y1+fmx_y2+fmx_z1+fmx_z2;
       fy=fy_x1+fy_x2+fy_z1+fy_z2;
       fmy=fmy_x1+fmy_x2+fmy_z1+fmy_z2;
       fz=fz_x1+fz_x2+fz_y1+fz_y2;
       fmz=fmz_x1+fmz_x2+fmz_y1+fmz_y2;
       Esita_far(i_an,j_an)=-1i*k*exp(-1i*k*r)/(4*pi*r)*(Z*(fx*cos(theta)*cos(fai)+fy*cos(theta)*sin(fai)-fz*sin(theta))...
                            -fmx*sin(fai)+fmy*cos(fai));
       Esita_far_analytical(i_an,j_an)=1i*(1i*2*pi*eps0*ds^3/dt)/(2*r*wave_length)*Z*sin(theta)*exp(-1i*k*r);
        x(i_an,j_an)=sin(theta)*cos(fai);
        y(i_an,j_an)=sin(theta)*sin(fai);
        z(i_an,j_an)=cos(theta);                
    end

end 
Esita_far=abs(Esita_far);
Esita_far_analytical=abs(Esita_far_analytical);
    
%    E=abs(Z*f_theta+f_m_phi).*abs(Z*f_theta+f_m_phi) +abs(-Z*f_phi+f_m_theta).*abs(-Z*f_phi+f_m_theta);
%    E=E/max(max(E));
  %***********************************************************************
  figure(1);            %E面方向性图
  E_e1=squeeze(Esita_far(1,:)); E_e2=squeeze(Esita_far(fix(length(azimuth)/2),:));
  polar([azimuth  -azimuth],[E_e1 E_e2],'r');
%   polar(theta,r_h1,'r');
   
  figure(2);             %H面方向性图
  E_h=squeeze(Esita_far(:,fix(length(elevation)/2)));
  polar(elevation,E_h','r');
  
  figure(3);
x1=x.*Esita_far_analytical;
y1=y.*Esita_far_analytical;
z1=z.*Esita_far_analytical;
surf(x1,y1,z1);  
  figure(4);
x2=x.*Esita_far;
y2=y.*Esita_far;
z2=z.*Esita_far;
surf(x2,y2,z2);