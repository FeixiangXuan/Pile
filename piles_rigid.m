% **********************************************************************
% Linear interpolation of input data along plie axis
% **********************************************************************

% assignment of figure size
% portatile
dim=[1 31 1366 662];
% fisso
%dim=[1 41 1600 784];


% reading of P-y curves
Pyread=load(filePy);
if size(Pyread,1)<2
    disp(' ')
    disp('ERROR: at least two P-y curves at different depth must be inserted')
    disp('ABORT calculation')
    disp(' ')
    return
end


% depth of pile nodes
z=[0:Lp/N:Lp]';

% discretization of w0 and u0
d1=[];
d1=interp1(displ(:,1),displ(:,2),[0:1:max(displ(:,1))]');
d1=[d1 interp1(displ(:,1),displ(:,3),[0:1:max(displ(:,1))]')];

% interpolation of soil properties
% unit weight of the soil [kN/m3]
gter  = interp1(Data_soil(:,1),Data_soil(:,2),z);
% friction angle of the soil [rad]
% Young modulus E of the soil [MPa]
Eter = interp1(Data_soil(:,1),Data_soil(:,3),z);
% Poisson coefficient of the soil [-]
nuter = interp1(Data_soil(:,1),Data_soil(:,4),z);
% Shear modulus of the soil [Mpa]
Gter = Eter./(2*(1+nuter));
% unit shaft resistance [kPa]
fs = interp1(Data_soil(:,1),Data_soil(:,5),z);
% Shear modulus of the soil at the base [Mpa]
Gb = Eb/(2*1+nub);


figure(1)
clf
set(gcf,'Position',dim)
set(gcf,'Color',[1 1 1])
set(gcf,'Name','Input Data')

t=subplot(3,3,1);
axis off
text(0,7/7,['lenght of the pile: L_p = ',num2str(Lp),' m'])
text(0,6/7,['diameter of the pile: D_p = ',num2str(Dp),' m'])
text(0,5/7,['water table depth: z_w = ',num2str(zw),' m'])
text(0,4/7,['water unit weight: g_w = ',num2str(gw),' kN/m^3'])
text(0,3/7,['top rotational stiffness: K_t = ',num2str(Ktet),' MNm/rad'])
text(0,2/7,['vertical bearing capacity factors: N_q = ',num2str(Nq),'; N_c = ',num2str(Nc)])
text(0,1/7,['base soil cohesion: c_b = ',num2str(cb),' kPa'])
text(0,0/7,['base elastic parameters: E_b = ',num2str(Eb),' MPa; n_b = ',num2str(nub)])

t=subplot(3,3,2);
hold on
box on
t=plot([0:1:max(displ(:,1))],d1(:,1),'-v',[0:1:max(displ(:,1))],d1(:,2),'->');
legend('w_0','u_0',2)
xlabel('Nstep')
ylabel('[m]')
title('imposed displacement history')
plot(0,0)

t=subplot(3,3,3);
set(t,'Ydir','reverse')
hold on
box on
t=plot(gter,z,'-b');
set(t,'marker','.')
xlabel('g_t_e_r [kN/m^3]')
ylabel('z [m]')
title('soil unit weight')
plot(0,0)

t=subplot(3,3,4);
set(t,'Ydir','reverse')
hold on
box on
t=plot(fs,z,'-b');
set(t,'marker','.')
xlabel('f_s [kPa]')
ylabel('z [m]')
title('unit shaft resistance')
plot(0,0)

t=subplot(3,3,5);
set(t,'Ydir','reverse')
hold on
box on
t=plot(Eter,z,'-',Gter,z,'-.');
set(t,'marker','.')
legend('E','G')
xlabel('E [MPa]')
ylabel('z [m]')
title('soil stiffness')
plot(0,0)

t=subplot(3,3,6);
hold on
box on
xlabel('y [m]')
ylabel('P [kN/m]')
title(['P-y curves, ',antype,' analysis'])
plot(0,0)
% assignment of P-y curves from file and plot
y=zeros((size(Pyread,2)-1)/2,size(Pyread,1));
P=zeros((size(Pyread,2)-1)/2,size(Pyread,1));
for i=1:size(Pyread,1)
    y(:,i)=Pyread(i,2:(size(Pyread,2)-1)/2+1)';
    P(:,i)=Pyread(i,(size(Pyread,2)-1)/2+2:end)';
    plot(y(:,i),P(:,i),'.-b');
    text(y(end,i),P(end,i),['  z=',num2str(Pyread(i,1)),'m']);
end
% addition of a last point at y=2*Dp
y=[y; 2*Dp*ones(1,size(Pyread,1))];
P=[P;P(end,:)];
XX=Pyread;
Pyread=[];
Pyread=XX(:,1);
Pyread=[Pyread y' P'];


% number of beam elements
nEtot = N;
% number of macroelements
nNtot = N+1;
% number of steps
nStep = max(displ(:,1)+1);

% definition of geometric properties of pile section
Rest=Dp/2;                               %Raggio esterno
Rint=Rest-tp;                            %Raggio interno
AreaSez=pi/4*(Dp^2-(Dp-2*tp)^2);         %Area della sezione trasversale
Jp=pi/64*(Dp^4-(Dp-2*tp)^4);             %Momento d'inerzia della sezione di palo

% definition of pile mechanical properties
EJ=Ep*Jp;                                %E*momento d'inerzia flessionale
EA=Ep*AreaSez;                           %E*AREA
EAI=zeros([nEtot,2]);
for ne=1:nEtot
    if antype=='NL'
        ze=Lp/N*(ne-1/2);
        EJ=interp1(Data_pile(:,1),Data_pile(:,2),ze);
        My(ne,1)=interp1(Data_pile(:,1),Data_pile(:,3),ze);
    end
    EAI(ne,:)=[EA,EJ];
end

t=subplot(3,3,9);
hold on
box on
title(['Bending stiffness, ',antype,' analysis'])
if antype=='LE'
    axis off
    text(0,6/7,['thickness of the pile: t_p = ',num2str(tp),' m'])
    text(0,5/7,['equivalent elastic modulus: E_p = ',num2str(Ep/1000),' GPa'])
    text(0,4/7,['sectional moment of inertia: J_p = ',num2str(Jp),' m^4'])
    text(0,3/7,['bending stiffness: EJ = ',num2str(EJ),' MNm'])
else
    xlabel('teta [1/m]')
    ylabel('M [MNm]')
    title(['M-teta curves, ',antype,' analysis'])
    plot(0,0)
    for ne=1:nEtot
       ze=Lp/N*(ne-1/2);
       plot([0,1,2]*My(ne,1)/EAI(ne,2),[0 My(ne,1) My(ne,1)],'--r')
       text(2*My(ne,1)/EAI(ne,2),My(ne,1),['  z_e = ',num2str(ze),' m'])
    end
end

%Calcolo lunghezze travi
dL=diff(z);

%Calcolo lunghezze macroelementi
dLequi=zeros([nNtot,1]);
for ii=2:nNtot-1
    dLequi(ii)=mean([z(ii+1) z(ii)])-mean([z(ii) z(ii-1)]);
end
dLequi(1)=dL(1)/2;
dLequi(nNtot)=dL(nEtot)/2;

%Costruzione matrice di rigidezza per struttura non assemblata
KK=[];
%Costruzione matrice di rigidezza per struttura assemblata
K=zeros([3*nNtot,3*nNtot]);
dKe=zeros([6,6,nEtot]);

%   dKe(:,:,ne)=[12*dEAI(ne,2)/dL(ne)^3,   6*dEAI(ne,2)/dL(ne)^2, -12*dEAI(ne,2)/dL(ne)^3, 6*dEAI(ne,2)/dL(ne)^2;
%                 6*dEAI(ne,2)/dL(ne)^2,   4*dEAI(ne,2)/dL(ne),    -6*dEAI(ne,2)/dL(ne)^2, 2*dEAI(ne,2)/dL(ne);
%               -12*dEAI(ne,2)/dL(ne)^3,  -6*dEAI(ne,2)/dL(ne)^2,  12*dEAI(ne,2)/dL(ne)^3,-6*dEAI(ne,2)/dL(ne)^2;
%                 6*dEAI(ne,2)/dL(ne)^2,   2*dEAI(ne,2)/dL(ne),    -6*dEAI(ne,2)/dL(ne)^2, 4*dEAI(ne,2)/dL(ne)];



for ne=1:nEtot
    dKe(:,:,ne)=[      EAI(ne,1)/dL(ne),                     0,                     0,     -EAI(ne,1)/dL(ne),                     0,                     0;
                                      0, 12*EAI(ne,2)/dL(ne)^3,  6*EAI(ne,2)/dL(ne)^2,                     0,-12*EAI(ne,2)/dL(ne)^3,  6*EAI(ne,2)/dL(ne)^2;
                                      0,  6*EAI(ne,2)/dL(ne)^2,    4*EAI(ne,2)/dL(ne),                     0, -6*EAI(ne,2)/dL(ne)^2,    2*EAI(ne,2)/dL(ne);
                      -EAI(ne,1)/dL(ne),                     0,                     0,      EAI(ne,1)/dL(ne),                     0,                     0;
                                      0,-12*EAI(ne,2)/dL(ne)^3, -6*EAI(ne,2)/dL(ne)^2,                     0, 12*EAI(ne,2)/dL(ne)^3, -6*EAI(ne,2)/dL(ne)^2;
                                      0,  6*EAI(ne,2)/dL(ne)^2,    2*EAI(ne,2)/dL(ne),                     0,-6*EAI(ne,2)/dL(ne)^2,     4*EAI(ne,2)/dL(ne);
                 ];

  % matrice assemblata
  nv=[3*ne-2:3*ne+3];
  K(nv,nv)=K(nv,nv)+dKe(:,:,ne);
  % matrice non assemblata
  nv=[6*(ne-1)+1:6*(ne-1)+6];
  KK(nv,nv)=dKe(:,:,ne);
end




% ***************************************************
% Definizione curve Q-w per ogni nodo (in verticale)
t=subplot(3,3,7);
hold on
box on
xlabel('w [m]')
ylabel('Q [kN/m]')
title(['Q-w curves, ',antype,' analysis'])
plot(0,0)
ww=[0:Dp/1000:2*Dp];
wbar=2*max(pi*Dp*fs./(pi*Gter/2*1e3));
for nn=1:nNtot
    kvv=pi*Gter(nn)/2*1e3; %[MPa-->kPa]
    Qss=pi*Dp*fs(nn); % [kN/m]
    Qwint(nn,:)=[z(nn) ww min(ww*kvv,Qss)];
    plot(ww,min(ww.*kvv,Qss),'--r')
    t=text(4*wbar/6,Qss,['z = ',num2str(z(nn)),' m']);
    if dim==[1 41 1600 784]
        set(t,'BackgroundColor',[1 1 1])
    end
end
v=axis;
v(2)=wbar;
axis(v)
clear v kvv Qss
% ***************************************************
% Definizione curve Q-w alla base del palo (in verticale)
kvv=4*Dp*Gb/(2*(1-nub))*1e3; % [conversione MPa-->kPa]
sv0=trapz(z,gter)-gw*(Lp-zw);
Qbs=(Nq*sv0 + Nc*cb)*pi/Dp^2/4; % [kN]
t=subplot(3,3,8);
hold on
box on
xlabel('w [m]')
ylabel('Q [kN]')
title(['Q-w curve at base, ',antype,' analysis'])
plot(ww,min(ww*kvv,Qbs),'--k')
% aggiungere la resistenza di base
Qwint(end,:)=Qwint(end,:)+[0 ww min(ww*kvv/dLequi(end),Qbs/dLequi(end))];
clear wbar kvv
% ***************************************************
% Definizione curve P-y per ogni nodo (in orizzontale)
yy=[0:Dp/1000:2*Dp];
xx=Pyread;
Pyread=[];
for i=1:size(xx,1);
    Pyread(i,:)=[xx(i,1) yy interp1(xx(i,2:(size(xx,2)-1)/2+1),xx(i,(size(xx,2)-1)/2+2:end),yy)];
end
ybar=2*max(xx(end,2:(size(xx,2)-1)/2));
clear xx

Pyint=z;
for i= 2:size(Pyread,2)
    AA=[];
    AA=interp1(Pyread(:,1),Pyread(:,i),z,'linear','extrap');
    AA=(AA+abs(AA))/2;
    Pyint=[Pyint AA];
    clear AA
end
% plot delle curve P-y estrapolate
t=subplot(3,3,6);
for i=1:size(Pyint,1)
    plot(Pyint(i,2:(size(Pyint,2)-1)/2),Pyint(i,(size(Pyint,2)-1)/2+2:end-1),'--r');
end
v=axis;
v(2)=ybar;
axis(v)
clear v
clear ybar
% ***************************************************

% ***************************************************
% definitions of intital tangent stiffness matrix for
% the soil reactions multiplied dy dLequi and converted
% [kN/m --> MN/m]
Ks=[];
for nn=1:nNtot
    % vertical stiffness
    QQ=Qwint(nn,(size(Qwint,2)-1)/2+2:end)*dLequi(nn)*1e-3;
    AA=diff(QQ)./diff(ww);
    Ks(3*(nn-1)+1,3*(nn-1)+1)=AA(1);
    clear AA QQ
    % horizontal stiffness
    PP=Pyint(nn,(size(Pyint,2)-1)/2+2:end)*dLequi(nn)*1e-3;
    AA=diff(PP)./diff(yy);
    Ks(3*(nn-1)+2,3*(nn-1)+2)=AA(1);
    clear AA PP
    % rotatinal stiffness
    Ks(3*(nn-1)+3,3*(nn-1)+3)=0;
end

% addition of the rotational stiffness of the foundation
Ks(3,3)=Ks(3,3)+Ktet;
% ***************************************************

% ***************************************************
% matrice con i valori delle curve Q-w e P-y ordinate per nodi
    for nn=1:nNtot        
        FF(3*(nn-1)+1,:)=Qwint(nn,(size(Qwint,2)-1)/2+2:end)*dLequi(nn)*1e-3;
        FF(3*(nn-1)+2,:)=Pyint(nn,(size(Pyint,2)-1)/2+2:end)*dLequi(nn)*1e-3;
        FF(3*(nn-1)+3,:)=0;        
    end    
% ***************************************************

% **************************
% calcolo delle sottomatrici
% **************************
% structure
K11=K(1:2,1:2);
K21=K(1:2,3:end);
K12=K(3:end,1:2);
K22=K(3:end,3:end);
% soil
Ks11=Ks(1:2,1:2);
Ks22=Ks(3:end,3:end);
% **************************