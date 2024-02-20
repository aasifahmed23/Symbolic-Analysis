% function [Ids,Vth,Qd, Qg, Qs, Qbulk,Qd_overlap, Qg_overlap, Qs_overlap,Cox] = N_18_MM_V2_C_symbolic(Vd_val,Vg_val,Vs_val,Vb_val,W_val,L_val)
%%
clear all;clc;
Vd_val=1.3;
Vg_val=0.8;
Vs_val=0;
Vb_val=0;
W_val=10e-6;
L_val=2e-6;
values={Vd_val,Vg_val,Vs_val,Vb_val,W_val,L_val};
syms Vd Vg Vs Vb real;
syms W L
variables={Vd,Vg,Vs,Vb,W,L};
Vgs=Vg-Vs;
Vds=Vd-Vs;
Vbs=Vb-Vs;
% Vds=1.8;
% W=1.8e-6;
% L=0.36E-6;
% Vgs=-1:0.01:3.8;
% Vbs=-0.5;%
%% NMOS model card
%% Refer pg# 161 of BSIM4.3.0 MOSFET Model here : https://cmosedu.com/cmos1/BSIM4_manual.pdf
TNOM    = 27 ;
TOX     = 4.1E-9;
XJ      = 1E-7; 
NCH     = 2.3549E17;
VTH0    = 0.3725327;
K1      = 0.5933684;
K2      = 2.050755E-3;
K3      = 1E-3;
K3B     = 4.5116437;
W0      = 1E-7;
NLX     = 1.870758E-7;
DVT0W   = 0;
DVT1W   = 0;
DVT2W   = 0;
DVT0    = 1.3621338;
DVT1    = 0.3845146;
DVT2    = 0.0577255;
U0      = 259.5304169;
UA      = -1.413292E-9;
UB      = 2.229959E-18;
UC      = 4.525942E-11;
VSAT    = 9.411671E4;
A0      = 1.7572867;
AGS     = 0.3740333;
B0      = -7.087476E-9;
B1      = -1E-7;
KETA    = -4.331915E-3; 
A1      = 0;
A2      = 1;
RDSW    = 111.886044;
PRWG    = 0.5;
PRWB    = -0.2;
WR      = 1;
WINT    = 0;    
LINT    = 1.701524E-8;
XL      = 0;
XW      = -1E-8;
DWG     = -1.365589E-8;
DWB     = 1.045599E-8;
VOFF    = -0.0927546;
NFACTOR = 2.4494296;
CIT     = 0; 
CDSC    = 2.4E-4;
CDSCD   = 0;
CDSCB   = 0;
ETA0    = 3.175457E-3;
ETAB    = 3.494694E-5;
DSUB    = 0.0175288;
PCLM    = 0.7273497;  
PDIBLC1 = 0.1886574;
PDIBLC2 = 2.617136E-3; 
PDIBLCB = -0.1; 
DROUT   = 0.7779462;
PSCBE1  = 3.488238E10;
PSCBE2  = 6.841553E-10;  
PVAG    = 0.0162206;
DELTA   = 0.01;
RSH     = 6.5; 
MOBMOD  = 1;
PRT     = 0;
UTE     = -1.5;
KT1     = -0.11;
KT1L    = 0;
KT2     = 0.022;
UA1     = 4.31E-9;
UB1     = -7.61E-18; 
UC1     = -5.6E-11;
AT      = 3.3E4;
WL      = 0;
WLN     = 1;  
WW      = 0;
WWN     = 1;
WWL     = 0;
LL      = 0;
LLN     = 1;
LW      = 0;
LWN     = 1;
LWL     = 0;
CAPMOD  = 2;
XPART   = 0.5;
CGDO    = 8.53E-10;
CGSO    = 8.53E-10;    
CGBO    = 1E-12;
CJ      = 9.513993E-4;   
PB      = 0.8;
MJ      = 0.3773625;
CJSW    = 2.600853E-10;
PBSW    = 0.8157101;
MJSW    = 0.1004233;
CJSWG   = 3.3E-10;
PBSWG   = 0.8157101;
MJSWG   = 0.1004233;
CF      = 0;     
PVTH0   = -8.863347E-4;   
PRDSW   = -3.6877287;
PK2     = 3.730349E-4; 
WKETA   = 6.284186E-3;
LKETA   = -0.0106193;
PU0     = 16.6114107;  
PUA     = 6.572846E-11; 
PUB     = 0;
PVSAT   = 1.112243E3;  
PETA0   = 1.002968E-4;    
PKETA   = -2.906037E-3;

% % %% PMOS Model Card
% % TNOM    = 27;
% % TOX     = 4.1E-9;
% % XJ      = 1E-7;
% % NCH     = 4.1589E17;
% % VTH0    = -0.3948389;
% % K1      = 0.5763529;
% % K2      = 0.0289236;
% % K3      = 0;
% % K3B     = 13.8420955;
% % W0      = 1E-6;
% % NLX     = 1.337719E-7;
% % DVT0W   = 0;
% % DVT1W   = 0;
% % DVT2W   = 0;
% % DVT0    = 0.5281977;
% % DVT1    = 0.2185978;
% % DVT2    = 0.1;
% % U0      = 109.9762536;
% % UA      = 1.325075E-9;
% % UB      = 1.577494E-21;
% % UC      = -1E-10;
% % VSAT    = 1.910164E5;
% % A0      = 1.7233027;
% % AGS     = 0.3631032;
% % B0      = 2.336565E-7;
% % B1      = 5.517259E-7;
% % KETA    = 0.0217218;
% % A1      = 0.3935816;
% % A2      = 0.401311;
% % RDSW    = 252.7123939;
% % PRWG    = 0.5;
% % PRWB    = 0.0158894;
% % WR      = 1;
% % WINT    = 0;
% % LINT    = 2.718137E-8;
% % XL      = 0;
% % XW      = -1E-8;
% % DWG     = -4.363993E-8;
% % DWB     = 8.876273E-10;
% % VOFF    = -0.0942201;
% % NFACTOR = 2;
% % CIT     = 0;
% % CDSC    = 2.4E-4;
% % CDSCD   = 0;
% % CDSCB   = 0;
% % ETA0    = 0.2091053;
% % ETAB    = -0.1097233;
% % DSUB    = 1.2513945;
% % PCLM    = 2.1999615;
% % PDIBLC1 = 1.238047E-3;
% % PDIBLC2 = 0.0402861;
% % PDIBLCB = -1E-3;
% % DROUT   = 0;
% % PSCBE1  = 1.034924E10;
% % PSCBE2  = 2.991339E-9;
% % PVAG    = 15;
% % DELTA   = 0.01;
% % RSH     = 7.5;
% % MOBMOD  = 1;
% % PRT     = 0;
% % UTE     = -1.5;
% % KT1     = -0.11;
% % KT1L    = 0;
% % KT2     = 0.022;
% % UA1     = 4.31E-9;
% % UB1     = -7.61E-18;
% % UC1     = -5.6E-11;
% % AT      = 3.3E4;
% % WL      = 0;
% % WLN     = 1;
% % WW      = 0;
% % WWN     = 1;
% % WWL     = 0;
% % LL      = 0;
% % LLN     = 1;
% % LW      = 0;
% % LWN     = 1;
% % LWL     = 0;
% % CAPMOD  = 2;
% % XPART   = 0.5;
% % CGDO    = 6.28E-10;
% % CGSO    = 6.28E-10;
% % CGBO    = 1E-12;
% % CJ      = 1.160855E-3;
% % PB      = 0.8484374;
% % MJ      = 0.4079216;
% % CJSW    = 2.306564E-10;
% % PBSW    = 0.842712;
% % MJSW    = 0.3673317;
% % CJSWG   = 4.22E-10;
% % PBSWG   = 0.842712;
% % MJSWG   = 0.3673317;
% % CF      = 0;
% % PVTH0   = 2.619929E-3;
% % PRDSW   = 1.0634509;
% % PK2     = 1.940657E-3;
% % WKETA   = 0.0355444;
% % LKETA   = -3.037019E-3;
% % PU0     = -1.0227548;
% % PUA     = -4.36707E-11;
% % PUB     = 1E-21;
% % PVSAT   = -50;
% % PETA0   = 1E-4;
% % PKETA   = -5.167295E-3;
% % %% Bias voltages
% % Vgs=0.55;
% % Vds=0.9;
% % Vbs=0.4;
% % %% Device sizes
% % W=(2e-6);
% % L=(0.72E-6);
Wdrawn=W;
Ldrawn=L;
%% Physical constants
kb=1.3806503e-23;
q=1.602176462e-19;
ephsilon_si=4.1;
ephsilon_knot=8.854187817e-12;
Vtmo=kb*TNOM/q;
Eg0=1.16-((7.02E-4)*TNOM^2/(TNOM+1108));
ni=(1.45e10)*((TNOM/300.15)^1.5)*exp(21.5565981-Eg0/(2*Vtmo));
Phi_s=2*Vtmo*log(NCH/ni);
%% Calc Threshold Voltage
T=TNOM;
Vt=Vtmo*T/TNOM;
Nds=1e20;
Vbi=Vt*log(NCH*Nds/ni^2);
Vbc=0.9*(Phi_s-(K1^2)/(4*K2^2));
delta1=0.001;
Vbs_eff=Vbc+0.5*(Vbs-Vbc-delta1+sqrt((Vbs-Vbc-delta1).^2 -4*delta1*Vbc)); %%%MATRIX
Xdep0=sqrt(2*ephsilon_si*Phi_s/(q*NCH));
Xdep=sqrt(2*ephsilon_si*(Phi_s-Vbs_eff)/(q*NCH));  %%%MATRIX
Cox=(ephsilon_si*ephsilon_knot)/TOX;
l_to=sqrt(ephsilon_si*Xdep0/Cox);
l_t_o=sqrt(ephsilon_si*Xdep/Cox); %%%MATRIX
l_tw=l_t_o.*(1+DVT2W.*Vbs_eff); %%%MATRIX
l_t=l_t_o.*(1+DVT2.*Vbs_eff); %%%MATRIX
TOXM=TOX;
K1_ox=K1*TOX/TOXM;
K2_ox=K2*TOX/TOXM;
Vth0_ox=VTH0-K1*sqrt(Phi_s);
%bifurcate Vth calc for Effective channel Length and Width (dW dL)  calc
dL=LINT+LL./L.^LLN+LW./W.^LWN+LWL./((L.^LLN).*(W.^LWN)); %%%MATRIX
dW_dash=WINT+WL./L.^WLN+WW./W.^WWN+WWL./((L.^WLN).*(W.^WWN)); %%%MATRIX
Leff=Ldrawn-2*dL; %%%MATRIX
Weff_dash=Wdrawn-2*dW_dash; %%%MATRIX
%bifurcate dL dW calc & Calc Vgst_eff
Cd=ephsilon_si./Xdep; %%%MATRIX
eta=1+NFACTOR*(Cd./Cox)+(CDSC+CDSCD.*Vds+CDSCB.*Vbs_eff).*(exp(-DVT1*Leff./(2*l_t))+2.*exp(-DVT1*Leff./l_t))./Cox+CIT./Cox; %%%MATRIX
%resume Vth calc
Vth=Vth0_ox+K1_ox*sqrt(Phi_s-Vbs_eff)-K2_ox*Vbs_eff; %%%MATRIX
Vth=Vth+K1_ox*(sqrt(1+NLX./Leff)-1)*sqrt(Phi_s)+(K3+K3B*Vbs_eff)*TOX*Phi_s./(Weff_dash+W0); %%%MATRIX
Vth=Vth-DVT0W*(exp(-DVT1W*Weff_dash.*Leff./(2*l_tw))+2*exp(-DVT1W*Weff_dash.*Leff./l_tw))*(Vbi-Phi_s);%%%MATRIX
Vth=Vth-DVT0*(exp(-DVT1*Leff./(2*l_t))+2*exp(-DVT1*Leff./l_t))*(Vbi-Phi_s);%%%MATRIX
Vth=Vth-(exp(-DSUB*Leff./(2.*l_to))+2*exp(-DSUB*Leff/l_to)).*(ETA0+ETAB*Vbs_eff).*Vds;%%%MATRIX
Vgst_eff=2*eta*Vt*log(1+exp((Vgs-Vth)/(2*eta*Vt)))./(1+2*eta*Cox*sqrt(2*Phi_s/(q*ephsilon_si*NCH))*exp(-(Vgs-Vth-2*VOFF)/(2*eta*Vt)));%%%MATRIX
%resume dL dW calc
dW=dW_dash+DWG*Vgst_eff+DWB*(sqrt(Phi_s-Vbs_eff)-sqrt(Phi_s));%%%MATRIX
% Leff=Ldrawn-2*dL;
Weff=Wdrawn-2*dW; %%%MATRIX
%% Mobility for Mobmod=1
mu_eff=U0./(1+(UA+UC*Vbs_eff).*((Vgst_eff+2*Vth)/TOX)+UB*((Vgst_eff+2*Vth)/TOX).^2); %%%MATRIX
%% Source/Drain Resistance Rds
Rds=RDSW*(1+PRWG*Vgst_eff+PRWB*(sqrt(Phi_s-Vbs_eff)-sqrt(Phi_s)))/(Weff_dash*1e6).^WR; %%%MATRIX
%% Drain Saturation Voltage Vdsat
Lambda=A1*Vgst_eff+A2; %%%MATRIX
Esat=2*VSAT./mu_eff; %%%MATRIX
% Find Abulk %%%MATRIX
Abulk=(1+(K1_ox./(2*sqrt(Phi_s-Vbs_eff))).*((A0*Leff./(Leff+2*sqrt(XJ*Xdep))).*(1-AGS*Vgst_eff.*(Leff./(Leff+2*sqrt(XJ*Xdep))).^2)+ B0./(B1+Weff_dash)))./(1+KETA*Vbs_eff);

if (subs(Rds,variables,values)>0) % does multiple substitution on symbolic expression of Rds(values,variables are cell arrays defined in line 10 and 12 resp.)
a=(Abulk.^2).*Weff*VSAT*Cox.*Rds+(1./Lambda -1).*Abulk;
b=-((Vgst_eff+2*Vt).*(2./Lambda -1)+Abulk.*Esat.*Leff+3*Abulk.*(Vgst_eff+2*Vt).*Weff.*VSAT.*Cox.*Rds);
c=(Vgst_eff+2*Vt).*Esat.*Leff+2*Weff.*VSAT*Cox.*Rds.*(Vgst_eff+2*Vt).^2;
Vdsat=(-b-sqrt(b.^2-4*a.*c))./(2*a);
else
   Vdsat=Esat*Leff*(Vgst_eff+2*Vt)/(Abulk*Esat*Leff+(Vgst_eff+2*Vt)); %%%MATRIX
end

%% Effective Vds
% DELTA=0.01;%effective Vds parameter
Vds_eff=Vdsat-0.5*(Vdsat-Vds-DELTA+sqrt((Vdsat-Vds-DELTA).^2+4*DELTA*Vdsat)); %%%MATRIX
%% Drain Current Expression
litl=sqrt(ephsilon_si*TOX*XJ/ephsilon_knot);
VAsat=(Esat.*Leff+Vdsat+2*Rds.*VSAT.*Cox.*Weff.*Vgst_eff.*(1-Abulk.*Vdsat*0.5./(Vgst_eff+2*Vt)))/(2./Lambda-1+Rds.*VSAT.*Cox.*Weff.*Abulk);
VASCBE=(Leff/PSCBE2).*exp(PSCBE1*litl./(Vds-Vds_eff));
THETA_ROUT=PDIBLC2+PDIBLC1*(exp(-DROUT*Leff/(2*l_to))+2*exp(-DROUT*Leff/l_to));
VADIBLC=((Vgst_eff+2*Vt)./(THETA_ROUT.*(1+PDIBLCB.*Vbs_eff))).*(1-Abulk.*Vdsat./(Abulk.*Vdsat+Vgst_eff+2*Vt));
VACLM=(Abulk.*Esat.*Leff+Vgst_eff).*(Vds-Vds_eff)./(PCLM*Abulk.*Esat*litl);
VA=VAsat+(1+PVAG*Vgst_eff./(Esat.*Leff))./(1./VACLM +1./VADIBLC);
Idso=Weff.*mu_eff.*Cox.*Vgst_eff.*Vds_eff.*(1-Abulk.*Vds_eff.*0.5./(Vgst_eff+2*Vt))./(Leff.*(1+Vds_eff./(Esat.*Leff)));
Ids=(Idso./(1+Rds.*Idso./Vds_eff)).*(1+(Vds-Vds_eff)./VA).*(1+(Vds-Vds_eff)./VASCBE);  
%% Substrate Current
alpha0=0; %first parameter of impact ionization current
alpha1=0.0;% Isub parameter for length scaling
beta0=30;% second parameter of impact ionization current
Isub=((alpha0+alpha1.*Leff).*(Vds-Vds_eff)./Leff).*(Idso./(1+Rds.*Idso./Vds_eff)).*(1+(Vds-Vds_eff)./VA).*exp(-beta0./(Vds-Vds_eff));
%% Temperature Effects
Vth_T=Vth+(KT1+KT1L./Leff+KT2*Vbs_eff).*(T/TNOM-1);
U0_T=U0*(T/TNOM)^UTE;
VSAT_T=VSAT-AT*(T/TNOM-1);
RDSW_T=RDSW+PRT*(T/TNOM-1);
UA_T=UA+UA1*(T/TNOM-1);
UB_T=UB+UB1*(T/TNOM-1);
UC_T=UC+UC1*(T/TNOM-1);

%% Intrinsic Capacitance Model (CAPMOD=2)

% Effective Channel Lenght & Width for Capacitance Model
DLC=LINT;
DWC=WINT;
LLC=LL;
LWC=LW;
LWLC=LWL;
WLC=WL;
WWC=WW;
WWLC=WWL;
delta_Lactive=DLC+LLC./L.^LLN+LWC./W.^LWN+LWLC./((L.^LLN).*(W.^LWN)); %%%MATRIX
delta_Wactive=DWC+WLC./L.^WLN+WWC./W.^WWN+WWLC./((L.^WLN).*(W.^WWN)); %%%MATRIX
Lactive=Ldrawn-2*delta_Lactive; 
Wactive=Wdrawn-2*delta_Wactive; 

% CAPMOD 2 COEFFS
Abulk0=(1+(K1_ox./(2*sqrt(Phi_s-Vbs_eff))).*(A0*Leff./(Leff+2*sqrt(XJ*Xdep))+B0./(B1+Weff_dash)))./(1+KETA*Vbs_eff);
CLC=1e-7;%const term for short channel model
CLE=0.6;%exponential term for short channel model
Abulk_dash=Abulk0.*(1+(CLC./Leff).^CLE);

noff=1.0;
voffcv=0;
Vgst_eff_cv=noff*eta*Vt*log(1+exp((Vgs-Vth-voffcv)/(noff*eta*Vt)));
Vdsat_cv=Vgst_eff_cv./Abulk_dash;

del4=0.02;
v4=Vdsat_cv-Vds-del4;
Vcv_eff=Vdsat_cv-0.5*(v4+sqrt(v4.^2+4*del4*Vdsat_cv));


%bias independent Vth is Vth_cv
Vth_cv=Vth0_ox+K1_ox*(sqrt(1+NLX./Leff)-1)*sqrt(Phi_s); %%%MATRIX
Vth_cv=Vth_cv-DVT0W*(exp(-DVT1W*Weff_dash.*Leff./(2*l_tw))+2*exp(-DVT1W*Weff_dash.*Leff./l_tw))*(Vbi-Phi_s);%%%MATRIX
Vth_cv=Vth_cv-DVT0*(exp(-DVT1*Leff./(2*l_t))+2*exp(-DVT1*Leff./l_t))*(Vbi-Phi_s);%%%MATRIX

vfb=Vth_cv-Phi_s-K1_ox*sqrt(Phi_s-Vbs_eff);
Vgb=Vgs-Vbs;
del3=0.02;
v3=vfb-Vgb-del3;
Vfb_eff=vfb-0.5*(v3+sqrt(v3.^2+4*del3*vfb));

%Bulk Terminal Charges
Qacc=-Wactive.*Lactive.*Cox.*(Vfb_eff-vfb);
Qsub0=-(0.5*Wactive.*Lactive.*Cox.*K1_ox^2).*(-1+sqrt(1+4*((Vgs-Vfb_eff-Vgst_eff_cv-Vbs_eff)/K1_ox^2)));
del_Qsub=Wactive.*Lactive.*Cox.*((1-Abulk_dash).*0.5.*Vcv_eff-((1-Abulk_dash).*Abulk_dash.*Vcv_eff.^2)./(12*Vgst_eff_cv-6*Abulk_dash.*Vcv_eff));
Qbulk=Qacc+Qsub0+del_Qsub;

%Src/Drain Terminal Charges
switch(XPART)
    case 0
        Qs=-(Wactive.*Lactive.*Cox).* ( Vgst_eff_cv.*0.5+ 0.25*Abulk_dash.*Vcv_eff -  ((Abulk_dash.*Vcv_eff).^2)/(24*(Vgst_eff_cv-0.5*Abulk.*Vcv_eff)) )  ;
        Qd=-(Wactive.*Lactive.*Cox).* ( Vgst_eff_cv.*0.5- 0.75*Abulk_dash.*Vcv_eff +  ((Abulk_dash.*Vcv_eff).^2)/(8*(Vgst_eff_cv-0.5*Abulk.*Vcv_eff)) )  ;
    case 0.5
        Qs=-(Wactive.*Lactive.*Cox*0.5).* ( Vgst_eff_cv-0.5*Abulk_dash.*Vcv_eff +     ((Abulk_dash.*Vcv_eff).^2)./(12*(Vgst_eff_cv-0.5*Abulk.*Vcv_eff)) )  ;
        Qd=Qs;
    case 1        
        Qs=-(Wactive.*Lactive.*Cox./(2*Vgst_eff_cv-Abulk_dash.*Vcv_eff)).*(Vgst_eff_cv.^3 -(4/3)*(Abulk_dash.*Vcv_eff).*Vgst_eff_cv.^2 +(2/3)*Vgst_eff_cv.*(Abulk_dash.*Vcv_eff).^2 -(2/15)*(Abulk_dash.*Vcv_eff).^3);
        Qd=-(Wactive.*Lactive.*Cox./(2*Vgst_eff_cv-Abulk_dash.*Vcv_eff)).*(Vgst_eff_cv.^3 -(5/3)*(Abulk_dash.*Vcv_eff).*Vgst_eff_cv.^2 +Vgst_eff_cv.*(Abulk_dash.*Vcv_eff).^2 -(1/5)*(Abulk_dash.*Vcv_eff).^3);
    otherwise
        Qs=-(Wactive.*Lactive*Cox*0.5)* ( Vgst_eff_cv-0.5*Abulk_dash*Vcv_eff +     ((Abulk_dash*Vcv_eff)^2)/(12*(Vgst_eff_cv-0.5*Abulk*Vcv_eff)) )  ;
        Qd=Qs;
end

% Gate terminal charges
Qinv=Qs+Qd;
Qg=-(Qinv+Qbulk);

%% Extrinsic Capacitance Model (CAPMOD=2)

%Source Overlap Charge
CGS1=0;%Light doped source-gate region overlap capacitance
CGD1=0;%Light doped drain-gate region overlap capacitance
CKAPPA=0.6;%Coefficient for lightly doped region overlap capacitance Fringing field capacitance
% NLDD=5e23;
% CKAPPA=2*ephsilon_si*q*NLDD/Cox^2;

del1=0.02;
v1=Vgs+del1;
Vgs_overlap=0.5*(v1-sqrt(v1.^2+4*del1));
Qs_overlap=Wactive.*(CGSO*Vgs+CGS1*(Vgs-Vgs_overlap-0.5*CKAPPA*(-1+sqrt(1-4*Vgs_overlap/CKAPPA))));

%Drain Overlap Charge
Vgd=Vgs-Vds;
v2=Vgd+del1;
Vgd_overlap=0.5*(v2-sqrt(v2.^2+4*del1));
Qd_overlap=Wactive.*(CGDO*Vgd+CGD1*(Vgd-Vgd_overlap-0.5*CKAPPA*(-1+sqrt(1-4*Vgd_overlap/CKAPPA))));

%Gate Overlap Charge
Qg_overlap=-(Qs_overlap+Qd_overlap+CGBO*Lactive*Vgb);

gm=diff(Ids,Vg);
% gm_val=subs(gm,variables,values);
gds=diff(Ids,Vd);
% gds_val=subs(gds,variables,values);
gmb=diff(Ids,Vb);
% gmb_val=subs(gmb,variables,values);
% 
Cgs=-diff(Qg,Vs);
% Cgs_val=subs(Cgs,variables,values);
Cgd=-diff(Qg,Vd);
% Cgd_val=subs(Cgd,variables,values);
Cgb=-diff(Qg,Vb);
% Cgb_val=subs(Cgb,variables,values);
% 
Csb=-diff(Qs,Vb);
% Csb_val=subs(Csb,variables,values);

Cdb=-diff(Qd,Vb);
% Cdb_val=subs(Cdb,variables,values);

% matlabFunction(gm,'file','test_gm.m');
matlabFunction(gds,'file','test_gds.m');
matlabFunction(gmb,'file','test_gmb.m');
% matlabFunction(Cgs,'file','test_Cgs.m');
% matlabFunction(Cgd,'file','test_Cgd.m');
% matlabFunction(Cgb,'file','test_Cgb.m');
% matlabFunction(Csb,'file','test_Csb.m');
% matlabFunction(Cdb,'file','test_Cdb.m');


