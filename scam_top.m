%% ######## Analysis of a Differential Amplifier##########
clc;clear all; close all;
warning off MATLAB:concatenation:integerInteraction

fid=fopen('diff_amp.cir');
Level=3;
[V,X,SymString] = scam_v5(fid,Level);
eval(SymString);
for i=1:length(V),
    eval([char(X(i)) '=' char(V(i)) ';']);
end

% %calculation of Ad
Vdd=0;
Vbias=0;
Vg1=0.5;
Vg2=-0.5;

% gm_MN2=gm_MN1;
% gm_MP4=gm_MP3;

Av=eval(v_3/(v_4-v_5))

% calculation of Zo
% Add this in line 3 of netlist 'diff_amp.cir'
% Vx  3 0 1
% Vdd=0;
% Vb=0;
% Vin=0;
% Zo=eval(Vx/I_Vx)
%%######## Analysis of a Differential Amplifier##########
% reset(symengine);
% %%
% clc;clear all; close all;
% warning off MATLAB:concatenation:integerInteraction
% 
% %######## Analysis of a Current Source Inverter##########
% fid=fopen('current_src_inv.cir');
% Level=1;
% [V,X,SymString] = scam_v5(fid,Level);
% eval(SymString);
% for i=1:length(V),
%     eval([char(X(i)) '=' char(V(i)) ';']);
% end
% 
% % %calculation of Av
% % Vdd=0;
% % Vb=0;
% % Vin=1;
% % Av=eval(v_3/v_4)
% 
% % % calculation of Zo
% % % Add this in line 3 of netlist 'current_src_inv.cir'
% % % Vx  3 0 1
% Vdd=0;
% Vb=0;
% Vin=0;
% Zo=eval(Vx/I_Vx)
% % ######## Analysis of a Current Source Inverter##########
