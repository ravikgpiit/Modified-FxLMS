
clear all;close all;clc;

% variation of P,R and V check

% [v_rms,x] = AllCombined_err(0);
% v(1) = v_rms(250);
% % visualization
% rms_style  = '-g'; rms_legend = ' \beta = 0'; plot_rms

% disp('Press a key !')  % Press a key here.
% pause;
% 
% [v_rms,x] = AllCombined_err(0.25);
% v(2) = v_rms(250);
% % visualization
% rms_style  = '-r'; rms_legend =  '\beta = 0.25'; plot_rms
% 
% disp('Press a key !')  % Press a key here.
% pause;
% 
% [v_rms,x] = AllCombined_err(0.5);
% v(3) = v_rms(250);
% % visualization
% rms_style  = '-b'; rms_legend =  '\beta = 0.5'; plot_rms
% 
% disp('Press a key !')  % Press a key here.
% pause;
% 
% [v_rms,x] = AllCombined_err(1);
% v(4) = v_rms(250);
% % visualization
% rms_style  = '-k'; rms_legend = '  \beta = 1'; plot_rms
% 
% disp('Press a key !')  % Press a key here.
% pause;
% 
% 
[v_rms,x] = AllCombined_err(2);
v(5) = v_rms(250);
% visualization
rms_style  = '-b'; rms_legend = '\beta = 2'; plot_rms
% 
% disp('Press a key !')  % Press a key here.
% pause;


