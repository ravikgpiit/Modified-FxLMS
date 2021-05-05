clear all;close all;clc;


[v_rms,x] = AllCombined_0(0);
v(1) = v_rms(250);
% visualization
rms_style  = '-g'; rms_legend = ' \beta = 0'; plot_rms


[v_rms,x] = AllCombined_0(0.5);
v(2) = v_rms(250);
% visualization
rms_style  = '-r'; rms_legend = ' \beta = 0.5'; plot_rms



[v_rms,x] = AllCombined_0(1);
v(3) = v_rms(250);
% visualization
rms_style  = '-b'; rms_legend = ' \beta = 1'; plot_rms


[v_rms,x] = AllCombined_0(1.5);
v(4) = v_rms(250);
% visualization
rms_style  = '-k'; rms_legend = ' \beta = 1.5'; plot_rms



[v_rms,x] = AllCombined_0(2);
v(5) = v_rms(250);
% visualization
rms_style  = '-m'; rms_legend = ' \beta = 2'; plot_rms



[v_rms,x] = AllCombined_0(2.5);
v(6) = v_rms(250);
% visualization
rms_style  = '-c'; rms_legend = ' \beta = 2.5'; plot_rms


[v_rms,x] = AllCombined_0(3);
v(7) = v_rms(250);
% visualization
rms_style  = '-b'; rms_legend = ' \beta = 3'; plot_rms


figure(111)
betaa = (0:0.5:3);
% betaa=[0 0.25 0.5 0.75 1 1.25 1.5 1.75 2.5 3.5 5];
betaaa = (0:.025:3);
vv = pchip(betaa,v,betaaa);
semilogy(betaa,v,'o',betaaa,vv)
xlabel('\beta'); ylabel('v_rms at x = 500');






