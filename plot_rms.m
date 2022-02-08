


hold on

semilogy(x,v_rms,rms_style,'Linewidth',2)

set(gca,'YScale','log')
axis([0 800 0.5*1e0 1e2/0.5])
xlabel('x'); ylabel('v_{rms}'); %grid on

legend_string = get(legend,'String');
legend_string{length(legend_string)+1} = rms_legend;
legend(legend_string,'Location','NW')

hold off