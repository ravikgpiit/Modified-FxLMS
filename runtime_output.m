
if abs(t(i)-round(t(i)/(Tend/10))*(Tend/10)) < dt/2
    fprintf('t(i) = %6.0f',t(i))
    if i == 1
        wall_time = 0;
        fprintf('\n')
    else
        wall_time = wall_time + toc;
        fprintf(' (ETF: %2.0f min)\n',ceil(wall_time/(i-1)*(nt-i-1)/60))
    end
    tic
end