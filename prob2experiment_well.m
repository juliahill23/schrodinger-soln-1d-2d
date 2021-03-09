more off;
format long;

% experiment 2 : rectangular well (vtype = 1)
idtype = 1;
vtype  = 1;
tmax   = 0.05;
lambda = 0.01;
level  = 9;
vpar   = [0.1, 0.3, 0, 1, -1000];
idpar  = [0.8, 0.5, 0.1, 0.1, 0.0, 0.0]; % boosted gaussian

% experiment 2: rectangular well (vtype = 1)
[x, y, t, psi, psire, psiim, psimod, v] = sch_2d_adi(tmax, level, ...
    lambda, idtype, idpar, vtype, vpar);
plotit = 1;

[nt, nx, ny] = size(psi);

% animation  -------------------------------------------
plotenable = 1;
pausesecs = 0.001;
avienable = 1;

if ~plotenable
   avienable = 0;
end

avifilename = 'sch_well.avi';
aviframerate = 100;

if avienable
   aviobj = VideoWriter(avifilename);
   open(aviobj);
end

if plotit
    for it = 1 : nt
        figure(level);
        psitemp = zeros(nx, ny);
        psitemp(:,:) = psimod(it, :,:);
        contourf(x,y,psitemp, 'LineColor', 'none');
    %    contourf(x, y, v, 'LineColor', 'none');
        colormap(summer)
        colorbar
        xlabel('x');
        ylabel('y') ;
        xlim([0, 1]);
        ylim([0, 1]);
        drawnow;
        if avienable
            if it == 1
                framecount = 5 * aviframerate ;
            else
                framecount = 1;
            end
            for iframe = 1 : framecount
                writeVideo(aviobj, getframe(gcf));
            end
        end
        
        pause(pausesecs);
    end
end

if avienable
   close(aviobj);
   fprintf('Created video file: %s\n', avifilename);
end
