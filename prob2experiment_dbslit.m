more off;

format long;

% experiment 3 : double slit (vtype = 2)
idtype = 1;
vtype  = 2;
tmax   = 0.03;
lambda = 0.01;
level  = 9;
vpar   = [0.33, 0.37, 0.63, 0.67, 4000];
idpar  = [0.5, 0.1, 0.9, 0.05, 0, 2]; % boosted gaussian

[x, y, t, psi, psire, psiim, psimod, v] = sch_2d_adi(tmax, level, ...
    lambda, idtype, idpar, vtype, vpar);
plotit = 1;

[nt, nx, ny] = size(psi);
% plot(v, 'Color', 'r');
% input('press enter for video')

% animation  -------------------------------------------
plotenable = 1;
pausesecs = 0.001;
avienable = 1;

if ~plotenable
   avienable = 0;
end

avifilename = 'sch_dbslit.avi';
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
   %     contourf(x, y, v, 'LineColor', 'none');
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
