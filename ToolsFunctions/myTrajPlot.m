function myTrajPlot(X_opt, xSimu ,U , x0 , xf , t, Drivers, n,animOpt )

%----------------- PCA -----------------------------------
%xOpenLoop= reshape(xOpenLoopbis , [n , nSamps])';
%xSelfDynamic= reshape(xSelfDynamicBis , [n , nSamps])';

[coeff,score,latent] = pca(X_opt);
pcOpt =  X_opt * coeff;
pcSimu = xSimu * coeff;
%pcSelf = xSelfDynamic * coeff;
%pcOpenLoop = xOpenLoop * coeff;
pcX0 = x0' * coeff;
pcXf = xf' * coeff;
%----------------- Plots -----------------------------------
%U = reshape(U_opt , [n , nSamps])';
figure 
plot(t,U(:,Drivers==1), 'LineWidth' , 3)
xlabel('t (s)' , 'FontWeight' , 'bold')
ylabel('U(t)' , 'FontWeight' , 'bold')
yline(0)
labels = {};
kLabel =0;
for k=1:n
    if Drivers(k)
         kLabel= kLabel+1;
        labels{kLabel} = sprintf('U %i',k);
    end
end
legend(labels, 'FontWeight' , 'bold')
title('Input signal')

% %plot variance of PCA
% figure
% plot(latent)
% title('Variance of Principal Components')
% ----------------- plot actual trajectories ---------------
if animOpt
    main_fixU(A' , xSimu', U' , 1 , length(t) );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
h2 = myPlot_dir(X_opt(:,1),X_opt(:,2) , 'd');
hold on;
% [h5, h6] = plot_dir(pcSimu(:,1),pcSimu(:,2) , 'g');
% hold on;
% [h3, h4] = plot_dir(xOpenLoop(:,1) ,xOpenLoop(:,2) , 'b');
% hold on;
h6 = myPlot_dir(xSimu(:,1),xSimu(:,2) , 'o');
% initial and target pointsxOpenLoop
hold on;
point1 = plot(x0(1) , x0(2) , 'Marker', 'x' , 'Color','k' ,'MarkerSize' ,10);
hold on;
point2 = plot(xf(1) ,xf(2) ,'Marker' , 'd' ,'Color', 'k', 'MarkerSize' , 10);
%plot target circle
[xTarget,yTarget] = genTargetCircle(xf,norm(xf)/100);
hold on;
plot(xTarget , yTarget , '.')

xlabel('X1' , 'FontWeight','bold')
ylabel('X2' , 'FontWeight' , 'bold')
legend([h2 ,h6 , point1 , point2] , 'optimal traj' , ...
    'simulated traj', 'start' , 'end', 'FontWeight' , 'bold')
title('optimal and simulated trajectories')
% -------------------- PCA plot --------------------------
figure
h1 = myPlot_dir(pcOpt(:,1),pcOpt(:,2) , 'd');
hold on;
% [h5, h6] = plot_dir(pcSimu(:,1),pcSimu(:,2) , 'g');
% hold on;
% [h3, h4] = plot_dir(pcOpenLoop(:,1) ,pcOpenLoop(:,2) , 'b');
% hold on;
h5 = myPlot_dir(pcSimu(:,1),pcSimu(:,2) , 'o');
% initial and target points;
hold on;
point1 = plot(pcX0(1),pcX0(2),'Marker','x' , 'Color' , 'k', 'MarkerSize',10);
hold on;
point2 = plot(pcXf(1),pcXf(2),'Marker','d' ,'Color', 'k', 'MarkerSize' , 10);

[xTarget,yTarget] = genTargetCircle(pcXf,norm(pcXf)/100);
hold on;
plot(xTarget , yTarget , '.')

xlabel('Z1' , 'FontWeight','bold')
ylabel('Z2' , 'FontWeight' , 'bold')
legend([h1, h5 , point1 , point2] , 'optimal traj' , ...
    'simulated traj', 'start' , 'end', 'FontWeight' , 'bold')
title('optimal and simulated PC trajectories')

% figure
% [X,Y,Z] = adjacency_plot_und(A,position) ; 
% plot3(X,Y,Z);
end