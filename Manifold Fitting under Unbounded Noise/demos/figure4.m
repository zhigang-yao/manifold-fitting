% xiayq @ 12/10/2022
% xiayq0121@zufe.edu.cn
% plot figure4 in
% Z. Yao and Y. Xia, Manifold Fitting under Unbounded Noise, arXiv:1909.10228

%% circle
rep = 9;
tau = 1;
n = 300;
sigma = 0.06;

fname =  sprintf('out/circle/Dist_n%d_s%.2f.mat',n, sigma);
D1 = load(fname);

fname = sprintf('out/circle/Manapprox_Dist_n%d_s%.2f.mat',n, sigma);
D2 = load(fname);

% ya21
for i = 1:2
        Mout = reshape(D2.Mouts(i,rep,:,:),2,n);
        proj_Mout = reshape(D2.proj_Mouts(i,rep,:,:),2,n);
        
        figure;
        
        plot(proj_Mout(1,:),proj_Mout(2,:),'r.','MarkerSize',6); hold on;
        plot(Mout(1,:), Mout(2,:),'k.','MarkerSize',6);
        
        a = gca;
        a.YAxis.FontSize = 14;
        a.XAxis.FontSize = 14;
        
        xlim([-1.5,1.5]);
        ylim([-1.5,1.5]);
        
        
        axis off;
        
        algo = sprintf('ya21d%d',i);
        sname = sprintf('figures/figure4/circle_n%d_%s_s%.2f.fig', n, algo, sigma);
        saveas(gcf,sname)
         
end

% ours, cf18, km17
num_algo = numel(D1.algos);
for i = 1 : num_algo
    
        moveflag = D1.infos{i,rep}.moveflag;
        Mout = D1.Mouts{i,rep};
        algo = D1.algos{i};
        dist = D1.Dist2_move{i,rep};
        proj_Mout = bsxfun(@times, Mout, tau./sqrt(sum(Mout.^2)));

        figure;
        plot(proj_Mout(1,:),proj_Mout(2,:),'r.','MarkerSize',6); hold on;
        plot(Mout(1,moveflag), Mout(2,moveflag),'k.','MarkerSize',6);
       
        a = gca;
        a.YAxis.FontSize = 14;
        a.XAxis.FontSize = 14;
        
        xlim([-1.5,1.5]);
        ylim([-1.5,1.5]);
        
        
        axis off;
        
        
        sname = sprintf('figures/figure4/circle_n%d_%s_s%.2f.fig', n, algo, sigma);
        saveas(gcf,sname)
        
end


%% sphere
rep = 19;
tau = 1;
n = 1000;
sigma = 0.06;

fname =  sprintf('out/sphere/Dist_n%d_s%.2f.mat',n, sigma);
D1 = load(fname);

fname = sprintf('out/sphere/Manapprox_Dist_n%d_s%.2f.mat',n, sigma);
D2 = load(fname);

% ya21
for i = 1:2
        Mout = reshape(D2.Mouts(i,rep,:,:),3,n);
        proj_Mout = reshape(D2.proj_Mouts(i,rep,:,:),3,n);
        
        figure;

        plot3(proj_Mout(1,:),proj_Mout(2,:), proj_Mout(3,:),'r.','MarkerSize',6); hold on;
        plot3(Mout(1,:), Mout(2,:), Mout(3,:),'k.','MarkerSize',6);
    
        xlim([-1.5,1.5]);
        ylim([-1.5,1.5]);
        zlim([-1.5,1.5]);
    
        view(3);
    
        axis off
        
        algo = sprintf('ya21d%d',i);
        sname = sprintf('figures/figure4/sphere_n%d_%s_s%.2f.fig', n, algo, sigma);
        saveas(gcf,sname)
         
         
end

% ours, cf18, km17
num_algo = numel(D1.algos);
for i = 1 : num_algo
    
        Mout = D1.Mouts{i,rep};
        proj_Mout = D1.proj_Mouts{i,rep};
        algo = D1.algos{i};

        figure;

        plot3(proj_Mout(1,:),proj_Mout(2,:), proj_Mout(3,:),'r.','MarkerSize',6); hold on;
        plot3(Mout(1,:), Mout(2,:), Mout(3,:),'k.','MarkerSize',6);
    
        xlim([-1.5,1.5]);
        ylim([-1.5,1.5]);
        zlim([-1.5,1.5]);
    
        view(3);
    
        axis off
        
        
        
        sname = sprintf('figures/figure4/sphere_n%d_%s_s%.2f.fig', n, algo, sigma);
        saveas(gcf,sname)
        
end

