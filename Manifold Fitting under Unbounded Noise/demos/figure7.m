% xiayq @ 12/11/2022
% xiayq0121@zufe.edu.cn
% plot figure7 in
% Z. Yao and Y. Xia, Manifold Fitting under Unbounded Noise, arXiv:1909.10228

% n = 500 or 800
rep = 8; % 15
a = 2/3; b = 1/3; tau = b;
n = 800; % 500
sigma = 0.04;

fname =  sprintf('out/torus/Dist_n%d_s%.2f.mat',n, sigma);
D1 = load(fname);

fname = sprintf('out/torus/Manapprox_Dist_n%d_s%.2f.mat',n, sigma);
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
    
        view(2);
    
        axis off
        
        algo = sprintf('ya21d%d',i);
        sname = sprintf('figures/figure5/torus_n%d_%s_s%.2f.fig', n, algo, sigma);
        saveas(gcf,sname)
           sname = sprintf('figures/figure5/torus_%s.eps', algo);
        saveas(gcf,sname,'psc')
         
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
    
        view(2);
    
        axis off
        
        
        
        sname = sprintf('figures/figure5/torus_n%d_%s_s%.2f.fig', n, algo, sigma);
        saveas(gcf,sname)
           sname = sprintf('figures/figure5/torus_%s.eps', algo);
        saveas(gcf,sname,'psc')
        
end


