% xiayq @ 12/11/2022
% xiayq0121@zufe.edu.cn
% refered to Z. Yao and Y. Xia, Manifold Fitting under Unbounded Noise, arXiv:1909.10228

clear; %clc

% parameters for data
D = 3; 
NumTrials = 20;
tau = 1; sigma = 0.04;%0.04

% method setup
algos = {'ours','cf18','km17'};


switch sigma
    case 0.09
        r1 = 2*sqrt(sigma);  
        r2 = 2.5*sqrt(sigma); 
        r3 = 1.2*sqrt(sigma); 
    case 0.08
        r1 = 2*sqrt(sigma);  
        r2 = 2.8*sqrt(sigma); 
        r3 = 1.2*sqrt(sigma); 
    case 0.07
        r1 = 2*sqrt(sigma);  
        r2 = 2*sqrt(sigma); 
        r3 = 1.2*sqrt(sigma); 
    case 0.06
        r1 = 2*sqrt(sigma);  
        r2 = 2.5*sqrt(sigma); 
        r3 = 1.2*sqrt(sigma); 
    case 0.05
        r1 = 2*sqrt(sigma);  
        r2 = 2.5*sqrt(sigma); 
        r3 = 1.2*sqrt(sigma); 
    case 0.04
        r1 = 2*sqrt(sigma); 
        r2 = 3*sqrt(sigma);
        r3 = 1.2*sqrt(sigma);
    case 0.03
        r1 = 2.2*sqrt(sigma); 
        r2 = 2.8*sqrt(sigma); 
        r3 = 1.4*sqrt(sigma); 
    case 0.02
        r1 = 2.5*sqrt(sigma); 
        r2 = 2.5*sqrt(sigma); 
        r3 = 1.6*sqrt(sigma); 
    case 0.01
        r1 = 3*sqrt(sigma);
        r2 = 3.6*sqrt(sigma);
        r3 = 2*sqrt(sigma);
end


num_algo = numel(algos);

avgdists = -ones(num_algo, NumTrials);
maxdists = -ones(num_algo, NumTrials);
ts = -ones(num_algo, NumTrials);

Mouts = cell(num_algo,NumTrials);
proj_Mouts = cell(num_algo,NumTrials);
infos = cell(num_algo,NumTrials);
Dist2 = cell(num_algo,NumTrials);
Dist2_move = cell(num_algo,NumTrials);


for rep = 1 : NumTrials

fprintf('------ Trial %d ------\n',rep);

fname = sprintf('simulations/sphere/t%d_s%.2f_trial%d.mat',tau, sigma, rep);

try 
    load(fname)
    NumSample = size(samples,2);
    NumIni = size(data_ini,2);
    fprintf('load data %d\n', rep);
    
catch
    
    % generate data
    NumSample = 1000;
    NumIni = 1000;
    rng(rep);
    samples = randn(3, NumSample);
    samples = samples*diag(1./sqrt(sum(samples.^2)))+ sigma*randn(3, NumSample);
    
    data_ini = randn(3, NumIni);
    data_ini = data_ini*diag(1./sqrt(sum(data_ini.^2)))+0.5*sqrt(sigma)/sqrt(D)*(2*rand(3,NumIni)-1);

    
    save(fname,'samples','data_ini');
    fprintf('generate data %d\n', rep);
end


%%

% parameters for algorithm
opts.epsilon = 1e-16;
opts.diff_tol = 1e-4;
opts.maxiter = 50;
opts.display = 0;
opts.initer = 10;

dim = 2;

for i = 1 : num_algo
    algo = algos{i};
    
    tic;
    switch algo
        case 'ours'
            [Mout, info] = manfit_our(samples, dim, r1, data_ini, opts);
        case 'cf18'
            [Mout, info] = manfit_cf18(samples, dim, r2, data_ini, opts);
        case 'km17'
            [Mout, info] = manfit_km17(samples, dim, r3, data_ini, opts);
        case 'uo11'
            Mout=pc_project_multidim(samples,data_ini,r,dim);
            Mout = Mout';
            info.moveflag = true(1,size(Mout,2));
    end
    t = toc;
    
    Mouts{i, rep} = Mout;
    infos{i, rep} = info;
    ts(i, rep) = t;
    
    fprintf('Trial %d with algo %s costs %.2f seconds \n', rep, algo, t);
     
 %   save(sprintf('out/sphere/Dist_s%.2f.mat', sigma),...
 %   'Mouts','infos', 'ts', 'Dist2', 'Dist2_move', 'avgdists','maxdists');
    
end
                        
%% calculate distances
for i = 1 : num_algo
    Mout = Mouts{i,rep};
    moveflag = infos{i,rep}.moveflag;
    
    proj_Mout = bsxfun(@times, Mout, tau./sqrt(sum(Mout.^2)));
    proj_Mouts{i,rep} = proj_Mout;
    Dist2{i,rep} = sqrt(sum((Mout-proj_Mout).^2));
    temp = Dist2{i,rep}(moveflag);
    Dist2_move{i,rep} = temp;
    
    fprintf('Trial %d, algo %s: Max = %.8f, Avg = %.8f, Num = %d \n',rep, algos{i}, max(temp), mean(temp), length(temp));
    
    avgdists(i,rep) = mean(temp);
    maxdists(i,rep) = max(temp);
end

end

 save(sprintf('out/sphere/Dist_s%.2f.mat', sigma),...
    'algos','Mouts','proj_Mouts','infos', 'ts', 'Dist2', 'Dist2_move', 'avgdists','maxdists','r1','r2','r3');

 
%% boxplot
MatName = sprintf('out/sphere/Dist_s%.2f.mat', sigma);
PyName = sprintf('out/sphere/Manapprox_Dist_s%.2f.mat',sigma);

D1 = load(MatName);
D2 = load(PyName);

algos = D1.algos;
algos{numel(algos)+1} = 'ya21(deg=1)';
algos{numel(algos)+1} = 'ya21(deg=2)';
num_algo = numel(algos);

maxdists = [D1.maxdists; D2.max_dist];
NumTrials = size(maxdists,2);

ordered = [5, 1, 4, 2, 3];


if NumTrials > 1
    figure;
    h = boxplot(maxdists(ordered,:)');
    a = gca;
    set(h, 'LineWidth', 1.5)
    for i = 1 : numel(algos)
        a.XTickLabel{i} = algos{ordered(i)};
    end
    a.XAxis.FontSize = 16;
    a.YAxis.FontSize = 16;
    %a.YScale = 'log';
    %ylim([0.001,0.1])
    title(['d = 2, \sigma = ', num2str(sigma)], 'FontSize',16)
   
    sname = sprintf('figures/sphere_max_s%.2f.fig', sigma);
    saveas(gcf,sname)

end


