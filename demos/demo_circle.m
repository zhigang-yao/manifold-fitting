% To generate Figure 5.
% The figure might be somewhat different because of the random noise
%
% xiayq @ 8/19/2019
%
% xiayq@zju.edu.cn
% refered to Yao, Z and Xia, Y. (2019). Manifold Fitting under Unbounded Noise, arXiv:1909.10228


clear; %clc

% parameters for data
NumTrials = 20;
tau = 1; sigma = 0.01; % or sigma = 0.04;
D = 2;
r = sqrt(sigma*tau);

% method setup
algos = {'ours','cf18','km17'};
num_algo = numel(algos);

avgdists = -ones(num_algo, NumTrials);
maxdists = -ones(num_algo, NumTrials);
ts = -ones(num_algo, NumTrials);

Mouts = cell(num_algo,NumTrials);
infos = cell(num_algo,NumTrials);
Dist2 = cell(num_algo,NumTrials);
Dist2_move = cell(num_algo,NumTrials);


for rep = 1 : NumTrials

fprintf('------ Trial %d ------\n',rep);

fname = sprintf('simulations/circle/circle_t%d_s%.2f_r%.2f_trial%d.mat',tau, sigma, r, rep);

try 
    load(fname)
    NumSample = size(samples,2);
    NumIni = size(data_ini,2);
    fprintf('load data %d\n', rep);
    
catch
    
    % generate data
    NumSample = 1000;
    NumIni = 1000;
    t = rand(1,NumSample)*2*pi;
    samples = [cos(t);sin(t)]+ sigma*randn(2, NumSample);
    
    t = rand(1,NumIni)*2*pi;
    data_ini = [cos(t);sin(t)] + 0.5*sqrt(sigma)/sqrt(D)*(2*rand(2,NumIni)-1);
 
%     samples = samples*tau;
%     data_ini = data_ini*tau;
    
    save(fname,'samples','data_ini');
    fprintf('generate data %d\n', rep);
end

%%

dim = 1;
% parameters for algorithm
opts.epsilon = 1e-8;
opts.diff_tol = 1e-4;
opts.maxiter = 50;
opts.display = 0;
opts.initer = 10;
opts.beta = dim+2;



for i = 1 : num_algo
    algo = algos{i};
    
    tic;
    switch algo
        case 'ours'
            [Mout, info] = manfit_our(samples, dim, r, data_ini, opts);
        case 'cf18'
            [Mout, info] = manfit_cf18(samples, dim, r, data_ini, opts);
        case 'km17'
            [Mout, info] = manfit_km17(samples, dim, r, data_ini, opts);
        case 'uo11'
            %Mout=pc_project_multidim(samples,data_ini,kernel_sigma,dim);
            Mout=pc_project_multidim(samples,data_ini,r,dim);
            Mout = Mout';
            info.moveflag = true(1,size(Mout,2));
    end
    t = toc;
    
    Mouts{i, rep} = Mout;
    infos{i, rep} = info;
    ts(i, rep) = t;
    
    fprintf('Trial %d with algo %s costs %.2f seconds \n', rep, algo, t);
    
%     save(sprintf('out/circle/Dist_t%d_s%.2f_r%.2f.mat',tau, sigma, r),...
%     'Mouts','infos', 'ts', 'Dist2', 'Dist2_move', 'avgdists','maxdists');
    
end
                        
%% calculate distances
for i = 1 : num_algo
    Mout = Mouts{i,rep};
    moveflag = infos{i,rep}.moveflag;
    
    proj_Mout = bsxfun(@times, Mout, tau./sqrt(sum(Mout.^2)));
    Dist2{i,rep} = sqrt(sum((Mout-proj_Mout).^2));
    temp = Dist2{i,rep}(moveflag);
    Dist2_move{i,rep} = temp;
    
    fprintf('Trial %d, algo %s: Max = %.8f, Avg = %.8f, Num = %d \n',rep, algos{i}, max(temp), mean(temp), length(temp));
    
    avgdists(i,rep) = mean(temp);
    maxdists(i,rep) = max(temp);
end

end

%  save(sprintf('out/circle/Dist_t%d_s%.2f_r%.2f.mat',tau, sigma, r),...
%     'Mouts','infos', 'ts', 'Dist2', 'Dist2_move', 'avgdists','maxdists');
 
%%
if NumTrials > 1
    
    figure;
    boxplot(maxdists');
    a = gca;
    for i = 1 : num_algo
        a.XTickLabel{i} = algos{i};
    end
    a.XAxis.FontSize = 18;
    a.YAxis.FontSize = 14;
    %a.YScale = 'log';
    %ylim([1e-4, 1e-2]);
    if sigma == 0.04
        title('d=1, \sigma = 0.04', 'FontSize',16)
    elseif sigma == 0.01
        title('d=1, \sigma = 0.01', 'FontSize',16)
    end
%     sname = sprintf('figures/circle_max_s%.2f.fig',sigma);
%     saveas(gcf,sname)

end


