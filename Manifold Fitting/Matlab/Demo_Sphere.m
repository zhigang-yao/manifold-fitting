% bjli @ Feb 11 2023.
% bjlistat@nus.edu.sg
% refered to Z. Yao, J. Su, B. Li, and S.T., Yau, Manifold Fitting

Jobs = {
    1, 'Visualization of the output';...
    2, 'Boxchart of Comparision';...
    3, 'Boxchart asymptotic property N Comparision';...
    };

job = Jobs{i_job,2};

switch  job
    case 'Visualization of the output'
        for fakeloop = 1

            rng(fakeloop);

            D = 3; dim = 2; tau = 1; sigma = 0.06;  r = 2*sqrt(tau*sigma); loadopt;
            algos = {'ysl22','yx19','cf18','km17'};  num_algo = numel(algos);

            NumSample = 100000;
            NumIni = 10000;
            samples = randn(3, NumSample);
            samples = samples*diag(1./sqrt(sum(samples.^2)))+ sigma*randn(3, NumSample);

            data_ini = randn(3, 2*NumIni);
            data_ini = data_ini*diag(1./sqrt(sum(data_ini.^2)))+2*sigma/sqrt(D)*(2*rand(3,2*NumIni)-1);
            proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
            norm_n2 = sum((data_ini - proj_data_ini).^2,1);
            [~, Index] = sort(norm_n2,'descend');
            data_ini = data_ini(:,Index(1:NumIni));

            for ii = 1%1:num_algo

                algo = algos{ii};

                switch algo
                    case 'ysl22'
                        Mout = manfit_ours(samples', sigma, data_ini',0); Mout = Mout';
                    case 'yx19'
                        [Mout, info] = manfit_yx23(samples, dim, r,  data_ini, opts);
                    case 'cf18'
                        [Mout, info] = manfit_cf18(samples, dim, r,  data_ini, opts);
                    case 'km17'
                        [Mout, info] = manfit_km17(samples, dim, r,  data_ini, opts);
                end

                subplot(1,2,2); hold on; box off; axis off;
                proj_Mout = bsxfun(@times, Mout, tau./sqrt(sum(Mout.^2)));
                plot(proj_Mout(1,:),proj_Mout(2,:),'k.'); hold on;
                plot(Mout(1,:),Mout(2,:),'r.');
                xlim([-1.2,1.2]);  ylim([-1.2,1.2]);
                pause(0.1),
            end
        end

    case 'Boxchart of Comparision'
        for fakeloop = 1
            % parameters for data
            NumTrials = 10;

            D = 3; dim = 2; tau = 1; sigma = 0.06; %[0.06;0.04;0.02]

            r = 2*sqrt(tau*sigma); loadopt;

            algos = {'ysl22','yx19','cf18','km17'};  num_algo = numel(algos);

            Mouts = cell(num_algo,NumTrials);
            Times = zeros(num_algo,NumTrials);
            avgdists = zeros(num_algo,NumTrials);
            maxdists = zeros(num_algo,NumTrials);

            for rep = 1:NumTrials

                fprintf('------ Trial %d ------\n',rep);

                rng(rep);

                NumSample = 1000;
                NumIni = 1000;
                samples = randn(3, NumSample);
                samples = samples*diag(1./sqrt(sum(samples.^2)))+ sigma*randn(3, NumSample);

                data_ini = randn(3, 2*NumIni);
                data_ini = data_ini*diag(1./sqrt(sum(data_ini.^2)))+2*sigma/sqrt(D)*(2*rand(3,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini = data_ini(:,Index(1:NumIni));

                for ii = 1:4

                    algo = algos{ii};

                    tic;

                    switch algo
                        case 'ysl22'
                            Mout = manfit_ours(samples', sigma, data_ini'); Mout = Mout';
                        case 'yx19'
                            [Mout, info] = manfit_yx23(samples, dim, r,  data_ini, opts);
                        case 'cf18'
                            [Mout, info] = manfit_cf18(samples, dim, r,  data_ini, opts);
                        case 'km17'
                            [Mout, info] = manfit_km17(samples, dim, r,  data_ini, opts);
                    end

                    t = toc;
                    Mouts{ii, rep} = Mout;
                    Times(ii,rep) = t;
                end

                for ii = 1:4
                    Mout = Mouts{ii, rep};

                    proj_Mout = bsxfun(@times, Mout, tau./sqrt(sum(Mout.^2)));
                    temp = sqrt(sum((Mout-proj_Mout).^2));

                    avgdists(ii,rep) = mean(temp);
                    maxdists(ii,rep) = max(temp);
                end
            end
            subplot(2,3,4);boxchart(maxdists'); box on;
            xticklabels({'ysl22','yx19','cf18','km17'});
            title('Sphere, Hauadorff Distance');
            subplot(2,3,5);boxchart(avgdists');  box on;
            xticklabels({'ysl22','yx19','cf18','km17'});
            title('Sphere, Average Distance');
            subplot(2,3,6);bar(mean(Times,2));  box on;
            xticklabels({'ysl22','yx19','cf18','km17'});
            title('Sphere, CPU Time (s)');
        end

    case 'Boxchart asymptotic property N Comparision'
        for fakeloop = 1
            NumTrials = 10;

            D = 3; dim = 2; tau = 1; sigma = 0.06;

            r = 2*sqrt(tau*sigma); loadopt;

            NumS = {'N = 1000','N = 2000','N = 3000'};  num_NumS = numel(NumS);

            Mouts1 = cell(num_NumS,NumTrials); Times1 = zeros(num_NumS,NumTrials);
            Mouts2 = cell(num_NumS,NumTrials); Times2 = zeros(num_NumS,NumTrials);

            avgdists1 = zeros(num_NumS,NumTrials); maxdists1 = zeros(num_NumS,NumTrials);
            avgdists2 = zeros(num_NumS,NumTrials); maxdists2 = zeros(num_NumS,NumTrials);

            for rep = 1:NumTrials

                fprintf('------ Trial %d ------\n',rep);

                rng(rep);

                % generate data
                NumSamples = [1000, 2000, 3000];
                NumIni = 1000;

                t1 = randn(3, NumSamples(1));
                samples1 = t1*diag(1./sqrt(sum(t1.^2)))+ sigma*randn(3, NumSamples(1));

                t2 = randn(3, NumSamples(2));
                samples2 = t2*diag(1./sqrt(sum(t2.^2)))+ sigma*randn(3, NumSamples(2));

                t3 = randn(3, NumSamples(3));
                samples3 = t3*diag(1./sqrt(sum(t3.^2)))+ sigma*randn(3, NumSamples(3));

                t = randn(3, 2*NumIni);
                data_ini = t*diag(1./sqrt(sum(t.^2)))+2*sigma/sqrt(D)*(2*rand(3,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini = data_ini(:,Index(1:NumIni));


                for ii = 1:3

                    num = NumS{ii};

                    switch num
                        case 'N = 1000'
                            tic; Mout1 = manfit_ours(samples1', sigma, data_ini'); Mout1 = Mout1'; t1 = toc;
                            tic; Mout2 = manfit_yx23(samples1, dim, r,  data_ini, opts); t2 = toc;
                        case 'N = 2000'
                            tic; Mout1 = manfit_ours(samples2', sigma, data_ini'); Mout1 = Mout1'; t1 = toc;
                            tic; Mout2 = manfit_yx23(samples2, dim, r,  data_ini, opts); t2 = toc;
                        case 'N = 3000'
                            tic; Mout1 = manfit_ours(samples3', sigma, data_ini'); Mout1 = Mout1'; t1 = toc;
                            tic; Mout2 = manfit_yx23(samples3, dim, r,  data_ini, opts); t2 = toc;
                    end

                    Mouts1{ii, rep} = Mout1; Times1(ii,rep) = t1;
                    Mouts2{ii, rep} = Mout2; Times2(ii,rep) = t2;
                end

                for ii = 1:3
                    Mout = Mouts1{ii, rep};
                    proj_Mout = bsxfun(@times, Mout, tau./sqrt(sum(Mout.^2)));
                    temp = sqrt(sum((Mout-proj_Mout).^2));
                    avgdists1(ii,rep) = mean(temp);
                    maxdists1(ii,rep) = max(temp);

                    Mout = Mouts2{ii, rep};
                    proj_Mout = bsxfun(@times, Mout, tau./sqrt(sum(Mout.^2)));
                    temp = sqrt(sum((Mout-proj_Mout).^2));
                    avgdists2(ii,rep) = mean(temp);
                    maxdists2(ii,rep) = max(temp);

                end
            end
            subplot(2,3,4);
            T = Res2Tab(maxdists1,maxdists2);
            boxchart(T.NumSample, T.Value, 'GroupByColor',T.Method);
            title('Shpere, Hauadorff Distance');
            xticks([1,2,3]);
            xticklabels(NumS);
            legend('ysl22','yx19');

            subplot(2,3,5);
            T = Res2Tab(avgdists1,avgdists2);
            boxchart(T.NumSample, T.Value, 'GroupByColor',T.Method);
            title('Sphere, Average Distance');
            xticks([1,2,3]);
            xticklabels(NumS);
            legend('ysl22','yx19');

            subplot(2,3,6);
            Time = [mean(Times1,2) mean(Times2,2)];
            b = bar(Time);
            b(1).FaceColor = '#CCE3F2';
            b(2).FaceColor = '#F7DDD1';
            xticks([1,2,3]);
            xticklabels(NumS);
            title('Sphere, CPU Time (s)');
            legend('ysl22','yx19');
        end


end
