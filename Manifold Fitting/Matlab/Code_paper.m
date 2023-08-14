clear; clc;
Colors = ['#F27970';'#BB9727';'#54B345';'#32B897';'#05B9E2';'#8983BF';'#C76DA2'];



Jobs = {
    1, 'View asymptotic property Circle N';...
    2, 'Boxchart asymptotic property Circle N';...
    3, 'Boxchart asymptotic property Cirlce S';...

    4, 'View asymptotic property Sphere N';...
    5, 'Boxchart asymptotic property Sphere N';...
    6, 'Boxchart asymptotic property Sphere S';...

    7, 'View asymptotic property Torus N';...
    8, 'Boxchart asymptotic property Torus N';...
    9, 'Boxchart asymptotic property Torus S';...

    };

i_job = 9;
job = Jobs{i_job,2};

switch  job

    case 'View asymptotic property Circle N'
        for fakeloop = 1

            NumTrials = 1;

            D = 2; dim = 1; tau = 1; sigma = 0.06;

            r = 2*sqrt(tau*sigma); loadopt;

            NumS = {'N = 300','N = 3000','N = 30000','N = 300000'};  num_NumS = numel(NumS);

            Mouts1 = cell(num_NumS,NumTrials); Times1 = zeros(num_NumS,NumTrials);

            avgdists1 = zeros(num_NumS,NumTrials); maxdists1 = zeros(num_NumS,NumTrials);

            for rep = 1:NumTrials

                fprintf('------ Trial %d ------\n',rep);

                rng(rep);

                NumSamples = [300, 3000, 30000 300000];
                NumIni = 300;

                t1 = rand(1,NumSamples(1))*2*pi;
                samples1 = [cos(t1);sin(t1)]+ sigma*randn(2, NumSamples(1));

                t2 = rand(1,NumSamples(2))*2*pi;
                samples2 = [cos(t2);sin(t2)]+ sigma*randn(2, NumSamples(2));

                t3 = rand(1,NumSamples(3))*2*pi;
                samples3 = [cos(t3);sin(t3)]+ sigma*randn(2, NumSamples(3));

                t4 = rand(1,NumSamples(4))*2*pi;
                samples4 = [cos(t4);sin(t4)]+ sigma*randn(2, NumSamples(4));

                t = rand(1,2*NumIni)*2*pi;
                data_ini = [cos(t);sin(t)] + 2*sigma/sqrt(D)*(2*rand(2,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini = data_ini(:,Index(1:NumIni));

                for ii = 1:4

                    num = NumS{ii};

                    switch num
                        case 'N = 300'
                            tic; Mout1 = manfit_ours(samples1', sigma, data_ini'); Mout1 = Mout1'; t1 = toc;
                        case 'N = 3000'
                            tic; Mout1 = manfit_ours(samples2', sigma, data_ini'); Mout1 = Mout1'; t1 = toc;
                        case 'N = 30000'
                            tic; Mout1 = manfit_ours(samples3', sigma, data_ini'); Mout1 = Mout1'; t1 = toc;
                        case 'N = 300000'
                            tic; Mout1 = manfit_ours(samples4', sigma, data_ini'); Mout1 = Mout1'; t1 = toc;

                    end

                    Mouts1{ii, rep} = Mout1; Times1(ii,rep) = t1;
                end

                for ii = 1:4
                    Mout = Mouts1{ii, rep};
                    proj_Mout = bsxfun(@times, Mout, tau./sqrt(sum(Mout.^2)));
                    subplot(3,4,ii);
                    hold on; box off; axis off;
                    plot(proj_Mout(1,:),proj_Mout(2,:),'k.'); hold on;
                    plot(Mout(1,:),Mout(2,:),'r.');
                    xlim([-1.2,1.2]);  ylim([-1.2,1.2]);
                    pause(0.1),
                end
            end
        end

    case 'Boxchart asymptotic property Circle N'
        for fakeloop = 1

            NumTrials = 50;

            D = 2; dim = 1; tau = 1; sigma = 0.06;

            r = 2*sqrt(tau*sigma); loadopt;

            NumS = {'N = 300','N = 3000','N = 30000','N = 300000'};
            
            num_NumS = numel(NumS);

            Mouts1 = cell(num_NumS,NumTrials); Times1 = zeros(num_NumS,NumTrials);

            avgdists1 = zeros(num_NumS,NumTrials); maxdists1 = zeros(num_NumS,NumTrials);

            for rep = 1:NumTrials

                fprintf('------ Trial %d ------\n',rep);

                rng(rep);

                NumSamples = [300, 3000, 30000 300000];
                NumIni = 300;

                t1 = rand(1,NumSamples(1))*2*pi;
                samples1 = [cos(t1);sin(t1)]+ sigma*randn(2, NumSamples(1));

                t2 = rand(1,NumSamples(2))*2*pi;
                samples2 = [cos(t2);sin(t2)]+ sigma*randn(2, NumSamples(2));

                t3 = rand(1,NumSamples(3))*2*pi;
                samples3 = [cos(t3);sin(t3)]+ sigma*randn(2, NumSamples(3));

                t4 = rand(1,NumSamples(4))*2*pi;
                samples4 = [cos(t4);sin(t4)]+ sigma*randn(2, NumSamples(4));

                t = rand(1,2*NumIni)*2*pi;
                data_ini = [cos(t);sin(t)] + 2*sigma/sqrt(D)*(2*rand(2,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini = data_ini(:,Index(1:NumIni));

                for ii = 1:4

                    num = NumS{ii};

                    switch num
                        case 'N = 300'
                            tic; Mout1 = manfit_ours(samples1', sigma, data_ini'); Mout1 = Mout1'; t1 = toc;
                        case 'N = 3000'
                            tic; Mout1 = manfit_ours(samples2', sigma, data_ini'); Mout1 = Mout1'; t1 = toc;
                        case 'N = 30000'
                            tic; Mout1 = manfit_ours(samples3', sigma, data_ini'); Mout1 = Mout1'; t1 = toc;
                        case 'N = 300000'
                            tic; Mout1 = manfit_ours(samples4', sigma, data_ini'); Mout1 = Mout1'; t1 = toc;

                    end

                    Mouts1{ii, rep} = Mout1; Times1(ii,rep) = t1;
                end

                for ii = 1:4
                    Mout = Mouts1{ii, rep};
                    proj_Mout = bsxfun(@times, Mout, tau./sqrt(sum(Mout.^2)));
                    temp = sqrt(sum((Mout-proj_Mout).^2));
                    avgdists1(ii,rep) = mean(temp);
                    maxdists1(ii,rep) = max(temp);

                end
            end
            subplot(2,2,3);
            T = Res2Tab_(maxdists1);
            boxchart(T.NumSample, T.Value);
            title('Hauadorff Distance (Circle, \sigma = 0.06)');
            xticks([1,2,3 4]);
            xticklabels({'3\times 10^2', '3\times 10^3', '3\times 10^4', '3\times 10^5'});


            subplot(2,2,4);
            T = Res2Tab_(avgdists1);
            boxchart(T.NumSample, T.Value);
            title('Average Distance (Circle, \sigma = 0.06)');
            xticks([1,2,3 4]);
            xticklabels({'3\times 10^2', '3\times 10^3', '3\times 10^4', '3\times 10^5'});
        end

    case 'Boxchart asymptotic property Cirlce S'
        for fakeloop = 1

            NumTrials = 50;

            D = 2; dim = 1; tau = 1;

            SigS = {'\sigma = 0.12','\sigma = 0.1','\sigma = 0.08','\sigma = 0.06','\sigma = 0.04','\sigma = 0.02'};  num_SigS = numel(SigS);

            Mouts1 = cell(num_SigS,NumTrials); Times1 = zeros(num_SigS,NumTrials);

            avgdists1 = zeros(num_SigS,NumTrials); maxdists1 = zeros(num_SigS,NumTrials);

            for rep = 1:NumTrials

                fprintf('------ Trial %d ------\n',rep);

                rng(rep);

                %generate data
                N = 30000;
                sigmas = [0.12 0.1 0.08 0.06 0.04 0.02];
                NumIni = 100;

                t1 = rand(1,N)*2*pi;
                samples1 = [cos(t1);sin(t1)]+ sigmas(1)*randn(2, N);

                t2 = rand(1,N)*2*pi;
                samples2 = [cos(t2);sin(t2)]+ sigmas(2)*randn(2, N);

                t3 = rand(1,N)*2*pi;
                samples3 = [cos(t3);sin(t3)]+ sigmas(3)*randn(2, N);

                t4 = rand(1,N)*2*pi;
                samples4 = [cos(t4);sin(t4)]+ sigmas(4)*randn(2, N);

                t5 = rand(1,N)*2*pi;
                samples5 = [cos(t5);sin(t5)]+ sigmas(5)*randn(2, N);

                t6 = rand(1,N)*2*pi;
                samples6 = [cos(t6);sin(t6)]+ sigmas(6)*randn(2, N);

                t = rand(1,2*NumIni)*2*pi;
                data_ini = [cos(t);sin(t)] + 2*sigmas(1)/sqrt(D)*(2*rand(2,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini1 = data_ini(:,Index(1:NumIni));

                t = rand(1,2*NumIni)*2*pi;
                data_ini = [cos(t);sin(t)] + 2*sigmas(2)/sqrt(D)*(2*rand(2,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini2 = data_ini(:,Index(1:NumIni));

                t = rand(1,2*NumIni)*2*pi;
                data_ini = [cos(t);sin(t)] + 2*sigmas(3)/sqrt(D)*(2*rand(2,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini3 = data_ini(:,Index(1:NumIni));

                t = rand(1,2*NumIni)*2*pi;
                data_ini = [cos(t);sin(t)] + 2*sigmas(4)/sqrt(D)*(2*rand(2,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini4 = data_ini(:,Index(1:NumIni));


                t = rand(1,2*NumIni)*2*pi;
                data_ini = [cos(t);sin(t)] + 2*sigmas(5)/sqrt(D)*(2*rand(2,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini5 = data_ini(:,Index(1:NumIni));


                t = rand(1,2*NumIni)*2*pi;
                data_ini = [cos(t);sin(t)] + 2*sigmas(6)/sqrt(D)*(2*rand(2,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini6 = data_ini(:,Index(1:NumIni));

                for ii = 1:6

                    i_sig = SigS{ii}; op = 0;

                    switch i_sig
                        case '\sigma = 0.12'
                            tic; Mout1 = manfit_ours(samples1', sigmas(1), data_ini1',op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.1'
                            tic; Mout1 = manfit_ours(samples2', sigmas(2), data_ini2',op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.08'
                            tic; Mout1 = manfit_ours(samples3', sigmas(3), data_ini3',op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.06'
                            tic; Mout1 = manfit_ours(samples4', sigmas(4), data_ini4',op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.04'
                            tic; Mout1 = manfit_ours(samples5', sigmas(5), data_ini5',op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.02'
                            tic; Mout1 = manfit_ours(samples6', sigmas(6), data_ini6',op); Mout1 = Mout1'; t1 = toc;

                    end

                    Mouts1{ii, rep} = Mout1; Times1(ii,rep) = t1;
                end

                for ii = 1:6
                    Mout = Mouts1{ii, rep};
                    proj_Mout = bsxfun(@times, Mout, tau./sqrt(sum(Mout.^2)));
                    temp = sqrt(sum((Mout-proj_Mout).^2));
                    avgdists1(ii,rep) = mean(temp);
                    maxdists1(ii,rep) = max(temp);

                end
            end

            subplot(2,2,1);
            T = Res2Tab_(maxdists1);
            boxchart(T.NumSample, T.Value);
            title('Hauadorff Distance (Circle, N = 3\times 10^4)');
            xticks([1,2,3,4,5,6]);
            xticklabels(SigS);

            subplot(2,2,2);
            T = Res2Tab_(avgdists1);
            boxchart(T.NumSample, T.Value);
            title('Average Distance (Circle, N = 3\times 10^4)');
            xticks([1,2,3,4,5,6]);
            xticklabels(SigS);

        end

    case 'View asymptotic property Sphere N'
        for fakeloop = 1

            NumTrials = 1;

            D = 3; dim = 2; tau = 1; sigma = 0.06;

            r = 2*sqrt(tau*sigma); loadopt;

            NumS = {'N = 1000','N = 5000','N = 25000','N = 125000'};  num_NumS = numel(NumS);

            Mouts1 = cell(num_NumS,NumTrials); Times1 = zeros(num_NumS,NumTrials);

            avgdists1 = zeros(num_NumS,NumTrials); maxdists1 = zeros(num_NumS,NumTrials);

            for rep = 1:NumTrials

                fprintf('------ Trial %d ------\n',rep);

                rng(rep);

                NumSamples = [1000, 5000, 25000 125000];
                NumIni = 1000;

                t1 = randn(3, NumSamples(1));
                samples1 = t1*diag(1./sqrt(sum(t1.^2)))+ sigma*randn(3, NumSamples(1));

                t2 = randn(3, NumSamples(2));
                samples2 = t2*diag(1./sqrt(sum(t2.^2)))+ sigma*randn(3, NumSamples(2));

                t3 = randn(3, NumSamples(3));
                samples3 = t3*diag(1./sqrt(sum(t3.^2)))+ sigma*randn(3, NumSamples(3));

                t4 = randn(NumSamples(4),3);
                samples4 = t4./sqrt(sum(t4.^2,2)) +  sigma*randn(NumSamples(4),3);
                samples4 = samples4';

                t = randn(3, 2*NumIni);
                data_ini = t*diag(1./sqrt(sum(t.^2)))+2*sigma/sqrt(D)*(2*rand(3,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini = data_ini(:,Index(1:NumIni));


                for ii = 1:4

                    num = NumS{ii};

                    switch num
                        case 'N = 1000'
                            tic; Mout1 = manfit_ours(samples1', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                        case 'N = 5000'
                            tic; Mout1 = manfit_ours(samples2', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                        case 'N = 25000'
                            tic; Mout1 = manfit_ours(samples3', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                        case 'N = 125000'
                            tic; Mout1 = manfit_ours(samples4', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                    end

                    Mouts1{ii, rep} = Mout1; Times1(ii,rep) = t1;
                end

                for ii = 1:4
                    Mout = Mouts1{ii, rep};
                    proj_Mout = bsxfun(@times, Mout, tau./sqrt(sum(Mout.^2)));
                    subplot(3,4,ii+4);
                    hold on; box off; axis off;
                    plot(proj_Mout(1,:),proj_Mout(2,:),'k.'); hold on;
                    plot(Mout(1,:),Mout(2,:),'r.');
                    xlim([-1.2,1.2]);  ylim([-1.2,1.2]);
                    pause(0.1),
                end
            end

        end

    case 'Boxchart asymptotic property Sphere N'
        for fakeloop = 1

            NumTrials = 50;

            D = 3; dim = 2; tau = 1; sigma = 0.06;

            r = 2*sqrt(tau*sigma); loadopt;

            NumS = {'N = 1000','N = 5000','N = 25000','N = 125000'};  num_NumS = numel(NumS);

            Mouts1 = cell(num_NumS,NumTrials); Times1 = zeros(num_NumS,NumTrials);

            avgdists1 = zeros(num_NumS,NumTrials); maxdists1 = zeros(num_NumS,NumTrials);

            for rep = 1:NumTrials

                fprintf('------ Trial %d ------\n',rep);

                rng(rep);

                NumSamples = [1000, 5000, 25000 125000];
                NumIni = 100;

                t1 = randn(3, NumSamples(1));
                samples1 = t1*diag(1./sqrt(sum(t1.^2)))+ sigma*randn(3, NumSamples(1));

                t2 = randn(3, NumSamples(2));
                samples2 = t2*diag(1./sqrt(sum(t2.^2)))+ sigma*randn(3, NumSamples(2));

                t3 = randn(3, NumSamples(3));
                samples3 = t3*diag(1./sqrt(sum(t3.^2)))+ sigma*randn(3, NumSamples(3));

                t4 = randn(NumSamples(4),3);
                samples4 = t4./sqrt(sum(t4.^2,2)) +  sigma*randn(NumSamples(4),3);
                samples4 = samples4';

                t = randn(3, 2*NumIni);
                data_ini = t*diag(1./sqrt(sum(t.^2)))+2*sigma/sqrt(D)*(2*rand(3,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini = data_ini(:,Index(1:NumIni));


                for ii = 1:4

                    num = NumS{ii};

                    switch num
                        case 'N = 1000'
                            tic; Mout1 = manfit_ours(samples1', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                        case 'N = 5000'
                            tic; Mout1 = manfit_ours(samples2', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                        case 'N = 25000'
                            tic; Mout1 = manfit_ours(samples3', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                        case 'N = 125000'
                            tic; Mout1 = manfit_ours(samples4', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                    end

                    Mouts1{ii, rep} = Mout1; Times1(ii,rep) = t1;
                end

                for ii = 1:4
                    Mout = Mouts1{ii, rep};
                    proj_Mout = bsxfun(@times, Mout, tau./sqrt(sum(Mout.^2)));
                    temp = sqrt(sum((Mout-proj_Mout).^2));
                    avgdists1(ii,rep) = mean(temp);
                    maxdists1(ii,rep) = max(temp);
                end
            end

            NumS = {'N = 1000','N = 5000','N = 25000','N = 125000'};            

            subplot(2,2,3);
            T = Res2Tab_(maxdists1);
            boxchart(T.NumSample, T.Value);
            title('Hauadorff Distance (Sphere, \sigma = 0.06)');
            xticks([1,2,3,4]);
            xticklabels({'1\times 10^3', '5\times 10^3', '2.5\times 10^4', '1.25\times 10^5'});

            subplot(2,2,4);
            T = Res2Tab_(avgdists1);
            boxchart(T.NumSample, T.Value);
            title('Average Distance (Sphere, \sigma = 0.06)');
            xticks([1,2,3,4]);
            xticklabels({'1\times 10^3', '5\times 10^3', '2.5\times 10^4', '1.25\times 10^5'});
        end

    case 'Boxchart asymptotic property Sphere S'
        for fakeloop = 1

            NumTrials = 50;

            D = 3; dim = 2; tau = 1;

            SigS = {'\sigma = 0.12','\sigma = 0.10','\sigma = 0.08','\sigma = 0.06','\sigma = 0.04','\sigma = 0.02'};  num_SigS = numel(SigS);

            Mouts1 = cell(num_SigS,NumTrials); Times1 = zeros(num_SigS,NumTrials);

            avgdists1 = zeros(num_SigS,NumTrials); maxdists1 = zeros(num_SigS,NumTrials);

            for rep = 1:NumTrials

                fprintf('------ Trial %d ------\n',rep);

                rng(rep);

                N = 25000;
                sigmas = [0.12 0.1 0.08 0.06 0.04 0.02];
                NumIni = 100;

                t1 = randn(3, N);
                samples1 = t1*diag(1./sqrt(sum(t1.^2)))+ sigmas(1)*randn(3, N);

                t2 = randn(3, N);
                samples2 = t2*diag(1./sqrt(sum(t2.^2)))+ sigmas(2)*randn(3, N);


                t3 = randn(3, N);
                samples3 = t3*diag(1./sqrt(sum(t3.^2)))+ sigmas(3)*randn(3, N);


                t4 = randn(3, N);
                samples4 = t4*diag(1./sqrt(sum(t4.^2)))+ sigmas(4)*randn(3, N);


                t5 = randn(3, N);
                samples5 = t5*diag(1./sqrt(sum(t5.^2)))+ sigmas(5)*randn(3, N);


                t6 = randn(3, N);
                samples6 = t6*diag(1./sqrt(sum(t6.^2)))+ sigmas(6)*randn(3, N);


                t = randn(3, 2*NumIni);
                data_ini = t*diag(1./sqrt(sum(t.^2)))+2*sigmas(1)/sqrt(D)*(2*rand(3,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini1 = data_ini(:,Index(1:NumIni));


                t = randn(3, 2*NumIni);
                data_ini = t*diag(1./sqrt(sum(t.^2)))+2*sigmas(2)/sqrt(D)*(2*rand(3,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini2 = data_ini(:,Index(1:NumIni));


                t = randn(3, 2*NumIni);
                data_ini = t*diag(1./sqrt(sum(t.^2)))+2*sigmas(3)/sqrt(D)*(2*rand(3,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini3 = data_ini(:,Index(1:NumIni));


                t = randn(3, 2*NumIni);
                data_ini = t*diag(1./sqrt(sum(t.^2)))+2*sigmas(4)/sqrt(D)*(2*rand(3,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini4 = data_ini(:,Index(1:NumIni));


                t = randn(3, 2*NumIni);
                data_ini = t*diag(1./sqrt(sum(t.^2)))+2*sigmas(5)/sqrt(D)*(2*rand(3,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini5 = data_ini(:,Index(1:NumIni));


                t = randn(3, 2*NumIni);
                data_ini = t*diag(1./sqrt(sum(t.^2)))+2*sigmas(6)/sqrt(D)*(2*rand(3,2*NumIni)-1);
                proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini6 = data_ini(:,Index(1:NumIni));


                for ii = 1:6

                    i_sig = SigS{ii}; op = 0;

                    switch i_sig
                        case '\sigma = 0.12'
                            tic; Mout1 = manfit_ours(samples1', sigmas(1), data_ini1',op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.10'
                            tic; Mout1 = manfit_ours(samples2', sigmas(2), data_ini2',op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.08'
                            tic; Mout1 = manfit_ours(samples3', sigmas(3), data_ini3',op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.06'
                            tic; Mout1 = manfit_ours(samples4', sigmas(4), data_ini4',op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.04'
                            tic; Mout1 = manfit_ours(samples5', sigmas(5), data_ini5',op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.02'
                            tic; Mout1 = manfit_ours(samples6', sigmas(6), data_ini6',op); Mout1 = Mout1'; t1 = toc;

                    end

                    Mouts1{ii, rep} = Mout1; Times1(ii,rep) = t1;
                end

                for ii = 1:6
                    Mout = Mouts1{ii, rep};
                    proj_Mout = bsxfun(@times, Mout, tau./sqrt(sum(Mout.^2)));
                    temp = sqrt(sum((Mout-proj_Mout).^2));
                    avgdists1(ii,rep) = mean(temp);
                    maxdists1(ii,rep) = max(temp);

                end
            end

            subplot(2,2,1);
            T = Res2Tab_(maxdists1);
            boxchart(T.NumSample, T.Value);
            title('Hauadorff Distance (Sphere, N = 2.5\times 10^4)');
            xticks([1,2,3,4,5,6]);
            xticklabels(SigS);

            subplot(2,2,2);
            T = Res2Tab_(avgdists1);
            boxchart(T.NumSample, T.Value);
            title('Average Distance (Sphere, N = 2.5\times 10^4)');
            xticks([1,2,3,4,5,6]);
            xticklabels(SigS);
        end

    case 'View asymptotic property Torus N'
        for fakeloop = 1

            NumTrials = 1;

            D = 3; dim = 2; a = 2/3; b = 1/3; tau = b; sigma = 0.06;

            NumS = {'N = 1000','N = 5000','N = 25000','N = 125000'};  num_NumS = numel(NumS);

            Mouts1 = cell(num_NumS,NumTrials); Times1 = zeros(num_NumS,NumTrials);

            avgdists1 = zeros(num_NumS,NumTrials); maxdists1 = zeros(num_NumS,NumTrials);

            for rep = 1:NumTrials

                fprintf('------ Trial %d ------\n',rep);

                rng(rep);

                NumSamples = [1000, 5000, 25000 125000];
                NumIni = 1000;

                samples = torusUnif(NumSamples(1), a, b);
                samples1 = samples + sigma*randn(D,NumSamples(1));

                samples = torusUnif(NumSamples(2), a, b);
                samples2 = samples + sigma*randn(D,NumSamples(2));

                samples = torusUnif(NumSamples(3), a, b);
                samples3 = samples + sigma*randn(D,NumSamples(3));

                samples = torusUnif(NumSamples(4), a, b);
                samples4 = samples + sigma*randn(D,NumSamples(4));


                data_ini = torusUnif(2*NumIni, a, b);
                data_ini = data_ini+ 2*sigma/sqrt(D)*(2*rand(D,2*NumIni)-1);
                proj_data_ini = pdtorus(a, b, data_ini);
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini = data_ini(:,Index(1:NumIni));


                for ii = 1:4

                    num = NumS{ii};

                    switch num
                        case 'N = 1000'
                            tic; Mout1 = manfit_ours(samples1', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                        case 'N = 5000'
                            tic; Mout1 = manfit_ours(samples2', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                        case 'N = 25000'
                            tic; Mout1 = manfit_ours(samples3', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                        case 'N = 125000'
                            tic; Mout1 = manfit_ours(samples4', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                    end

                    Mouts1{ii, rep} = Mout1; Times1(ii,rep) = t1;
                end

                for ii = 1:4
                    Mout = Mouts1{ii, rep};
                    [proj_Mout, tempdist] = pdtorus(a, b, Mout);
                    subplot(3,4,ii+8);
                    hold on; box off; axis off;
                    plot(proj_Mout(1,:),proj_Mout(2,:),'k.'); hold on;
                    plot(Mout(1,:),Mout(2,:),'r.');
                    xlim([-1.2,1.2]);  ylim([-1.2,1.2]);
                    pause(0.1),                    

                end
            end
        end

    case 'Boxchart asymptotic property Torus N'
        for fakeloop = 1

            NumTrials = 50;

            D = 3; dim = 2; a = 2/3; b = 1/3; tau = b; sigma = 0.06;

            NumS = {'N = 1000','N = 5000','N = 25000','N = 125000'};  num_NumS = numel(NumS);

            Mouts1 = cell(num_NumS,NumTrials); Times1 = zeros(num_NumS,NumTrials);

            avgdists1 = zeros(num_NumS,NumTrials); maxdists1 = zeros(num_NumS,NumTrials);

            for rep = 1:NumTrials

                fprintf('------ Trial %d ------\n',rep);

                rng(rep);

                NumSamples = [1000, 5000, 25000 125000];
                NumIni = 100;

                samples = torusUnif(NumSamples(1), a, b);
                samples1 = samples + sigma*randn(D,NumSamples(1));

                samples = torusUnif(NumSamples(2), a, b);
                samples2 = samples + sigma*randn(D,NumSamples(2));

                samples = torusUnif(NumSamples(3), a, b);
                samples3 = samples + sigma*randn(D,NumSamples(3));

                samples = torusUnif(NumSamples(4), a, b);
                samples4 = samples + sigma*randn(D,NumSamples(4));


                data_ini = torusUnif(2*NumIni, a, b);
                data_ini = data_ini+ 2*sigma/sqrt(D)*(2*rand(D,2*NumIni)-1);
                proj_data_ini = pdtorus(a, b, data_ini);
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini = data_ini(:,Index(1:NumIni));


                for ii = 1:4

                    num = NumS{ii};

                    switch num
                        case 'N = 1000'
                            tic; Mout1 = manfit_ours(samples1', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                        case 'N = 5000'
                            tic; Mout1 = manfit_ours(samples2', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                        case 'N = 25000'
                            tic; Mout1 = manfit_ours(samples3', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                        case 'N = 125000'
                            tic; Mout1 = manfit_ours(samples4', sigma, data_ini',0); Mout1 = Mout1'; t1 = toc;
                    end

                    Mouts1{ii, rep} = Mout1; Times1(ii,rep) = t1;
                end

                for ii = 1:4
                    Mout = Mouts1{ii, rep};
                    [proj_Mout, tempdist] = pdtorus(a, b, Mout);
                    temp = sqrt(sum((Mout-proj_Mout).^2));
                    avgdists1(ii,rep) = mean(temp);
                    maxdists1(ii,rep) = max(temp);
                end
            end
            subplot(2,2,3);
            T = Res2Tab_(maxdists1);
            boxchart(T.NumSample, T.Value);
            title('Hauadorff Distance (Torus, \sigma = 0.06)');
            xticks([1,2,3,4]);
            xticklabels({'1\times 10^3', '5\times 10^3', '2.5\times 10^4', '1.25\times 10^5'});

            subplot(2,2,4);
            T = Res2Tab_(avgdists1);
            boxchart(T.NumSample, T.Value);
            title('Average Distance (Torus, \sigma = 0.06)');
            xticks([1,2,3,4]);
            xticklabels({'1\times 10^3', '5\times 10^3', '2.5\times 10^4', '1.25\times 10^5'});

        end

    case 'Boxchart asymptotic property Torus S'
        for fakeloop = 1

            NumTrials = 50;

            D = 3; dim = 2; a = 2/3; b = 1/3; tau = b;

            SigS = {'\sigma = 0.12','\sigma = 0.10','\sigma = 0.08','\sigma = 0.06','\sigma = 0.04','\sigma = 0.02'};  num_SigS = numel(SigS);

            Mouts1 = cell(num_SigS,NumTrials); Times1 = zeros(num_SigS,NumTrials);

            avgdists1 = zeros(num_SigS,NumTrials); maxdists1 = zeros(num_SigS,NumTrials);

            for rep = 1:NumTrials

                fprintf('------ Trial %d ------\n',rep);

                rng(rep);

                N = 25000;
                sigmas = [0.12 0.1 0.08 0.06 0.04 0.02];
                NumIni = 100;

                samples = torusUnif(N, a, b);
                samples1 = samples + sigmas(1)*randn(D,N);

                samples = torusUnif(N, a, b);
                samples2 = samples + sigmas(2)*randn(D,N);

                samples = torusUnif(N, a, b);
                samples3 = samples + sigmas(3)*randn(D,N);

                samples = torusUnif(N, a, b);
                samples4 = samples + sigmas(4)*randn(D,N);

                samples = torusUnif(N, a, b);
                samples5 = samples + sigmas(5)*randn(D,N);

                samples = torusUnif(N, a, b);
                samples6 = samples + sigmas(6)*randn(D,N);

                data_ini = torusUnif(2*NumIni, a, b);
                data_ini = data_ini+ 2*sigmas(1)/sqrt(D)*(2*rand(D,2*NumIni)-1);
                proj_data_ini = pdtorus(a, b, data_ini);
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini1 = data_ini(:,Index(1:NumIni));

                data_ini = torusUnif(2*NumIni, a, b);
                data_ini = data_ini+ 2*sigmas(2)/sqrt(D)*(2*rand(D,2*NumIni)-1);
                proj_data_ini = pdtorus(a, b, data_ini);
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini2 = data_ini(:,Index(1:NumIni));

                data_ini = torusUnif(2*NumIni, a, b);
                data_ini = data_ini+ 2*sigmas(3)/sqrt(D)*(2*rand(D,2*NumIni)-1);
                proj_data_ini = pdtorus(a, b, data_ini);
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini3 = data_ini(:,Index(1:NumIni));


                data_ini = torusUnif(2*NumIni, a, b);
                data_ini = data_ini+ 2*sigmas(4)/sqrt(D)*(2*rand(D,2*NumIni)-1);
                proj_data_ini = pdtorus(a, b, data_ini);
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini4 = data_ini(:,Index(1:NumIni));

                data_ini = torusUnif(2*NumIni, a, b);
                data_ini = data_ini+ 2*sigmas(5)/sqrt(D)*(2*rand(D,2*NumIni)-1);
                proj_data_ini = pdtorus(a, b, data_ini);
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini5 = data_ini(:,Index(1:NumIni));


                data_ini = torusUnif(2*NumIni, a, b);
                data_ini = data_ini+ 2*sigmas(6)/sqrt(D)*(2*rand(D,2*NumIni)-1);
                proj_data_ini = pdtorus(a, b, data_ini);
                norm_n2 = sum((data_ini - proj_data_ini).^2,1);
                [~, Index] = sort(norm_n2,'descend');
                data_ini6 = data_ini(:,Index(1:NumIni));

                for ii = 1:6

                    i_sig = SigS{ii}; op = 0;

                    switch i_sig
                        case '\sigma = 0.12'
                            tic; Mout1 = manfit_ours(samples1', sigmas(1), data_ini1',op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.10'
                            tic; Mout1 = manfit_ours(samples2', sigmas(2), data_ini2',op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.08'
                            tic; Mout1 = manfit_ours(samples3', sigmas(3), data_ini3',op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.06'
                            tic; Mout1 = manfit_ours(samples4', sigmas(4), data_ini4',op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.04'
                            tic; Mout1 = manfit_ours(samples5', sigmas(5), data_ini5',op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.02'
                            tic; Mout1 = manfit_ours(samples6', sigmas(6), data_ini6',op); Mout1 = Mout1'; t1 = toc;

                    end

                    Mouts1{ii, rep} = Mout1; Times1(ii,rep) = t1;
                end

                for ii = 1:6
                    Mout = Mouts1{ii, rep};
                    [proj_Mout, tempdist] = pdtorus(a, b, Mout);
                    temp = sqrt(sum((Mout-proj_Mout).^2));
                    avgdists1(ii,rep) = mean(temp);
                    maxdists1(ii,rep) = max(temp);
                end
            end

            subplot(2,2,1);
            T = Res2Tab_(maxdists1);
            boxchart(T.NumSample, T.Value);
            title('Hauadorff Distance (Torus, N = 2.5\times 10^4)');
            xticks([1,2,3,4,5,6]);
            xticklabels(SigS);

            subplot(2,2,2);
            T = Res2Tab_(avgdists1);
            boxchart(T.NumSample, T.Value);
            title('Average Distance (Torus, N = 2.5\times 10^4)');
            xticks([1,2,3,4,5,6]);
            xticklabels(SigS);
        end

end
