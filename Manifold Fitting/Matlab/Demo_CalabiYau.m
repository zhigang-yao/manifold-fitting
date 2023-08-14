% bjli @ Feb 11 2023.
% bjlistat@nus.edu.sg
% refered to Z. Yao, J. Su, and B. Li, Manifold Fitting


Jobs = {
    1, 'Visualization of the Calabi-Yau';...
    2, 'Visualization of the Performance';...
    3, 'Boxchart asymptotic property S';...
    };

i_job = 3;

job = Jobs{i_job,2};


switch  job
    case 'Visualization of the Calabi-Yau'
        drawcy;

    case 'Visualization of the Performance'
        subplot(1,3,1);
        drawcy;
        subplot(1,3,2);
        drawcy; hold on; rng(1);
        D = 4; d = 2; s = 0.01; sig = 0.06; numb = 200; [samples,data_ini] = con_CalabiYau(s,sig,numb);
        scatter3(data_ini(:,1),data_ini(:,2), cos(1)*data_ini(:,3) + sin(1)*data_ini(:,4),20,'filled','MarkerEdgeColor','g');
        subplot(1,3,3);
        drawcy; hold on;
        Mout = manfit_ours(samples, sig, data_ini,0);
        scatter3(Mout(:,1),Mout(:,2), cos(1)*Mout(:,3) + sin(1)*Mout(:,4),20,'filled','MarkerEdgeColor','g');

    case 'Boxchart asymptotic property S'
        for fakeloop = 1

            NumTrials = 50;

             D = 4; d = 2; 

            SigS = {'\sigma = 0.03','\sigma = 0.025','\sigma = 0.02','\sigma = 0.015','\sigma = 0.01','\sigma = 0.005'};  num_SigS = numel(SigS);

            Mouts1 = cell(num_SigS,NumTrials); Times1 = zeros(num_SigS,NumTrials);

            avgdists1 = zeros(num_SigS,NumTrials); maxdists1 = zeros(num_SigS,NumTrials);

            for rep = 1:NumTrials

                fprintf('------ Trial %d ------\n',rep);

                rng(rep);

                s = 0.05;
                sigmas = 0.5*[0.06 0.05 0.04 0.03 0.02 0.01];
                NumIni = 100;
                [samples1,data_ini1] = con_CalabiYau(s,sigmas(1),NumIni);                
                [samples2,data_ini2] = con_CalabiYau(s,sigmas(2),NumIni);                
                [samples3,data_ini3] = con_CalabiYau(s,sigmas(3),NumIni);                
                [samples4,data_ini4] = con_CalabiYau(s,sigmas(4),NumIni);                
                [samples5,data_ini5] = con_CalabiYau(s,sigmas(5),NumIni);                
                [samples6,data_ini6] = con_CalabiYau(s,sigmas(6),NumIni);                


                for ii = 1:6

                    i_sig = SigS{ii}; op = 0;

                    switch i_sig
                        case '\sigma = 0.03'
                            tic; Mout1 = manfit_ours(samples1, sigmas(1), data_ini1,op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.025'
                            tic; Mout1 = manfit_ours(samples2, sigmas(2), data_ini2,op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.02'
                            tic; Mout1 = manfit_ours(samples3, sigmas(3), data_ini3,op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.015'
                            tic; Mout1 = manfit_ours(samples4, sigmas(4), data_ini4,op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.01'
                            tic; Mout1 = manfit_ours(samples5, sigmas(5), data_ini5,op); Mout1 = Mout1'; t1 = toc;
                        case '\sigma = 0.005'
                            tic; Mout1 = manfit_ours(samples6, sigmas(6), data_ini6,op); Mout1 = Mout1'; t1 = toc;

                    end

                    Mouts1{ii, rep} = Mout1; Times1(ii,rep) = t1;
                end

                for ii = 1:6
                    Mout = Mouts1{ii, rep};
                    temp = (real(1-((Mout(1,:)+Mout(3,:)*1i).^4+(Mout(2,:)+Mout(4,:)*1i).^4))).^2;
                    avgdists1(ii,rep) = mean(temp);
                    maxdists1(ii,rep) = max(temp);
                end
            end

            subplot(1,2,1);
            T = Res2Tab_(maxdists1);
            boxchart(T.NumSample, T.Value);
            title('Hauadorff Distance (Calabi-Yau, N = 313296)');
            xticks([1,2,3,4,5,6]);
            xticklabels(SigS);
            %set(gca,'Yscale','log');

            subplot(1,2,2);
            T = Res2Tab_(avgdists1);
            boxchart(T.NumSample, T.Value);
            title('Average Distance (Calabi-Yau, N = 313296)');
            xticks([1,2,3,4,5,6]);
            xticklabels(SigS);
            %set(gca,'Yscale','log');


        end

end