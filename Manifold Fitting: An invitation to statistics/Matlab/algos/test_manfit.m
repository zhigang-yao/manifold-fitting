D = 2; dim = 1; tau = 1; sigma = 0.06; sig = sigma;rng(1);

NumSample = 30000;
NumIni = 300;
t = rand(1,NumSample)*2*pi;
sample_ = [cos(t);sin(t)];
sample = sample_ + sigma*randn(2, NumSample);

t = rand(1,2*NumIni)*2*pi;
data_ini = [cos(t);sin(t)] + 2*sigma/sqrt(D)*(2*rand(2,2*NumIni)-1);
proj_data_ini = bsxfun(@times, data_ini, tau./sqrt(sum(data_ini.^2)));
norm_n2 = sum((data_ini - proj_data_ini).^2,1);
[~, Index] = sort(norm_n2,'descend');
data_ini = data_ini(:,Index(1:NumIni));


sample_ = sample_'; sample = sample'; data_ini = data_ini';

for ii = 1:300

clf;
box on; hold on;
scatter(sample(:,1),sample(:,2),2,'MarkerEdgeColor',[200,200,200]/256);
scatter(data_ini(:,1),data_ini(:,2),2,'MarkerEdgeColor',[100,100,100]/256);
scatter(sample_(:,1),sample_(:,2),2,'MarkerEdgeColor','k');

x = data_ini(ii,:);

dists = pdist2(x,sample); BNbr = sample(dists<10*sig/log10(NumSample),:); 

xbar = mean(BNbr,1);

a = plot(x(1),x(2),'r.','MarkerSize',20);

b = plot(xbar(:,1), xbar(:,2),'b.','MarkerSize',20);


dx = x - xbar; dx = dx/norm(dx); ndx = null(dx);
Q = [dx' ndx];
sample_s = sample - xbar;
sample_s = sample_s*Q;
CNbr = (abs(sample_s(:,1)) < 0.9 & sum(sample_s(:,2:end).^2,2)< 0.01^2); %Cylinder.
B = sample(CNbr,:);
C = plot(B(:,1), B(:,2),'r.','MarkerSize',5);

upoo = xbar + 6*sig*(dx/norm(dx));
dooo = xbar - 6*sig*(dx/norm(dx));
plot([upoo(1) dooo(1)],[upoo(2) dooo(2)],'k--','linewidth',1.5);

pause,
end