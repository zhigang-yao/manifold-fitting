function Mout = manfit_ours(sample, sig, sample_init,op_average)
if nargin <= 3
   op_average = 1;
end

Mout  = sample_init;  N = size(sample_init,1); 
N0 = size(sample,1);  ns = 1:N0;

r =  5*sig/log10(N0); R =  10*sig*sqrt(log(1/sig))/log10(N0);

parfor ii = 1:N

    x = sample_init(ii,:);

    dists = pdist2(x,sample);

    IDX1 = dists<2*r; IDX1 = ns(IDX1);

    IDX2 = knnsearch(sample,x,'K',5);

    IDX = union(IDX1,IDX2);

    BNbr = sample(IDX,:);

    xbar = mean(BNbr,1) + eps;

    dx = x - xbar; dx = dx/norm(dx);
    
    Q = [dx' null(dx)];
    
    sample_s = sample - x;
    
    sample_s = sample_s*Q;
    
    CNbr = (abs(sample_s(:,1)) < R & sum(sample_s(:,2:end).^2,2)< r^2); %Cylinder.

    if sum(CNbr)>10
    
       Mout(ii,:) = mean(sample(CNbr,:),1);

    else

       Mout(ii,:) = xbar;     
    
    end

end

if op_average == 1
    Mout_new = Mout;


    parfor ii = 1:N
        x = Mout(ii,:);

        dists = pdist2(x,Mout);

        IDX1 = dists<r/4; IDX1 = ns(IDX1);

        IDX2 = knnsearch(Mout,x,'K',5);

        IDX = union(IDX1,IDX2);

        BNbr = Mout(IDX,:);

        Mout_new(ii,:) = mean(BNbr,1);
    end

    Mout = Mout_new;
end

