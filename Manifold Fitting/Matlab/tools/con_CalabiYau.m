function [sample, sample_init] = con_CalabiYau(s,sig,numb)
n = 4;  [theta,xi]=meshgrid(-1.5:s:1.5,1*(pi/2)*(0:s:16)/16);
z=theta+xi*1i; [row_z,col_z] = size(z);
Z1 = zeros(row_z,col_z,16); Z2 = zeros(row_z,col_z,16);  tt = 0;
for k1 = 0:3
    for k2 = 0:3
        tt = tt + 1;
        Z1(:,:,tt) = exp(k1*2*pi*1i/n)*cosh(z).^(2/n);
        Z2(:,:,tt) = exp(k2*2*pi*1i/n)*(sinh(z)/1i).^(2/n);
    end
end
sample_ = [real(Z1(:)) real(Z2(:)) imag(Z1(:)) imag(Z2(:))]; N0 = size(sample_,1);

noise = sig*randn(N0,4);

sample  = sample_ + noise;

init_seed = randperm(size(sample,1));

sample_init = sample(init_seed(1:numb),:);