n=4; s=0.1; alp=1; ca=cos(alp); sa=sin(alp); % Projection

[theta,xi]=meshgrid(-1.5:s:1.5,1*(pi/2)*(0:1:16)/16);

z=theta+xi*1i;

% Color scheme
%tt=2*pi*(1:200)'/200; co=.5+.5*[cos(tt) cos(tt+1) cos(tt+2)];
%colormap(co)
colormap summer;

%% Plot
for k1=0:(n-1)
    for k2=0:(n-1)
        z1=exp(k1*2*pi*1i/n)*cosh(z).^(2/n);
        z2=exp(k2*2*pi*1i/n)*(sinh(z)/1i).^(2/n);
        sample_=real(z1);
        Y=real(z2);
        Z=ca*imag(z1)+sa*imag(z2);
        h=surf(sample_,Y,Z,z1*0+(k1+k2*n));hold on; %pause,
        set(h,'EdgeAlpha',0.4)
    end
end

view([2 3 1])
camlight
h=camlight('left');
set(h,'Color',[1 1 1]*.5)
axis equal
axis vis3d
axis off