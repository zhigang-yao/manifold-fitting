load('pts.mat'); pts = pts./(max(pts)); pts = pts((pts(:,1)==0|pts(:,1)==1),1:3); 

% figure(1); subplot(1,2,1);
% 
% plot3(pts(:,1),pts(:,2),pts(:,3),'r.');
% 
% hold on; title('Before manifold fitting');
% 
% Mout = manfit_ours(pts, 0.002, pts); 
% 
% subplot(1,2,2); plot3(Mout(:,1),Mout(:,2),Mout(:,3),'r.'); 
% 
% title('After manifold fitting');
% 
% hold on;

figure(1); subplot(2,2,1);

plot(pts(pts(:,1)==0,2),pts(pts(:,1)==0,3),'r.');

hold on; title('The 1st slice, Before manifold fitting');

subplot(2,2,2);

plot(pts(pts(:,1)==1,2),pts(pts(:,1)==1,3),'r.');

hold on; title('The last slice, Before manifold fitting');

Mout = manfit_ours(pts, 0.002, pts); 

subplot(2,2,3);

plot(Mout(pts(:,1)==0,2),Mout(pts(:,1)==0,3),'r.');

hold on; title('The 1st slice, After manifold fitting');

subplot(2,2,4);

plot(Mout(pts(:,1)==1,2),Mout(pts(:,1)==1,3),'r.');

hold on; title('The last slice, After manifold fitting');

