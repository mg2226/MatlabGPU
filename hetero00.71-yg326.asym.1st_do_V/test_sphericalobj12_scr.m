minindex1=[-100 -100 -100];
maxindex1=[100 100 100];
delta=2.74;
rho0=1.0;
R2=200.0;
projimage1=test_sphericalobj1(minindex1,maxindex1,delta,rho0,R2);

figure;
imagesc(projimage1);
axis equal
colorbar

minindex2=minindex1(1:2);
maxindex2=maxindex1(1:2);
projimage2=test_sphericalobj2(minindex2,maxindex2,delta,rho0,R2);

figure;
imagesc(projimage1);
axis equal
colorbar
