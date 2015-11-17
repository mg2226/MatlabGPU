rho0=10.0;
Ro=200.0;

N123=200;
deltax123=2.74;

k1=0;
k2=1/deltax123; %so the symbolic and FFT formulas have approx the same range
deltak=(k2-k1)/5000.0;

[vecPbysymbolic,ksymbolic,vecPbyFFT,kFFT]=test_sphericalobjrecipspace(rho0,Ro,k1,deltak,k2,N123,deltax123);

if isreal(vecPbysymbolic)
  fprintf(1,'test_sphericalobjrecipspace_scr: vecPbysymbolic is real\n');
else
  fprintf(1,'test_sphericalobjrecipspace_scr: vecPbysymbolic is complex\n');
end

%delete points that are very small relative to the peak so that the log plot has less dynamic range
ii=find(abs(vecPbysymbolic)<(1.0e-6)*max(abs(vecPbysymbolic(:))));
vecPbysymbolic(ii)=NaN;

figure;
semilogy(ksymbolic,abs(vecPbysymbolic));
xlabel('k');
ylabel('P(k) (symbolic)');

figure;
semilogy(ksymbolic(1:500),abs(vecPbysymbolic(1:500)));
xlabel('k');
ylabel('P(k) (symbolic)');

figure;
semilogy(kFFT,abs(vecPbyFFT));
xlabel('k');
ylabel('P(k) (FFT)');

figure;
semilogy(kFFT(1:30),abs(vecPbyFFT(1:30)));
xlabel('k');
ylabel('P(k) (FFT)');

figure;
semilogy(ksymbolic,abs(vecPbysymbolic),kFFT,abs(vecPbyFFT));
xlabel('k');
ylabel('P(k)');
legend('symbolic','FFT');

figure;
semilogy(ksymbolic(1:500),abs(vecPbysymbolic(1:500)),kFFT(1:30),abs(vecPbyFFT(1:30)));
xlabel('k');
ylabel('P(k)');
legend('symbolic','FFT');
