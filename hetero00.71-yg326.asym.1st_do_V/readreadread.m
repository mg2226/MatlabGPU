function imagestack=readreadread(fn_imagestack)
%function imagestack=readreadread(fn_imagestack)
%Fake function that returns fake images.
%fn_imagestack argument is not used.

state4rng=rng; %save state of pseudorandom number generator
pseudorandomnumberseed=383267;
pseudorandomnumberseed=abs(fix(pseudorandomnumberseed));
rng(pseudorandomnumberseed);

Na=10;
Nb=11;
dimensions=[Na Nb];
delta=[2.57 2.57];
radius=0.2*max(delta)*max(dimensions);

[X1,X2]=ndgrid(delta(1)*([0:dimensions(1)-1]-(dimensions(1)-1)/2), ...
               delta(2)*([0:dimensions(2)-1]-(dimensions(2)-1)/2));
R=sqrt(X1.^2 + X2.^2);
mask=(R<=radius);

virusvalue=3.0;
pixelnoisestddev=1.0;
meanvalue=virusvalue;
scaling=1.5;

Nv=3;
imagestack=cell(Nv,1);
for ii=1:Nv
  imagestack{ii}=pixelnoisestddev.*randn(Na,Nb);
  fprintf(1,'readreadread: ii %d mean %g variance %g\n',ii,mean(imagestack{ii}(:)),cov(imagestack{ii}(:)));
end

for ii=1:Nv
  imagestack{ii}(mask)=imagestack{ii}(mask)+virusvalue;
end

imagestack{1}(1,1)=virusvalue;
imagestack{2}=meanvalue + scaling.*imagestack{2};
imagestack{3}(end,end)=virusvalue/2.0;

%fprintf(1,'readreadread: imagestack{1}:\n');
%disp(imagestack{1});
%fprintf(1,'readreadread: imagestack{2}:\n');
%disp(imagestack{2});
%fprintf(1,'readreadread: imagestack{3}:\n');
%disp(imagestack{3});

rng(state4rng); %restore state of pseudorandom number generator
