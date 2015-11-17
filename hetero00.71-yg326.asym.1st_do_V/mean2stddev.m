function stddev=mean2stddev(cbar,FractionOfMeanForMinimum,FractionOfMean)
%function stddev=mean2stddev(cbar,FractionOfMeanForMinimum,FractionOfMean)
%cbar(1:Nc) - mean of coefficients for one model
%FractionOfMeanForMinimum - fraction of largest element of mean to use as minimum standard deviation for perturbation from mean
%FractionOfMean - fraction of mean to use as standard deviation

if FractionOfMeanForMinimum<0.0
  error('mean2stddev: FractionOfMeanForMinimum %g < 0.0\n',FractionOfMeanForMinimum);
end
if FractionOfMean<0.0
  fprintf(1,'mean2stddev: FractionOfMean %g < 0.0 -- flag for using in calling programs, use abs(FractionOfMean) here\n',FractionOfMean);
end

%Zhye Yin JSB 2003 paper: compute minimum from l=0,n=0,p=1 weight, assume that l=0,n=0,p=1 weight is cbar(1).
%minstddev=FractionOfMeanForMinimum*abs(cbar(1));
%Could search for the l=0,n=0,p=1 element -- have not written code.
%Here I use the largest element of cbar.
minstddev=FractionOfMeanForMinimum*max(abs(cbar));
fprintf('mean2stddev: FractionOfMeanForMinimum is %d, minimum stddev is %d\n',FractionOfMeanForMinimum,minstddev);
stddev=abs(FractionOfMean)*abs(cbar);
%display(stddev);
last=max(stddev);
%display(last);
for nNc=1:length(stddev)

  if stddev(nNc)==0.0
    %fprintf(1,'mean2stddev: stddev(nNc=%d) is zero\n',nNc);
    stddev(nNc)=last;
  end

  if stddev(nNc)<minstddev
      %fprintf(1,'mean2stddev: stddev(nNC=%d) is less than minimum %d\n',nNc,minstddev);
    stddev(nNc)=minstddev;
  end

  last=stddev(nNc);
end
