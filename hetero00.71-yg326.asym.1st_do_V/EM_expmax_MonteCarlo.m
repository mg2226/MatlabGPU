function [vobj_best,pixelnoisevar_best,loglikelihood_best]=EM_expmax_MonteCarlo(vk,y,vobj,pixelnoisevar,EM_MC,EM_iter,all_tilde_b,absthres)
%function [vobj_best,pixelnoisevar_best,loglikelihood_best]=EM_expmax_MonteCarlo(vk,y,vobj,pixelnoisevar,EM_MC,EM_iter,all_tilde_b)

if EM_MC.Nrandic<=0
  error('EM_expmax_MonteCarlo: EM_MC.Nrandic %d < 0\n',EM_MC.Nrandic);
end

fprintf(1,'EM_expmax_MonteCarlo: Nrandic %d\n',EM_MC.Nrandic);

[vobj_best,pixelnoisevar_best,loglikelihood_best,toomanyiterations]=EM_expmax_iterate2convergence(vk,y,vobj,pixelnoisevar,EM_iter,all_tilde_b,absthres);

assert_UpdateBothCopiesOfcbar(vobj);

if EM_MC.Nrandic==1
  return;
end

Neta=length(vobj);

%EM_expmax_iterate2convergence computes the log likelihood only for homogeneous problems.
%So Monte Carlo applied to the initial conditions is only permitted for homogeneous problems.
for eta=1:Neta
  if ~isempty(vobj{eta}.nu)
    error('EM_expmax_MonteCarlo: eta %d The case of multiple initial conditions for heterogeneous problems is not implemented\n',eta);
  end
end

%Nc(eta) can vary with eta so use a cell array, not a 2-D matrix
cbar_ic_mean=cell(Neta,1);
cbar_ic_stddev=cell(Neta,1);
for eta=1:Neta
  cbar_ic_mean{eta}=vobj{eta}.cbar;
  cbar_ic_stddev{eta}=mean2stddev(cbar_ic_mean{eta},EM_MC.FractionOfMeanForMinimum,EM_MC.FractionOfMean);
  fprintf(1,'EM_expmax_MonteCarlo: eta %d cbar_ic_mean(eta):\n',eta);
  disp(cbar_ic_mean{eta}');
  fprintf(1,'EM_expmax_MonteCarlo: eta %d cbar_ic_stddev(eta):\n',eta);
  disp(cbar_ic_stddev{eta}');
end

for nMC=2:EM_MC.Nrandic

  fprintf(1,'EM_expmax_MonteCarlo: nMC %d\n',nMC);

  for eta=1:Neta
    vobj{eta}.cbar=cbar_ic_mean{eta} + cbar_ic_stddev{eta}.*randn(size(cbar_ic_stddev{eta}));
    vobj{eta}.clnp.c=vobj{eta}.cbar;
  end

  [vobj,pixelnoisevar,loglikelihood,toomanyiterations]=EM_expmax_iterate2convergence(vk,y,vobj,pixelnoisevar,EM_iter,all_tilde_b,absthres);

  assert_UpdateBothCopiesOfcbar(vobj);

  if loglikelihood>loglikelihood_best
    fprintf(1,'EM_expmax_MonteCarlo: nMC %d current answer is best so far loglikelihood %.12g loglikelihood_best %.12g\n', ...
            nMC,loglikelihood,loglikelihood_best);
    loglikelihood_best=loglikelihood;
    vobj_best=vobj;
    pixelnoisevar_best=pixelnoisevar;
  else
    fprintf(1,'EM_expmax_MonteCarlo: nMC %d current answer is not best so far loglikelihood %g loglikelihood_best %g\n', ...
            nMC,loglikelihood,loglikelihood_best);
  end

end %close for nMC=2:EM_MC.Nrandic



