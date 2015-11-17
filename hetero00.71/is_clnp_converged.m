function [rconverged,fconverged]=is_clnp_converged(cbarnew,cbar,ftol,rtol,ftol_dividebyNc)
%function [rconverged,fconverged]=is_clnp_converged(cbarnew,cbar,ftol,rtol,ftol_dividebyNc)

%This function applies the convergence criteria to each class separately and requires that all classes converge.
%Alternatively, the cbar values for all classes can be combined into one vector and the convergence criteria can be applied to that vector.

rconverged=true;
fconverged=true;

for eta=1:length(cbar)
  errorsignal=norm(cbarnew{eta}-cbar{eta},1);
  averagesignal=0.5*(norm(cbarnew{eta},1)+norm(cbar{eta},1));
  if ftol_dividebyNc
    %Has no effect on the rtol criteria.
    errorsignal=errorsignal./length(cbar{eta});
    averagesignal=averagesignal./length(cbar{eta});
  end
  fconverged=(fconverged && (errorsignal<=ftol));
  rconverged=(rconverged && ((errorsignal/averagesignal)<=rtol));
end
