function [vobj,pixelnoisevar,loglikefinal,toomanyiterations]=EM_expmax_iterate2convergence(vk,y,vobj,pixelnoisevar,EM_iter,all_tilde_b,absthres)
%function [vobj,pixelnoisevar,loglikefinal,toomanyiterations]=EM_expmax_iterate2convergence(vk,y,vobj,pixelnoisevar,EM_iter,all_tilde_b)

%The "old" values of cbar, nu, and q are in vobj.
%The "new" values of cbar, nu, and q are in cbarNEW, nuNEW, and qNEW.
%The "old" value of the pixel noise variance Q is pixelnoisevar.
%The "new" value of the pixel noise variance Q is pixelnoisevarNEW.
%cbar occurs twice in vobj: vobj{eta}.clnp.c and vobj{eta}.cbar and
%these are kept equal to one another.
%The "state" is really cbar, nu, and q in vobj and pixelnoisevar.

%Issues:
%(1) Make sure that the old values of the parameters being estimated are always in vobj
%and are appropriately passed to the various functions.

%Future:
%(1) Want to use the same software for homogeneous and for heterogeneous problems.
%(a) Save computation in the homogeneous problem by not updating nu.  Item (a)
%is fully implemented in this version of the code.
%(b) Should be able to save additional computation in the homogeneous problem by
%avoiding computation rather than just having nu=0 in the heterogeneous
%formulas.  Item (b) is not fully implemented.
%Furthermore, it makes sense to have some classes be heterogeneous while other
%classes are homogeneous.  Nothing has been done for a mixed
%homogeneous-heterogeneous problem.  The sign that a model is homogeneous is the
%vobj{eta}.nu is empty.  Having vobj{eta}.nu empty means that the heterogeneous
%formulas will fail.  So, save the true situation in is_hetero(1:Neta) and
%is_all_hetero and then set vobj{eta}.nu equal to a vector of 0.  Replace the
%vector with an empty matrix [] before returning from this function.
%(2) Evaluate the log likelihood for heterogeneous problems.  Currently, it is
%evaluated only for homogeneous problems.  The log likelihood is essential for
%comparing different initial conditions in the Monte Carlo process.  I am
%willing either to return ZhengWangDoerschuk JOSA-A 2012 Eq. 19 (the true log
%likelihood) or Eq. 21 (the conditional log likelihood of the expectation
%maximization algorithm) since either should be a fair comparison of different
%answers.  The old C and C/C++/OpenMP/MPI codes for the homogeneous problem
%could do this and I believe the current code reproduces that calculation.  For
%a heterogeneous problem, all three terms in Eq. 21 should be present in the
%returned value.
 
p_theta_eta=[]; %2015-07-13 pd83: p_theta_eta was not found at the final save

fprintf(1,'EM_expmax_iterate2convergence: size(vk) %d %d size(y) %d %d\n',size(vk),size(y));

Neta=length(vobj);

cbarNEW=cell(Neta,1);
nuNEW=cell(Neta,1);
qNEW=zeros(Neta,1);

cbarTMP=cell(Neta,1); %Software workaround so that is_clnp_converged(.) can be used.

%Determine if the problem is all heterogeneous versus all homogeneous versus mixed.
is_hetero=false(Neta,1);
for eta=1:Neta
  if isempty(vobj{eta}.nu)
    is_hetero(eta)=false;
  else
    is_hetero(eta)=true;
  end
end
if all(is_hetero==true)
  is_all_hetero=true;
  fprintf(1,'EM_expmax_iterate2convergence: all classes are heterogeneous\n');
elseif all(is_hetero==false)
  is_all_hetero=false;
  fprintf(1,'EM_expmax_iterate2convergence: all classes are homogeneous\n');
else
  error('EM_expmax_iterate2convergence: mixed heterogeneous and homogeneous problem\n');
end

%If the problem is all homogeneous, set all of the empty nu's to zero vectors.
if ~is_all_hetero
  fprintf(1,'EM_expmax_iterate2convergence: problem has only homogeneous classes, setting vobj{eta}.nu to 0 vector\n');
  for eta=1:Neta
    vobj{eta}.nu=zeros(size(vobj{eta}.cbar));
  end
end

%Print the initial conditions.
fprintf(1,'EM_expmax_iterate2convergence: Initial pixelnoisevar: %g\n',pixelnoisevar);
for eta=1:Neta
  fprintf(1,'EM_expmax_iterate2convergence: Initial cbar{%d} from vobj:\n',eta);
  disp_v(EM_iter.verbosity,1,vobj{eta}.cbar')
  fprintf(1,'EM_expmax_iterate2convergence: Initial nu{%d} from vobj:\n',eta);
  disp_v(EM_iter.verbosity,1,vobj{eta}.nu');
  fprintf(1,'EM_expmax_iterate2convergence: Initial q(%d) from vobj: %g\n',eta,vobj{eta}.q);
end

%Allocate and initialize history variables.  Except for loglike_history, the history
%variables preserve the state (cbar, nu, and q in vobj and pixelnoisevar) at the start
%of each iteration (so the values that were used to call this function are preserved)
%with the final entries being the state at the close of the terminal iteration.
cbar_history=cell(Neta,EM_iter.maxiter+1);
nu_history=cell(Neta,EM_iter.maxiter+1);
q_history=zeros(Neta,EM_iter.maxiter+1);
pixelnoisevar_history=zeros(EM_iter.maxiter+1,1);
loglike_history=-inf(EM_iter.maxiter+1,1); %set to -Inf

convergedtwice=false; %require convergence to be achieved for two iterations in a row
toomanyiterations=false; %exited the loop not because convergence was achieved but because too many iterations were used.
if is_all_hetero %all classes are heterogenous models.
  %Previously assumed that the nu initial condition was not set and therefore executed havesetnuic=false
  %Now test if the initial condition on nu is already set <==> nu>0.0 for all classes and all elements.
  %If desire to set nu initial conditions here, can always read nu=0.0 from a file or set nu=0.0 when changing from a homogeneous to a heterogeneous model.
  havesetnuic=true;
  for eta=1:Neta
    if any(vobj{eta}.nu<=0.0)
      havesetnuic=false;
      break;
    end
  end
else %all classes are homogenous models.
  havesetnuic=false;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main loop over iterations starts here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter=1;
% profile -history
while true
  
  fprintf(1,'EM_expmax_iterate2convergence: Iteration %d\n',iter);

  if get(0,'Diary')=='on' %test if the diary is on and if yes then flush to disk
    diary off; diary on
  end
    
  %Copy the current value of the state into the history variables.
  for eta=1:Neta
    cbar_history{eta,iter}=vobj{eta}.cbar;
    nu_history{eta,iter}=vobj{eta}.nu;
    q_history(eta,iter)=vobj{eta}.q;
  end
  pixelnoisevar_history(iter)=pixelnoisevar;
  if iter==1
    loglike_history(iter)=-inf;
  else
    loglike_history(iter)=loglikeNEW
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Allow special actions if this is the first attempt to compute a heterogeneous model from a homogeneous model initial condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%new members of EM_iter which control the special actions
%EM_iter.homo2hetero_1st_V_MaxIter: integer >= 0
%EM_iter.homo2hetero_estpixelnoisevar: used exclusively in "if ..." test
%EM_iter.homo2hetero_estV: used exclusively in "if ..." test

  if iter==1
    V_MaxIterSave=EM_iter.V_MaxIter; %save original value
    EM_iter.V_MaxIter=EM_iter.homo2hetero_1st_V_MaxIter; %change to special value
    [vobjNEW,pixelnoisevarNEW]=homo2hetero_estVandpixelnoisevar(vk,y,vobj,pixelnoisevar,EM_iter,all_tilde_b,absthres,EM_iter.homo2hetero_estpixelnoisevar,EM_iter.homo2hetero_estV);
    EM_iter.V_MaxIter=V_MaxIterSave; %restore original value
    clear V_MaxIterSave; %delete temporary variable
    if ~isempty(vobjNEW)

      for etap=1:length(vobj)
        fprintf(1,'EM_expmax_iterate2convergence: homo2hetero calculation: vobj{%d}.cbar, sqrt(vobj{%d}.nu), and sqrt(vobjNEW{%d}.nu):\n',etap,etap,etap);
        for lnp=1:length(vobj{etap}.nu)
          fprintf(1,'%g %g %g\n',vobj{etap}.cbar(lnp),sqrt(vobj{etap}.nu(lnp)),sqrt(vobjNEW{etap}.nu(lnp)));
        end
      end

      vobj=vobjNEW;
    end
    if ~isempty(pixelnoisevarNEW)

      fprintf(1,'EM_expmax_iterate2convergence: homo2hetero calculation: pixelnoisevar %g pixelnoisevarNEW %g\n',pixelnoisevar,pixelnoisevarNEW);

      pixelnoisevar=pixelnoisevarNEW;
    end
    loglikeNEW=-inf; %not computed
    %treat these special actions as an iteration
    iter=iter+1;
    continue;
  end %if iter==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Beginning of joint update of pixelnoisevar, q, and cbar.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %Want to do a joint update of pixelnoisevar, q, and cbar.  So do not want to copy the new
  %value of pixelnoisevar (which is contained in pixelnoisevarNEW) into pixelnoisevar until
  %the new values of q and cbar are computed (which are contained in qNEW and cbarNEW).

  pixelnoisevarNEW=pixelnoisevar;
  if iter<=EM_iter.maxiter4pixelnoisevarupdate && ...
    (is_all_hetero || (~is_all_hetero && EM_iter.estimate_noise_var_in_homogeneous_problem))
    %Update pixelnoisevar.

    fprintf(1,'EM_expmax_iterate2convergence: iter %d compute p(nuisance|data,parameters) for update of pixelnoisevar\n',iter);
    tstart=tic;
    %Compute $p(\theta_i,\eta_i|y_i,\bar{c}, V, q)$ by ZhengWangDoerschuk JOSA-A 2012 Eq. 22.
    [p_theta_eta,loglike2]=set_p_theta_eta_loglike_diagV(EM_iter.rule,vk,y,pixelnoisevar,vobj,all_tilde_b);
    im2use=Image2use(p_theta_eta,EM_iter.rule,absthres);
%     save('im2use.mat','im2use');
%     fprintf('saved one im2use!!!');
%     im2use=cell(1,5000);
%     vect=ones(1,1200);
%     for j=1:1200
%         vect(j)=j;
%     end
%     for i=1:5000
%         im2use{i}=vect;
%     end
    fprintf(1,'EM_expmax_iterate2convergence: iter %d (%g) done p(nuisance|data,parameters) for update of pixelnoisevar\n',iter,toc(tstart));

    fprintf(1,'EM_expmax_iterate2convergence: iter %d compute pixelnoisevarNEW (pixelnoisevar %g)\n',iter,pixelnoisevar);
    for eta=1:Neta
      if ~isfinite(p_theta_eta{eta})
        error('EM_expmax_iterate2convergence: ~isfinite(p_theta_eta{eta=%d}) is true\n',eta);
      end
    end
    tstart=tic;
    %ZhengWangDoerschuk JOSA-A 2012 Section 4.B.3 p. 964.
    pixelnoisevarNEW=set_Q_noise_solve(EM_iter.rule,vk,y,vobj,p_theta_eta,EM_iter.Na,im2use,all_tilde_b);
    fprintf(1,'EM_expmax_iterate2convergence: iter %d (%g) done pixelnoisevarNEW %g\n',iter,toc(tstart),pixelnoisevarNEW);
    if ~isfinite(pixelnoisevarNEW)
      error('EM_expmax_iterate2convergence: ~isfinite(pixelnoisevarNEW) is true\n');
    end
  else
    fprintf(1,'EM_expmax_iterate2convergence: iter %d DID NOT compute p(nuisance|data,parameters) for update of pixelnoisevar\n',iter);

  end

  %Update cbar and q.
  fprintf(1,'EM_expmax_iterate2convergence: iter %d update cbar and q\n',iter)
  tstart=tic;
  %ZhengWangDoerschuk JOSA-A 2012 Eqs. 26 and 33.
  %Updating cbar and q require $p(\theta_i,\eta_i|y_i,\bar{c}, V, q)$ and that is computed
  %within EMstep_cbar_dV_large (by set_p_theta_eta_loglike_diagV) using the new value of the
  %pixel noise covariance.
  fprintf(1,'EM_expmax_iterate2convergence: iter %d compute new cbar and q and loglike via EMstep_cbar_dV_large\n',iter)
  %profile on -history
  [cbarNEW,qNEW,loglikeNEW,solvedlinearsystem]=EMstep_cbar_dV_large(EM_iter.rule,vk,y,pixelnoisevar,vobj,all_tilde_b,EM_iter.MinimumClassProb,absthres);
%   stats=profile('info');
% profsave(stats,'profile_3');
% profile off
% fprintf('TESTTTTTTTTTTTT STOP!');
  fprintf(1,'EM_expmax_iterate2convergence: iter %d (%g) done new cbar and q and loglike via EMstep_cbar_dV_large\n',iter,toc(tstart));
  fprintf(1,'EM_expmax_iterate2convergence: iter %d done new loglikelihood %.16g\n',iter,loglikeNEW);
  if all(solvedlinearsystem==0)
    fprintf(1,'EM_expmax_iterate2convergence: EMstep_cbar_dV_large: did not solve F\g for any class\n');
  end
  for eta=1:Neta
    fprintf(1,'EM_expmax_iterate2convergence: iter %d done cbarNEW{%d}:\n',iter,eta);
    disp_v(EM_iter.verbosity,2,cbarNEW{eta}');
    if any(~isfinite(cbarNEW{eta}))
      error('EM_expmax_iterate2convergence: cbarNEW{eta=%d} fails isfinite\n',eta);
    end
  end
  fprintf(1,'EM_expmax_iterate2convergence: iter %d done qNEW for all classes:\n',iter);
  disp_v(EM_iter.verbosity,1,qNEW');

  %Copy the old cbar from vobj so that is_clnp_converged(.) can be used.
  for eta=1:Neta
    cbarTMP{eta}=vobj{eta}.cbar;
  end

  %Test if the loglikelihood has converged.
  absoluteloglike=iter>1 && ...
    (EM_iter.loglikeftol<0.0 || abs(loglikeNEW-loglike_history(iter))<=EM_iter.loglikeftol);
  relativeloglike=iter>1 && ...
    (EM_iter.loglikertol<0.0 || abs(loglikeNEW-loglike_history(iter))/(0.5*(abs(loglikeNEW)+abs(loglike_history(iter))))<=loglikertol);
  loglike_converge=absoluteloglike && relativeloglike;

  %Test if cbar has converged.
  [relativecbar,absolutecbar]=is_clnp_converged(cbarNEW,cbarTMP,EM_iter.cbarftol,EM_iter.cbarrtol,EM_iter.cbarftol_dividebyNc);
  fprintf(1,'EM_expmax_iterate2convergence: iter %d relativecbar %g absolutecbar %g\n',iter,relativecbar,absolutecbar);
  cbar_converge=(EM_iter.cbarftol<0.0 || absolutecbar) && (EM_iter.cbarrtol<0.0 || relativecbar);

  %Update part of the state=(virusobj,pixelnoisevar).
  pixelnoisevar=pixelnoisevarNEW; %This may do nothing if, in fact, pixelnoisevar was not updated.
  for eta=1:Neta
    vobj{eta}.clnp.c=cbarNEW{eta};
    vobj{eta}.cbar=cbarNEW{eta};
    vobj{eta}.q=qNEW(eta);
  end

  %Combine the tests, and make sure that convergence has occurred for two updates in a row.
  %Note that these two updates will be separated by an update of nu and q.
  %Qiu Wang code uses || not && and uses an extremely strict criteria for loglike that is essentially never satisfied.
  if cbar_converge && loglike_converge
    fprintf(1,'EM_expmax_iterate2convergence: The cbar estimation converged at least once (cbar_converge %d loglike_converge %d).\n', ...
            cbar_converge,loglike_converge);
    if convergedtwice
      fprintf(1,'EM_expmax_iterate2convergence: The cbar estimation has converged a second time in a row.\n');
      loglikefinal=loglikeNEW;
      toomanyiterations=false;
      break; %break out of "while true" loop, only exit for convergence
    else
      fprintf(1,'EM_expmax_iterate2convergence: failed convergedtwice test\n');
      convergedtwice=true; %save the fact that convergence occured once (but not twice)
    end
  else %neither cbar nor the loglikelihood converged.
    fprintf(1,'EM_expmax_iterate2convergence: failed cbar_converge AND loglike_converge test\n');
    convergedtwice=false;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of joint update of pixelnoisevar, q, and cbar.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if convergedtwice && is_all_hetero

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Beginning of joint update of q and nu.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %While updating cbar and q, have converged once on cbar and/or loglikelihood, now update nu and q.
    %Only do this if the problem has only heterogeneous classes.

    %Always start the nu optimization from an initial condition that is
    %EM_iter.nu_ic_FractionOfcbar times the current value of cbar.  The Qiu
    %Wang code sets this initial condition from the current cbar every time
    %nu is optimized.  In this code, that is made optional: either do as in
    %the Qiu Wang code or in the first optimization of nu the initial
    %condition is set from the current cbar but in later optimizations of
    %nu the initial condition is the result from the immediately preceding
    %optimization of nu.

    fprintf(1,'EM_expmax_iterate2convergence: iter %d about to update nu and q\n');

    if ~havesetnuic || ~isempty(EM_iter.nu_ic_always_proportional2cbar)
      fprintf(1,'EM_expmax_iterate2convergence: iter %d about to set nu initial conditional proportional to the current cbar\n');
      for eta=1:Neta
        if ~havesetnuic || (~isempty(EM_iter.nu_ic_always_proportional2cbar) && EM_iter.nu_ic_always_proportional2cbar(eta))
          if EM_iter.nu_ic_FractionOfcbar>=0.0
            vobj{eta}.nu=mean2stddev(vobj{eta}.cbar,EM_iter.nu_ic_FractionOfcbarForMinimum,EM_iter.nu_ic_FractionOfcbar).^2;
          else
            vobj{eta}.nu=mean2stddev(vobj{eta}.cbar,EM_iter.nu_ic_FractionOfcbarForMinimum,EM_iter.nu_ic_FractionOfcbar);
          end
        end
      end
      havesetnuic=true;
    else
      fprintf(1,'EM_expmax_iterate2convergence: iter %d nu initial conditional not changed\n');
    end

    for eta=1:Neta

      fprintf(1,'EM_expmax_iterate2convergence: iter %d compute p(nuisance|data,parameters) for update of nu and q\n',iter);
      tstart=tic;
      %Compute $p(\theta_i,\eta_i|y_i,\bar{c}, V, q)$ by ZhengWangDoerschuk JOSA-A 2012 Eq. 22 or equivalently
      %p223.yzheng.qiuwang.heterogeneous.GaussianMixture.tex 2011-07-31 Eq. 25.
      %Do this inside of the eta loop because vobj is updated at the bottom -- eta=1 uses all older information, eta=2 uses new information for Class 1 and old for Class 2 and higher classes, etc
      [p_theta_eta,loglike2]=set_p_theta_eta_loglike_diagV(EM_iter.rule,vk,y,pixelnoisevar,vobj,all_tilde_b);
      im2use=Image2use(p_theta_eta,EM_iter.rule,absthres);
      fprintf(1,'EM_expmax_iterate2convergence: iter %d (%g) done p(nuisance|data,parameters) for update of nu and q\n', ...
              iter,toc(tstart));

      %Update q.
      WEIGHT_INDEX=6; %the index of the abscissa weights in the rule data structure
      fprintf(1,'EM_expmax_iterate2convergence: iter %d compute qNEW\n',iter);
      tstart=tic;
      %ZhengWangDoerschuk JOSA-A 2012 Eq. 26.
      Nv=size(y,2);
      qNEW(eta)=sum(p_theta_eta{eta}*EM_iter.rule(:,WEIGHT_INDEX))/Nv;
      fprintf(1,'EM_expmax_iterate2convergence: iter %d (%g) done qNEW(eta=%d) %g:\n',iter,toc(tstart),eta,qNEW(eta));

      %Set lower and upper bounds for the nu optimization problem.
      lb=zeros(length(vobj{eta}.cbar),1);
      ub=ones(length(vobj{eta}.cbar),1)*Inf;

      %Update nu.
      fprintf(1,'EM_expmax_iterate2convergence: iter %d compute nuNEW using fmincon\n',iter);
      tstart=tic;
%       im2use=cell(1,5000);
%       vect=ones(1,1200);
%       for j=1:1200
%           vect(j)=j;
%       end
%       for i=1:5000
%           im2use{i}=vect;
%       end
      %ZhengWangDoerschuk JOSA-A 2012 from Eq. 35 using Eqs. 43 and 44.
      nuNEW{eta}=fmincon(@(nu) set_Q_dV_unc(EM_iter.rule,vk,y,pixelnoisevar,vobj,p_theta_eta,nu,im2use,all_tilde_b,eta), ...
        vobj{eta}.nu,[],[],[],[],lb,ub,[],optimset('Algorithm','trust-region-reflective','GradObj','on', ...
        'Hessian','on','TolX',EM_iter.V_TolX,'MaxIter',EM_iter.V_MaxIter,'Display','iter'));
      fprintf(1,'EM_expmax_iterate2convergence: iter %d (%g) done nuNEW(eta):\n',iter,toc(tstart),eta);
      disp_v(EM_iter.verbosity,2,nuNEW{eta}');

      %Update part of the state=(virusobj,pixelnoisevar).
      vobj{eta}.nu=nuNEW{eta};
      vobj{eta}.q=qNEW(eta);

    end %close 'for eta=1:Neta' loop

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %End of joint update of q and nu.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  else
    fprintf(1,'EM_expmax_iterate2convergence: failed convergedtwice AND is_all_hetero test\n');

  end %close 'if convergedtwice && is_all_hetero' test

  %Test if the maximum number of iterations has been reached.  If yes, break.  If no, increment iter and continue.
  if iter>=EM_iter.maxiter
    toomanyiterations=true;
    break; %break out of "while true" loop, only exit for NON convergence
  end
%   if iter==1
%       save('pthetaeta1.mat','p_theta_eta');
%       stats=profile('info');
%       profsave(stats,'profile_results');
%       profile off
%   end
  iter=iter+1;
  
  if mod(iter,20)==0
    save(['pthetaeta' num2str(iter) '.mat'],'p_theta_eta');
  end

end %close "while true" loop
save('pthetaeta.final.mat','p_theta_eta'); %2015-07-13 pd83: p_theta_eta was not found at the final save
%Something is always done during an iteration of the loop and the state is immediately updated.
%So retain the last group of events by copying the state to the history variables one more time.
%It is possible to leave the loop via two 'break' statements: "convergence" or "too many iterations".
%In both cases, the iter-th value of the history variables has already been set at the beginning of
%the loop and the last group of events should be placed in the (iter+1)-th value.
for eta=1:Neta
  cbar_history{eta,iter+1}=vobj{eta}.cbar;
  nu_history{eta,iter+1}=vobj{eta}.nu;
  q_history(eta,iter+1)=vobj{eta}.q;
end
pixelnoisevar_history(iter+1)=pixelnoisevar;


%set loglike_history(iter+1) to -Inf above.

if toomanyiterations
  fprintf(1,'EM_expmax_iterate2convergence: no convergence in maxiter %d iterations\n',EM_iter.maxiter);
  loglikefinal=loglike_history(iter);
else
  fprintf(1,'EM_expmax_iterate2convergence: failed toomanyiterations test\n');
end

if ~is_all_hetero
  %All classes are homogeneous classes.  Therefore, return to vobj{eta}.nu=[].
  fprintf(1,'EM_expmax_iterate2convergence: problem has only homogeneous classes, returning vobj{eta}.nu to [] vector\n');
  for eta=1:Neta
    if any(vobj{eta}.nu~=0.0)
      error('EM_expmax_iterate2convergence: eta %d homogeneous problem but vobj{eta}.nu ~= 0 vector\n',eta);
    end
    vobj{eta}.nu=[];
  end
else
  fprintf(1,'EM_expmax_iterate2convergence: failed ~is_all_hetero test\n');
end

if ~isempty(EM_iter.fn_savehistory)
  fprintf(1,'EM_expmax_iterate2convergence: about to save history and return\n');
  save(EM_iter.fn_savehistory,'cbar_history','nu_history','q_history','pixelnoisevar_history','loglike_history');
else
  fprintf(1,'EM_expmax_iterate2convergence: about to not save history and return\n');
end
