function [vobjNEW,pixelnoisevarNEW]=homo2hetero_estVandpixelnoisevar(vk,y,vobj,pixelnoisevar,EM_iter,all_tilde_b,absthres,homo2hetero_estpixelnoisevar,homo2hetero_estV)
%function [vobjNEW,pixelnoisevarNEW]=homo2hetero_estVandpixelnoisevar(vk,y,vobj,pixelnoisevar,EM_iter,all_tilde_b,absthres,homo2hetero_estpixelnoisevar,homo2hetero_estV)

fprintf(1,'homo2hetero_estVandpixelnoisevar: at entry\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set default return values.
vobjNEW=[];
pixelnoisevarNEW=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if homo2hetero_estpixelnoisevar %estimate pixelnoisevar ZhengWangDoerschuk JOSA-A 2012 Q

  fprintf(1,'homo2hetero_estVandpixelnoisevar: will update pixelnoisevarNEW\n');

  %return to the homogeneous model
  vobj_homo=vobj;
  for eta=1:length(vobj_homo)
    vobj_homo{eta}.nu=zeros(size(vobj{eta}.nu));
  end

  %estimate pixelnoisevar Q using the homogeneous model vobj_homo
  %Compute $p(\theta_i,\eta_i|y_i,\bar{c}, V, q)$ by ZhengWangDoerschuk JOSA-A 2012 Eq. 22.
  [p_theta_eta,loglike2]=set_p_theta_eta_loglike_diagV(EM_iter.rule,vk,y,pixelnoisevar,vobj_homo,all_tilde_b);
  im2use=Image2use(p_theta_eta,EM_iter.rule,absthres);
  for eta=1:length(vobj)
    if ~isfinite(p_theta_eta{eta})
      error('homo2hetero_estVandpixelnoisevar: ~isfinite(p_theta_eta{eta=%d}) is true\n',eta);
    end
  end
  %ZhengWangDoerschuk JOSA-A 2012 Section 4.B.3 p. 964, using the homogeneous model vobj_homo
  pixelnoisevarNEW=set_Q_noise_solve(EM_iter.rule,vk,y,vobj_homo,p_theta_eta,EM_iter.Na,im2use,all_tilde_b);
  fprintf(1,'homo2hetero_estVandpixelnoisevar: pixelnoisevarNEW %g\n',pixelnoisevarNEW);
  if ~isfinite(pixelnoisevarNEW)
    error('homo2hetero_estVandpixelnoisevar: ~isfinite(pixelnoisevarNEW) is true\n');
  end

  %done with the homogeneous model
  clear vobj_homo

end %if homo2hetero_estpixelnoisevar %estimate pixelnoisevar \nu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if homo2hetero_estV %estimate nu ZhengWangDoerschuk JOSA-A 2012 \nu = diag(V_\eta)

  fprintf(1,'homo2hetero_estVandpixelnoisevar: will update vobjNEW\n');

  vobjNEW=vobj; %will change vobj{eta}.nu below

  %an initial condition on nu should already exist
  for eta=1:length(vobj)
    if any(vobj{eta}.nu<=0.0)
      fprintf(1,'homo2hetero_estVandpixelnoisevar: nu initial condition not set for eta %d\n',eta);
      error('homo2hetero_estVandpixelnoisevar: nu initial condition not set\n');
    end
  end

  pixelnoisevar2use=pixelnoisevar;
  if pixelnoisevarNEW~=[] %this means that pixelnoisevarNEW was computed above
    pixelnoisevar2use=pixelnoisevarNEW;
  end

  nuNEW=cell(length(vobj),1);

  for eta=1:length(vobj)

    %Compute $p(\theta_i,\eta_i|y_i,\bar{c}, V, q)$ by ZhengWangDoerschuk JOSA-A 2012 Eq. 22 or equivalently
    %p223.yzheng.qiuwang.heterogeneous.GaussianMixture.tex 2011-07-31 Eq. 25.
    %Notice the use of pixelnoisevar2use.
    %Notice the use of vobj, which contains the initial condition for nu
    [p_theta_eta,loglike2]=set_p_theta_eta_loglike_diagV(EM_iter.rule,vk,y,pixelnoisevar2use,vobj,all_tilde_b);
    im2use=Image2use(p_theta_eta,EM_iter.rule,absthres);

    %Set lower and upper bounds for the nu optimization problem.
    lb=zeros(length(vobj{eta}.cbar),1);
    ub=ones(length(vobj{eta}.cbar),1)*Inf;

    %Update nu.
    %ZhengWangDoerschuk JOSA-A 2012 from Eq. 35 using Eqs. 43 and 44.
    %Notice the use of pixelnoisevar2use.
    %Notice the use of vobj, which contains the initial condition for nu
    nuNEW{eta}=fmincon(@(nu) set_Q_dV_unc(EM_iter.rule,vk,y,pixelnoisevar2use,vobj,p_theta_eta,nu,im2use,all_tilde_b,eta), ...
        vobj{eta}.nu,[],[],[],[],lb,ub,[],optimset('Algorithm','trust-region-reflective','GradObj','on', ...
        'Hessian','on','TolX',EM_iter.V_TolX,'MaxIter',EM_iter.V_MaxIter,'Display','iter'));
    disp_v(EM_iter.verbosity,2,nuNEW{eta}');
    vobjNEW{eta}.nu=nuNEW{eta};

  end %close 'for eta=1:length(vobj)' loop

end %if homo2hetero_estV %estimate nu ZhengWangDoerschuk JOSA-A 2012 \nu = diag(V_\eta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
