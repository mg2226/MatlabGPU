function [Qem gQem HQem]=set_Q_dV_unc(rule, vk, y, noise_var, virusobjs, p_theta_eta, nu, im2use, all_tilde_b, eta)
%function [Qem gQem HQem]=set_Q_dV_unc(rule, vk, y, noise_var, virusobjs, p_theta_eta, nu, im2use, all_tilde_b, eta)
% function to compute Q, its gradient and Hessian for fimunc

% read tilde_b in only once

% Qiu Wang
% 06/06/2010
%
%modified by pd83 on 08-07-2015:
%(1) pass tilde_b to set_Q_dV_unc rather than read tilde_b within set_Q_dV_unc.
%(2) make set_Q_dV_unc work for multiclass problems by passing the value of eta to set_Q_dV_unc rather than assuming that eta=1.

t_Q=clock;

WEIGHT_INDEX = 6; % the index of the abscissa weights in the rule data structure

vo = virusobjs; % use a shorter name inside the function

%% initialize variables
%Neta = size(vo,1);
Ny = size(y, 1);
Nv = size(y, 2);
Nzi = size(rule, 1);
Nc = size(vo{1}.cbar,1);


%% compute the integral $$Q$$ Eq. (178) its gradient and Hessian
Qem = 0;
gQem=zeros(Nc,1);
HQem=zeros(Nc,Nc);

% compute the expectation with numerical integrations
%gQem_eta=0;
%detV0=prod(vo{eta}.nu);
lndetV=sum(log(nu));

parfor (n=1:Nzi)
%     maxNumCompThreads(4);
    Rabc = euler2R(rule(n, 1:3));
    L=setL_nord(Rabc, vo{eta}.clnp.il, vo{eta}.clnp.in, vk, vo{eta}.Htable, vo{eta}.map_unique2lp, all_tilde_b{eta}.tilde_b, vo{eta}.Ftable, vo{eta}.BasisFunctionType);
    %     Sigma0=L*diag(vo{eta}.nu)*L';
    %     Sigma=L*diag(nu)*L';
    %     % add noise_var on the diagonal of Sigma; this equals to adding a Q
    %     % matrix that is noise_var * eye
    %     for ii=1:size(Sigma, 1)
    %         Sigma0(ii,ii) = Sigma0(ii,ii) + noise_var;
    %         Sigma(ii,ii) = Sigma(ii,ii) + noise_var;
    %     end
    %     lndet_Sigma = lndet(Sigma);

    % Use fast lndet algorithm
    [lndet_Sigma,Mid]=lndet_fast(noise_var,L,nu,lndetV);

    Lc = L*vo{eta}.cbar; % \mu = L*c;
    y_Lc = bsxfun(@minus, y(:,im2use{eta,n}), Lc);
    wri = rule(n,WEIGHT_INDEX)*p_theta_eta{eta}(im2use{eta,n},n);
    sum_wri=sum(wri);
    %D=L'*(Sigma0\y_Lc);
    [D,M]= LT_SigmaInv_fast(noise_var,L,vo{eta}.nu,y_Lc);
    repwri=repmat(sqrt(wri'),Nc,1);
    D=D.*repwri;
    D=D*D';

    y_Lc_Sigma = y_Lc_Sigma_fast(noise_var,L,nu,y_Lc,Mid);

    %     % Use fast inversion algorithm
    %     Sigma0_inv=SigmaInv_fast(noise_var,L,vo{eta}.nu);
    %     Sigma_inv=SigmaInv_fast(noise_var,L,nu);
    %
    %     %LT_Sigma0_inv_L = L'*(Sigma0\L);  % L'*Sigma_inv_L; compute \Xi
    %     LT_Sigma0_inv_L = L'*Sigma0_inv*L;  % L'*Sigma_inv_L; compute \Xi
    %     M=LT_Sigma0_inv_L;
    %     Lc = L*vo{eta}.cbar; % \mu = L*c;
    %     y_Lc = bsxfun(@minus, y, Lc);
    %     wri = rule(n,WEIGHT_INDEX)*p_theta_eta{eta}(:,n);
    %     sum_wri=sum(wri);
    %     %D=L'*(Sigma0\y_Lc);
    %     D=L'*(Sigma0_inv*y_Lc);
    %     repwri=repmat(sqrt(wri'),Nc,1);
    %     D=D.*repwri;
    %     D=D*D';
    %
    %     %y_Lc_Sigma = sum( y_Lc .* (Sigma\y_Lc) );
    %     y_Lc_Sigma = sum( y_Lc .* (Sigma_inv*y_Lc) );
    Qem = Qem + (-0.5 * (lndet_Sigma+y_Lc_Sigma)*wri) ;
    %Qem = Qem + (-0.5 * ((lndet_Sigma+y_Lc_Sigma)*p_theta_eta{eta}(:,n)) * rule(n, WEIGHT_INDEX));
    gQem = gQem + diag(D - M*sum_wri);
    HQem = HQem + (-2*D + M*sum_wri).*M;

end
%end
% Qem = Qem - 0.5*Ny*log(2*pi)*Nv;

%fprintf(1,'gQem size = %d\n',size(gQem));
%fprintf(1,'HQem size = %d\n',size(HQem));
Qem=-Qem;
gQem=-gQem;
HQem=-HQem;


fprintf(1,'Q = %.16g\n',-Qem);
fprintf(1,'set_Q_dV_unc time: %d\n',etime(clock, t_Q));

