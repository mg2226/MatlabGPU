function noise_var=set_Q_noise_solve(rule, vk, y,virusobjs, p_theta_eta, Na, im2use, all_tilde_b)
%function noise_var=set_Q_noise_solve(rule, vk, y,virusobjs, p_theta_eta, Na, im2use, all_tilde_b)
%OLD, LOOKS LIKE ML: function [Qem gQem HQem]=set_Q_noise(rule, vk, y, noise_var, virusobjs, p_theta_eta)
% function to compute Q as a function of noise variance

% read tilde_b in only once

% Qiu Wang
% 06/06/2010
%
%modified by pd83 08-07-2015:
%(1) pass all_tilde_b to set_Q_noise_solve rather than read all_tilde_b in set_Q_noise_solve.
%(2) added loop over eta to make the algorithm appropriate for multi-class problems.

t_Q=clock;

WEIGHT_INDEX = 6; % the index of the abscissa weights in the rule data structure

vo = virusobjs; % use a shorter name inside the function

%% initialize variables
Neta = size(vo,1);
Ny = size(y, 1);
Nv = size(y, 2);
Nzi = size(rule, 1);
Nc = size(vo{1}.cbar,1);


%% compute the integral $$Q$$ Eq. (178) its gradient and Hessian
Qem = 0;
A=0;
B=0;

% compute the expectation with numerical integrations
for eta=1:Neta
%gQem_eta=0;
%detV0=prod(vo{eta}.nu);
%nu=vo{eta}.nu;
%unneeded since the lndet_fast function is commented out detV=prod(nu);

tilde_b=all_tilde_b{eta}.tilde_b;

parfor (n=1:Nzi)
%     maxNumCompThreads(4);
    Rabc = euler2R(rule(n, 1:3));
    L=setL_nord(Rabc, vo{eta}.clnp.il, vo{eta}.clnp.in, vk, vo{eta}.Htable, vo{eta}.map_unique2lp, tilde_b, vo{eta}.Ftable, vo{eta}.BasisFunctionType);
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
%     lndet_Sigma=lndet_fast(noise_var,L,nu,detV);

    Lc = L*vo{eta}.cbar; % \mu = L*c;
%     vect=zeros(1,1200);
%     for i=1:1200
%         vect(i)=i;
%     end
    y_Lc = bsxfun(@minus, y(:,im2use{eta,n}), Lc);
    %y_Lc = bsxfun(@minus, y(:,vect), Lc);
    wri = rule(n,WEIGHT_INDEX)*p_theta_eta{eta}(im2use{eta,n},n);
%     sum_wri=sum(wri);
    %D=L'*(Sigma0\y_Lc);
%     [D M]= LT_SigmaInv_fast(noise_var,L,vo{eta}.nu,y_Lc);
%     repwri=repmat(sqrt(wri'),Nc,1);
%     D=D.*repwri;
%     D=D*D';

%     y_Lc_Sigma = y_Lc_Sigma_fast(noise_var,L,nu,y_Lc);

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
%     Qem = Qem + (-0.5 * (lndet_Sigma+y_Lc_Sigma)*wri) ;
    %Qem = Qem + (-0.5 * ((lndet_Sigma+y_Lc_Sigma)*p_theta_eta{eta}(:,n)) * rule(n, WEIGHT_INDEX));
%     gQem = gQem + diag(D - M*sum_wri);
%     HQem = HQem + (-2*D + M*sum_wri).*M;

    A=A+sum(wri);
    B=B+ sum(y_Lc.^2)*wri;

end %parfor (n=1:Nzi)

A=A*Na*Na;
noise_var=B/A;
% noise_var=noise_var/8;

end %for eta=1:Neta
% Qem = Qem - 0.5*Ny*log(2*pi)*Nv;

%fprintf(1,'gQem size = %d\n',size(gQem));
%fprintf(1,'HQem size = %d\n',size(HQem));
% Qem=-Qem;
% gQem=-gQem;
% HQem=-HQem;


% fprintf(1,'Q = %.16g\n',-Qem);
fprintf(1,'set_Q_noise_solve time: %d\n',etime(clock, t_Q));

