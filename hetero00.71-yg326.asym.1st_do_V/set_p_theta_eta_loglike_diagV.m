function [p, loglike]=set_p_theta_eta_loglike_diagV(rule, vk, y, noise_var, virusobjs, all_tilde_b)
%function [p, loglike]=set_p_theta_eta_loglike_diagV(rule, vk, y, noise_var, virusobjs, all_tilde_b)
% compute $$p(\theta_i,\eta_i|y_i,\bar{c}, V, q)$$ by Eq. (25) and the log
% likelihood $$p(y|\bar{c}, V, q)$$ by Eq. (320)
%
% Input parameters:
%   rule: integration rule object
%   Lall: precomputed data structures, please set_Lall()
%   y: image data
%   Qnoise: covariance matrix of the noise
%   cbar, V, q: current estimation of mean $$\bar{c}$$, variance $$V$ and
%   prior class probability $$q$, respectively.
%
% Output parameters:
%  p: $$p(\theta_i,\eta_i|y_i,\bar{c}, V, q)$$ from Eq. (25)
%
% This is the diagonal V version of set_p_theta_eta
% Use fast algorithm to compute lndet(Sigma) and Sigma_inverse
% read tilde_b once
%
% Qiu Wang
% 05/22/2010
%
%modified by pd83 08-07-2015 so that tilde_b is passed not read

fprintf(1, 'set_p_theta_eta_loglike_diagV: entering\n');

t_p=clock;

WEIGHT_INDEX = 6; % the index of the abscissa weights in the rule data structure

vo = virusobjs; % use a shorter name inside the function

%% initialize variables
Neta = size(vo, 1);
%Ny = size(y, 1);
Ni = size(y, 2);
Nzi = size(rule, 1);
exponents = cell(Neta, 1);


%% compute the scaling factor of the exponents
%fprintf(1,'computing the scaling factor of the exponents...\n');
for eta=1:Neta
    Nc=length(vo{eta}.nu);
    exponents_eta = zeros(Ni, Nzi);
    %for n=1:Nzi
    lndetV=sum(log(vo{eta}.nu));

    % read tilde b
    lmax = max(vo{eta}.clnp.il);
    tilde_b = all_tilde_b{eta}.tilde_b;
    parfor (n=1:Nzi)
%     parfor (n=1:1)
        %       if mod(n,1000)==0
        %           fprintf(1,'n = %d\n',n);
        %       end
        %maxNumCompThreads(4);
        Rabc = euler2R(rule(n,1:3));
        L =setL_nord(Rabc, vo{eta}.clnp.il, vo{eta}.clnp.in, vk, vo{eta}.Htable, vo{eta}.map_unique2lp, tilde_b, vo{eta}.Ftable, vo{eta}.BasisFunctionType);
        
        %fprintf(1,'elapsed time %d\n',etime(clock,t_p));
        % replace lndet by lndet_fast
        [lndet_Sigma,M]=lndet_fast(noise_var,L,vo{eta}.nu,lndetV);
        %fprintf(1,'elapsed time %d\n',etime(clock,t_p));

        Lc = L * vo{eta}.cbar;
        y_Lc = bsxfun(@minus, y, Lc);
        y_Lc_Sigma=y_Lc_Sigma_fast(noise_var,L,vo{eta}.nu,y_Lc,M); % use fast algorithm
        assert(all(isfinite(y_Lc_Sigma)),'set_p_theta_eta_loglike_diagV: assert(all(isfinite(y_Lc_Sigma))): eta %d n %d\n',eta,n)
        %fprintf(1,'elapsed time %d\n',etime(clock,t_p));
        % use fast inverse algorithm
        %     [SigmaInv nuzero midmx]=SigmaInv_fast(noise_var,L,vo{eta}.nu);
        %
        % %     Sigma=L*diag(vo{eta}.nu)*L';
        % %     % add noise_var on the diagonal of Sigma; this equals to adding a Q
        % %     % matrix that is noise_var * eye
        % %     for ii=1:size(Sigma, 1);
        % %       Sigma(ii,ii) = Sigma(ii,ii) + noise_var;
        % %     end
        % %     lndet_Sigma = lndet(Sigma);
        %     %fprintf(1,'elapsed time %d\n',etime(clock,t_p));
        %     Lc = L * vo{eta}.cbar;
        %     y_Lc = bsxfun(@minus, y, Lc);
        %     fprintf(1,'elapsed time %d\n',etime(clock,t_p));
        %     %y_Lc_Sigma = sum( y_Lc .* (Sigma\y_Lc) );
        %     y_Lc_sq=1/noise_var*sum(y_Lc.*y_Lc);
        %     if nuzero
        %         y_Lc_Sigma1=y_Lc_sq;
        %     else
        %         LTy=L'*y_Lc;
        %         y_Lc_Sigma1=y_Lc_sq-noise_var^(-2)*sum(LTy.*(midmx*LTy));toc;
        %     end
        %     fprintf(1,'elapsed time %d\n',etime(clock,t_p));
        %     y_Lc_Sigma = sum( y_Lc .* (SigmaInv*y_Lc) );
        %     fprintf(1,'elapsed time %d\n',etime(clock,t_p));
        %fprintf(1,'size of y_Lc_Sigma=%d\n',size(y_Lc_Sigma));
        %fprintf(1,'size of lndet=%d\n',size(lndet_Sigma));
        exponents_eta(:,n) = -0.5*(y_Lc_Sigma(:)+lndet_Sigma);
        %fprintf(1,'elapsed time %d\n',etime(clock,t_p));
    end

    exponents{eta} = exponents_eta;
    if (isreal(exponents{eta}) == 0)
        error('The exponents matrix is not real! eta %d\n', eta)
    end
    if (eta==1)
        maxexp = max(exponents{eta},[],2); % column vector of size Ni
        minexp = min(exponents{eta},[],2); % column vector of size Ni
    else
        maxexp = max(maxexp, max(exponents{eta},[],2)); % column vector of size Ni
        minexp = min(minexp, min(exponents{eta},[],2)); % column vector of size Ni
    end
end

% compute the scaling factor for the exponents in order to prevent overflow
% maxexp and minexp are column vectors with size Ni.
if (size(maxexp,1) ~= Ni || size(minexp,1) ~= Ni)
    error('size of maxexp %d, size of minexp %d, shhould equal to Ni %d!\n', size(maxexp,1), size(minexp,1), Ni)
end
expshift = zeros(Ni,1); % vector of size Ni
for img=1:Ni
    expshift(img)=scaling_rule(minexp(img), maxexp(img));
end

%% compute $$p(\theta_i,\eta_i | y_i, \bar{c}, V, q}$$
Zeta=cell(Neta, 1); % numerator $$p(y_i|\theta_i,\eta_i,\bar{c}^{\eta_i},V_{\eta_i})q_{\eta_i}$$
Zetabar=zeros(Ni,1); % denominator
p=cell(Neta, 1);

for eta=1:Neta
    Zeta{eta} = zeros(Ni, Nzi);
    for n=1:Nzi
        Zeta{eta}(:,n) = exp(exponents{eta}(:,n)+expshift(:)) * vo{eta}.q;
        Zetabar(:) = Zetabar(:) + Zeta{eta}(:,n)*rule(n,WEIGHT_INDEX);
    end
end
assert(all(Zetabar(:)~=0.0),'set_p_theta_eta_loglike_diagV: assert(all(Zetabar(:)~=0.0))');

for eta=1:Neta
    p{eta} = bsxfun(@rdivide, Zeta{eta}, Zetabar);
end

%% compute the log likelihood
loglike = sum(log(Zetabar)) - sum(expshift); % need to undo the effect of shifting exponents
%fprintf(1, 'set_p_theta_eta_nu_noL: log likelihood %.16g.\n', loglike);

fprintf(1,'set_p_theta_eta_loglike_diagV time: %d\n',etime(clock, t_p));
