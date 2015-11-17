function [cbarnew, qnew, loglike, solvedlinearsystem]=EMstep_cbar_dV_large(rule, vk, y, noise_var, virusobjs, all_tilde_b, MinimumClassProb, absthres)
%function [cbarnew, qnew, loglike, solvedlinearsystem]=EMstep_cbar_dV_large(rule, vk, y, noise_var, virusobjs, all_tilde_b, MinimumClassProb, absthres)
%
% find cbar for the diagonal V problem
% tilde_b is passed to this function

% Qiu Wang
% 05/25/2010

%t_p=clock;

fprintf('EMstep_cbar_dV_large: entering\n');
WEIGHT_INDEX = 6; % the index of the abscissa weights in the rule data structure

vo = virusobjs; % use a shorter name inside the function

%% initialize variables
Neta = size(vo, 1);
%Ny = size(y, 1);
Nv = size(y, 2);
Nzi = size(rule, 1);

%solvedlinearsystem(1)=1 if solved F\g, =0 if used previous answer
solvedlinearsystem=ones(Neta,1);

% allocating storage
cbarnew = cell(Neta, 1);
qnew = zeros(Neta, 1);

%% compute $$p(\theta_i,\eta_i|y_i,\bar{c}, V, q)$$ by Eq. (25)
%p_theta_eta{1:Neta}(1:Nv,1:Nzi)
[p_theta_eta, loglike]=set_p_theta_eta_loglike_diagV(rule, vk, y, noise_var, vo, all_tilde_b);
im2use=Image2use(p_theta_eta,rule,absthres);
% im2use=cell(1,5000);
% vect=ones(1,1200);
% for j=1:1200
%     vect(j)=j;
% end
% for i=1:5000
%     im2use{i}=vect;
% end
n_nonzero_p_theta_eta=zeros(Neta,1);
n_total_p_theta_eta=zeros(Neta,1);
for eta=1:Neta
  n_nonzero_p_theta_eta(eta)=sum(sum(p_theta_eta{eta}>0.0));
  n_total_p_theta_eta(eta)=prod(size(p_theta_eta{eta}));
end
fprintf(1,'EMstep_cbar_dV_large: n_total_p_theta_eta:');
fprintf(1,' %d',n_total_p_theta_eta);
fprintf(1,' n_nonzero_p_theta_eta:');
fprintf(1,' %d',n_nonzero_p_theta_eta);
fprintf(1,'\n');

%fprintf(1,'initial loglike = %f\n',loglike);

%% compute the expectation and maximize it to get new $$\bar{c}$$ and $$V$$
for eta=1:Neta
    Nc = size(vo{eta}.cbar,1);
    F = zeros(Nc, Nc);
    g = zeros(Nc, 1);
    qnew_eta = 0;
    % compute the expectation with numerical integrations
    %for n=1:Nzi

    %will use tilde_b values all_tilde_b{eta}.tilde_b
    
    parfor (n=1:Nzi)
        
        
        %       if mod(n,1000)==0
        %           fprintf(1,'n = %d\n',n);
        %       end
        %maxNumCompThreads(4);
        % for computing the log likelihood Qem
        Rabc = euler2R(rule(n,1:3));
        L=setL_nord(Rabc, vo{eta}.clnp.il, vo{eta}.clnp.in, vk, vo{eta}.Htable, vo{eta}.map_unique2lp, all_tilde_b{eta}.tilde_b, vo{eta}.Ftable, vo{eta}.BasisFunctionType);

        [LT_SigmaInv_y LT_SigmaInv_L]=LT_SigmaInv_fast(noise_var,L,vo{eta}.nu,y(:,im2use{eta,n}));

        wri = rule(n,WEIGHT_INDEX)*p_theta_eta{eta}(im2use{eta,n},n);
%         pa=p_theta_eta{eta}(im2use{eta,n},n);
%         pb=p_theta_eta{eta}(:,n);
%         if pa~=pb
%             fprintf('n=%d, pa~=pb',n);
%         end
%         ya=y(:,im2use{eta,n});
%         yb=y;
%         if ya~=yb
%             fprintf('n=%d,ya~=yb',n);
%         end
        
        % update F and g
        %wri = rule(n,WEIGHT_INDEX)*sum(p_theta_eta{eta}(:,n)); % vector
        g = g + LT_SigmaInv_y*wri; % vector of size Nc
        sum_wri = sum(wri);
        F = F + LT_SigmaInv_L*sum_wri; % matrix of size Nc-by-Nc
        qnew_eta = qnew_eta + sum_wri;

        %     % fast compute SigmaInv
        %     SigmaInv=SigmaInv_fast(noise_var,L,vo{eta}.nu);
        %
        % %     Sigma=L*diag(vo{eta}.nu)*L';
        % %
        % %     % add noise_var on the diagonal of Sigma; this equals to adding a Q
        % %     % matrix that is noise_var * eye
        % %     for ii=1:size(Sigma, 1);
        % %       Sigma(ii,ii) = Sigma(ii,ii) + noise_var;
        % %     end
        %
        %     %LT_Sigma_inv = L'/Sigma;  % L'*Sigma_inv; compute the right inverse
        %     LT_Sigma_inv = L'*SigmaInv;  % L'*Sigma_inv; compute the right inverse
        %     wri = rule(n,WEIGHT_INDEX)*p_theta_eta{eta}(:,n);
        %     % update F and g
        %     %wri = rule(n,WEIGHT_INDEX)*sum(p_theta_eta{eta}(:,n)); % vector
        %     g = g + LT_Sigma_inv*y*wri; % vector of size Nc
        %     sum_wri = sum(wri);
        %     F = F + LT_Sigma_inv*L*sum_wri; % matrix of size Nc-by-Nc
        %     qnew_eta = qnew_eta + sum_wri;
        %fprintf(1,'elapsed time %d\n',etime(clock,t_p));

    end
    qnew(eta) = qnew_eta;

    % maximization step
    MaxCondNumber4F=1.0e8;
    if cond(F)<MaxCondNumber4F
      %do the update
      cbarnew{eta}=F\g;
    else
      %do not do the update
      if any(vo{eta}.clnp.c-vo{eta}.cbar~=0.0) %see assert_UpdateBothCopiesOfcbar.m
        error('EMstep_cbar_dV_large: cond(F) large and .clnp.c~=.cbar\n');
      end
      cbarnew{eta}=vo{eta}.clnp.c;
      solvedlinearsystem(eta)=0;
    end

end % for eta=1:Neta

% update the class probabilites
qnew = qnew / Nv;
qnew=applyclassproblimit(qnew,MinimumClassProb);
