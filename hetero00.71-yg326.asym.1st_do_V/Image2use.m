function [im2use] = Image2use(p_theta_eta,rule,thres)
%Select image to use for each abscissa 
% Yunye Gong
% 11/06/2013

f=thres;
Neta=size(p_theta_eta,1);
% if isempty(p_theta_eta{1})
%     error('empty p_theta_eta');
% end
Nv=size(p_theta_eta{1},1);
Nzi=size(p_theta_eta{1},2);

ab_im=zeros(Neta,Nzi,Nv);
WEIGHT_INDEX=6;

for eta=1:Neta
    p_eta=p_theta_eta{eta};
    for ii=1:Nv
        %pw=bsxfun(@times,p_eta(ii,:),rule(:,WEIGHT_INDEX)');
        pw=p_eta(ii,:).*rule(:,WEIGHT_INDEX)';
        pstar=max(pw);
        pmin=min(pw);
        if pmin<0.0
            fprintf('minimum probability is lower than 0');
        end
%         fprintf('minimum prob is %d',pmin);
        ab_indices=find(pw>=f*pstar);
        
        for j=1:length(ab_indices)
            ab_im(eta,ab_indices(j),ii)=1;
        end
    end
end

im2use=cell(Neta,Nzi);
for eta=1:Neta
    for ii=1:Nzi
        imlist=find(ab_im(eta,ii,:)==1);
%         imlist=cumsum(ones(1,1200));
        im2use{eta,ii}=imlist;
%         if length(imlist)~=1200
%             fprintf('not selecting fullset');
%         end
    end
end
%clear index data

%save('im2use.mat','im2use');
clear ab_im;
