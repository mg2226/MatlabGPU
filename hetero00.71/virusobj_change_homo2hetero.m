function vobj=virusobj_change_homo2hetero(vobj,homo2hetero)
%function vobj=virusobj_change_homo2hetero(vobj,homo2hetero)
%vobj: virus object
%homo2hetero{1:Neta}.action -- 1 means heterogeneous, 0 means no change, -1 means homogeneous
%homo2hetero{1:Neta}.FractionOfMeanForMinimum
%homo2hetero{1:Neta}.FractionOfMean

%additional case options could include additional rules for the heterogeneous nu.

for eta=1:length(vobj)

  switch homo2hetero{eta}.action

    case 1 %make this class heterogeneous
      fprintf(1,'virusobj_change_homo2hetero: eta %d: make heterogeneous\n',eta);
      if homo2hetero{eta}.FractionOfMean>=0.0
        vobj{eta}.nu=mean2stddev(vobj{eta}.cbar,homo2hetero{eta}.FractionOfMeanForMinimum,homo2hetero{eta}.FractionOfMean).^2;
      else
        vobj{eta}.nu=mean2stddev(vobj{eta}.cbar,homo2hetero{eta}.FractionOfMeanForMinimum,homo2hetero{eta}.FractionOfMean);
      end

    case 0 %no change
      fprintf(1,'virusobj_change_homo2hetero: eta %d: no change\n',eta);

    case -1 %make this class homogeneous
      fprintf(1,'virusobj_change_homo2hetero: eta %d: make homogeneous\n',eta);
      vobj{eta}.nu=[];

    otherwise
      error('virusobj_change_homo2hetero: eta %d unknown operator %s\n',ii,homo2hetero{ii}.action);

    end

end
