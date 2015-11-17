function d=rot_metric(R1,R2,type)
% given two 3x3 rotation matices, 5 definitions of the metric
% d1 -norm of the differences of two quaternions
if type==1
    q1=rotation2q(R1);
    q2=rotation2q(R2);
    d=min(norm(q1-q2,2),norm(q1+q2,2));
    if (d1<0)||(d1>sqrt(2))
        fprintf('d1 is not in range %d',d1);
    end
    d=d1;
% % d2 -inner product of unit quaternions
elseif type==2
    d2=acos(abs(q1*(q2.')));
    if (d2<0)||(d2>pi/2)
        fprintf('d2 is not in range %d',d2);
    end
    d=d2;
% % d3 - inner product without arccos calculation
elseif type==3
    d3=1-abs(q1*(q2.'));
    if (d3<0)||(d3>1)
        fprintf('d3 is not in range %d',d3);
    end
    d=d3;
% % d4 - deviation from the identity matrix
elseif type==4
    I=[1 0 0; 0 1 0; 0 0 1];
    d4=norm(I-R1*(R2.'),2);
    if (d4<0)||(d4>2.001)
        fprintf('d4 is not is range %d',d4);
    end
    d=d4;
% d5 - geodesic on the unit sphere
elseif type==5
    d5=norm(logm(R1*(R2.')),2);
    if (d5<0)||(d5>pi)
        fprintf('d5 is not in range %d',d5);
    end
    d=d5;
else
    error('no such metric type');
end
