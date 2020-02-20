function e_out = EulerParams(R)

if size(R) == 2
    R_new = eye(3);
    R_new(1:2,1:2) = R;
    R = R_new;
    clear R_new;
end
e_out = zeros(4,1);

e_out(1) = sqrt(trace(R)+1)/2;

if e_out(1) ~= 0
    e_out(2) = (R(2,3)-R(3,2))/(4*e_out(1));
    e_out(3) = (R(3,1)-R(1,3))/(4*e_out(1));
    e_out(4) = (R(1,2)-R(2,1))/(4*e_out(1));
else
    e22 = (1+2*R(1,1)-trace(R))/4;
    e32 = (1+2*R(2,2)-trace(R))/4;
    e42 = (1+2*R(3,3)-trace(R))/4;
    if e22 ~= 0
        e_out(2) = sqrt(e22);
        e_out(3) = (R(2,1)+R(1,2))/(4*e_out(2));
        e_out(4) = (R(3,1)+R(1,3))/(4*e_out(2));
    elseif e32 ~= 0
        e_out(3) = sqrt(e32);
        e_out(2) = (R(2,1)+R(1,2))/(4*e_out(3));
        e_out(4) = (R(3,2)+R(2,3))/(4*e_out(3));
    else
        e_out(4) = sqrt(e42);
        e_out(2) = (R(3,1)+R(1,3))/(4*e_out(4));
        e_out(3) = (R(2,3)+R(3,2))/(4*e_out(4));
    end
end
