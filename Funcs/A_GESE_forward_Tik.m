function [ rec ] = A_GESE_forward_Tik(a, lambda,T_for,P_for,S_for,F_for,M_for)

rec = M_for(F_for(S_for(P_for(T_for(a)))));
rec = [rec; lambda*a(:)];

end
