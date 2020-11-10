function [ rec ] = A_GESE_adjoint_Tik(y, Nsize, lambda, T_adj,P_adj,S_adj,F_adj,M_adj)


rec = T_adj(P_adj(S_adj(F_adj(M_adj( y(1:Nsize) )))));
rec = rec+lambda*y(Nsize+1:end);

end
