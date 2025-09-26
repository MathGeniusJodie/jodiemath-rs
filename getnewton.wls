


(* Define your function *)
g[x_] := CubeRoot[x]

inv[x_] := FullSimplify[InverseFunction[g][x], Assumptions -> x > 0]

f[s_] := 1/(inv[s]-x)

Print["       Newton: ", FullSimplify[f[s]/f'[s]]]
Print["       Halley: ", FullSimplify[2*f'[s]/f''[s]]]
Print["Householder 3: ", FullSimplify[3*f''[s]/f'''[s]]]
Print["Householder 4: ", FullSimplify[4*f'''[s]/f''''[s]]]
Print["Householder 5: ", FullSimplify[5*f''''[s]/f'''''[s]]]
Print["Householder 6: ", FullSimplify[6*f'''''[s]/f''''''[s]]]