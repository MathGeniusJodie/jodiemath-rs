(* Define your function eg. CubeRoot[x] *)
g[x_] := CubeRoot[x]

inv[x_] := FullSimplify[InverseFunction[g][x], Assumptions -> x > 0]

f[s_] := 1/(inv[s]-x)

Print["       Newton: ", FullSimplify[s+f[s]/f'[s]]]
Print["       Halley: ", FullSimplify[s+2*f'[s]/f''[s]]]
Print["Householder 3: ", FullSimplify[s+3*f''[s]/f'''[s]]]
Print["Householder 4: ", FullSimplify[s+4*f'''[s]/f''''[s]]]
Print["Householder 5: ", FullSimplify[s+5*f''''[s]/f'''''[s]]]
Print["Householder 6: ", FullSimplify[s+6*f'''''[s]/f''''''[s]]]