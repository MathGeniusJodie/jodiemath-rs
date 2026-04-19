(* Define your function eg. CubeRoot[x] *)
g[x_] := CubeRoot[x]

inv[x_] := FullSimplify[InverseFunction[g][x], Assumptions -> x > 0]

f[y_] := 1/(inv[y]-x)

Print["       Newton: ", FullSimplify[y+f[y]/f'[y]]]
Print["       Halley: ", FullSimplify[y+2*f'[y]/f''[y]]]
Print["Householder 3: ", FullSimplify[y+3*f''[y]/f'''[y]]]
Print["Householder 4: ", FullSimplify[y+4*f'''[y]/f''''[y]]]
Print["Householder 5: ", FullSimplify[y+5*f''''[y]/f'''''[y]]]
Print["Householder 6: ", FullSimplify[y+6*f'''''[y]/f''''''[y]]]