using BEAST
using CompScienceMeshes



h = 2π / 4; Γ = meshcircle(1.0, h)
X = lagrangecxd0(Γ)

elements, ad, nr = assemblydata(X)
test_elements = elements
trial_elements = elements


e = test_elements[1]

dom_e = domain(e)

I = [2,1]
ichart_e = CompScienceMeshes.permute_vertices(dom_e, I)






dom1 = domain(chart1)
dom2 = domain(chart2)


ichart2 = CompScienceMeshes.permute_vertices(dom2, J)