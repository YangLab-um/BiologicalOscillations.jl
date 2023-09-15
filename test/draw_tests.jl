using BiologicalOscillations, Images, FileIO

test_base = dirname((Base.source_path()))

image = draw_connectivity([0 0 1; 1 0 0; 0 -1 0])
ref = load(joinpath(test_base, "draw_connectivity", "01.png"))
@test Gray.(image) == Gray.(ref)

image = draw_connectivity([1 0 1; 1 0 0; 0 -1 0])
ref = load(joinpath(test_base, "draw_connectivity", "02.png"))
@test Gray.(image) == Gray.(ref)

image = draw_connectivity([1 -1 1; 1 0 0; 0 -1 0])
ref = load(joinpath(test_base, "draw_connectivity", "03.png"))
@test Gray.(image) == Gray.(ref)

image = draw_connectivity([0 0 0 0 1; -1 0 0 0 0; 0 -1 -1 0 0 ; -1 0 1 0 0; 0 0 0 -1 0])
ref = load(joinpath(test_base, "draw_connectivity", "04.png"))
@test Gray.(image) == Gray.(ref)

image = draw_connectivity([0 0 0 0 1; -1 0 0 0 0; 0 -1 -1 0 0 ; -1 0 1 0 0; 0 0 0 -1 0], coherent_nodes=[3], incoherent_nodes=[4], positive_cycles=[[1, 4, 5]], negative_cycles=[[3]])
ref = load(joinpath(test_base, "draw_connectivity", "05.png"))
@test Gray.(image) == Gray.(ref)

image = draw_connectivity([1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1; 1 1 1 1 1 1 1])
ref = load(joinpath(test_base, "draw_connectivity", "06.png"))
@test Gray.(image) == Gray.(ref)
