using StructuresKit

#Ix Iy Ixy J Cw
section_properties = [(3.230E6,449530,-865760, 397.09, 3.4104E9)]

#E  ν
material_properties = [(200,0.30)]

#kx kϕ
spring_stiffness = [(0.0,300/1000)]

#ay_kx
spring_location = [(101.6)]

#qx qy start end
loads = [(0.0001, 0.00002, 0.0, 7620.0)]

#ax ay
load_locations = [(27.826,101.6)]

#end boundary conditions
#type=1 u''=v''=ϕ''=0 (simply supported), type=2 u'=v'=ϕ'=0  (fixed), type=3 u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free end, e.g., a cantilever)
end_boundary_conditions = [1 1]

#supports
#location where u=v=ϕ=0
supports = [0.0 7620]

#member information
#L(1) dL(2) section_properties(3) material_properties(4) spring_stiffness(5) spring_location(6) load(7) load_location(8) 
member_definitions = [(7620, 7620/12, 1, 1, 1, 1, 1, 1)]

#Discretize beam.
z, m = ThinWalledBeam.discretize(member_definitions)

#Apply loads.
num_nodes = length(z)
qx = 0.001 * ones(num_nodes)
qy = 0.001 * ones(num_nodes)
ax = 27.8 * ones(num_nodes)
ay = 101.6 * ones(num_nodes)

#Define springs.
kx = 0.0 * ones(num_nodes)
kϕ = 0.3 * ones(num_nodes)
ay_kx = 101.6 * ones(num_nodes)


#Define model.
model = ThinWalledBeam.define(z, m, member_definitions, section_properties, material_properties, kx, kϕ, ay_kx, qx, qy, ax, ay, end_boundary_conditions, supports)

#Solve model.
model = ThinWalledBeam.solve(model)