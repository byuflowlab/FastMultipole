using WriteVTK
using HDF5
using BenchmarkTools

struct One
    bodies
end

n_bodies = 1000000
n_position = 3
n_strength=4
n_velocity = 3
i_position = 1:n_position
i_strength_s = n_position+1
i_strength_v = n_position+2:n_position+n_strength
i_velocity = n_position+n_strength+1:n_position+n_strength+n_velocity
n_total = n_position+n_strength+n_velocity
one_bodies = rand(n_total,n_bodies,1,1)
one = One(one_bodies)

function xmf_file(tag, h5fname, n_bodies)
    # xml file for paraview
    # Generates XDMF file specifying fields for paraview
    xmf = open("test_one"*string(tag)*".xmf", "w")

    # Open xmf block
    print(xmf, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n")
    print(xmf, "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"3.0\">\n")
        print(xmf, "\t<Domain>\n")
        print(xmf, "\t\t<Grid Name=\"particles\" GridType=\"Uniform\">\n")

            # print(xmf, "\t\t\t\t<Time Value=\"", time, "\" />\n")

            # Nodes: particle positions
            print(xmf, "\t\t\t\t<Geometry Type=\"XYZ\">\n")
                print(xmf, "\t\t\t\t\t<DataItem DataType=\"Float\"",
                            " Dimensions=\"", n_bodies, " ", 3,
                            "\" Format=\"HDF\" Precision=\"8\">",
                            h5fname, ":position</DataItem>\n")
            print(xmf, "\t\t\t\t</Geometry>\n")

            # Topology: every particle as a point cell
            print(xmf, "\t\t\t\t<Topology Dimensions=\"", n_bodies, "\" Type=\"Polyvertex\"/>\n")

            # Attribute: Gamma
            print(xmf, "\t\t\t\t<Attribute Center=\"Node\" Name=\"Gamma\" Type=\"Vector\">\n")
                print(xmf, "\t\t\t\t\t<DataItem DataType=\"Float\"",
                            " Dimensions=\"", n_bodies, " ", 3, "\" Format=\"HDF\" Precision=\"8\">",
                            h5fname, ":vector_strength</DataItem>\n")
            print(xmf, "\t\t\t\t</Attribute>\n")

            # Attribute: Gamma
            print(xmf, "\t\t\t\t<Attribute Center=\"Node\" Name=\"Velocity\" Type=\"Vector\">\n")
            print(xmf, "\t\t\t\t\t<DataItem DataType=\"Float\"",
                        " Dimensions=\"", n_bodies, " ", 3, "\" Format=\"HDF\" Precision=\"8\">",
                        h5fname, ":velocity</DataItem>\n")
            print(xmf, "\t\t\t\t</Attribute>\n")

            # Attribute: sigma
            print(xmf, "\t\t\t\t<Attribute Center=\"Node\" Name=\"scalar_strength\" Type=\"Scalar\">\n")
                print(xmf, "\t\t\t\t\t<DataItem DataType=\"Float\"",
                            " Dimensions=\"", n_bodies, "\" Format=\"HDF\" Precision=\"8\">",
                            h5fname, ":scalar_strength</DataItem>\n")
            print(xmf, "\t\t\t\t</Attribute>\n")

            

        print(xmf, "\t\t</Grid>\n")
        print(xmf, "\t</Domain>\n")
    print(xmf, "</Xdmf>\n")

    close(xmf)
end

function vtk_noview(one::One; wvtk=true, hdf5=false)
    x = reshape(one.bodies[1:n_position, :, 1, 1], n_position, n_bodies,1,1)
    ss = reshape(one.bodies[i_strength_s,:, 1, 1], n_bodies,1,1)
    vs = reshape(one.bodies[i_strength_v,:, 1, 1], n_position, n_bodies,1,1)
    v = reshape(one.bodies[i_velocity,:,1,1], n_position, n_bodies,1,1)
    tag = rand(Int64)
    if wvtk
        vtk_grid("test_one"*string(tag), x; compress = false) do vtk
            vtk["scalar_strength"] = ss
            vtk["vector_strength"] = vs
            vtk["velocity"] = v
        end
    elseif hdf5
        # binary data to file
        h5fname = "test_one"*string(tag)*".h5"
        h5 = HDF5.h5open(h5fname, "w")
        h5["position"] = x
        h5["scalar_strength"] = ss
        h5["vector_strength"] = vs
        h5["velocity"] = v
        close(h5)

        xmf_file(tag, h5fname, n_bodies)
    end
end

struct Many
    position
    scalar_strength
    vector_strength
    velocity
end

positions = rand(3,n_bodies,1,1)
scalar_strength = rand(n_bodies,1,1)
vector_strength = rand(3,n_bodies,1,1)
velocity = rand(3,n_bodies,1,1)
many = Many(positions, scalar_strength, vector_strength, velocity)

function vtk_noview(many::Many; wvtk=true, hdf5=false)
    tag = rand(Int64)
    if wvtk
        vtk_grid("test_many"*string(tag), many.position; compress = false) do vtk
            vtk["scalar_strength"] = many.scalar_strength
            vtk["vector_strength"] = many.vector_strength
            vtk["velocity"] = many.velocity
        end
    elseif hdf5
        h5fname = "test_many"*string(tag)*".h5"
        h5 = HDF5.h5open(h5fname, "w")
        h5["position"] = many.position
        h5["scalar_strength"] = many.scalar_strength
        h5["vector_strength"] = many.vector_strength
        h5["velocity"] = many.velocity
        close(h5)

        xmf_file(tag, h5fname, n_bodies)
    end
end

@btime vtk_noview(one)
@btime vtk_noview(one)
@btime vtk_noview(one; wvtk=false, hdf5=true)
@btime vtk_noview(one; wvtk=false, hdf5=true)

@btime vtk_noview(many)
@btime vtk_noview(many)
@btime vtk_noview(many; wvtk=false, hdf5=true)
@btime vtk_noview(many; wvtk=false, hdf5=true)

# check structure
# Ps = [
#     rand(3),
#     rand(3),
#     rand(3),
#     rand(3)
# ]

# X = [P[i] for i in 1:3, P in Ps]

# function hdf5(one::One)
#     # Creates/overwrites HDF5 file
#     h5 = HDF5.h5open("test.hdf5", "w")

#     # Writes fields
#     # NOTE: It is very inefficient to convert the data structure to a matrices
#     # like this. This could help to make it more efficient: https://stackoverflow.com/questions/58983994/save-array-of-arrays-hdf5-julia
#     # UPDATE 2021/11: I tried multiple ways of pre-allocating memory in the disk
#     #   through HDF5 and then dumping data into it from pfield through
#     #   iterators, but for some reason HDF5 always re-allocates memory
#     #   when trying to write anything but arrays.
#     h5["X"] = [P.X[i] for i in 1:3, P in iterate(self; include_static=true)]
#     h5["Gamma"] = [P.Gamma[i] for i in 1:3, P in iterate(self; include_static=true)]
#     h5["sigma"] = [P.sigma[1] for P in iterate(self; include_static=true)]
#     h5["circulation"] = [P.circulation[1] for P in iterate(self; include_static=true)]
#     h5["vol"] = [P.vol[1] for P in iterate(self; include_static=true)]
#     h5["static"] = Int[P.static[1] for P in iterate(self; include_static=true)]
#     h5["i"] = [P.index[1] for P in iterate(self; include_static=true)]

#     if isLES(self)
#         h5["C"] = [P.C[i] for i in 1:3, P in iterate(self; include_static=true)]
#     end

#     # # Connectivity information
#     # h5["connectivity"] = [i%3!=0 ? 1 : Int(i/3)-1 for i in 1:3*np]

#     # # Write fields
#     # dtype = HDF5.datatype(T)
#     #
#     # for (field, dim) in [("X", 3), ("Gamma", 3), ("sigma", 1)] # Iterate over fields
#     #
#     #     dims = dim==1 && false ? HDF5.dataspace(np) : HDF5.dataspace(dim, np)
#     #     chunk = dim==1 && false ? (np,) : (1, np)
#     #     dset = HDF5.d_create(h5, field, dtype, dims, "chunk", chunk)
#     #
#     #     for (pi, P) in enumerate(iterator(self; include_static=true))
#     #         dset[:, pi] .= getproperty(P, Symbol(field))
#     #     end
#     #
#     # end

#     # # Connectivity information
#     # dtype = HDF5.datatype(Int)
#     # dims = HDF5.dataspace(3*np, 1)
#     # chunk = (np, 1)
#     # dset = HDF5.d_create(h5, "connectivity", dtype, dims, "chunk", chunk)
#     # for i in 1:np
#     #     dset[3*(i-1)+1, 1] = 1
#     #     dset[3*(i-1)+2, 1] = 1
#     #     dset[3*(i-1)+3, 1] = i-1
#     # end

#     close(h5)