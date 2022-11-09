import FLOWFMM
fmm = FLOWFMM
import WriteVTK

function save_vtk(filename, elements; compress=false)
    n_bodies = size(elements.bodies)[2]
    WriteVTK.vtk_grid(filename, reshape(elements.bodies[fmm.i_POSITION,:], 3, n_bodies,1,1); compress) do vtk
        vtk["scalar strength"] = reshape(elements.bodies[fmm.i_STRENGTH[1],:],n_bodies,1,1)
        vtk["vector strength"] = reshape(elements.bodies[fmm.i_STRENGTH[2:4],:],3,n_bodies,1,1)
        vtk["scalar potential"] = reshape(elements.potential[fmm.i_POTENTIAL[1],:],n_bodies,1,1)
        vtk["vector potential"] = reshape(elements.potential[fmm.i_POTENTIAL[2:4],:],3,n_bodies,1,1)
        vtk["velocity"] = reshape(elements.velocity,3,n_bodies,1,1)
    end
end
