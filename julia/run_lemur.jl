
module run_lemur

using PyPlot
using Cxx
using Libdl
include("/Users/gr_1/Documents/IsoFlex/IsoFlex.jl")
pw = pwd() * "/../src/"
Libdl.dlopen(pw * "liblemur", Libdl.RTLD_GLOBAL)
cxxinclude(pw * "lemur.h")

export run_lemur
function set_lemur(maker, nm, var)
    if length(var) > 1
        var = vec(var)
        x = convert(cxxt"double *",pointer(var));
        @cxx maker -> set(pointer(nm), x, length(var))
        
    else
        x = convert(cxxt"double",var);
        try
            @cxx maker -> set(pointer(nm), x)
        catch
            print(string("invalid param", nm))
        end
    end
    return

    end
function get_lemur(maker, nm, ny, nx)
    y = @cxx maker -> get(pointer(nm))
    z = unsafe_wrap(Array,pointer(y),ny*nx)
    z = reshape(z,ny,nx)
    return z
end
function run(lemur_params)
    model = @cxxnew lemur(lemur_params.ny,lemur_params.nx)

    for nm = fieldnames(typeof(lemur_params))
        print(string(nm,' '))
        try
            set_lemur(model, string(nm), getfield(lemur_params,nm))
        catch
            print(" error ")
        end
    end


    z = get_lemur(model, "z", lemur_params.ny, lemur_params.nx)
    zi = zeros(size(z))
    flex = zeros(size(z))
    ero = zeros(size(z))

    for t = 0:lemur_params.dt:lemur_params.t

        @time begin

        @cxx model -> lakefill()
        @cxx model -> erosion_fluvial()
        z =  copy(get_lemur(model, "z", lemur_params.ny, lemur_params.nx))
        ero = copy(z .- zi)

        zi = copy(z)

        z[lemur_params.bcx .== 0] .+= 10
        #u = IsoFlex.flexural(ero,dx=lemur_params.dx,dy=lemur_params.dy,Te=lemur_params.flex,buffer=200)

        #u = IsoFlex.flexural(ero,dx=lemur_params.dx,dt=lemur_params.dt, tt = lemur_params.t, t = t,
         #           dy=lemur_params.dy,Te=lemur_params.flex,buffer=200)
       # z[lemur_params.bcx .== 0] .+= u[lemur_params.bcx .== 0]
        set_lemur(model,"z", vec(z))
        end

        subplots(figsize=[5,5])
        imshow(z)
        colorbar()
        println("1")

    end
end
end