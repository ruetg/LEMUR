mutable struct lemur_obj
    k::Array{Float64,2}
    undercapacity::Array{Float64,2}
    z::Array{Float64,2}
    dx::Float64
    dy::Float64
    t::Float64
    dt::Float64
    display::Float64
    Udt::Float64
    flex::Float64
    kd::Float64
    m::Float64
    n::Float64
    bc::Array{Float64,1}
    wdt::Float64
    firstcall::Float64
    u::Array{Float64,2}
    udt::Float64
    ny::Float64
    nx::Float64
    deposit::Float64
    srho::Float64
    drawdt::Float64
    l::Float64
    ks::Float64
    sinkfill::Float64
    uselandsed::Float64
    maxareasinkfill::Float64
    precip::Float64
    evaprate::Float64
    bcx::Array{Float64,2}
    k_sed::Array{Float64,2}

    function lemur_obj()
        obj = new()
        obj.uselandsed = 0;
        obj.sinkfill = 1;
        obj.ny = 500
        obj.nx = 500
        obj.z = rand(Int(obj.ny), Int(obj.nx)) .* 10 .+ 100;
        #obj.z ./=maximum(vec(obj.z));
        obj.display = 1
        obj.udt = 100000
        obj.deposit = false
        bc = zeros(Int(obj.ny), Int(obj.nx))
        bc[:, 1] .= 0
        bc[:, end] .= 0
        bc[1, :] .= 1
        bc[end, :] .= 0
        obj.bcx = zeros(Int(obj.ny),Int(obj.nx))
        obj.bcx[:] .= bc[:]
        obj.z[bc .== 1] .= 0
        bc = vec(bc)
        obj.srho = 2400
        obj.udt = 1e5
        obj.bc = float(findall(bc .== 1))
        obj.k = zeros(Int(obj.ny),Int(obj.nx)) .+ 1e-6
        obj.firstcall = 1
        obj.l = 1
        obj.evaprate = 1e-30
        obj.maxareasinkfill = 1000000
        obj.precip = 1.0
        obj.kd = 1.0
        obj.ks = 0
        obj.dt = 1e5
        obj.t = 3e5
        obj.undercapacity = zeros(Int(obj.ny), Int(obj.nx)) 
        obj.k_sed = obj.k
        obj.flex = 10e3
        obj.dx = 1000
        obj.dy = 1000
        obj.m = .5
        obj.n = 1
        obj.u = zeros(Int(obj.ny), Int(obj.nx)) .+ .000000000000001
        return obj
    end
    lemur_obj()
    
end
module run_lemur

using PyPlot
using Cxx
using Libdl
include("/Users/gr_1/Documents/IsoFlex/IsoFlex.jl")



pw = pwd() * "/../src/"
Libdl.dlopen(pw * "liblemur", Libdl.RTLD_GLOBAL)
cxxinclude(pw * "lemur.h")

export run_lemur, lemur_obj, set_lemur, get_lemur

function set_lemur(maker, nm, var::Array{Float64})
    var = vec(var)
    x = convert(cxxt"double *",pointer(var));
    @cxx maker -> set(pointer(nm), x, length(var))
    return

end
function set_lemur(maker, nm, var::Float64)
    print("here2")
    x = convert(cxxt"double",var);
    @cxx maker -> set(pointer(nm), x)
end
function get_lemur(maker, nm, ny, nx)
    
    y = @cxx maker -> get(pointer(nm))
    z = unsafe_wrap(Array,pointer(y),Int(ny*nx))
    z = reshape(z,Int(ny),Int(nx))
    return z
    
end
function run(lemur_params)
    model = @cxxnew lemur(lemur_params.ny,lemur_params.nx)

    for nm = fieldnames(typeof(lemur_params))
        print(string(nm,' '))
        if length(getfield(lemur_params,nm)) <= 1
            println(getfield(lemur_params,nm))
        end

        set_lemur(model, string(nm), getfield(lemur_params,nm))
    end
    print("here")


    z = get_lemur(model, "z", lemur_params.ny, lemur_params.nx)
    zi = zeros(size(z))
    flex = zeros(size(z))
    ero = zeros(size(z))
    z[vec(lemur_params.bcx .== 1)] .= 0
    set_lemur(model,"z", vec(z))
    @cxx model -> lakefill()
    z = copy(get_lemur(model, "z", lemur_params.ny, lemur_params.nx))

    zi[:] .= z[:]

    for t = 0:lemur_params.dt:lemur_params.t    
        @cxx model -> erosion_fluvial()
        z =  copy(get_lemur(model, "z", lemur_params.ny, lemur_params.nx))
        z[lemur_params.bcx .== 1] .= 0
        ero .= copy(zi .- z)
        zi .= z

        z[lemur_params.bcx .== 0] .+= lemur_params.u[lemur_params.bcx .== 0]
        u = IsoFlex.flexural(ero,dx=lemur_params.dx,dy=lemur_params.dy,Te=lemur_params.flex,buffer=200)

       # u = IsoFlex.flexural(ero,dx=lemur_params.dx,dt=lemur_params.dt, tt = lemur_params.t, t = t,
       # dy=lemur_params.dy,Te=lemur_params.flex,buffer=200)
        
        z[lemur_params.bcx .== 0] .+= u[lemur_params.bcx .== 0]

        set_lemur(model,"z", vec(z))
        @cxx model -> lakefill()
        
        subplots(figsize=[5,5])
        imshow(z)
        colorbar()
        println("1")
    end
    model = nothing

    return ero
end
end