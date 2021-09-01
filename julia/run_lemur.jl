
export set_param
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
    t_c::Float64
    sinkfill::Float64
    uselandsed::Float64
    maxareasinkfill::Float64
    precip::Float64
    evaprate::Float64
    bcx::Array{Float64,2}
    k_sed::Array{Float64,2}
    ufunc::Function
    function lemur_obj()
        function tempfunc(t)
            return nothing
        end
        obj = new()
        obj.uselandsed = 0;
        obj.sinkfill = 1;
        obj.ny = 500
        obj.nx = 500
        obj.z = rand(Int(obj.ny), Int(obj.nx)) .* 10 .+ 2000;
        #obj.z ./=maximum(vec(obj.z));
        obj.display = 1
        obj.udt = 100000
        obj.deposit = false
        bc = zeros(Int(obj.ny), Int(obj.nx))
        bc[:, 1] .= 1
        bc[:, end] .= 1
        bc[1, :] .= 1
        bc[end, :] .= 1
        obj.bcx = zeros(Int(obj.ny),Int(obj.nx))
        obj.bcx[:] .= bc[:]
        obj.z[bc .== 1] .= 0
        bc = vec(bc)
        obj.srho = 2400
        obj.udt = 1e5
        obj.bc = float(findall(bc .== 1))
        obj.k = zeros(Int(obj.ny),Int(obj.nx)) .+ 1e-6
        
        obj.ufunc = tempfunc 
        obj.firstcall = 1
        obj.l = 1
        obj.evaprate = 1e-30
        obj.maxareasinkfill = 1000000
        obj.precip = 1.0
        obj.kd = 1.0
        obj.ks = 0
        obj.dt = 1e6
        obj.t = 3e6
        obj.undercapacity = zeros(Int(obj.ny), Int(obj.nx)) 
        obj.k_sed = obj.k
        obj.flex = 20e3
        obj.dx = 1000
        obj.dy = 1000
        obj.m = .5
        obj.t_c = 100e6
        obj.n = 1
        obj.u = zeros(Int(obj.ny), Int(obj.nx)) .+ .00000000000000
        return obj
    end
    lemur_obj()


end

function set_param(obj::lemur_obj, nm::AbstractString, val::Any)
   # print(obj.k)
    setfield!(obj,Symbol(nm),val)
    m,n = size(val)
    if (m != obj.nx) | (n!= obj.ny)
        println("Warn: inconsistent size")
    end
    if all(size(obj.undercapacity) .!= size(val))
        println("Warn: inconsistent size in undercapacity")
    end
    if all(size(obj.k) .!= size(val))
        println("Warn: inconsistent size in undercapacity")
    end
    if all(size(obj.k) .!= size(val))
        println("Warn: inconsistent size in k")
    end
    if all(size(obj.bcx) .!= size(val))
        println("Warn: inconsistent size in bc")
    end
    if all(size(obj.z) .!= size(val))
        println("Warn: inconsistent size in bc")
    end
    return obj
end

module run_lemur

using Plots
using Cxx
using Libdl
include("/Users/gr_1/Documents/IsoFlex/IsoFlex.jl")



pw = "//Users/gr_1/Documents/lemur/src/"
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


mutable struct datas
    u::Array{Float64,3}
    z::Array{Float64,3}
    u2::Array{Float64,3}
    chi::Array{Float64,3}
    a::Array{Float64,3}
    function datas()
        obj = new()
    end
    datas()

end
using Statistics

function run(lemur_params; compute_sedflux = false, calc_chi = true)
    model = @cxxnew lemur(lemur_params.ny, lemur_params.nx)

    for nm = fieldnames(typeof(lemur_params))
        if (typeof(getfield(lemur_params,nm)) == Float64) || (typeof(getfield(lemur_params,nm)) == Array{Float64,1}) || (typeof(getfield(lemur_params,nm)) == Array{Float64,2})
            print(string(nm,' '))

            if length(getfield(lemur_params,nm)) <= 1
                println(getfield(lemur_params,nm))
            end
            set_lemur(model, string(nm), getfield(lemur_params,nm))
        end
    end
        
    
    print(lemur_params.k[1])
    z = get_lemur(model, "z", lemur_params.ny, lemur_params.nx)
    zi = zeros(size(z))
    flex = zeros(size(z))
    ero = zeros(size(z))
    z[vec(lemur_params.bcx .== 1)] .= 0
    set_lemur(model,"z", vec(z))
    @cxx model -> lakefill()
    z = copy(get_lemur(model, "z", lemur_params.ny, lemur_params.nx))

    zi[:] .= z[:]
    sedflux = nothing 
    data = datas()
    data.z = zeros(floor(Int16,lemur_params.ny),floor(Int16,lemur_params.nx),ceil(Int16,lemur_params.t/lemur_params.dt+1))
    data.u = zeros(floor(Int16,lemur_params.ny),floor(Int16,lemur_params.nx),ceil(Int16,lemur_params.t/lemur_params.dt+1))
    data.u2 = zeros(floor(Int16,lemur_params.ny),floor(Int16,lemur_params.nx),ceil(Int16,lemur_params.t/lemur_params.dt+1))
    data.chi = zeros(floor(Int16,lemur_params.ny),floor(Int16,lemur_params.nx),ceil(Int16,lemur_params.t/lemur_params.dt+1))
    data.a = zeros(floor(Int16,lemur_params.ny),floor(Int16,lemur_params.nx),ceil(Int16,lemur_params.t/lemur_params.dt+1))
    
    
    for t = 0:lemur_params.dt:lemur_params.t
        t2 = @time begin
            @cxx model -> erosion_fluvial()

            z =  copy(get_lemur(model, "z", lemur_params.ny, lemur_params.nx))
            z[lemur_params.bcx .== 1] .= 0
            ero .= copy(zi .- z)
            println(mean(vec(ero)))
            zi .= z


            z[lemur_params.bcx .== 0] .+= lemur_params.u[lemur_params.bcx .== 0]
            if lemur_params.t_c == 0
                u = IsoFlex.flexural(ero,dx=lemur_params.dx,dy=lemur_params.dy,Te=lemur_params.flex,buffer=200)
            else
                u,u2=IsoFlex.viscoelastic_lithos(ero; dx=lemur_params.dx,dt=lemur_params.dt, tt = lemur_params.t, t = t,
                dy=lemur_params.dy,Te=lemur_params.flex,buffer = 200, t_c = lemur_params.t_c)
            end
            #u2 = u;
            z[lemur_params.bcx .== 0] .+= lemur_params.u[lemur_params.bcx .== 0] *lemur_params.dt
            if !isnothing(lemur_params.ufunc(0))
                z .+= lemur_params.ufunc(t) * lemur_params.dt
            end

            set_lemur(model,"z", vec(z))
            

            @cxx model -> lakefill()
        end

        i = floor(Int16,t/lemur_params.dt+1)

        
        data.z[:,:,i] = z;
        data.u[:,:,i] = u;
        acc = get_lemur(model,"acc", lemur_params.ny, lemur_params.nx)
        data.a[:,:,i] = acc
            
        if lemur_params.t_c > 0
            data.u2[:,:,i] = u2;
        end
        if (calc_chi == true)
            chi =  zeros(Int(lemur_params.ny), Int(lemur_params.nx)) 
            I = get_lemur(model, "stack", lemur_params.ny, lemur_params.nx)
            R = get_lemur(model, "rec", lemur_params.ny, lemur_params.nx)
            println(I[end])

            for j =1:length(I)
                if I[j]!=R[j]
                    chi[Int(I[j])] = 1/acc[Int(I[j])] + chi[Int(R[Int(I[j])])]
                end
            end
            chi[acc .< 0] .= NaN
            
            data.chi[:,:,i] .= chi[:,:];
        end
        if compute_sedflux
            sed = copy(ero)
            I = get_lemur(model, "stack", lemur_params.ny, lemur_params.nx)
            R = get_lemur(model, "rec", lemur_params.ny, lemur_params.nx)
            println(I[end])

            for j =length(I):-1:1
                if I[j]!=R[j]
                    sed[Int(R[Int(I[j])])] += sed[Int(I[j])]
                end
            end
            
            data.u2[:,:,i] = sed

        end
    end
    println("here")
    return data,model;
end
end

