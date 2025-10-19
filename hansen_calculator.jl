using CSV, DataFrames, SpecialFunctions

const NEWCOMB_FILE = "n_cache.csv"
const HENSEN_FILE = "h_cache.csv"


# === Load CSV into memory ===
function load_newcomb_cache(path::String)
    if isfile(path)
        df = CSV.read(path, DataFrame)
        return Dict{NTuple{4, Int}, Float64}((row.n, row.m, row.rho, row.sigma) => row.value for row in eachrow(df))
    else
        return Dict{NTuple{4, Int}, Float64}()
    end
end

# Load existing Hansen coefficients from file
function load_hansen_csv(path::String)
    if !isfile(path)
        return Dict{Tuple{Int, Int, Int, Float64}, Float64}()
    end
    df = CSV.read(path, DataFrame)
    return Dict((row.n, row.m, row.k, row.e) => row.value for row in eachrow(df))
end

newcomb_cache = load_newcomb_cache(NEWCOMB_FILE)
newcomb_new_entries = Dict{NTuple{4, Int}, Float64}()

hansen_cache = load_hansen_csv(HENSEN_FILE)
hansen_new_entries = Dict{Tuple{Int, Int, Int, Float64}, Float64}()


# === Save to file ===
function append_newcomb_to_csv(path::String, new_data::Dict{NTuple{4, Int}, Float64})
    write_header = !isfile(path)
    open(path, "a") do io
        if write_header
            println(io, "n,m,rho,sigma,value")
        end
        for ((n, m, ρ, σ), val) in new_data
            println(io, "$n,$m,$ρ,$σ,$val")
        end
    end
end


# Append only missing entries from `new_data` to file
function append_hansen_to_csv(path::String, new_data::Dict{Tuple{Int, Int, Int, Float64}, Float64})
    write_header = !isfile(path)
    open(path, "a") do io
        if write_header
            println(io, "n,m,k,e,value")
        end
        for ((n, m, k, e), val) in new_data
            println(io, "$n,$m,$k,$e,$val")
        end
    end
end



# === Main recursive + persistent function ===
function newcomb(n::Int, m::Int, ρ::Int, σ::Int;
                 cache=newcomb_cache, buffer=newcomb_new_entries)

    key = (n, m, ρ, σ)
    if haskey(cache, key)
        return cache[key]
    end

    # Base cases (define explicitly!)
    if ρ < 0 || σ < 0
        return 0.0
    elseif ρ == 0 && σ == 0
        return 1.0
    elseif σ>ρ
        return newcomb(n, -m, σ, ρ)
    elseif  ρ == 1 && σ == 0
        return m-n/2
    elseif σ == 0
        val = (2*(2*m-n)*newcomb(n,m+1,ρ-1,0) + 
                (m-n)*newcomb(n,m+2,ρ-2,0))  /4/ρ
    else
        val = -2*(2*m+n)*newcomb(n,m-1,ρ,σ-1) -
                (m+n)*newcomb(n,m-2,ρ,σ-2) -
                (ρ-5*σ+4+4*m+n)*newcomb(n,m,ρ-1,σ-1)
        for j in 2:min(ρ, σ)
            val += 2*(ρ-σ+m)*(-1)^j*gamma(5/2)/gamma(j+1)/gamma(5/2-j)*
                    newcomb(n,m,ρ-j,σ-j)
        end
        val = val/4/σ    
    end
    # Cache result
    cache[key] = val
    buffer[key] = val
    return val
end

function hensen(n::Int, m::Int, k::Int, e::Float64; 
    hcache=hansen_cache, hbuffer=hansen_new_entries,
    ncache=newcomb_cache, nbuffer=newcomb_new_entries,
    max_terms=50)

    key = (n, m, k, e)
    # If already computed, return cached value
    if haskey(hcache, key) # for some reason was cache
        return hcache[key]
    end

    if k == 0
        if n == 0 && m == 0
            return 1.0
        elseif n == 1 && m == 0
            hensen = 1+e^2/2
        elseif n == 1 && m == 1
            hensen = -3*e/2
        elseif n == 2 && m == 0
            hensen = 1 + 3/2*e^2 
        elseif n == 2 && m == 1
            hensen = -2*e - 1/2*e^3
        elseif m < 0
            hensen = hensen(n,-m,k,e)
        elseif m > n
            hensen = 1/(n-m+1)/e*(e*(n+m+1)*hensen(n,m-2,k,e)+2*m*hensen(n,m-1,k,e))
        else
            hensen = (2*n+3)/(n+2)*hensen(n-1,m,k,e)-(n+1-m)*(n+1+m)/(n+2)/(n+1)*(1-e^2)*hensen(n-2,m,k,e)
        end
    else
        hensen = 0.0
        α = max(0, k - m)
        β = max(0, m - k)
        for σ in 0:max_terms
            hensen += newcomb(n,m,σ+α,σ+β) * e^(2*σ)
        end
        hensen = hensen * e^abs(m-k)
        
    end
    hcache[key] = hensen
    hbuffer[key] = hensen
    return hensen
end




# 
for ℓ in 0:20
    for m in -ℓ:ℓ
        println("l: $ℓ, m: $m")
        for e in 0.05:0.05:0.9
            hensen(ℓ,m,0, e)
        end
    end
end
# for j in -25:25
#     if j == 0
#         continue
#     end
#     for l in 0:20
#         println("j: $j, l: $l")
#         for m in -l:l
#             # if m == 0
#             #     continue
#             # end
#             for q in 0:90
#                 for e in 0.05:0.05:0.9
#                     hensen(2*q+l,m, j, e)
#                     hensen(l,m,q, e)
#                     hensen(l,m,-q, e)
#                     hensen(-l-1,m,q, e)
#                     hensen(-l-1,m,-q, e)
#                 end
#             end
#         end
#     end
# end


append_newcomb_to_csv(NEWCOMB_FILE, newcomb_new_entries)
append_hansen_to_csv(HENSEN_FILE, hansen_new_entries)

println(hensen(-3,6,2,0.24905))
println(newcomb(2, 2, 0, 0))