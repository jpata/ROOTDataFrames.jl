module ROOTDataFrames

const MAX_SIZE = 500
using ROOT, DataFrames, DataArrays
import Base.length, Base.getindex

import DataFrames.nrow, DataFrames.size, DataFrames.index

type TreeDataFrame <: AbstractDataFrame
    tf::TObjectA
    tt::TTreeA
    bvars::Vector{Any}
    index::DataFrames.Index
    types::Vector{Type}
    leafsizes::Vector{Any}
end

import DataFrames.showcols
function showcols(io::IO, df::TreeDataFrame)
    println(io, summary(df))
    metadata = DataFrame(Name = names(df),
        Eltype = eltypes(df),
        #Missing = colmissing(df)
    )
    showall(io, metadata, true, symbol("Col #"), false)
    return
end
import DataFrames.show
function show(io::IO, df::TreeDataFrame)
    showcols(io, df)
    return
end


function TreeDataFrame(fns::AbstractVector, treename="dataframe")

    tch = TChain(treename)

    for fn in fns
        #println("adding file $fn")
        AddFile(tch, convert(ASCIIString, fn), -1)
    end
    #println("TChain $tch created")
    #tf = TFile(fn)
    #tt = root_cast(TTree, Get(root_cast(ROOT.TDirectory, tf), "dataframe"))

    brs = GetListOfBranches(tch);
    brs = [root_cast(TBranch, brs[i]) for i=1:length(brs)];
    #println("TTree has $(length(brs)) branches")
    bd = Dict()
    bd_isna = Dict()
    for b in brs
        #println("creating branch $b")
        bn = GetName(root_cast(TObject, b)) |> bytestring;
        if contains(bn, "ISNA")
            bd_isna[join(split(bn, "_")[1:end-1], "_")] = b
        else
            bd[bn] = b
        end
    end

    bridx = 1
    bvars = Any[]
    bidx = Dict{Symbol,Union(Real,AbstractArray{Real,1})}()
    types = Type[]
    leafsizes = Any[]
    #println(collect(keys(bd)))
    for k in keys(bd)
        #println("branch $k")
        leaves = GetListOfLeaves(bd[k])
        if length(leaves)!=1
            #warn("$k, nleaf=$(length(leaves))")
            continue
        end
        leaf = root_cast(TLeaf, leaves[1])
        leaf_staticsize = ROOT.GetLenStatic(leaf)
        leafsize = 1
        if leaf_staticsize!=1
            leafsize = leaf_staticsize 
            #warn("$k, nleaf size == $(leaf_staticsize)")
            #continue
        end

        lc = ROOT.GetLeafCount(leaf)
        if lc != C_NULL
            lc = TLeaf(lc)
            bname = lc |> GetName |> bytestring
            leafsize = symbol(bname)
            #warn("$k, nleaf size dynamic $bname")
            #continue
        end

        leaf = root_cast(TLeaf, leaves[1])
        t = GetTypeName(leaf)|>bytestring|>parse
        if !haskey(ROOT.type_replacement, t)
            warn("branch $k with type $(bytestring(GetTypeName(leaf))) does not have a julia typle replacement") 
            continue
        end
        t = eval(ROOT.type_replacement[t])::Type
        push!(types, t)

        push!(leafsizes, leafsize)

        _leafsize = isa(leafsize, Real) ? leafsize : MAX_SIZE 

        bvar = (zeros(t, _leafsize), Bool[true for i=1:_leafsize])
        push!(bvars, bvar)

        #does not work reliably for TChain
        #br = GetBranch(tch, k)
        #SetAddress(br, convert(Ptr{Void}, bvar[1]))
        #br = GetBranch(tch, "$(k)_ISNA")
        #SetAddress(br, convert(Ptr{Void}, bvar[2]))

        #println(bvar[1])
        SetBranchAddress(tch, k, convert(Ptr{Void}, pointer(bvar[1])))
        #SetBranchAddress(tch, k, bvar[1])
        if haskey(bd_isna, k)
            #println("NA branch activated for $k")
            SetBranchAddress(tch, "$(k)_ISNA", convert(Ptr{Void}, bvar[2]))
            bvar[2][:] = true
        else
            #println("NA branch de-activated for $k")
            bvar[2][:] = false
        end
        bidx[symbol(k)] = bridx

        bridx += 1
    end

    idx = DataFrames.Index(bidx, collect(keys(bidx)))

    TreeDataFrame(TObject(C_NULL), tch, bvars, idx, types, leafsizes)
end

TreeDataFrame(fn::String) = TreeDataFrame([fn])

function Base.length(t::TreeDataFrame)
    if t.tt != C_NULL
        return GetEntries(t.tt)
    else
        return 0
    end
end
Base.size(df::TreeDataFrame) = (nrow(df), ncol(df))
Base.size(df::TreeDataFrame, n) = size(df)[n]

function load_row(df::TreeDataFrame, i::Integer)
    return GetEvent(df.tt, i-1)
end

function Base.getindex{T <: Any}(df::TreeDataFrame, i::Int64, s::Symbol, t::Type{T}=Any, get_entry::Bool=false)
    if get_entry
        load_row(df, i)
    end
    v::Vector{T}, na::Vector{Bool} = df.bvars[df.index[s]]
    ret = DataArray(v, na)
    _size = df.leafsizes[df.index[s]]
    if isa(_size, Symbol)
        _size = df[i, _size, df.types[df.index[_size]]]
    end
    return length(ret)>1 ? ret[1:_size] : first(ret)
end

import DataFrames.nrow, DataFrames.ncol
DataFrames.nrow(df::TreeDataFrame) = length(df)
DataFrames.ncol(df::TreeDataFrame) = length(df.bvars)

import Base.names
Base.names(x::TreeDataFrame) = names(x.index)

import DataFrames.eltypes
DataFrames.eltypes(x::TreeDataFrame) = x.types

set_branch_status!(df, pat, status) = SetBranchStatus(
    df.tt, pat, status, convert(Ptr{Cuint}, 0)
)

function enable_branches(df, brs)
    #println("enabling branches $brs")
    #SetCacheSize(df.tt, 0)
    #SetCacheSize(df.tt, 256 * 1024 * 1024)
    set_branch_status!(df, "*", false)
    for b in brs
        set_branch_status!(df, "$b", true)
        #AddBranchToCache(df.tt, "$b") #cache does not work properly, causes crash
    end
end

function Base.getindex(df::TreeDataFrame, mask::AbstractVector, ss::AbstractVector{Symbol})
    length(mask) == nrow(df) || error("mask=$(length(mask)), nrow=$(nrow(df))")
    enable_branches(df, ["$(s)*" for s in ss])
    names_types = Dict{Symbol, DataType}([n=>df.types[df.index[n]] for n in names(df)])
    sizes = Dict{Symbol, Any}([n=>df.leafsizes[df.index[n]] for n in names(df)])

    const n = sum(mask)

    #need to disable GC to prevent silent corruption and crash
    #still un-understood 
    #gc_disable()
   
    arrlist = Any[]
    for s in ss
        _size = sizes[s]
        _type = names_types[s]
        if typeof(_size) <: Symbol
            _size = MAX_SIZE
        end
        if _size > 1
            _type = DataVector{_type}
        end
        a = DataFrames.DataArray(_type, n)
        push!(arrlist, a)
    end
    const ret = DataFrame(
        arrlist,
        DataFrames.Index(ss)
    )
    j = 1


    t0 = time()
    for i=1:nrow(df)
        (!isna(mask[i]) && mask[i]) || continue
        nloaded = load_row(df, i)
        for nn in ss
            x = df[j, nn, names_types[nn]]
            ret[j, nn] = x
        end
        j += 1
    end
    t1 = time()
    #gc_enable()
    return ret
end

Base.getindex(df::TreeDataFrame, s::Symbol) = df[[s]][s]
Base.getindex(df::TreeDataFrame, ss::AbstractVector{Symbol}) = df[[true for i=1:nrow(df)], ss]
Base.getindex(df::TreeDataFrame, inds::AbstractVector{Int64}, ss::AbstractVector{Symbol}) = df[[i in inds for i=1:nrow(df)], ss]
Base.getindex(df::TreeDataFrame, i::Int64, ss::AbstractVector{Symbol}) = df[[i], ss]
Base.getindex(df::TreeDataFrame, mask::AbstractVector, s::Symbol) = df[mask, [s]][s]

function TreeDataFrame(fn, ns::AbstractVector, types::AbstractVector; treename="dataframe")
    tf = TFile(fn, "RECREATE")
    const tree = TTree(
        treename, treename
    )

    bnames = Symbol[]
    btypes = Type[]
    leafsizes = Any[]
    bvars = Any[]
    bidx = Dict{Symbol,Union(Real,AbstractArray{Real,1})}()

    nb = 1
    for (cn, ct) in zip(ns, types)
        push!(bnames, convert(Symbol, cn))
        push!(btypes, convert(Type, ct))

        bv = ct[0]
        push!(bvars, bv)

        br = Branch(
            tree, string(cn),
            convert(Ptr{Void}, pointer(bv)),
            "$cn/$(SHORT_TYPEMAP[ct])"
        )
        bidx[symbol(cn)] = nb

        cn_na = symbol("$(cn)_ISNA")
        bv = Bool[true]
        push!(bvars, bv)
        push!(leafsizes, 1)

        br = Branch(
            tree, string(cn_na),
            convert(Ptr{Void}, pointer(bv)),
            "$cn_na/O"
        )
        bidx[cn_na] = nb

        nb += 2
    end

    dtf = TreeDataFrame(
        tf, tree, bvars,
        DataFrames.Index(bidx, collect(keys(bidx))),
        btypes,
        leafsizes
    )
end

function writetree(fn, df::AbstractDataFrame;progress=true, treename="dataframe")
    dtf = TreeDataFrame(fn, names(df), eltypes(df);treename=treename)

    for i=1:nrow(df)
        for j=1:ncol(df)
            const nc = 2 * j - 1
            const nc_isna = 2 * j
            if !isna(df[i, j])
                dtf.bvars[nc_isna][1] = false
                dtf.bvars[nc][1] = df[i, j]
            else
                dtf.bvars[nc][1] = convert(dtf.types[j], 0)
                dtf.bvars[nc_isna][1] = true
            end
        end

        Fill(dtf.tt)
    end

    Write(dtf.tt)
    Close(dtf.tf)
end

index(df::TreeDataFrame) = df.index

function writetree_temp(outfile, df::DataFrame)
    tempf = mktemp()[1]
    print("writing to $tempf...");writetree(tempf, df);println("done")
    
    for i=1:5
        try
            println("cleaning $outfile...");isfile(outfile) && rm(outfile)
            println("copying...");cp(tempf, outfile)
            s = stat(outfile)
            s.size == 0 && error("file corrupted")
            break
        catch err
            warn("$err: retrying after sleep")
            sleep(5)
        end
    end
end

export writetree, TreeDataFrame
export writetree_temp
export load_row
export enable_branches
end # module
