module ROOTDataFrames

using ROOT, DataFrames, DataArrays
import Base.length, Base.getindex

import DataFrames.nrow, DataFrames.size, DataFrames.index

type TreeDataFrame <: AbstractDataFrame
    tf::TObjectA
    tt::TTreeA
    bvars::Vector{Any}
    index::DataFrames.Index
    types::Vector{Type}
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
    #println(collect(keys(bd)))
    for k in keys(bd)
        #println("branch $k")
        leaves = GetListOfLeaves(bd[k])
        if length(leaves)!=1
            warn("$k, Nleaf=$(length(leaves)), skipping")
            continue
        end
        leaf = root_cast(TLeaf, leaves[1])
        t = GetTypeName(leaf)|>bytestring|>parse
        if !haskey(ROOT.type_replacement, t)
            warn("branch $k with type $(bytestring(GetTypeName(leaf))) does not have a julia typle replacement") 
            continue
        end
        #println("type $t")
        t = eval(ROOT.type_replacement[t])::Type
        push!(types, t)

        bvar = (t[0.0], Bool[true]);
        push!(bvars, bvar)

        #does not work reliably for TChain
        #br = GetBranch(tch, k)
        #SetAddress(br, convert(Ptr{Void}, bvar[1]))
        #br = GetBranch(tch, "$(k)_ISNA")
        #SetAddress(br, convert(Ptr{Void}, bvar[2]))

        #print(bvar[1])
        SetBranchAddress(tch, k, convert(Ptr{Void}, pointer(bvar[1])))
        #SetBranchAddress(tch, k, bvar[1])
        if haskey(bd_isna, k)
            #println("NA branch activated for $k")
            SetBranchAddress(tch, "$(k)_ISNA", convert(Ptr{Void}, bvar[2]))
            bvar[2][1] = true
        else
            #println("NA branch de-activated for $k")
            bvar[2][1] = false
        end
        bidx[symbol(k)] = bridx

        bridx += 1
    end

    idx = DataFrames.Index(bidx, collect(keys(bidx)))

    TreeDataFrame(TObject(C_NULL), tch, bvars, idx, types)
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

function Base.getindex(df::TreeDataFrame, i::Int64, s::Symbol, get_entry::Bool=false)
    if get_entry
        load_row(df, i)
    end
    v::Vector{Any}, na::Vector{Bool} = df.bvars[df.index[s]]
    return na[1] ? NA : deepcopy(v[1])
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
    names_types = Dict{Any,Any}([n=>df.types[df.index[n]] for n in names(df)])

    const n = sum(mask)

    #need to disable GC to prevent silent corruption and crash
    #still un-understood 
    #gc_disable()
    
    const ret = DataFrame(
        Any[DataFrames.DataArray(names_types[s], n) for s in ss],
        DataFrames.Index(ss)
    )
    j = 1


    t0 = time()
    for i=1:nrow(df)
        (!isna(mask[i]) && mask[i]) || continue
        nloaded = load_row(df, i)
        #println(i, " ", nloaded)
        i%50000 == 0 && print(".")
        for nn in ss
            ret[j, nn] = df[j, nn]
        end
        j += 1
    end
    t1 = time()
    println()
    println(nrow(df)/(t1-t0))
    #gc_enable()
    return ret
end

Base.getindex(df::TreeDataFrame, s::Symbol) = df[[s]][s]
Base.getindex(df::TreeDataFrame, ss::AbstractVector{Symbol}) = df[[true for i=1:nrow(df)], ss]
Base.getindex(df::TreeDataFrame, mask::AbstractVector, s::Symbol) = df[mask, [s]][s]

function TreeDataFrame(fn, ns::AbstractVector, types::AbstractVector; treename="dataframe")
    tf = TFile(fn, "RECREATE")
    const tree = TTree(
        treename, treename
    )

    bnames = Symbol[]
    btypes = Type[]
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
        btypes
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
