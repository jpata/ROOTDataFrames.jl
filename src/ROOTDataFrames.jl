module ROOTDataFrames

using ROOT, DataFrames, DataArrays
import Base.length, Base.getindex

import DataFrames.nrow, DataFrames.size

immutable TreeDataFrame <: AbstractDataFrame
    tf::TObjectA
    tt::TTreeA
    bvars::Vector{Any}
    index::DataFrames.Index
    types::Vector{Type}
end

function TreeDataFrame(fns::AbstractVector, treename="dataframe")

    tch = TChain(treename)

    for fn in fns
        AddFile(tch, convert(ASCIIString, fn), -1)
    end
    #tf = TFile(fn)
    #tt = root_cast(TTree, Get(root_cast(ROOT.TDirectory, tf), "dataframe"))

    brs = GetListOfBranches(tch);
    brs = [root_cast(TBranch, brs[i]) for i=1:length(brs)];

    bd = Dict()
    bd_isna = Dict()
    for b in brs
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
    println(collect(keys(bd)))
    for k in keys(bd)
        leaves = GetListOfLeaves(bd[k])
        length(leaves)==1 || error("$k, Nleaf=$(length(leaves))")
        leaf = root_cast(TLeaf, leaves[1])
        t = GetTypeName(leaf)|>bytestring|>parse
        t = eval(ROOT.type_replacement[t])::Type
        push!(types, t)

        bvar = (t[0.0], Bool[true]);
        push!(bvars, bvar)

        #does not work reliably for TChain
        #br = GetBranch(tch, k)
        #SetAddress(br, convert(Ptr{Void}, bvar[1]))
        #br = GetBranch(tch, "$(k)_ISNA")
        #SetAddress(br, convert(Ptr{Void}, bvar[2]))


        SetBranchAddress(tch, k, convert(Ptr{Void}, bvar[1]))
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

Base.length(t::TreeDataFrame) = GetEntries(t.tt)
Base.size(df::TreeDataFrame) = (nrow(df), ncol(df))
Base.size(df::TreeDataFrame, n) = size(df)[n]

function load_row(df::TreeDataFrame, i::Integer)
    return GetEvent(df.tt, i-1)
end

function Base.getindex(df::TreeDataFrame, i::Int64, s::Symbol, get_entry=false)
    if get_entry
        load_row(df, i)
    end
    v, na = df.bvars[df.index[s]]
    return na[1] ? NA : v[1]
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
    SetCacheSize(df.tt, 0)
    SetCacheSize(df.tt, 256 * 1024 * 1024)
    set_branch_status!(df, "*", false)
    for b in brs
        set_branch_status!(df, "$b", true)
        AddBranchToCache(df.tt, "$b")
    end
end

function Base.getindex(df::TreeDataFrame, mask::AbstractVector, ss::AbstractVector{Symbol})
    length(mask) == nrow(df) || error("mask=$(length(mask)), nrow=$(nrow(df))")
    enable_branches(df, ["$(s)*" for s in ss])
    names_types = {
        n => df.types[df.index[n]] for n in names(df)
    }

    const n = sum(mask)

    const ret = DataFrame(
        { DataArray(names_types[s], n) for s in ss},
        DataFrames.Index(ss)
    )
    j = 1
    for i=1:nrow(df)
        (!isna(mask[i]) && mask[i]) || continue
        load_row(df, i)
        for nn in ss
            ret[j, nn] = df[j, nn]
        end
        j += 1
    end
    return ret
end

Base.getindex(df::TreeDataFrame, s::Symbol) = df[[s]][s]
Base.getindex(df::TreeDataFrame, ss::AbstractVector{Symbol}) = df[[true for i=1:nrow(df)], ss]
Base.getindex(df::TreeDataFrame, mask::AbstractVector, s::Symbol) = df[mask, [s]][s]

function writetree(fn, df::AbstractDataFrame;progress=true, treename="dataframe")
    tf = TFile(fn, "RECREATE")
    const tree = TTree(
        treename, treename
    )

    bnames = Symbol[]
    btypes = Type[]
    bvars = Any[]
    bidx = Dict{Symbol,Union(Real,AbstractArray{Real,1})}()

    nb = 1
    for (cn, ct) in zip(names(df), eltypes(df))
        push!(bnames, cn)
        push!(btypes, ct)

        bv = ct[0]
        push!(bvars, bv)

        br = Branch(
            tree, string(cn),
            convert(Ptr{None}, pointer(bv)),
            "$cn/$(SHORT_TYPEMAP[ct])"
        )
        bidx[symbol(cn)] = nb

        cn_na = symbol("$(cn)_ISNA")
        bv = Bool[true]
        push!(bvars, bv)

        br = Branch(
            tree, string(cn_na),
            convert(Ptr{None}, pointer(bv)),
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

    Write(tree)
    Close(tf)
end

export writetree, TreeDataFrame
export load_row
end # module
