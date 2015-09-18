module ROOTDataFrames

include("log.jl")

const MAX_SIZE = 50
const DEBUG = false

using ROOT, DataFrames, DataArrays
import Base.length, Base.getindex

import DataFrames.nrow, DataFrames.size, DataFrames.index

const LIBROOT = "/Users/joosep/.julia/v0.4/ROOT/libroot"

function SetBranchAddress{T <: Real}(__obj::TTreeA, bname::ASCIIString, add::Vector{T}, p::Ptr{Void}=C_NULL)
    ccall(("TTree_SetBranchAddress1",LIBROOT), Int32, (Ptr{Void}, Ptr{UInt8}, Ref{T}, Ptr{Ptr{TBranch}}), __obj.p, bname, add, p)
end

type BranchValue{T, N}
    data::Vector{T}
    size::Symbol
    branch::TBranch
end

BranchValue{T,N <: Symbol}(::Type{T}, s::N, br::TBranch=TBranch(C_NULL)) =
    BranchValue{T,N}(zeros(T, MAX_SIZE), s, br)
BranchValue{T,N}(::Type{T}, s::Type{Val{N}}, br::TBranch=TBranch(C_NULL)) =
    BranchValue{T,N}(zeros(T, N), :nothing, br::TBranch)

getsize{T, N}(bval::BranchValue{T, N}) = N
getsize{T, N <: Symbol}(bval::BranchValue{T, N}) = bval.size

import Base.eltype
# eltype{T, N <: Symbol}(bval::BranchValue{T, N}) = Vector{T}
# eltype{T, N}(bval::BranchValue{T, N}) = T
# eltype{T}(bval::BranchValue{T, 1}) = T
eltype{T, N}(bval::BranchValue{T, N}) = T

datatype{T, N <: Symbol}(bval::BranchValue{T, N}) = Vector{T}
datatype{T, N}(bval::BranchValue{T, N}) = Vector{T}
datatype{T}(bval::BranchValue{T, 1}) = T

value{T, N <: Symbol}(bval::BranchValue{T, N}) = bval.data::Vector{T}
value{T, N}(bval::BranchValue{T, N}) = bval.data::Vector{T}
value{T}(bval::BranchValue{T, 1}) = bval.data[1]::T

import Base.call
call{T}(bval::BranchValue{T, 1}) = bval.data[1]::T
call{T, N}(bval::BranchValue{T, N}) = bval.data[1:N]::Vector{T}
call{T, N <: Symbol}(bval::BranchValue{T, N}) = bval.data::Vector{T}

#export value
function makeclass(
    names::Vector{Symbol},
    types::Vector{BranchValue},
    )
    typename = gensym("RowData")
    ret = :(
        type $typename
        
        $typename() = new()
        end
    )
    #members
    for (name, typ) in zip(names, types)
        push!(ret.args[3].args, Expr(:(::), name, typeof(typ)))
    end
    
    #constructor
    for (name, typ) in zip(names, types)
        leafsize = getsize(typ)
        if typeof(leafsize) <: Real #make a scalar
            push!(ret.args[3].args[2].args[2].args[2].args,
                Expr(:call, :BranchValue, eltype(typ), Val{leafsize})
            )
        else #variable length branch
            push!(ret.args[3].args[2].args[2].args[2].args,
                Expr(:call, :BranchValue, eltype(typ), QuoteNode(leafsize))
            )
        end
    end
    eval(ret)
    
    return eval(typename)
end

type TreeDataFrame{T} <: AbstractDataFrame
    tf::TObjectA
    tt::TTreeA
    row::T
    rowindex::Int64
    treeindex::Int64
end

import DataFrames.showcols
function showcols(io::IO, df::TreeDataFrame)
    println(io, summary(df))
    metadata = DataFrame(
        Name = names(df),
        Eltype = eltypes(df),
        sizes = [getsize(getfield(df.row, n)) for n in names(df)],
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


function update_branches!{T}(df::TreeDataFrame{T})
    for field in fieldnames(df.row)
        br = GetBranch(df.tt, string(field))
        getfield(df.row, field).branch = br
    end
end

function TreeDataFrame(fns::AbstractVector; treename="dataframe")

    @debugprint "opening $fns:$treename"
    #DEBUG && println("[TreeDataFrame] opening $fns:$treename")
    tch = TChain(treename)

    for fn in fns
        @debugprint "adding file $fn"
        #-1 -> load number of entries from file
        AddFile(tch, convert(ASCIIString, fn), -1)
        nentries = GetEntries(tch)
        @debugprint "N=$nentries"
    end
    @debugprint "TChain $tch created"
    #tf = TFile(fn)
    #tt = root_cast(TTree, Get(root_cast(ROOT.TDirectory, tf), "dataframe"))

    branches = GetListOfBranches(tch);
    branches = [root_cast(TBranch, branches[i]) for i=1:length(branches)];
    @debugprint "TTree has $(length(branches)) branches"
    branches_d = Dict()
    branch_names = Symbol[]
    branch_types = BranchValue[]
    for b::TBranch in branches
        #println("creating branch $b")
        @assert !is_null(b)
        name = GetName(root_cast(TObject, b)) |> bytestring |> symbol;
        branches_d[name] = b
        push!(branch_names, name)
    end

    row_types = Dict{Symbol, BranchValue}()
    for name in branch_names
        #println("branch $k")
        leaves = GetListOfLeaves(branches_d[name])
        if length(leaves)!=1
            warn("branch=$k, nleaf=$(length(leaves)), skipping")
            continue
        end
        
        leaf = root_cast(TLeaf, leaves[1])
        #Leaf data type in ROOT types (Double_t etc)
        leaftype = GetTypeName(leaf)|>bytestring|>parse

        leaf_staticsize = ROOT.GetLenStatic(leaf)
        leafsize = 1
        if leaf_staticsize!=1
            leafsize = leaf_staticsize 
            #warn("$name, nleaf size == $(leaf_staticsize)")
            #continue
        end

        #Branch is variable size, size is determined by another branch
        leafcount = ROOT.GetLeafCount(leaf)
        if leafcount != C_NULL
            leafcount = TLeaf(leafcount)
            leafname = leafcount |> GetName |> bytestring
            leafsize = symbol(leafname)
            #warn("$name, nleaf size dynamic $bname")
            #continue
        end
        
        if !haskey(ROOT.type_replacement, leaftype)
            warn("branch=$name with type=$leaftype does not have a julia type replacement, skipping") 
            continue
        end

        #convert leaf type to julia type
        leaftype = eval(ROOT.type_replacement[leaftype])::Type

        if typeof(leafsize) <: Integer
            leafsize = Val{leafsize}
        end

        push!(
            branch_types,
            BranchValue(
                leaftype, leafsize, branches_d[name]
            )
        )
    end

    rowtype = makeclass(branch_names, branch_types)
    rowdata = rowtype()
    
    for field in fieldnames(rowdata)
        typ = fieldtype(rowtype, field)
        #getfield(rowdata, field).branch = branches_d[field]
        SetBranchAddress(tch, string(field), getfield(rowdata, field).data)
    end

    df = TreeDataFrame{rowtype}(
        TObject(C_NULL),
        tch,
        rowdata,
        1,
        1
    )
    update_branches!(df)
    return df
end

TreeDataFrame(fn::String) = TreeDataFrame([fn])

function Base.length(t::TreeDataFrame)
    if t.tt != C_NULL
        @assert !is_null(t.tt)
        return GetEntries(t.tt)
    else
        return 0
    end
end
Base.size(df::TreeDataFrame) = (nrow(df), ncol(df))
Base.size(df::TreeDataFrame, n) = size(df)[n]

function load_row(df::TreeDataFrame, i::Integer)
    localentry = LoadTree(df.tt, i-1)
    localentry >= 0 || error("could not get entry $i")
    return GetEntry(df.tt, localentry)
end

function update!(df::TreeDataFrame, i::Int)
    treeindex = GetTreeNumber(df.tt) + 1
    if treeindex != df.treeindex
        update_branches!(df)
        df.treeindex = treeindex
    end
    df.rowindex = i
end

function load_row(df::TreeDataFrame, i::Integer, cols::Vector{Symbol})
    @assert !is_null(df.tt)

    localentry = LoadTree(df.tt, i-1)
    localentry >= 0 || error("could not get entry $i")
    update!(df, i)

    ntot = 0
    for col in cols
        br = getfield(df.row, col).branch
        @assert !is_null(br)
        ntot += GetEntry(br, localentry)
    end
    return ntot
end

function load_row(df::TreeDataFrame, i::Integer)
    localentry = LoadTree(df.tt, i-1)
    localentry >= 0 || error("could not get entry $i")
    update!(df, i)
    return GetEntry(df.tt, localentry)
end

function Base.getindex{T <: Any, X <: Any}(
    df::TreeDataFrame{T},
    i::Int64,
    s::Symbol,
    getentry=false,
    ::Type{BranchValue{X, 1}}=fieldtype(T, s),
    )
    if getentry
        load_row(df, i, [s])>0 || error("could not load row $i")
    end
    ret = getfield(df.row, s)
    return value(ret)::X
end

function Base.getindex{T <: Any, R <: Symbol, X <: Any}(
    df::TreeDataFrame{T},
    i::Int64,
    s::Symbol,
    _size::Int64=-1,
    getentry::Bool=false,
    ::Type{BranchValue{X, R}}=fieldtype(T, s),
    )
    if _size == -1
        sizebranch = getsize(getfield(df.row, s))
        if getentry
            load_row(df, i) > 0 || error("could not load row $i")
        end
        _size = value(getfield(df.row, sizebranch))
    end
    ret = getfield(df.row, s)

    return value(ret)[1:_size]
end

import DataFrames.nrow, DataFrames.ncol
DataFrames.nrow(df::TreeDataFrame) = length(df)
DataFrames.ncol(df::TreeDataFrame) = length(names(df))

import Base.names
Base.names(df::TreeDataFrame) = fieldnames(df.row)

import DataFrames.eltypes
DataFrames.eltypes(df::TreeDataFrame) = Type[datatype(getfield(df.row, n)) for n in names(df)]

set_branch_status!(df, pattern, status) = SetBranchStatus(
    df.tt, pattern, status, convert(Ptr{Cuint}, 0)
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

function Base.getindex{T}(
    df::TreeDataFrame{T},
    mask::AbstractVector,
    cols::AbstractVector{Symbol}
    )
    length(mask) == nrow(df) || error("mask=$(length(mask)), nrow=$(nrow(df))")
    enable_branches(df, ["$(s)*" for s in cols])

    #prepare output array
    names_types = Dict{Symbol, DataType}([n => eltype(getfield(df.row, n)) for n in names(df)])
    sizes = Dict{Symbol, Any}([n => getsize(getfield(df.row, n)) for n in names(df)])

    #output vector length
    const n = sum(mask)

    #need to disable GC to prevent silent corruption and crash
    #still un-understood 
    #gc_disable()
   
    arrlist = Any[]
    for col in cols
        _size = sizes[col]
        _type = names_types[col]
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
        DataFrames.Index(cols)
    )
    j = 1

    t0 = time()
    for i=1:nrow(df)
        (!isna(mask[i]) && mask[i]) || continue
        nloaded = load_row(df, i, cols)
        for col in cols
            ret[j, col] = value(getfield(df.row, col))
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

TreeDataFrame(
    fn::ASCIIString,
    names::Vector{Symbol},
    types::Vector{Type}
    ;kwargs...
    ) = TreeDataFrame(
        fn, names,
        BranchValue[BranchValue(dt, Val{1}) for dt in types]; kwargs...
)

function TreeDataFrame(
    fn::ASCIIString,
    names::Vector{Symbol},
    types::Vector{BranchValue}
    ;treename="dataframe"
    )

    rowdatatype = makeclass(names, types)
    rowdata = rowdatatype()

    tf = TFile(fn, "RECREATE")

    #make a TTree pointer
    const tree = TTree(
        treename, treename
    )
    
    for (ibranch, (name, typ)) in enumerate(
        zip(names, types)
        )
        
        eltyp = eltype(typ)
        br = Branch(
            tree, string(name),
            convert(Ptr{Void}, pointer(getfield(rowdata, name).data)),
            "$name/$(SHORT_TYPEMAP[eltyp])"
        )
        getfield(rowdata, name).branch = br
    end

    dtf = TreeDataFrame(
        tf,
        tree,
        rowdata,
        1,
        1
    )
end

import Base.setindex!
function Base.setindex!{T, K <: Real}(
    df::TreeDataFrame{T}, val::K, i::Integer, col::Symbol
    )
    getfield(df.row, col).data[1] = val
end

function Base.setindex!{T, K <: Real}(
    df::TreeDataFrame{T}, val::Vector{K}, i::Integer, col::Symbol
    )
    getfield(df.row, col).data[:] = val[:]
end

function writetree(fn::ASCIIString, df::AbstractDataFrame
    ;progress=true, treename="dataframe")
    dtf = TreeDataFrame(fn, names(df), eltypes(df);treename=treename)

    colnames = names(df)
    for i=1:nrow(df)
        for j=1:ncol(df)
            dtf[i, colnames[j]] = df[i, colnames[j]]
        end

        Fill(dtf.tt)
    end

    Write(dtf.tt)
    Close(dtf.tf)
end

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

is_selected(i, row, selector::Function) = selector(row)
is_selected(i, row, selector::AbstractArray) = selector[i]

function with{T, R <: Real}(df::TreeDataFrame{T}, func::Function, selector, branches::Vector{Symbol}, rng=1:length(df), ::Type{R}=Float64)
    enable_branches(df, branches)
    ntot = 0
    tic()
    ret = zeros(R, length(rng))
    bitmask = BitArray(length(rng))
    enable_branches(df, branches)
    for i=rng
        ntot += load_row(df, i, branches)

        @inbounds const sel = is_selected(i, df.row, selector)::Bool
        @inbounds bitmask[i] = sel
        if sel
            @inbounds const res = func(df.row)::R
            @inbounds ret[i] = res
        end
    end
    nmb = float(ntot) / 1024.0 / 1024.0
    dt = toq()
    speed = ntot/dt / 1024.0 / 1024.0
    println("Read ", round(nmb), " Mb in ", dt, " s speed=", round(speed,2), " Mb/s")
    return ret[bitmask]
end

function loop{T}(df::TreeDataFrame{T}, f1::Function, f2::Function, branches::Vector{Symbol}, rng=1:length(df))
    enable_branches(df, branches)
    ntot = 0
    tic()
    enable_branches(df, branches)
    for i=rng
        ntot += load_row(df, i, branches)
        @inbounds const sel = f2(df.row)::Bool
        if sel
            @inbounds f1(df.row)
        end
    end
    nmb = float(ntot) / 1024.0 / 1024.0
    dt = toq()
    speed = ntot/dt / 1024.0 / 1024.0
    println("Read ", round(nmb), " Mb in ", dt, " s speed=", round(speed,2), " Mb/s")
end


export with, loop
export BranchValue
export writetree, TreeDataFrame
export writetree_temp
export load_row
export enable_branches
end # module
