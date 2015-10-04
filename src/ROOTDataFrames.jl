module ROOTDataFrames

include("log.jl")

#maximum size of preallocation for vector branches
const MAX_SIZE = 50

#print debugging info
const DEBUG = false

using ROOT, DataFrames, DataArrays
import Base.length, Base.getindex

import DataFrames.nrow, DataFrames.size, DataFrames.index

#holds data from TBranch
type BranchValue{T, N}
    data::Vector{T} #buffer
    size::Symbol #branch size (number or other branch name)
    branch::TBranch
end

BranchValue{T<:Real,N <: Symbol}(::Type{T}, s::N, br::TBranch=TBranch(C_NULL)) =
    BranchValue{T,N}(zeros(T, MAX_SIZE), s, br)

BranchValue{T<:AbstractVector, N <: Symbol}(::Type{T}, s::N, br::TBranch=TBranch(C_NULL)) =
    BranchValue{eltype(T),N}(zeros(eltype(T), MAX_SIZE), s, br)

#static array
BranchValue{T<:Real, N}(::Type{T}, s::Type{Val{N}}, br::TBranch=TBranch(C_NULL)) =
    BranchValue{T,N}(zeros(T, N), :nothing, br::TBranch)

#static array from vector
BranchValue{T<:AbstractVector, N}(::Type{T}, s::Type{Val{N}}, br::TBranch=TBranch(C_NULL)) =
    BranchValue{eltype(T),N}(zeros(eltype(T), N), :nothing, br::TBranch)


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

#return value held in branch buffer
value{T, N <: Symbol}(bval::BranchValue{T, N}) = bval.data::Vector{T}
value{T, N}(bval::BranchValue{T, N}) = bval.data::Vector{T}
value{T}(bval::BranchValue{T, 1}) = bval.data[1]::T

import Base.call
call{T}(bval::BranchValue{T, 1}) = bval.data[1]::T
call{T, N}(bval::BranchValue{T, N}) = bval.data[1:N]::Vector{T}
call{T, N <: Symbol}(bval::BranchValue{T, N}) = bval.data::Vector{T}

#makes a class that represents the structure of a tree
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


#return value of scalar branch
@generated function Base.getindex{T <: Any, X <: Any, N}(
    df::TreeDataFrame{T},
    i::Int64,
    s::Symbol,
    getentry::Bool=true,
    ft::Type{BranchValue{X, N}}=fieldtype(T, s),
    )

    if N == Symbol
        ex = quote
            #get dynamic branch size
            sizebranch = getsize(getfield(df.row, s))
            if getentry
                load_row(df, i, [s, sizebranch]) > 0 || error("could not load row $i")
            end
            _size = value(getfield(df.row, sizebranch))
            ret = getfield(df.row, s)
            return value(ret)[1:_size]::Vector{X}
        end
    else
        if N == 1
            ret  = :(value(ret)::X)
        elseif N > 1
            ret = :(value(ret)[1:$N]::Vector{X})
        end

        ex = quote
            if getentry
                load_row(df, i, [s])>0 || error("could not load row $i")
            end
            ret = getfield(df.row, s)
            return $ret
        end
    end

    return ex

end
# 
# #return value of scalar branch
# function Base.getindex{T <: Any, X <: Any}(
#     df::TreeDataFrame{T},
#     i::Int64,
#     s::Symbol,
#     getentry::Bool=true,
#     ::Type{BranchValue{X, N}}=fieldtype(T, s),
#     )
#     println("scalar called ", s)
#     if getentry
#         load_row(df, i, [s])>0 || error("could not load row $i")
#     end
#     ret = getfield(df.row, s)
#     @show ret
#     return value(ret)[1:_size]
# end

#return value of dynamic branch
# function Base.getindex{T <: Any, R <: Symbol, X <: Any}(
#     df::TreeDataFrame{T},
#     i::Int64,
#     s::Symbol,
#     getentry::Bool=true,
#     _size::Int64=-1,
#     ft::Type{BranchValue{X, R}}=fieldtype(T, s),
#     )
#     println(ft)
#     #get dynamic branch size
#     if _size == -1
#         sizebranch = getsize(getfield(df.row, s))
#         if getentry
#             load_row(df, i) > 0 || error("could not load row $i")
#         end
#         _size = value(getfield(df.row, sizebranch))
#     end
#     ret = getfield(df.row, s)
#     @show ret, _size
#     return value(ret)[1:_size]::Vector{X}
# end

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
    coltypes = Type[]
    for col in cols
        _size = sizes[col]
        _type = names_types[col]
        if typeof(_size) <: Symbol
            _size = MAX_SIZE
        end
        if _size > 1
            _type = DataVector{_type}
        end
        push!(coltypes, _type)
        a = DataFrames.DataArray(_type, n)
        push!(arrlist, a)
    end

    #construct output dataframe
    const ret = DataFrame(
        arrlist,
        DataFrames.Index(cols)
    )
    j = 1

    t0 = time()
    nloaded = 0
    for i=1:nrow(df)
        (!isna(mask[i]) && mask[i]) || continue
        nloaded += load_row(df, i, cols)
        for (icol, col) in enumerate(cols)
            ret[j, col] = df[i, col]
        end
        j += 1
    end
    t1 = time()
    enable_branches(df, ["*"])
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
    types::Vector{Type},
    colsizes::Vector{Any}
    ;kwargs...
    ) = TreeDataFrame(
        fn, names,
        BranchValue[BranchValue(dt, cs) for (dt, cs) in
        zip(types, colsizes)]; kwargs...
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

        sz = getsize(typ)
        rtype = SHORT_TYPEMAP[eltyp]


        
        brstring = ""

        #variable size branch
        if isa(sz, Symbol)
            brstring = "$name[$sz]/$rtype"
        #fixed size branch
        elseif sz > 1
            brstring = "$name[$sz]/$rtype"
        # scalar branch
        else
            brstring = "$name/$rtype"
        end

        br = Branch(
            tree, string(name),
            convert(Ptr{Void}, pointer(getfield(rowdata, name).data)),
            brstring
        )
        assert(root_pointer(br) != C_NULL)
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
    v = val[:]
    lv = length(v)
    assert(lv < MAX_SIZE)
    #fill with zeros
    getfield(df.row, col).data[:] = zero(eltype(getfield(df.row, col)))

    #fill with correct length
    getfield(df.row, col).data[1:lv] = v[1:lv] 
end

function TreeDataFrame(df::AbstractDataFrame, filename; size_branches::Associative=Dict(), treename="dataframe")
    colsizes = Any[]

    _eltypes = eltypes(df)
    _names = names(df)
    _sizes = Array(Any, ncol(df))
    for i=1:ncol(df)
        _el = _eltypes[i]
        _cname = _names[i]
        
        size_branch = get(size_branches, _cname, symbol(:n, _cname))

        idx_size = findfirst(_names, size_branch)
        if idx_size > 0
            idx_size > i && error("found size branch $size_branch for $_cname, but must be declared before.")
            _size = size_branch

        #constant size
        elseif _el <: AbstractVector && eltype(_el) <: Real
            maxsize = maximum(map(x->length(df[x, _cname]), 1:nrow(df)))
            assert(maxsize < MAX_SIZE)
            _size = Val{maxsize}
        #scalar
        else
            _size = Val{1}
        end
        _sizes[i] = _size
    end
    dtf = TreeDataFrame(filename, _names, _eltypes, _sizes;treename=treename)

    for i=1:nrow(df)
        for j=1:ncol(df)
            dtf[i, _names[j]] = df[i, _names[j]]
        end

        Fill(dtf.tt)
    end

    Write(dtf.tt)
    Close(dtf.tf)
    return TreeDataFrame([filename]; treename=treename)
end

function writetree(
    fn::ASCIIString,
    df::AbstractDataFrame;
    size_branches=Dict(),
    treename="dataframe"
    )
    TreeDataFrame(df, fn; size_branches=size_branches, treename=treename)
end

is_selected(i, row, selector::Function) = selector(row)
is_selected(i, row, selector::AbstractArray) = selector[i]

function with{T, R <: Real}(
    df::TreeDataFrame{T}, #input data
    func::Function, #process func ev->R
    selector, #selector func ev->Bool
    branches::Vector{Symbol}, #branches to enable
    rng=1:length(df), #rows to process
    ::Type{R}=Float64 #return type
    )
    enable_branches(df, branches)
    ntot = 0
    tic()
    ret = zeros(R, length(rng))
    bitmask = BitArray(length(rng))
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
    #println("Read ", round(nmb), " Mb in ", dt, " s speed=", round(speed,2), " Mb/s")
    #re-enable branches
    enable_branches(df, ["*"])
    return ret[bitmask]
end

function loop{T}(df::TreeDataFrame{T}, f1::Function, f2::Function, branches::Vector{Symbol}, rng=1:length(df))
    enable_branches(df, branches)
    ntot = 0
    tic()
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
    enable_branches(df, ["*"])
end

export value
export with, loop
export BranchValue
export writetree, TreeDataFrame
export writetree_temp
export load_row
export enable_branches
end # module
