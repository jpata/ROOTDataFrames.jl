using DataFrames, ROOTDataFrames

N = 10000

#pre-allocate input data
bigdf = similar(

    #here we define TTree structure
    DataFrame(
        evt=Int64[],
        njet=Int64[],
        jet_pt=Vector{Float64}[],
        jet_eta=Vector{Float64}[],
        jet_id=Vector{Int64}[],

        nlep=Int64[],
        lep_pt=Vector{Float64}[],
        lep_eta=Vector{Float64}[],
        lep_id=Vector{Int64}[]

    ),
    N
)

#fill event by event
for iev=1:N
    nj = round(Int, rand()*10)
    bigdf[iev, :evt] = iev
    bigdf[iev, :njet] = nj

    #create jets with random vectors
    bigdf[iev, :jet_pt] = [rand()*200 for i=1:nj]
    bigdf[iev, :jet_eta] = [randn()*2 for i=1:nj]
    bigdf[iev, :jet_id] = [i for i=1:nj]

    nl = round(Int, rand()*2)
    bigdf[iev, :nlep] = nl
    bigdf[iev, :lep_pt] = [rand()*200 for i=1:nl]
    bigdf[iev, :lep_eta] = [randn()*2 for i=1:nl]
    bigdf[iev, :lep_id] = [i for i=1:nl]

end

#write output tree
dtf = TreeDataFrame(bigdf, "dyntree.root";
    #need to specify which branches correspond to dynamic branch lengths
    size_branches=Dict(
        :jet_pt=>:njet, :jet_eta=>:njet, :jet_id=>:njet,
        :lep_pt=>:nlep, :lep_eta=>:nlep, :lep_id=>:nlep,
    )
)
