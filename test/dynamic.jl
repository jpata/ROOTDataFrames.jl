using DataFrames, ROOT, ROOTDataFrames
using Base.Test



dyndf = DataFrame(
    njet=Int64[3, 1],
    jet_pt=Vector{Float64}[Float64[1.0, 2.0, 3.0], Float64[4.0]],
    jet_eta=Vector{Float32}[Float32[3.0, 4.0, 5.0], Float32[6.0]],
    jet_id=Vector{Int64}[Int64[1,2,3], Float32[-4]],
    arr=Vector{Float64}[rand(5), rand(5)],
)

dtf = TreeDataFrame(dyndf, "dyn.root"; size_branches=Dict(:jet_pt=>:njet, :jet_eta=>:njet, :jet_id=>:njet))
@test dtf[:jet_pt] == dyndf[:jet_pt]
@test dtf[:jet_pt] == dyndf[:jet_pt]
@test dtf[:jet_id] == dyndf[:jet_id]
@test dtf[:arr] == dyndf[:arr]
print_with_color(:green, "small dynamic tree OK\n")

N = 10000
bigdf = similar(
    DataFrame(
        njet=Int64[],
        jet_pt=Vector{Float64}[],
        jet_eta=Vector{Float64}[]
    ),
    N
)

for iev=1:N
    nj = round(Int, rand()*10)
    bigdf[iev, :njet] = nj
    bigdf[iev, :jet_pt] = [rand()*200 for i=1:nj]
    bigdf[iev, :jet_eta] = [randn()*5 for i=1:nj]
end

dtf2 = TreeDataFrame(bigdf, "dyn2.root"; size_branches=Dict(:jet_pt=>:njet, :jet_eta=>:njet))
@test mean(map(sum, bigdf[:jet_pt])) == mean(map(sum, dtf2[:jet_pt]))
print_with_color(:green, "big dynamic tree OK\n")

df = TreeDataFrame([string(Pkg.dir(), "/ROOTDataFrames/dynamic.root")];
    treename="tthNtupleAnalyzer/events",
)

sumpt = 0.0
sumjet = 0
for i=1:nrow(df)
    load_row(df, i)
    n__jet = df.row.n__jet()
    jet__pt = df.row.jet__pt()[1:n__jet]
    sumjet += n__jet
    sumpt += sum(jet__pt)
end

@test_approx_eq_eps sumpt/sumjet 56.964209222875205 0.01
print_with_color(:green, "realistic dynamic tree OK\n")


r = with(df,
    row->Float64(sum(row.jet__pt()[1:row.n__jet()])),
    row->row.n__jet()>3, [:jet__pt, :n__jet], 1:length(df),
    Float64
)
@test_approx_eq_eps mean(r) 688.4442511628289 0.01
print_with_color(:green, "with() OK\n")

outdf = DataFrame()
v = Array(Vector{Float64}, 3)
v[1] = [1.0,2.0,3.0]
v[2] = [4.0,5.0]
v[3] = [6.0]

v2 = Array(Vector{Float64}, 3)
v2[1] = round(rand(2),2)
v2[2] = round(rand(4),2)
v2[3] = round(rand(6),2)

outdf[:b] = v2
outdf[:na] = [3,2,1]
outdf[:a] = v
writetree("dyntree.root", outdf)

df2 = TreeDataFrame(["dyntree.root"])
@test all(df2[:na] .== outdf[:na])

for i=1:length(df2)
    load_row(df2, i)
    @test all(df2[i, :a] .== outdf[i, :a])
end
print_with_color(:green, "dynamic tree with auto branch size OK\n")
