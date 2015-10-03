using DataFrames, ROOT, ROOTDataFrames
using Base.Test


df = TreeDataFrame([string(Pkg.dir(), "/ROOTDataFrames/dynamic.root")]; treename="tthNtupleAnalyzer/events")

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


r = with(df,
    row->Float64(sum(row.jet__pt()[1:row.n__jet()])),
    row->row.n__jet()>3, [:jet__pt, :n__jet], 1:length(df),
    Float64
)
println(mean(r))
@test_approx_eq_eps mean(r) 688.4442511628289 0.01

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
