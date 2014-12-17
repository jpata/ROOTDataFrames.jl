using DataFrames, ROOT, ROOTDataFrames
using Base.Test

df = TreeDataFrame([string(Pkg.dir(), "/ROOTDataFrames/dynamic.root")], "tthNtupleAnalyzer/events")

sumpt = 0.0
sumjet = 0
for i=1:nrow(df)
    load_row(df, i)
    n__jet = df[i, :n__jet, Int64]
    jet__pt = df[i, :jet__pt, Float64]
    sumjet += n__jet
    sumpt += sum(jet__pt)
end

@test_approx_eq_eps sumpt/sumjet 56.964209222875205 0.01
