using DataFrames, ROOT, ROOTDataFrames
using Base.Test

df = TreeDataFrame([string(Pkg.dir(), "/ROOTDataFrames/dynamic.root")], "tthNtupleAnalyzer/events")

for i=1:nrow(df)
    nloaded = load_row(df, i, Symbol[:n__lep, :lep__pt])
    nl = df[i, :n__lep]
    pts = df[i, :lep__pt, Float32][1:nl]
    nloaded2 = load_row(df, i)
    @test nl == df[i, :n__lep]
    @test pts == df[i, :lep__pt, Float32][1:nl]
end
