using ROOTDataFrames, DataFrames
using Base.Test

df = DataFrame(a=Float64[1.0, 2.0, 3.0], b=Int64[0, 1, 2], e=Int32[1, 1, 1])

writetree("test.root", df)
@test isfile("test.root")

tdf = TreeDataFrame("test.root")
@test size(tdf) == size(df)
#@test all(names(tdf) .== [:a, :b, :e])

@test all(tdf[:a] .== df[:a])

@test all(tdf[[:a, :b]][:a] .== df[[:a, :b]][:a])
@test all(tdf[[:a, :b]][:b] .== df[[:a, :b]][:b])

@test all(tdf[:b] .== df[:b])
@test all(tdf[:e][1:2] .== df[:e][1:2])

@test all(tdf[[true, true, false], [:a]][:a] .== [1.0, 2.0])

tdf = TreeDataFrame(["test.root", "test.root"])
@test nrow(tdf) == 6
@test all(tdf[:a] .== vcat(df[:a], df[:a]))
#rm("test.root")

println("bigdf")
df = DataFrame(x = Int64[i for i=1:10^7])

writetree("big.root", df)
tdf = TreeDataFrame("big.root")
@test nrow(tdf) == nrow(df)

tdf = TreeDataFrame(["big.root", "big.root"])
@test nrow(tdf) == 2 * nrow(df)
#rm("big.root")
