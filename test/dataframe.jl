using ROOTDataFrames, DataFrames
using Base.Test

df = DataFrame(a=Float64[1.0, 2.0, 3.0], b=Int64[0, 1, 2], c=['a', 'c', 'e'], e=[1, 1, 1])
df[3, :e] = NA

writetree("test.root", df)
@test isfile("test.root")

tdf = TreeDataFrame("test.root")
@test size(tdf) == size(df)

@test all(tdf[:a] .== df[:a])
@test all(tdf[:b] .== df[:b])
@test all(tdf[:c] .== df[:c])
@test isna(tdf[3, :e])
@test all(tdf[:e][1:2] .== df[:e][1:2])

rm("test.root")