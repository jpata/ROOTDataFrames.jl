using ROOTDataFrames, DataFrames
using Base.Test

df = DataFrame(a=Float64[1,2,3])

writetree("test.root", df)

tdf = TreeDataFrame("test.root")


@test all(tdf[:a] .== df[:a])
