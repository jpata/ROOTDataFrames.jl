using ROOTDataFrames, DataFrames
using Base.Test

df = DataFrame(a=Float64[1.0, 2.0, 3.0], b=Int64[0, 1, 2], e=Int32[1, 1, 1])

writetree("test.root", df)
@test isfile("test.root")
print_with_color(:green, "writetree OK\n")

tdf = TreeDataFrame("test.root")
print_with_color(:green, "TreeDataFrame from file OK\n")
@test size(tdf) == size(df)
print_with_color(:green, "size OK\n")

#@test all(names(tdf) .== [:a, :b, :e])

@test all(tdf[:a] .== df[:a])
print_with_color(:green, "tdf[:a] .== df[:a] OK\n")

@test all(tdf[[:a, :b]][:a] .== df[[:a, :b]][:a])
print_with_color(:green, "tdf[[:a, :b]][:a] .== df[[:a, :b]][:a] OK\n")
@test all(tdf[[:a, :b]][:b] .== df[[:a, :b]][:b])

@test all(tdf[:b] .== df[:b])
@test all(tdf[:e][1:2] .== df[:e][1:2])

@test all(tdf[[true, true, false], [:a]][:a] .== [1.0, 2.0])

print_with_color(:green, "masked index OK\n")


tdf = TreeDataFrame(["test.root", "test.root"])
print_with_color(:green, "two-file chain opened OK\n")
@test nrow(tdf) == 6
print_with_color(:green, "two-file chain size OK\n")
@test all(tdf[:a] .== vcat(df[:a], df[:a]))
print_with_color(:green, "two-file :a OK\n")

#rm("test.root")

println("bigdf")
df = DataFrame(x = Int64[i for i=1:10^7])
writetree("big.root", df)
print_with_color(:green, "bigdf written OK\n")

tdf = TreeDataFrame("big.root")
@test nrow(tdf) == nrow(df)
print_with_color(:green, "bigdf opened OK\n")

tdf = TreeDataFrame(["big.root", "big.root"])
@test nrow(tdf) == 2 * nrow(df)
print_with_color(:green, "2x bigdf opened OK\n")

#rm("big.root")
