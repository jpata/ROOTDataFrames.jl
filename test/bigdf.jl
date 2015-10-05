using ROOTDataFrames, DataFrames
using Base.Test

N = 500000
m = 25

df = DataFrame()

for i=1:m
    df[symbol("col$i")] = Float64[j for j=1:N]
end

for i=m+1:m+10
    df[symbol("dyncol$i")] = Vector{Float64}[rand(10) for j=1:N]
end
# 
# for i=m+1:m+10
#     df[symbol("dyncol$i")] = Vector{Float64}[rand(10) for j=1:N]
# end

println("writing DataFrame")
writetree("test.root", df)
