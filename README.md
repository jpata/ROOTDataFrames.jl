# ROOTDataFrames

Wraps ROOT TTrees through the AbstractDataFrame interface. Supports on-demand access of columns and rows which are based either on disk or on the network.

#Installation

1. Install ROOT.jl: http://github.com/jpata/ROOT.jl
2. Install this package using `Pkg.clone("https://github.com/jpata/ROOTDataFrames.jl.git")`
3. Test using `Pkg.test("ROOTDataFrames")`

#Usage

~~~
using DataFrames, ROOT, ROOTDataFrames
data = TreeDataFrame(ASCIIString["file1.root", "file2.root"], "my_event_tree")
N = nrow(data)
c1 = data[[:col1, :col2]]

for i=1:N
  load_row(data, i)
  x = data[i, :x]
  y = data[i, :y]
end
~~~
