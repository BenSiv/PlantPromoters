"""
Script layout
Combines the tables of all the plants
"""

# activate project enviroment
using Pkg
Pkg.activate("/home/bensiv/Dropbox/202207-RemilkProject/")

using CSV, DataFrames, DataFramesMeta

directories = [dir for dir in readdir() if !occursin(".", dir)]

AllTables = DataFrame()
for dir in directories
    append!(AllTables, CSV.read("$dir/OutputTable.csv", DataFrame))
    println(dir)
end

CSV.write("OutputTable.csv",AllTables)