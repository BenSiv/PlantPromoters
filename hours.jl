
cd("/home/bensiv/Dropbox/202207-RemilkProject/Localized_output/")

# activate project enviroment
using Pkg
Pkg.activate("/home/bensiv/Dropbox/202207-RemilkProject/")

using CSV, DataFrames, DataFramesMeta
using Dates

WorkLog = DataFrame(Start = DateTime[], End = DateTime[])

push!(WorkLog, [Start, End])

WorkLog.Total = Minute.(WorkLog.End .- WorkLog.Start)

using Statistics

Hour(sum(WorkLog.Total))

CSV.write("WorkHouts.csv", WorkLog)

# ===========================================================

WorkLog = CSV.read("WorkHouts.csv", DataFrame)

Start = DateTime(2022,08,08,15,15)
End = DateTime(2022,08,08,16,30)
push!(WorkLog, [Start, End, Minute(End-Start)])
