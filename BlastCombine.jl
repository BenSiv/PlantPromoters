"""
Combine adjacent sequences for downstream analysis
"""

using Pkg
Pkg.activate("/home/bensiv/Dropbox/202207-RemilkProject/")


using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
            help = "blastn results (outfmt 6) as a tsv file"
            required = true
        "--output", "-o"
            help = "outpput name as a tsv file"
            default = "output.tsv"
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end

    return parsed_args
end

parsed_args = main()


using CSV, DataFrames, DataFramesMeta

Blast = CSV.read(parsed_args["input"], DataFrame, delim = "\t")

Blast.SubGroup = fill(1, nrow(Blast))

BlastGrouped = @linq Blast |> groupby([:Query, :Reference, :Strand])

BlastCombined = DataFrame()
for group in BlastGrouped
    # group = BlastGrouped[1]
    if nrow(group) > 1
        sort!(group, :Ref_Start)

        SubGroup = [1]
        for i in 2:nrow(group)
            if (group.Ref_Start[i] - group.Ref_Start[i-1] < group.Alignment_length[i-1]*10) && (group.Ref_End[i] - group.Ref_End[i-1] < group.Alignment_length[i-1]*10)
                push!(SubGroup, last(SubGroup))
            else
                push!(SubGroup, last(SubGroup)+1)
            end
        end


        group.SubGroup = SubGroup
        group_subgrouped = @linq group |> groupby([:SubGroup])

        for subgroup in group_subgrouped
            # subgroup = group_subgrouped[1]
            if nrow(subgroup) > 1
                combined = first(subgroup)
                if subgroup.Strand == "+"
                    combined.Ref_Start, combined.Que_Start = minimum(subgroup.Ref_Start), minimum(subgroup.Que_Start)
                    combined.Ref_End, combined.Que_End = maximum(subgroup.Ref_End), maximum(subgroup.Que_End)
                    combined.Alignment_length = combined.Ref_End - combined.Ref_Start
                else
                    combined.Ref_Start, combined.Que_Start = maximum(subgroup.Ref_Start), maximum(subgroup.Que_Start)
                    combined.Ref_End, combined.Que_End = minimum(subgroup.Ref_End), minimum(subgroup.Que_End)
                    combined.Alignment_length = combined.Ref_Start - combined.Ref_End
                end
                push!(BlastCombined, combined)
            else
                push!(BlastCombined, first(subgroup))
            end
        end
    else
        push!(BlastCombined, first(group))
    end
end


CSV.write(parsed_args["output"], BlastCombined, delim = "\t")