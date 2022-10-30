"""
This is a script that takes as input tsv as outputted from blastn (outfmt 6) and returns the same file after preliminary processing.
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

using CSV, DataFrames

Blast = CSV.read(parsed_args["input"], DataFrame, delim = "\t", header = false)

ColumnNames = ["Query", "Reference", "Precent_identity", "Alignment_length", "Mismach", "Gap", "Que_Start", "Que_End", "Ref_Start", "Ref_End", "Pvalue", "Escore"]

rename!(Blast, names(Blast) .=> ColumnNames)

Strand = [row.Ref_Start < row.Ref_End ? "+" : "-" for row in eachrow(Blast)]

Blast.Strand = Strand

CSV.write(parsed_args["output"], Blast, delim = "\t")