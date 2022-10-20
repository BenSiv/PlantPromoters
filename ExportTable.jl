"""
Script layout
Exports the final output table
"""

# activate project enviroment
using Pkg
Pkg.activate("/home/bensiv/Dropbox/202207-RemilkProject/")

# get command line arguments from the user
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--output_dir", "-o"
            help = "outpput directory"
            default = pwd()
        "--accessions", "-a"
            help = "input accessions list"
            required = true
        "--promoter", "-p"
            help = "input promoter file as fasta"
            required = true
        "--minimal_promoter", "-m"
            help = "input minimal promoter file as fasta"
            required = true
        "--features", "-f"
            help = "input features per promoter"
            required = true
        "--blast", "-b"
            help = "output of blast on promoter after proccesing"
            required = true
        "--interactive", "-i"
            help = "save the commandline arguments as a yaml file for interactive running"
            action = :store_true
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

using YAML
include("stop.jl")

if parsed_args["interactive"]
    YAML.write_file("parsed_args.yml", parsed_args)
    stop("stopped and saved the arguments as parsed_args.yml")
    # run `parsed_args = YAML.load_file("parsed_args.yml")` to retrive it in the script
end

# precompile packages
using CSV, DataFrames, DataFramesMeta
using BioSequences, FASTX


Accessions = CSV.read(parsed_args["accessions"], DataFrame)

promoter_df = DataFrame(accession = String[], sequence = LongDNA[])
open(FASTA.Reader, parsed_args["promoter"]) do reader
    for record in reader
        push!(promoter_df, [FASTA.identifier(record), FASTA.sequence(record)])
    end
end

minimal_promoter_df = DataFrame(accession = String[], minimal_sequence = LongDNA[])
open(FASTA.Reader, parsed_args["minimal_promoter"]) do reader
    for record in reader
        push!(minimal_promoter_df, [FASTA.identifier(record), FASTA.sequence(record)])
    end
end

OutputTable = innerjoin(promoter_df, minimal_promoter_df, on = "accession")
OutputTable = outerjoin(OutputTable, Accessions, on = "accession")
OutputTable.length = length.(OutputTable.sequence)
OutputTable.minimal_length = length.(OutputTable.minimal_sequence)

features = CSV.read(parsed_args["features"], DataFrame)
feature_groups = @linq features |> groupby(:Promoter)
features_df = DataFrame(accession = vcat([unique(group.Promoter) for group in feature_groups]...), features = [unique(group.Motif) for group in feature_groups])

OutputTable = innerjoin(OutputTable, features_df, on = "accession")

RYrepeats = ["RY-repeat", "RYREPEAT4", "RYREPEATLEGUMINBOX", "RYREPEATVFLEB4", "RYREPEATGMGY2", "RYREPEATBNNAPA"]
OutputTable.contains_RYrepeat = [any(in(f).(RYrepeats)) for f in OutputTable.features]

PromBlast = CSV.read(parsed_args["blast"], DataFrame, delim = "\t")

using Pipe
BlastCount = @pipe PromBlast |> groupby(_, :Query) |> combine(_, nrow => :genome_count) |> rename(_, :Query => :accession)

OutputTable = outerjoin(OutputTable, BlastCount, on = "accession")

select!(OutputTable, :accession => :name, :description => :protein_annotation, :sequence, :length, :minimal_sequence, :minimal_length, :features, :contains_RYrepeat, :genome_count)

CSV.write(join([parsed_args["output_dir"], "OutputTable.csv"], "/"), OutputTable)

