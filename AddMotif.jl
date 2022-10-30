"""
Script layout
Takes a fasta file and inserts a sequence in desirable locations
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
        "--promoter", "-p"
            help = "input promoter file as fasta"
            required = true
        "--regulatory_elements", "-r"
            help = "input regulatory_elements file as tsv"
            default = "PLACE_plus_audited.tsv"
            required = true
        "--motif", "-m"
            help = "name of the motif to insert"
            required = true
        "--locations", "-l"
            help = "locations to insert motif seperated by ','"
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


motifs = CSV.read(parsed_args["regulatory_elements"], DataFrame)

motif = @linq motifs |> where(:Name .== parsed_args["motif"]) |> select(:Sequence) |> first |> collect |> first
motif = LongDNA{4}(motif)

locations = parse.(Int64, string.(split(parsed_args["locations"],",")))

promoters = open(FASTA.Reader, parsed_args["promoter"])
promoters_plus_motif = open(FASTA.Writer, join([parsed_args["output_dir"],"promoters.fasta"],"/"))
for record in promoters
    sequence = FASTA.sequence(record)
    for location in locations
        sequence = join([sequence[1:location],motif,sequence[location+1:end]],"")
    end
    new_promoter = FASTA.Record(FASTA.identifier(record), sequence)
    write(promoters_plus_motif, new_promoter)
end
close(promoters_plus_motif)
