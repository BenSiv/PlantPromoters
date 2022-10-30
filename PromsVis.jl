"""
Script layout
Visualizing the promoter motifs
"""

cd("/home/bensiv/Dropbox/202207-RemilkProject/Localized_output/")

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

# load promoters
promoters = open(FASTA.Reader, parsed_args["promoter"])

# access the promotor and store them as a dictionary
promotersDict = Dict()
for record in promoters
    push!(promotersDict, "$(FASTA.identifier(record))" => FASTA.sequence(record))
end

# load regulatory elements as a dataframe
regulatory_elements = CSV.read(parsed_args["regulatory_elements"], DataFrame)
regulatory_elements = @linq regulatory_elements |> transform(:Sequence = LongDNA{4}.(:Sequence))


# find all the occurences of regulatory elements in the promoters
function find_features(promoters_dict, reg_df)
    FeaturesTable = DataFrame(Promoter = [], Motif = [], Occurences = [], Direction = [])

    for key in keys(promoters_dict)
        for motif in eachrow(reg_df)
            
            query = ExactSearchQuery(motif.Sequence, iscompatible)
            rev_query = ExactSearchQuery(reverse_complement(motif.Sequence), iscompatible)
    
            sense = findall(query, promoters_dict[key])
            anti_sense = findall(rev_query, promoters_dict[key])
            
            occurences = [sense;anti_sense]
            # occurences = sense
            directions = [fill("+", length(sense)); fill("-", length(anti_sense))]
            # directions = fill("+", length(sense))
            n = length(occurences)
            if !isempty(occurences)
                FeaturesTable = [FeaturesTable; DataFrame(Promoter = fill(key, n), Motif = fill(motif.Name, n), Occurences = occurences, Direction = directions)]
            end
        end
    end
    
    return FeaturesTable
end

FeaturesTable = find_features(promotersDict, regulatory_elements)
CSV.write(join([parsed_args["output_dir"], "FeaturesTable.tsv"],"/"), FeaturesTable, delim = "\t")

function motif_locations(features_table, motif, direction)
    temp = @linq features_table |> where(:Motif .== motif .&& :Direction .== direction)
        
    return first.(temp.Occurences)
end

RYrepeat = motif_locations(FeaturesTable, "RYREPEATLEGUMINBOX", "+")
TATAbox = motif_locations(FeaturesTable, "TATABOX1", "+")
TATAbox = motif_locations(FeaturesTable, "TATAbox-1", "+")

using Plots, StatsPlots

density(TATAbox)
density(RYrepeat)
density(TATAbox, label = "TATAbox-1", linewidth = 3, grid = false, ylabel = "motif density", xlabel = "location (bp)", title = "motif density across multiple promoters")

density([TATAbox, RYrepeat], label = ["TATA-box" "RYrepeat"], linewidth = 3, grid = false, ylabel = "motif density", xlabel = "location (bp)", title = "motif density across multiple promoters")

savefig(join([parsed_args["output_dir"], "motif_density.png"],"/"))