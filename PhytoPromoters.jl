"""
Script layout
Input: promoters.fasta, regulatory_elements.tsv
Output: promoters.gb, minimal_promoters.fasta, FinalOutput.tsv
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
# regulatory_elements = @linq regulatory_elements |> unique(:Sequence) |> transform(:Sequence = LongDNA{4}.(:Sequence))
regulatory_elements = @linq regulatory_elements |> transform(:Sequence = LongDNA{4}.(:Sequence))

# find all the occurences of regulatory elements in the promoters
function find_features(promoters_dict, reg_df)
    FeaturesTable = DataFrame(Promoter = [], Motif = [], Occurences = [], Direction = [])

    for key in keys(promoters_dict)
        for motif in eachrow(reg_df)
            
            query = ExactSearchQuery(motif.Sequence, iscompatible)
            # rev_query = ExactSearchQuery(reverse_complement(motif.Sequence), iscompatible)
    
            sense = findall(query, promoters_dict[key])
            # anti_sense = findall(rev_query, promoters_dict[key])
            
            # occurences = [sense;anti_sense]
            occurences = sense
            # directions = [fill("+", length(sense)); fill("-", length(anti_sense))]
            directions = fill("+", length(sense))
            n = length(occurences)
            if !isempty(occurences)
                FeaturesTable = [FeaturesTable; DataFrame(Promoter = fill(key, n), Motif = fill(motif.Name, n), Occurences = occurences, Direction = directions)]
                # push!(Features, ["$key", "$(motif.Name)", occurences, Directions])
            end
        end
    end
    
    return FeaturesTable
end

FeaturesTable = find_features(promotersDict, regulatory_elements)


# # save the output as a jld file
# using JLD
# save("Features.jld", "Features", Features)

CSV.write(join([parsed_args["output_dir"], "Features.csv"],"/"), FeaturesTable)

# include functions from the GenBank file
include("GenBank.jl")

# generate a genbank records
GBRecords = write_gb_records(promotersDict, FeaturesTable)
# GBRecords = [GBRecords[1:285];GBRecords[287:299]]

# write the records to a GenBank file
gb = open(join([parsed_args["output_dir"], "promoters.gb"],"/"), "w")
for record in GBRecords
    write(gb, writeGB(record))
end
close(gb)

# stop("stop here")

# include functions from the minimal_promoter file
include("minimal_promoter.jl")

# load("Features.jld")["Features"]

FeatureLoc = Dict()
for prom in unique(FeaturesTable.Promoter)
    tempdf = @linq FeaturesTable |> where(:Promoter .== prom)
    FeatureRange = sort(unique(vcat(tempdf.Occurences...)))
    push!(FeatureLoc, "$prom" => FeatureRange)
end 

FeatureLocDilu = Dict()
for key in keys(FeatureLoc)
    push!(FeatureLocDilu, "$key" => Dilution(FeatureLoc[key], 8, length(promotersDict[key])))
end 

MinimalProm = Dict()
for key in keys(promotersDict)
    MinProm = GetLocSeq(promotersDict[key], FeatureLocDilu[key])
    push!(MinimalProm, "$key" => MinProm)
end

minimal_promoters = open(FASTA.Writer, join([parsed_args["output_dir"], "minimal_promoters.fasta"],"/"))
for seq in keys(MinimalProm)
    rec = FASTA.Record(seq, MinimalProm[seq])
    write(minimal_promoters, rec)
end
close(minimal_promoters)

# ==================================================================

# load promoters
minimal_promoters = open(FASTA.Reader, join([parsed_args["output_dir"], "minimal_promoters.fasta"],"/"))

# access the promotor and store them as a dictionary
minimal_promotersDict = Dict()
for record in minimal_promoters
    push!(minimal_promotersDict, "$(FASTA.identifier(record))" => FASTA.sequence(record))
end

# find all the occurences of regulatory elements in the promoters
FeaturesTable = find_features(minimal_promotersDict, regulatory_elements)

# include functions from the GenBank file
include("GenBank.jl")

# generate a genbank records
GBRecords = write_gb_records(minimal_promotersDict, FeaturesTable)


# write the records to a GenBank file
gb = open(join([parsed_args["output_dir"], "minimal_promoters.gb"],"/"), "w")
for record in GBRecords
    write(gb, writeGB(record))
end
close(gb)
