"""
Script layout
Returns a sequence of gegulatory region based on accession
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
        "--gtf_file", "-t"
            help = "input gtf_file"
            default = "/home/bensiv/Documents/Remilk_plants/Genomes/Glycine_max/ncbi_dataset/data/GCF_000004515.6/genomic.gtf"
        "--sequence_report", "-s"
            help = "input sequence_report"
            default = "/home/bensiv/Documents/Remilk_plants/Genomes/Glycine_max/ncbi_dataset/data/GCF_000004515.6/sequence_report.jsonl"
        "--genome", "-g"
            help = "input genome"
            default = "/home/bensiv/Documents/Remilk_plants/Genomes/Glycine_max/GenomeBlast/Gmax.fna"
        "--upstream", "-u"
            help = "number of bp upsteam to transctiprion"
            default = 1000
            arg_type = Int
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

# load accessions
accessions_df = CSV.read(parsed_args["accessions"], DataFrame, header = false)
rename!(accessions_df, names(accessions_df) .=> ["accession"])

# load GTF data 
GTF = CSV.read(parsed_args["gtf_file"], DataFrame, header = false, comment = "#")

# name the columns
colnames = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
rename!(GTF, names(GTF) .=> colnames)

# split attribute column to subcolumns

# split attribute column to subcolumns
transcript_ids = []
descriptions = []
for row in eachrow(GTF)
    attributes = string.(split(row.attribute, ";"))
    attributes_kv = split.(attributes, "\"")
    description_length = length(descriptions)
    transcript_id_length = length(transcript_ids)

    for (i,x) in enumerate(attributes_kv)
        # attributes_kv[i] = replace.(x, " " => "")[findall(isempty.(x) .== false)]
        kv = replace.(x, " " => "")[findall(isempty.(x) .== false)]
        if length(kv) > 1
            if "transcript_id" in kv
                push!(transcript_ids, kv[2])
            elseif "product" in kv
                push!(descriptions, kv[2])
            end
        end
    end
    if length(descriptions) == description_length
        push!(descriptions, missing)
    end
    if length(transcript_ids) == transcript_id_length
        push!(transcript_ids, missing)
    end
end

GTF.transcript_id = transcript_ids
GTF.description = descriptions

dropmissing!(GTF)

# transcript_ids = []
# for attribute in GTF.attribute
#     idloc = findfirst("transcript_id", attribute)
#     push!(transcript_ids, attribute[idloc[end]+3:idloc[end]+16])
# end

# filter for accession
GTF = @linq GTF |> where(in(accessions_df.accession).(:transcript_id))
# GTF = GTF[in(accessions_df.accession).(transcript_ids),:]

# filter for transcript location only
GTF_transcripts = @linq GTF |> where(:feature .== "transcript")

# load genome
Genome = open(FASTA.Reader, parsed_args["genome"])

# get chromosomes length
ChrLengths = Dict()
for record in Genome
    push!(ChrLengths, "$(FASTA.identifier(record))" => length(FASTA.sequence(record)))
end

using JSON

# uplaod the squence report to convert from genbankAccession in the genome to refseqAccession in the transcriptome
sequence_report_file = open(parsed_args["sequence_report"], "r")
sequence_report = DataFrame()
for line in eachline(sequence_report_file)
    line_data = JSON.parse(line)
    line_df = DataFrame(collect(keys(line_data)) .=> collect(values(line_data)))
    select!(line_df, [:genbankAccession, :refseqAccession])
    push!(sequence_report, first(line_df))
end
close(sequence_report_file)

# filter for only chromosomes in genome
sequence_report = @linq sequence_report |> where(in(collect(keys(ChrLengths))).(:genbankAccession))

# filter for chromosomes with accession
GTF_transcripts = @linq GTF_transcripts |> where(in(sequence_report.refseqAccession).(:seqname))

# add genbankAccession column to GTF table
GTF_transcripts = innerjoin(GTF_transcripts, rename(sequence_report, "refseqAccession" => "seqname"), on = :seqname)

# add a column of the transcript name 
GTF_transcripts.name = [element[2] for element in split.([element[4] for element in split.(GTF_transcripts.attribute)], '"')]

# add a column of the range upstream 
function RegRange(GTFs, sequence_report, chr_lengths, bps)
    # reg_range = []
    reg_range = DataFrame(name = [], from = [], to = [])
    for row in eachrow(GTF_transcripts)
        if row.strand == "+"
            if row.start-bps < 1
                push!(reg_range, [row.name, 1, row.start])
            else
                push!(reg_range, [row.name, row.start-bps, row.start])
            end
        else
            max_length = ChrLengths[row.genbankAccession]
            if row.end+bps > max_length
                push!(reg_range, [row.name, row.end, max_length])
            else
                push!(reg_range, [row.name, row.end, row.end+bps])
            end
        end
    end    

    return reg_range
end

reg_range = RegRange(GTF_transcripts, sequence_report, ChrLengths, parsed_args["upstream"])

# add the results as a new colmuns
RegRegion = @linq GTF_transcripts |> leftjoin(reg_range, on = :name) |> select(:seqname, :name, :from, :to, :strand) |> rename(:seqname => :chr)
# GTF_transcripts.RegRange = reg_range

# store the locations in new dataframe
# RegRegion = DataFrame(chr = GTF_transcripts.seqname, name = GTF_transcripts.name, from = first.(GTF_transcripts.RegRange), to = last.(GTF_transcripts.RegRange), strand = )

CSV.write(join([parsed_args["output_dir"], "RegRegion.tsv"],"/"), RegRegion, delim = "\t")

# extract sequences and save in fasta format
Genome = open(FASTA.Reader, parsed_args["genome"])
RegSeq = open(FASTA.Writer, join([parsed_args["output_dir"], "RegSeq.fasta"],"/"))
for record in Genome
    Acc = @linq sequence_report |> where(:genbankAccession .== FASTA.identifier(record))
    ChrDF = @linq RegRegion |> where(:chr .== Acc.refseqAccession[1])
    for row in eachrow(ChrDF)
        if row.strand == "+"
            rec = FASTA.Record(row.name, FASTA.sequence(record)[row.from:row.to])
            write(RegSeq, rec)
        else
            rec = FASTA.Record(row.name, reverse_complement(FASTA.sequence(record)[row.from:row.to]))
            write(RegSeq, rec)
        end
    end
end
close(RegSeq)
