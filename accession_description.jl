"""
Extracts the description of accessions
"""

# activate parent enviroment
using Pkg
Pkg.activate("/home/bensiv/Dropbox/202207-RemilkProject/")

# load dataframe packages
using CSV, DataFrames, DataFramesMeta

# load GTF data of soy
GTF = CSV.read("/home/bensiv/Documents/Remilk_plants/Genomes/Arabidopsis_thaliana/ncbi_dataset/data/GCF_000001735.4/genomic.gtf", DataFrame, header = false, comment = "#")

# name the columns
colnames = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
rename!(GTF, names(GTF) .=> colnames)

# filter for transcript location only
GTF_transcripts = @linq GTF |> where(:feature .== "transcript")

# split attribute column to subcolumns
function SplitAttributes(gtf_table)
    transcript_ids = []
    descriptions = []
    for row in eachrow(gtf_table)
        attributes = string.(split(row.attribute, ";"))
        attributes_kv = split.(attributes, "\"")
        description_length = length(descriptions)
        transcript_id_length = length(transcript_ids)

        for (i,x) in enumerate(attributes_kv)
            attributes_kv[i] = replace.(x, " " => "")[findall(isempty.(x) .== false)]
            kv = replace.(x, " " => "")[findall(isempty.(x) .== false)]
            if "transcript_id" in kv
                push!(transcript_ids, kv[2])
            elseif "product" in kv
                if length(kv) > 1
                    push!(descriptions, kv[2])
                else
                    push!(descriptions, missing)
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

    gtf_table.transcript_id = transcript_ids
    gtf_table.description = descriptions

    return dropmissing!(gtf_table)
end

GTF_transcripts = SplitAttributes(GTF_transcripts)

Accessions = CSV.read("/home/bensiv/Dropbox/202207-RemilkProject/Localized_output/SeedProteins/Arabidopsis/accessions.csv", DataFrame, header = false)
rename!(Accessions, names(Accessions) .=> "transcript_id")
Accessions = innerjoin(Accessions, GTF_transcripts, on = :transcript_id)

select!(Accessions, :transcript_id => :accession, :description)

CSV.write("/home/bensiv/Dropbox/202207-RemilkProject/Localized_output/SeedProteins/Bnap/Accessions_description.csv", Accessions_used)


# if there is no feature transcript, need to find it by the start and stop codons

GTF = SplitAttributes(GTF)

GTF_group = @linq GTF |> groupby(:transcript_id)

GTF_transcripts_only = DataFrame() 
for group in GTF_group
    if all(in(group.feature).(["start_codon", "stop_codon"]))
        temp_df = first(group)
        temp_df.feature = "transcript"
        temp_df.start = group.start[findfirst(group.feature .== "start_codon")]
        temp_df.end = group.end[findfirst(group.feature .== "stop_codon")]
        push!(GTF_transcripts_only, temp_df)
    end
end

Accessions = innerjoin(Accessions, GTF_transcripts_only, on = :transcript_id)

select!(Accessions, :transcript_id => :accession, :description)

CSV.write("/home/bensiv/Dropbox/202207-RemilkProject/Localized_output/SeedProteins/Arabidopsis/Accessions_description.csv", Accessions)
CSV.write("/home/bensiv/Documents/Remilk_plants/Genomes/Arabidopsis_thaliana/ncbi_dataset/data/GCF_000001735.4/genomic_transcripts.csv", GTF_transcripts_only[:,1:end-1], header = false)
