using GeneticsMakie, CairoMakie, CSV, DataFrames, SnpArrays, Arrow
const GM = GeneticsMakie

CairoMakie.activate!(type = "png")
set_theme!(font = "Arial")

@info "Loading GENCODE annotation"
ispath(joinpath(@__DIR__, "../data")) || mkpath(joinpath(@__DIR__, "../data"))
if !isfile(joinpath(@__DIR__, "../data/gencode.v19.annotation.parsed.gtf.arrow"))
    @info "Downloading GENCODE"
    run(`curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh37_mapping/gencode.v39lift37.annotation.gtf.gz --output ./data/gencode.v39lift37.annotation.gtf.gz`)
    @time gencode = CSV.read(joinpath(@__DIR__, "../data/gencode.v39lift37.annotation.gtf.gz"), DataFrame,
        header = ["seqnames", "source", "feature", "start", "end", "score", "strand", "phase", "info"],
        delim = "\t", skipto = 6)
    GM.parsegtf!(gencode)
    select!(gencode, :seqnames, :feature, :start, :end, :strand, :gene_id, :gene_name, :gene_type, :transcript_id)
    @time Arrow.write(joinpath(@__DIR__, "../data/gencode.v39lift37.annotation.parsed.gtf.arrow"), gencode)
    run(`rm ./data/gencode.v39lift37.annotation.gtf.gz`)
end
@time gencode = Arrow.Table(joinpath(@__DIR__, "../data/gencode.v39lift37.annotation.parsed.gtf.arrow"))|> DataFrame
@assert (3_247_110, 9) == size(gencode)

@info "Loading 1000 Genomes reference panel"
@time kgp = SnpData(joinpath(@__DIR__, "../data/kgp"))

function subsetref(ref::SnpData, chr::AbstractString, range1::Real, range2::Real, path::AbstractString)
    SnpArrays.filter(ref, trues(size(ref)[1]), GM.findlocus(ref, chr, range1, range2); des = path)
    SnpData(path)
end

@info "Loading GWAS results"
phenotypes = ["scz", "bd", "asd", "adhd", "neuroticism", "alz", "menopause", "weight", "height"]
gwas = []
for p in phenotypes
    if !isfile(joinpath(@__DIR__, "../data/", p * ".gwas.arrow"))
        @time GM.downloadgwas(joinpath(@__DIR__, "../data/"), pheno = p)
        @time sumstat = CSV.read(joinpath(@__DIR__, "../data", GM.gwas[p].file), DataFrame, comment = "##", missingstring = ["NA"])
        @time GeneticsMakie.mungesumstats!(sumstat)
        @time Arrow.write(joinpath(@__DIR__, "../data", p * ".gwas.arrow"), sumstat)
        rm(joinpath(@__DIR__, "../data", GM.gwas[p].file))
    end
    push!(gwas, DataFrame(Arrow.Table(joinpath(@__DIR__, "../data", p * ".gwas.arrow"))))
end

function subsetgwas(gwas, chr::AbstractString, range1::Real, range2::Real)
    gwas_subset = Vector{DataFrame}(undef, length(gwas))
    for i in 1:length(gwas)
        gwas_subset[i] = gwas[i][GM.findlocus(gwas[i], chr, range1, range2), :]
    end
    gwas_subset
end

issig(P::AbstractVector; p = 5e-8) = any(P .< p)
issig(df::DataFrame; p = 5e-8) = issig(df.P; p = p)

@info "Working on visualization"
ispath(joinpath(@__DIR__, "../figs")) || mkpath(joinpath(@__DIR__, "../figs"))
function locuszoom(genes)
    for gene in genes
        @info "Working on $gene gene"
        window = 1e6
        chr, start, stop = GM.findgene(gene, gencode)
        range1, range2 = start - window, stop + window
        @info "Subsetting 1000 Genomes"
        @time kgp_subset = subsetref(kgp, chr, range1, range2, "./data/kgp.filtered")
        @info "Subsetting GWAS results."
        @time gwas_subset = subsetgwas(gwas, chr, range1, range2)
        titles = [GM.gwas[p].title for p in phenotypes]
        @info "Plotting phenome-wide LocusZoom."
        n = length(titles)
        f = Figure(resolution = (306, 1500))
        axs = [Axis(f[i, 1]) for i in 1:(n + 1)]
        for i in 1:n
            if issig(gwas_subset[i])
                GM.plotlocus!(axs[i], chr, range1, range2, gwas_subset[i]; ld = kgp_subset)
                if gwas_subset[i].BP[argmin(gwas_subset[i].P)] < (range1 + range2) / 2
                    Label(f[i, 1, Top()], "$(titles[i])", textsize = 6, halign = :right, padding = (0, 7.5, -5, 0))
                else
                    Label(f[i, 1, Top()], "$(titles[i])", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
                end    
            else
                GM.plotlocus!(axs[i], chr, range1, range2, gwas_subset[i])
                Label(f[i, 1, Top()], "$(titles[i])", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
            end
            rowsize!(f.layout, i, 30)
        end
        rs = GM.plotgenes!(axs[n + 1], chr, range1, range2, gencode; height = 0.1)
        rowsize!(f.layout, n + 1, rs)
        GM.labelgenome(f[n + 1, 1, Bottom()], chr, range1, range2)
        Colorbar(f[1:n, 2], limits = (0, 1), ticks = 0:1:1, height = 20,
            colormap = (:gray60, :red2), label = "LD", ticksize = 0, tickwidth = 0,
            tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,
            labelsize = 6, width = 5, spinewidth = 0.5)
        Label(f[1:n, 0], text = "-log[p]", textsize = 6, rotation = pi / 2)
        for i in 1:(n + 1)
            vlines!(axs[i], start, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs[i], stop, color = (:gold, 0.5), linewidth = 0.5)
        end
        for i in 1:n
            lines!(axs[i], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
        end
        rowgap!(f.layout, 5)
        colgap!(f.layout, 5)
        resize_to_layout!(f)
        save(joinpath(@__DIR__, "../figs/$(gene)-locuszoom.png"), f, px_per_unit = 4)
        for ext in ["bed", "bim", "fam"]
            rm(joinpath(@__DIR__, "../data/kgp.filtered." * ext))
        end
    end
end

genes = ["CHRNA5", "XRN2"]
locuszoom(genes)