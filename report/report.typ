// Project report. Numbers are read from the workflow outputs, so the text follows the data.
// Build with: python report/build_report.py

#let s = json("/results/summary.json")
#let qc = csv("/results/qc/sample_qc.tsv", delimiter: "\t").slice(1)
#let de = s.differential_expression
#let saha = de.contrasts.SAHA_vs_DMSO
#let rmd = de.contrasts.RMD_vs_DMSO
#let direct = de.contrasts.RMD_vs_SAHA
#let agree = s.drug_concordance
#let ercc = s.ercc
#let pub = s.published_comparison

// ---------------------------------------------------------------- helpers

#let commas(n) = {
  let digits = str(calc.round(n)).clusters()
  let out = ""
  for (i, ch) in digits.enumerate() {
    if i > 0 and calc.rem(digits.len() - i, 3) == 0 { out += "," }
    out += ch
  }
  out
}
#let num(x, digits: 2) = {
  // fixed number of decimals, with a true minus sign
  let text = str(calc.round(calc.abs(x), digits: digits))
  if digits > 0 {
    if not text.contains(".") { text += "." }
    text += "0" * (digits - text.split(".").at(1).len())
  }
  if x < 0 { "−" + text } else { text }
}
#let pct(x, digits: 0) = num(100 * x, digits: digits) + "%"
#let millions(n) = str(calc.round(n / 1e6, digits: 1)) + " million"
#let row(table, key, value) = table.filter(r => r.at(key) == value).first()

#let accent = rgb("#2a78d6")
#let ink2 = rgb("#52514e")
#let hair = rgb("#c3c2b7")
#let repo = "github.com/soumyadeepreddu04-web/hdaci-rnaseq"

// ---------------------------------------------------------------- page

#set document(
  title: "Vorinostat and romidepsin in primary human CD4+ T cells: an RNA-seq re-analysis from raw reads",
  author: "Soumyadeep",
)
#set page(
  paper: "a4",
  margin: (x: 2.1cm, top: 2.3cm, bottom: 2.2cm),
  header: context {
    if counter(page).get().first() > 1 {
      set text(size: 8pt, fill: ink2)
      [Vorinostat and romidepsin in primary CD4+ T cells: RNA-seq re-analysis #h(1fr) Soumyadeep]
    }
  },
  footer: context {
    set text(size: 8pt, fill: ink2)
    [#link("https://" + repo)[#repo] #h(1fr) #counter(page).display("1 / 1", both: true)]
  },
)
#set text(font: ("Helvetica Neue", "Arial"), size: 9.5pt, lang: "en", region: "GB")
#set par(justify: true, leading: 0.62em, spacing: 0.95em)
#set heading(numbering: "1.1")
#show heading.where(level: 1): it => block(above: 1.5em, below: 0.75em, text(size: 12.5pt, weight: "bold", it))
#show heading.where(level: 2): it => block(above: 1.2em, below: 0.6em, text(size: 10pt, weight: "bold", it))
#show link: set text(fill: accent)
#set figure(gap: 0.7em)
#show figure.where(kind: image): set figure(placement: auto)
#show raw: set text(font: ("Menlo", "DejaVu Sans Mono"))
#show figure.caption: it => {
  set text(size: 8.3pt, fill: ink2)
  set par(justify: false, leading: 0.55em)
  set align(left)
  it
}
#show figure: set block(above: 1.2em, below: 1.2em)
#set table(stroke: none, inset: (x: 5pt, y: 3.6pt))
#show table: set text(size: 8.5pt)
#show table: set par(justify: false)
#let rule(weight) = table.hline(stroke: weight + hair)
#let head(..cells) = table.header(..cells.pos().map(c => text(weight: "bold", c)))

#let step(title, tool) = rect(radius: 3pt, stroke: 0.6pt + hair, inset: (x: 4pt, y: 6pt), width: 100%)[
  #set par(justify: false, leading: 0.45em)
  #set align(center)
  #text(weight: "bold", size: 8pt, title) \
  #text(size: 6.8pt, fill: ink2, tool)
]
#let arrow = text(fill: ink2, size: 9pt, sym.arrow.r)

// ---------------------------------------------------------------- title

#text(size: 17pt, weight: "bold")[Vorinostat and romidepsin in primary human CD4+ T cells]
#v(-0.5em)
#text(size: 11.5pt, fill: ink2)[An RNA-seq re-analysis from raw reads, with two independent checks]
#v(0.2em)
#text(size: 9pt)[Soumyadeep · MSc Bioinformatics, Queen Mary University of London · November 2026 \
Code, results and this report: #link("https://" + repo)[#repo]]
#v(0.3em)
#line(length: 100%, stroke: 0.8pt + accent)

#let saha_n = saha.up + saha.down
#let rmd_n = rmd.up + rmd.down
#let gene(contrast, name) = contrast.highlight.filter(r => r.gene_name == name).first()
#let shift(name) = ercc.global_shift.filter(r => r.contrast == name).first()
#let spike(contrast, group) = ercc.by_subgroup.filter(r => r.contrast == contrast and r.subgroup == group).first().observed_median
#let strand(setting) = s.strandedness.filter(r => r.setting == setting).first().assigned_fraction
#let pub_saha = pub.contrasts.SAHA_vs_DMSO
#let pub_rmd = pub.contrasts.RMD_vs_DMSO
#let loo = s.leave_one_donor_out
#let flagged = s.qc.flagged_libraries
// The text below describes one flagged library by name; stop if a new run flags something else.
#assert(flagged == ("D7_SAHA",), message: "QC flags changed: review the text of sections 4.1 and 5")
#let odd = qc.filter(r => r.at(0) == "D7_SAHA").first()
#let rest = qc.filter(r => r.at(0) != "D7_SAHA").map(r => float(r.at(11)))
#let shift_both = (shift("SAHA_vs_DMSO").median_shift_log2 + shift("RMD_vs_DMSO").median_shift_log2) / 2
#let less_rna = pct(1 - calc.pow(2, shift_both))

#block(fill: rgb("#f1f5fb"), inset: (x: 12pt, y: 10pt), radius: 4pt, width: 100%)[
  #set par(justify: false)
  #text(weight: "bold", size: 10.5pt)[Summary]

  *Question.* Vorinostat and romidepsin are HDAC inhibitors tested as HIV latency-reversing drugs. What do they do to gene expression in uninfected primary human CD4+ T cells, and do they do the same thing?

  *Approach.* Twelve public RNA-seq libraries (four donors, each treated with DMSO, vorinostat and romidepsin) were re-processed from the raw reads with a Snakemake workflow: fastp, Subread, featureCounts, DESeq2 and fgsea.

  *Findings.*
  - Both drugs change more than half of the expressed genes: #commas(saha_n) of #commas(de.genes_tested) (#pct(saha_n / de.genes_tested)) for vorinostat and #commas(rmd_n) (#pct(rmd_n / de.genes_tested)) for romidepsin, at a false discovery rate of #pct(de.alpha).
  - They drive largely the same programme. Fold changes correlate at r = #num(agree.pearson_r), and #pct(agree.significant_both_same_direction / agree.significant_both) of the #commas(agree.significant_both) genes changed by both drugs move in the same direction. Romidepsin, at a 67-fold lower concentration, pushes the changes about #num(agree.major_axis_slope, digits: 1) times as far.
  - Immune gene sets go down (interferon responses, allograft rejection, inflammatory response). Cholesterol homeostasis, glycolysis and gene sets belonging to other tissues go up.

  *How far it can be trusted.*
  - Fold changes agree with those from the authors' own full-depth counts (r = #num(pub_saha.lfc_pearson_r) and #num(pub_rmd.lfc_pearson_r)), although this analysis used a quarter of the reads, a different aligner and a newer annotation.
  - Spike-in controls with known ratios are recovered accurately (slope #num(ercc.accuracy_slope)), and they show no sign of a large genome-wide shift that would break the normalisation.
  - One library fails a QC threshold. The conclusions hold when that donor is left out.
]

= Background and question

Antiretroviral therapy suppresses HIV but does not remove the latent provirus that persists in resting CD4+ T cells. One proposed route to a cure is to force the latent virus back into expression with a latency-reversing agent, so that infected cells can be cleared. Histone deacetylase (HDAC) inhibitors were among the first classes tested in people. Two of them are licensed cancer drugs: vorinostat (SAHA), which inhibits HDACs broadly, and romidepsin (RMD), which is more potent and more selective for class I HDACs. Both reactivate HIV in the laboratory, but neither has reduced the size of the reservoir in clinical trials [1].

HDACs regulate host genes as well as the provirus, so a drug given to reverse latency also changes the cell it acts in. Beliakova-Bethell and colleagues [1] measured this by RNA-seq in primary CD4+ T cells from four donors, treated for 24 hours with vorinostat (1 µM), romidepsin (15 nM) or the solvent DMSO. This project re-analyses the 12 uninfected libraries of that study from the raw reads, with a different toolchain, and asks three questions.

+ How many genes respond to each drug, and which pathways do they belong to?
+ Do the two drugs drive the same programme at different strengths, or different programmes?
+ How far can the answer be trusted? Two things are tested: whether a smaller, differently processed version of the data agrees with the authors' own counts, and whether the usual normalisation assumption holds for drugs that act on chromatin across the whole genome.

= Data and design

The data are public: GEO series GSE114883, ENA project PRJNA472973. Each library is total RNA sequenced as 51-base read pairs on an Illumina HiSeq 2500, with 77 to 89 million pairs per library. The study holds 24 libraries; this analysis uses the 12 from mock-infected cells, which isolate the effect of the drug on the host cell.

The design is paired. Each of the four donors (D1, D2, D4, D7) contributes one DMSO, one vorinostat and one romidepsin library, so every drug-treated sample has a control from the same person. The authors also added ERCC spike-in RNA to each lysate in proportion to the number of cells: Mix 1 to DMSO samples and Mix 2 to drug-treated samples. The two mixes contain the same 92 synthetic transcripts at known, different ratios, which is what makes the checks in section 5 possible.

= Methods

#figure(
  grid(
    columns: (1fr, 9pt, 1fr, 9pt, 1fr, 9pt, 1fr, 9pt, 1fr, 9pt, 1fr),
    align: horizon + center,
    step("Fetch reads", "ENA, first 20 M pairs"), arrow,
    step("Trim and QC", "fastp, MultiQC"), arrow,
    step("Align", "Subread, GRCh38 + ERCC"), arrow,
    step("Count", "featureCounts, GENCODE 50"), arrow,
    step("Test", "DESeq2, paired model"), arrow,
    step("Interpret", "fgsea, Hallmark sets"),
  ),
  caption: [The workflow. Snakemake [12] runs every step and rebuilds only what is out of date. Two checks branch off after counting: the spike-in analysis and the comparison with the authors' deposited counts.],
  placement: none,
)

*Reads.* Each library is streamed from ENA and the first 20 million read pairs are kept, about a quarter of each file. This keeps the whole run on a laptop. Where a read sits in the file reflects its position on the flow cell, not the gene it came from, so the subset is a fair sample; section 5 measures what the reduced depth costs. fastp [4] removes adapters and low-quality reads, and MultiQC [13] collects the per-library reports.

*Alignment and counting.* Read pairs are aligned to the GRCh38 primary assembly plus the 92 ERCC sequences with the Subread aligner [2], keeping only pairs with a single best location. featureCounts [3] then counts the pairs that overlap the exons of exactly one GENCODE release 50 gene on the expected strand. The genome was chosen over a transcript-only reference because these are total-RNA libraries: many pairs come from introns, and genome alignment accounts for them explicitly instead of forcing them onto spliced transcripts.

*Differential expression.* DESeq2 [5] fits a negative binomial model with the design `~ donor + treatment`, so each drug is compared with DMSO within donors and differences between people do not count as noise. A gene is tested if it has at least 10 read pairs in at least 4 libraries. Significance is called at a false discovery rate of #pct(de.alpha). Fold changes shown in plots and rankings are shrunk with `ashr` [6], which pulls the noisy estimates of low-count genes towards zero.

*Pathways.* Every tested gene is ranked by its Wald statistic and gene set enrichment is run with fgsea [7] on the 50 MSigDB Hallmark sets [8]. Ranking the whole list avoids choosing a cut-off, which matters when thousands of genes change.

*Checks.* The spike-in analysis and the comparison with the published counts are described with their results in section 5.

= Results

== Data quality

All 12 libraries pass the basic checks (Table 1). Between #pct(s.qc.q30_rate.min) and #pct(s.qc.q30_rate.max) of bases have a quality score of 30 or more, and trimming removes less than #calc.ceil(100 * (1 - s.qc.pairs_after_trimming_fraction.min))% of read pairs. A median of #pct(s.qc.uniquely_mapped_fraction.median) of pairs map to a single location and #pct(s.qc.assigned_fraction.median) are assigned to a gene, which leaves about #millions(s.qc.assigned_pairs.median) counted pairs per library. Most unassigned pairs map outside annotated exons, as expected for total RNA.

#figure(
  table(
    columns: (auto, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr),
    align: (left, right, right, right, right, right, right),
    rule(0.8pt),
    head([Library], [Pairs after trimming], [Bases ≥ Q30], [Uniquely mapped], [Assigned to a gene], [Spike-in share], [Genes detected]),
    rule(0.4pt),
    ..qc.map(r => (
      r.at(0).replace("_", " ") + if flagged.contains(r.at(0)) { " †" } else { "" },
      commas(float(r.at(3))),
      pct(float(r.at(2)), digits: 1),
      pct(float(r.at(5)), digits: 1),
      pct(float(r.at(8)), digits: 1),
      pct(float(r.at(9)), digits: 1),
      commas(float(r.at(10))),
    )).flatten(),
    rule(0.8pt),
  ),
  caption: [Per-library quality metrics. Each library starts from 20 million read pairs. "Assigned to a gene" is the share of all pairs counted by featureCounts; "spike-in share" is the share of counted pairs that are ERCC; a gene is "detected" with 10 or more pairs. † flagged: mapping or assignment rate more than five percentage points below the median.],
) <qc>

*Library type.* Counting one library under each strand setting assigns #pct(strand("reverse"), digits: 1) of its pairs as reverse-stranded, #pct(strand("forward"), digits: 1) as forward-stranded and #pct(strand("unstranded"), digits: 1) as unstranded. The libraries are therefore reverse-stranded, as the authors state, and the wrong setting would have discarded #pct(1 - strand("forward") / strand("reverse")) of the counts.

*One library stands out.* In the vorinostat library of donor D7 only #pct(float(odd.at(8))) of pairs are assigned to a gene, because #pct(float(odd.at(11))) map outside annotated exons, against #pct(calc.min(..rest)) to #pct(calc.max(..rest)) in the other eleven. It still gives #millions(s.qc.assigned_pairs.min) counted pairs and sits with the other vorinostat libraries in Figure 2, so it was kept. Section 5.3 shows that the conclusions do not depend on it.

== The drug is the main source of variation

#figure(
  image("/results/figures/pca.png", width: 88%),
  caption: [Principal component analysis of the 500 most variable genes (variance-stabilised counts). Left: all variation. Right: after subtracting each donor's own mean for every gene, which is the comparison the paired model makes. Labels are donors.],
)

The first principal component orders the libraries DMSO, vorinostat, romidepsin and carries most of the variance; the second separates the two drugs (Figure 2). Differences between donors are small beside the drug effect, and once each donor's own mean is removed the four donors lie almost on top of each other within each treatment.

== More than half of the expressed genes respond

Of the #commas(de.genes_in_annotation) annotated genes, #commas(de.genes_tested) have enough reads to test. Vorinostat changes #commas(saha_n) of them (#pct(saha_n / de.genes_tested)): #commas(saha.up) go up and #commas(saha.down) go down. Romidepsin changes #commas(rmd_n) (#pct(rmd_n / de.genes_tested)): #commas(rmd.up) up and #commas(rmd.down) down (Figure 3). Counting only changes of at least two-fold after shrinkage leaves #commas(saha.up_strong + saha.down_strong) genes for vorinostat and #commas(rmd.up_strong + rmd.down_strong) for romidepsin.

Two known responses to HDAC inhibition serve as positive controls. The most significant gene for both drugs is _H1-0_, which encodes the linker histone H1.0 and is induced by HDAC inhibitors in lymphocytes [14]; it rises #num(calc.pow(2, gene(saha, "H1-0").log2FoldChange), digits: 0)-fold with vorinostat and #num(calc.pow(2, gene(rmd, "H1-0").log2FoldChange), digits: 0)-fold with romidepsin. _CDKN1A_ (p21), the classic vorinostat-induced gene [15], rises #num(calc.pow(2, gene(saha, "CDKN1A").log2FoldChange), digits: 1)-fold and #num(calc.pow(2, gene(rmd, "CDKN1A").log2FoldChange), digits: 1)-fold. The two genes that the original authors followed up also change in the same direction with both drugs here: _SMARCB1_ falls (log2 fold change #num(gene(saha, "SMARCB1").log2FoldChange) with vorinostat and #num(gene(rmd, "SMARCB1").log2FoldChange) with romidepsin) and _PARP1_ rises (+#num(gene(saha, "PARP1").log2FoldChange) and +#num(gene(rmd, "PARP1").log2FoldChange)).

#figure(
  image("/results/figures/volcano.png", width: 94%),
  placement: top,
  caption: [Each point is a gene. The vertical axis is on a square-root scale and values above 150 are drawn at 150. Labels mark the largest changes among protein-coding genes with an adjusted p-value below $10^(-20)$.],
)

== The same programme, pushed further by romidepsin

Across all tested genes the two sets of fold changes correlate at r = #num(agree.pearson_r) (Figure 4, left). Of the #commas(agree.significant_saha) genes changed by vorinostat, #pct(agree.saha_genes_also_significant_with_rmd) are also changed by romidepsin, and of the #commas(agree.significant_both) genes changed by both, #commas(agree.significant_both_same_direction) (#pct(agree.significant_both_same_direction / agree.significant_both, digits: 1)) move in the same direction.

The drugs differ mainly in strength. The slope of the relationship is #num(agree.major_axis_slope): a gene that vorinostat doubles is raised about #num(calc.pow(2, agree.major_axis_slope), digits: 1)-fold by romidepsin, on average. The concentrations were 1 µM and 15 nM, so romidepsin has the larger effect at a 67-fold lower concentration, in line with its higher potency against class I HDACs [1]. The same gradient shows gene by gene: in most of the genes in Figure 4 (right), expression moves from DMSO through vorinostat to romidepsin in every donor.

#figure(
  grid(
    columns: (1.2fr, 1fr),
    gutter: 8pt,
    align: horizon,
    image("/results/figures/drug_concordance.png"),
    image("/results/figures/top_genes_heatmap.png"),
  ),
  placement: top,
  caption: [Left: fold changes for the two drugs, all #commas(agree.genes) tested genes. The slope is from a major-axis fit, which treats the two axes alike. Right: the 20 protein-coding genes most consistently raised and the 20 most consistently lowered by both drugs, ranked by the weaker of their two test statistics.],
)

Strength is not the whole story. The direct comparison of romidepsin with vorinostat finds #commas(direct.up_strong + direct.down_strong) genes that differ at least two-fold between the drugs, and the next section shows which gene sets they fall in.

== Pathways: immune programmes down, metabolic and out-of-place programmes up

#figure(
  image("/results/figures/hallmark_gsea.png", width: 84%),
  caption: [Gene set enrichment on the MSigDB Hallmark collection. Sets shown are significant (FDR below 0.05) for at least one drug. A positive score means that the set's genes are concentrated among genes that go up with the drug.],
)

Of the #s.gsea.sets_tested Hallmark sets large enough to test, #s.gsea.significant_saha are significant for vorinostat and #s.gsea.significant_rmd for romidepsin. Of these, #s.gsea.significant_both are significant for both drugs and #s.gsea.significant_both_same_direction of them point the same way (Figure 5).

- *Down with both drugs:* the interferon alpha and gamma responses, allograft rejection (a set rich in T-cell and antigen-presentation genes), complement and the inflammatory response.
- *Up with both drugs:* cholesterol homeostasis, glycolysis and hypoxia, together with sets that describe other tissues, such as myogenesis and adipogenesis. The second group is what would be expected if HDAC inhibition releases genes that are normally silent in T cells, although this analysis does not test that explanation.
- *Different between the drugs:* in the direct comparison, cell-cycle sets (E2F targets, G2M checkpoint), MYC targets and the unfolded protein response are higher with romidepsin, and the interferon responses are lower. The unfolded protein response is the one set the two drugs move in opposite directions.

For latency reversal the relevant point is that both drugs, at the concentrations used, suppress the cell's own interferon and inflammatory programmes and change thousands of genes that have nothing to do with HIV. This is the off-target burden the original authors describe [1], and here it grows with the potency of the drug.

= How far can the result be trusted?

== Agreement with the authors' own counts

The authors deposited a count table made with TopHat and HTSeq from every read, against GENCODE release 21. Fitting the same DESeq2 model to those counts gives fold changes that correlate with this workflow's at r = #num(pub_saha.lfc_pearson_r, digits: 3) for vorinostat and #num(pub_rmd.lfc_pearson_r, digits: 3) for romidepsin, across the #commas(pub_saha.genes_compared) genes tested in both (Figure 6). Of the genes called significant here, #pct(pub_saha.confirmed_fraction, digits: 1) and #pct(pub_rmd.confirmed_fraction, digits: 1) are also significant in the same direction in the authors' counts. Only #pub_saha.opposite_direction and #pub_rmd.opposite_direction genes are significant in both with opposite signs.

#figure(
  image("/results/figures/published_comparison.png", width: 88%),
  caption: [Fold changes from this workflow (20 million pairs, Subread, GENCODE 50) against fold changes from the authors' deposited counts (all reads, TopHat and HTSeq, GENCODE 21), same statistical model. The line marks equality.],
)

The full-depth counts give more significant genes (#commas(pub_saha.published_counts_up + pub_saha.published_counts_down) and #commas(pub_rmd.published_counts_up + pub_rmd.published_counts_down), against #commas(saha_n) and #commas(rmd_n) here), which is the cost of using a quarter of the reads. Expression levels agree less closely than fold changes do: the median correlation of log counts per million is #num(pub.per_sample_log_cpm_pearson.median) per library. A likely reason is that the two annotation releases define many genes differently, which shifts a gene's level in every library alike and therefore cancels in a fold change.

The paper itself reports fewer genes, #commas(pub.reported_by_authors.SAHA_vs_DMSO.up + pub.reported_by_authors.SAHA_vs_DMSO.down) for vorinostat and #commas(pub.reported_by_authors.RMD_vs_DMSO.up + pub.reported_by_authors.RMD_vs_DMSO.down) for romidepsin out of #commas(pub.reported_by_authors.genes_after_filtering), from a different statistical model (edgeR with RUVSeq factors). Both analyses find that romidepsin changes more genes than vorinostat.

== Spike-in controls

*Fold-change accuracy.* #ercc.spike_ins_detected of the 92 spike-ins have enough reads to analyse, #ercc.subgroup_b_used of them in subgroup B, which is present at the same concentration in both mixes. With subgroup B as the reference, the other three subgroups should differ between drug and DMSO libraries by log2 fold changes of −2, +0.58 and +1. The observed medians are #num(spike("SAHA_vs_DMSO", "A")), +#num(spike("SAHA_vs_DMSO", "C")) and +#num(spike("SAHA_vs_DMSO", "D")) for vorinostat and #num(spike("RMD_vs_DMSO", "A")), +#num(spike("RMD_vs_DMSO", "C")) and +#num(spike("RMD_vs_DMSO", "D")) for romidepsin, and a line through all the test spike-ins has a slope of #num(ercc.accuracy_slope) (Figure 7, left). The workflow measures fold changes of known size without compressing them.

#figure(
  image("/results/figures/ercc.png", width: 100%),
  caption: [Left: observed against expected fold change for each detected spike-in; bars are subgroup medians and the line marks equality. Right: the amount by which spike-in normalisation would move every host fold change, estimated separately in each donor; bars are medians.],
)

*The normalisation assumption.* DESeq2 scales libraries on the assumption that most genes do not change, which fails if a treatment raises or lowers transcription across the board [10]. Scaling by subgroup B instead moves every host fold change by a constant. The median of that shift across donors is #num(shift("SAHA_vs_DMSO").median_shift_log2) log2 units for vorinostat and #num(shift("RMD_vs_DMSO").median_shift_log2) for romidepsin, which corresponds to about #less_rna less host RNA per spike-in in treated cells than the default assumes (Figure 7, right). Three things follow.

+ The shift is small beside the effects reported. Among genes changed by both drugs the median absolute log2 fold change is #num(agree.median_abs_lfc_saha_in_shared) for vorinostat and #num(agree.median_abs_lfc_rmd_in_shared) for romidepsin, and because the shift is the same for every gene it does not change which genes respond most.
+ It cannot be attributed to the drugs. Treated libraries received a different spike-in mix from controls, so a #less_rna difference in how the two mixes were prepared or added would look identical. That both drugs show the same shift is what this explanation predicts.
+ Spike-in scaling would not be a better normalisation here. The estimates range from #num(shift("RMD_vs_DMSO").min_shift_log2) to +#num(shift("SAHA_vs_DMSO").max_shift_log2) between libraries, the extreme being the library flagged in QC, and refitting the model with spike-in size factors cuts the number of significant genes for vorinostat from #commas(shift("SAHA_vs_DMSO").significant_default) to #commas(shift("SAHA_vs_DMSO").significant_spike_in_normalised), because the noise of #ercc.subgroup_b_used spike-ins is passed on to every gene [11].

The spike-ins are therefore used as a check, not as the normalisation, and the check finds no evidence of the large global shift that would invalidate the default.

== Leaving out one donor at a time

#figure(
  table(
    columns: (auto, 1fr, 1fr, 1fr, 1fr),
    align: (left, right, right, right, right),
    rule(0.8pt),
    table.cell(rowspan: 2, align: bottom, text(weight: "bold")[Donor left out]),
    table.cell(colspan: 2, align: center, text(weight: "bold")[Vorinostat vs DMSO]),
    table.cell(colspan: 2, align: center, text(weight: "bold")[Romidepsin vs DMSO]),
    text(weight: "bold")[Fold-change r], text(weight: "bold")[Still significant],
    text(weight: "bold")[Fold-change r], text(weight: "bold")[Still significant],
    rule(0.4pt),
    ..("D1", "D2", "D4", "D7").map(d => {
      let a = loo.filter(r => r.donor_left_out == d and r.contrast == "SAHA_vs_DMSO").first()
      let b = loo.filter(r => r.donor_left_out == d and r.contrast == "RMD_vs_DMSO").first()
      (d, num(a.lfc_pearson_r, digits: 3), pct(a.still_significant_fraction), num(b.lfc_pearson_r, digits: 3), pct(b.still_significant_fraction))
    }).flatten(),
    rule(0.8pt),
  ),
  caption: [The model refitted with three donors. "Fold-change r" is the correlation with the four-donor fold changes; "still significant" is the share of four-donor significant genes that remain significant.],
)

Removing any one donor leaves the fold changes almost unchanged (r of #num(calc.min(..loo.map(r => r.lfc_pearson_r)), digits: 2) or more), and no gene that stays significant changes direction. Between #pct(calc.min(..loo.map(r => r.still_significant_fraction))) and #pct(calc.max(..loo.map(r => r.still_significant_fraction))) of the significant genes remain significant with three donors, a loss that reflects the smaller sample. Leaving out D7, the donor with the flagged library, loses the fewest significant genes, although it moves the vorinostat fold changes slightly more than leaving out any other donor.

= Limitations

- *Four donors, one dose, one time point.* The comparison of strength holds for 1 µM vorinostat against 15 nM romidepsin at 24 hours and says nothing about other doses or times.
- *Uninfected cells.* The analysis describes the host response and does not measure HIV reactivation. The study's 12 libraries from a latency model were not analysed.
- *Bulk RNA from mixed CD4+ T-cell subsets.* A change in expression cannot be told apart from a change in the mix of cell types, for example if a drug is more toxic to one subset.
- *A quarter of the reads, uniquely mapped pairs only.* Genes in repeated or multi-copy families are under-counted, and weakly expressed genes have less power than in the original study.
- *"Significant" is a low bar here.* When more than half of all genes pass, the count of genes says little; effect sizes and gene-set results carry the information.
- *No experimental validation.* This is a computational re-analysis of public data.

The natural next steps are to add the latency-model libraries and quantify HIV transcripts alongside host genes, and to quantify endogenous retroviruses, which the same group has shown respond to vorinostat in CD4+ T cells [16].

= Reproducing this analysis

Everything in this report, including the PDF itself, is produced by one `snakemake` command from a clean checkout. The workflow downloads the reads, the reference genome and the gene sets itself, and the numbers in this PDF are read from `results/summary.json`, which the workflow writes.

#block(fill: rgb("#f4f4f2"), inset: 8pt, radius: 3pt, width: 100%, breakable: false)[
  #set text(size: 7.6pt)
  ```
  git clone https://github.com/soumyadeepreddu04-web/hdaci-rnaseq.git && cd hdaci-rnaseq
  python3 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt
  Rscript workflow/envs/install_r_packages.R
  snakemake --cores 12 --resources ena_download=1
  ```
]

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, left, left),
    rule(0.8pt),
    head([Component], [Version], [Used for]),
    rule(0.4pt),
    [Snakemake], [#s.software.snakemake], [workflow management],
    [fastp (Rfastp)], [#s.software.Rfastp], [trimming and read QC],
    [Rsubread], [#s.software.Rsubread], [alignment (Subread) and counting (featureCounts)],
    [MultiQC], [#s.software.multiqc.replace("multiqc, version ", "")], [combined QC report],
    [DESeq2], [#s.software.DESeq2], [differential expression],
    [ashr], [#s.software.ashr], [fold-change shrinkage],
    [fgsea], [#s.software.fgsea], [gene set enrichment],
    [R / Bioconductor], [#s.software.R.replace("R version ", "").split(" ").first() / #s.software.Bioconductor], [runtime],
    [Reference], [GENCODE #s.gencode_release], [GRCh38 primary assembly and gene annotation, plus ERCC92],
    rule(0.8pt),
  ),
  caption: [Software and reference versions used for the results in this report.],
)

= References

#set par(justify: false, leading: 0.5em, spacing: 0.6em)
#set text(size: 8.3pt)
+ Beliakova-Bethell N, Mukim A, White CH, et al. Histone deacetylase inhibitors induce complex host responses that contribute to differential potencies of these compounds in HIV reactivation. _J Biol Chem_ 2019;294(14):5576–5589. doi:10.1074/jbc.RA118.005185. Data: GEO GSE114883.
+ Liao Y, Smyth GK, Shi W. The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. _Nucleic Acids Res_ 2019;47(8):e47. doi:10.1093/nar/gkz114
+ Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. _Bioinformatics_ 2014;30(7):923–930. doi:10.1093/bioinformatics/btt656
+ Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. _Bioinformatics_ 2018;34(17):i884–i890. doi:10.1093/bioinformatics/bty560
+ Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. _Genome Biol_ 2014;15:550. doi:10.1186/s13059-014-0550-8
+ Stephens M. False discovery rates: a new deal. _Biostatistics_ 2017;18(2):275–294. doi:10.1093/biostatistics/kxw041
+ Korotkevich G, Sukhov V, Budin N, et al. Fast gene set enrichment analysis. bioRxiv 060012. doi:10.1101/060012
+ Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database Hallmark gene set collection. _Cell Syst_ 2015;1(6):417–425. doi:10.1016/j.cels.2015.12.004
+ Jiang L, Schlesinger F, Davis CA, et al. Synthetic spike-in standards for RNA-seq experiments. _Genome Res_ 2011;21(9):1543–1551. doi:10.1101/gr.121095.111
+ Lovén J, Orlando DA, Sigova AA, et al. Revisiting global gene expression analysis. _Cell_ 2012;151(3):476–482. doi:10.1016/j.cell.2012.10.012
+ Risso D, Ngai J, Speed TP, Dudoit S. Normalization of RNA-seq data using factor analysis of control genes or samples. _Nat Biotechnol_ 2014;32(9):896–902. doi:10.1038/nbt.2931
+ Mölder F, Jablonski KP, Letcher B, et al. Sustainable data analysis with Snakemake. _F1000Research_ 2021;10:33. doi:10.12688/f1000research.29032.2
+ Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. _Bioinformatics_ 2016;32(19):3047–3048. doi:10.1093/bioinformatics/btw354
+ Sourlingas TG, Tsapali DS, Kaldis AD, Sekeri-Pataryas KE. Histone deacetylase inhibitors induce apoptosis in peripheral blood lymphocytes along with histone H4 acetylation and the expression of the linker histone variant, H1°. _Eur J Cell Biol_ 2001;80(11):726–732. doi:10.1078/0171-9335-00201
+ Richon VM, Sandhoff TW, Rifkind RA, Marks PA. Histone deacetylase inhibitor selectively induces p21WAF1 expression and gene-associated histone acetylation. _Proc Natl Acad Sci USA_ 2000;97(18):10014–10019. doi:10.1073/pnas.180316197
+ White CH, Beliakova-Bethell N, Lada SM, et al. Transcriptional modulation of human endogenous retroviruses in primary CD4+ T cells following vorinostat treatment. _Front Immunol_ 2018;9:603. doi:10.3389/fimmu.2018.00603
