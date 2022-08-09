# Long-read-sequencing-README

## Brief summary

Long-read Sequencing
To date, most large-scale genetic sequencing efforts for Alzheimer's disease and related dementias have been performed using short-read DNA sequencing. Although these approaches can identify single nucleotide changes and small indels, they are not optimized to identify large structural variations or repeat expansions. Furthermore, many areas of the genome cannot be accurately sequenced with this technology, like homologous elements, highly GC-rich regions, centromeric regions, and telomeres.

Long-read sequencing enables us to generate accurate genetic sequencing data for challenging genomic regions to identify structural variants driving Alzheimer's disease and related dementias pathology. A greater understanding of the genetic architecture of the Alzheimer's disease and related dementias genome will lead to further insight into the disease and pathway mechanisms underlying them and new potential therapeutic targets for these diseases.

With this research, we will build a public resource consisting of long-read genome sequencing data from a large number of confirmed people with Alzheimer's disease and related dementias and healthy individuals. We will make both the raw and processed data publicly available to the community, along with our analysis pipeline, algorithms, and optimized DNA isolation protocols.

## Guppy Analysis

**TL/DR:** MIG2 mode optimized on GCP with persistent use discount is
roughly **$47.25 ** for average RNA samples at **1.2 TB** and
**$130.73** for average DNA samples at **1.2 TB**.

## Introdunction

Here we present the benchamrking results for guppy basecaller on A100
GPUs. We benchmarked guppy on a single A100 GPU in multi-instance GPU
(MIG) mode as well as on 1, 2 and 4 A100 GPU’s in full mode. All our
estimates are with persistent use mode (30% reduced cost) on GCP.
Following table shows general specification of A100 (we have used the
only available (upto our knowledge) A100, 40 GB SMX, 108 streaming
multiprocessors (SMs) and a toal of 6912 cuda cores).


    ------------------------------------------------------------------------------------
        GPU      Memory   Cores   SMs   Cores_SM   GPCs   Cost_hr   vCPUs   RAM    Disk 
    ----------- -------- ------- ----- ---------- ------ --------- ------- ------ ------
       A100       40GB    6912    108      64       7      $2.79     12     85GB   10GB 

       A100       40GB    6912    108      64       7      $3.75     12     85GB   4TB  

     A100-MIG2    20GB    2688    42       64       -      $3.75     12     85GB   4TB  

     A100-MIG3    10GB    1792    28       64       -      $3.75     12     85GB   4TB  
    ------------------------------------------------------------------------------------

Table 1. A100 (in full and MIG MODE) configurations. Rows 1-2. A100 in
**full** mode with 10GB default and 4TB Persistent disk respectively.
Rows 3-4. A100 in MIG2 and MIG3 modes used for current benchmarking.

For more details on A100 GPU, see
<https://www.techpowerup.com/gpu-specs/a100-pcie.c3623>

<https://www.nvidia.com/content/dam/en-zz/Solutions/Data-Center/a100/pdf/nvidia-a100-datasheet-us-nvidia-1758950-r4-web.pdf>

> Please note that Guppy base caller crashes if run on data from gbucket
> (due to disk latency issues), all these tests were done by attaching a
> 4TB persistent SSD disk with a2-highgpu-1g VM (default VM for A100
> GPU). To run more than one samples simultaneously (say 3 with 3
> multi-instances, we need \~ 7TB of persistent disk as some samples are
> \> 3 TB)

### A100 in MIG mode

Single A100 GPU (40 GB) can be configured in 5 different GPU instances
(MIG) for best performance (see Table below). We have tested A100 in 7,
3 and 2 GPU instances both for RNA and DNA samples. **Results for MIG7
mode are slower than MIG3 and MIG2 modes and hence are not reported
here.**

![Table 2. A100 MIG Configurations](A100MIG.jpg)

Source
<https://docs.nvidia.com/datacenter/tesla/mig-user-guide/#:~:text=The%20new%20Multi%2DInstance%20GPU,resources%20for%20optimal%20GPU%20utilization.>

## Summary

A100 GPU in MIG mode is **almost 0.54** time cheaper than that of single
A100 in Full mode for RNA samples. In the case of DNA samples, MIG mode
is **almost 0.68** times cheaper than that of 4 A100 GPUs.

## A100 Performance

We benchmarked guppy base caller in 3 different MIG modes on a single
A100 GPU as well as single guppy base caller was run on 1, 2 and 4 A100
GPUs to identify best A100 utilization. See Figure 1 & 2 for details.

> Results for MIG7 mode are slower than MIG3 and MIG2 mode and hence are
> not reported here.

#### Results

Here we discuss our results for both RNA and DNA samples in MIG and full
mode.

> Please note that in MIG2 mode, an RNA and DNA sample was split into 2
> different folders (equal number of FASTA5 files) and 2 guppy
> basecaller’s were run on these sub-samples. Also note that basecalling
> rates (Figure 1A) in mig modes are shown as sum of rates for all
> instances and total time (Figure 1B) is the average of all instances.

All results in Figure 1 & 2 are based on a sample size of 1.2 TB for
both RNA and DNA samples.

> Please note that calculations for DNA in full A100 mode (A100-1) were
> not completed (due to long run time, we terminated the job).

Figure 1A shows guppy base calling rate for a single A100 in MIG3 and
MIG2 modes as well as A100 in full mode with guppy running on 1
(A100-1), 2 (A100-2) and 4 (A100-4) A100 GPUs for a single RNA/DNA
sample. As can be seen from Figure 1A, guppy basecalling rate in MIG2
mode is **\~2X** faster than in single A100 in full mode and gave
comparable performance on 2 A100 GPUs in full mode for both RNA and DNA
samples. Furthermore, total base calling completion time for both RNA
and DNA samples in MIG2 mode is compareable to that on 2 A100 GPUs in
full mode (Figure 1B).

![](A100Benchmarks2_files/figure-markdown_github/unnamed-chunk-5-1.png)

Figure 1. A100 GPU Perfomance in MIG and full mode.

Consequently, total cost for RNA and DNA replicates is **\~ 1/2** in
MIG2 mode than on 2 A100 GPUs in full mode (Figure 1C). Figure 1D gave
estimate for cost per TB for both DNA and RNA samples.

> We conclude that MIG2 mode gave best cost for both RNA and DNA
> samples.

> Major driver for higer cost in DNA replicates as compared to RNA
> replicates is due to larger number of samples (raw signals per TB) in
> DNA (3X as many as RNA samples for same size (in TB) of FASTA5 files)
> replicates.

> Note that in the case of DNA, MIG2 mode jobs finishes with 3 to 4 hr
> time difference causing little higher cost. A workarond for this would
> be submit two parallel jobs, one on each MIG machine to maximimize
> A100 utilization.

### Cost Prediction for RNA and DNA samples

Next, based on above estimates, we have calculated cost per sample for
each RNA and DNA sample in various A100 GPU configurations.

#### RNA Samples cost analysis

Figure 2 and Table 1 presents cost estimation for 10 RNA samples. As
seen in Table 1 and Figure 3, cost estimate for each RNA sample is **\~
1/2** in MIG2 mode than in the case of single A100 in full mode.

<table>
<caption>
Predicted Cost
</caption>
<thead>
<tr>
<th style="text-align:left;">
MODE
</th>
<th style="text-align:right;">
Cost
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;font-weight: bold;font-style: italic;color: white !important;background-color: black !important;">
MIG2
</td>
<td style="text-align:right;font-weight: bold;font-style: italic;color: white !important;background-color: black !important;">
707
</td>
</tr>
<tr>
<td style="text-align:left;">
MIG3
</td>
<td style="text-align:right;">
811
</td>
</tr>
<tr>
<td style="text-align:left;">
a100-1
</td>
<td style="text-align:right;">
1300
</td>
</tr>
</tbody>
</table>

Table 1. Total cost for 10 RNA samples

![](A100Benchmarks2_files/figure-markdown_github/unnamed-chunk-7-1.png)

Figure 2. RNA samples cost analysis. A) RNA samples sizes, B) Price
comparison between MIG2 and single A100 in full mode.

#### DNA Samples cost analysis

Figure 3 and Table 2 show cost analysis for 97 DNA samples. In the case
of DNA samples, MIG2 mode is the cheapest mode among MIG2, 2 A100 and 4
A100 GPUs with A100-2 being almost **\~2X** as expensive as MIG2 mode.

<table>
<caption>
Predicted Cost
</caption>
<thead>
<tr>
<th style="text-align:left;">
MODE
</th>
<th style="text-align:right;">
Cost
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;font-weight: bold;font-style: italic;color: white !important;background-color: black !important;">
MIG2
</td>
<td style="text-align:right;font-weight: bold;font-style: italic;color: white !important;background-color: black !important;">
19146
</td>
</tr>
<tr>
<td style="text-align:left;">
A100-2
</td>
<td style="text-align:right;">
36425
</td>
</tr>
<tr>
<td style="text-align:left;">
A100-4
</td>
<td style="text-align:right;">
28034
</td>
</tr>
</tbody>
</table>

Table 2. Total cost for 97 DNA samples

![](A100Benchmarks2_files/figure-markdown_github/unnamed-chunk-9-1.png)

Figure 3. DNA samples cost analysis. A) DNA sample sizes, B) Price
comparison between MIG2, 2 A100 and 4 A100 GPUs in full mode.


## Publications

XXX

## Graphical overview

XXX

## Overview of this code repo

1. XXXX
    1. Clinvar analyses INDI lines
    2. Genetic risk score analysis of INDI lines
    3. Genomic alteration analysis of INDI lines
2. XXX
   1. XXX
   2. XXX
3. XXX
   1. XXX
   2. XXX

## Contact details

For more information or questions please email:



## Full long-read sequencing Team:

The long-read sequencing team:

Adam Phillippy, NHGRI
Benedict Paten, UCSC
Bryan Traynor, NIA
Cornelis Blauwendraat, NIA
Fritz Sedlazeck, Baylor College of Medicine
Sonja Scholz, NINDS
and more





