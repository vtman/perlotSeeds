# perlotSeeds
Periodic lossless ternary seeds of maximum weight

<nav>
  <ul>
    <li><a href="#link_comp">Compilation</a></li>
    <li><a href="#link_intro">Introduction</a></li>
    <li><a href="#link_binaryLevel">periodicBinaryBlockLevel: binary blocks of less than maximum weight</a></li>
    <li><a href="#link_ternaryBlock">periodicTernaryBlocks: ternary blocks of maximum weight</a></li>	  
 <li><a href="#link_bin2text">convertBin2Text: convert ternary blocks in binary format into a text</a></li>
    <li><a href="#link_maxWeight">ternarySeedMaxWeight: ternary seeds of maximum weight</a></li>
    <li><a href="#link_data">Linked data</a></li>
    <li><a href="#link_fna2acgt">fna2acgt: convert a reference FNA file into a binary one</a></li>
    <li><a href="#link_acgt2lib">acgt2lib: creating an unsorted library of records (signature, position)</a></li>
    <li><a href="#link_sortLib">sortLib: sorting the library of records</a></li>
    <li><a href="#link_fastq2bin">fastq2bin: convert FASTQ files (reads) into a binary ones</a></li>
    <li><a href="#link_alignReads">alignReads: simplistic alignment of paired-end reads</a></li>
  </ul>
  </nav>

<h2 id="link_comp">Compilation</h2>

It is possible to download Intel C compiler (can be obtained from <a href="https://software.intel.com/content/www/us/en/develop/tools/oneapi/all-toolkits.html">oneAPI</a>) and use the following string for compilation

<tt>icpc codeName.cpp -mssse3 -std=c++17 -qopenmp -o codeName.exe</tt>

In some (serial application) <tt>-qopenmp</tt> option can be omitted. For Windows, it is possible to install <a href="https://visualstudio.microsoft.com/vs/">Visual Studio</a> (before installation of the oneAPI toolkits). Some libraries are slightly different for Windows and Linux, so for Windows uncomment

<tt>#define WIN32 1</tt>

<h2 id="link_intro">Introduction</h2>
We consider sequences of symbols <tt>A</tt>, <tt>C</tt>, <tt>G</tt> and <tt>T</tt>. In practical applications, we have a long <i>reference</i> sequence (up to billions of symbols) and a short sequence (<i>read</i>, often 50-300 symbols) obtained experimentally. The read is a small chunk of an unknown sequence, similar to the reference sequence. Therefore, we expect these two long sequences to have similar subsequences. There will also be some deviations, like SNPs (single-nucleotide polymorphisms), when some symbols are replaced by other symbols or insertions/deletions when some symbols are added or removed. 

The most straightforward approach to form the unknown long sequence is to pre-align its chunks (reads) to the known reference sequence and then perform a full-scale comparison to consider possible SNPs and insertions/deletions. For the pre-aligning step, we usually create a library of pairs (position within the reference sequence, the sequence of symbols) and sort all records by "sequence of symbols" value. So, if a read contains one of the subsequences from the library's records, we can find the corresponding positions within the reference sequence and thus pre-align the read. Of course, this approach only works if all elements of the chunks are the same.

For example, we have two sequences <tt>TTGGAGATCG</tt> and <tt>TAGGTGCTCG</tt> (of length 10). We compare the corresponding elements of the sequences and get 1 if there is a match and 0 if there is a mismatch.

<table>
  <tr><th>Sequence 1</th><th><tt>TTGGAGATCG</tt></th></tr>
  <tr><th>Sequence 2</th><th><tt>TAGGTGCTCG</tt></th></tr>
  <tr><th>Match</th><th><tt>1011010111</tt></th></tr>
</table>

We see that if we know the position of the first sequence <tt>TTGGAGATCG</tt> (for example, position 16 to the reference <tt>ACGACAACCTTGTCGTTGGAGATCGGAAGAGCACACGTCTGAAC</tt>), then we may pre-align another string containing <tt>TTGGAGATCG</tt>. However, it is not possible to pre-align a string containing <tt>TAGGTGCTCG</tt> since the library of records does not contain this chunk.

Reads often SNPs, so it is crucial to deal with possible mismatches. One of the approaches is to use <i>seeds</i>, i.e. a sequence of 0 and 1 elements. Let two sequences of symbols and a seed of the same length. When an element of the seed is 1, we compare the corresponding elements of two symbol sequences; otherwise, we ignore possible deviations. 

If we use seed <tt>1010101010</tt>, then we only need to compare symbols at odd positions, five in total. Three symbols match in the sequences. 
<table>
  <tr><th>Sequence 1</th><th><tt>TTGGAGATCG</tt></th></tr>
  <tr><th>Sequence 2</th><th><tt>TAGGTGCTCG</tt></th></tr>
  <tr><th>Seed</th>      <th><tt>1010101010</tt></th></tr>
  <tr><th>Match</th>     <th><tt>1_1_0_0_1_</tt></th></tr>
</table>

If we use seed <tt>1010010101</tt> (length is 10, weight is 5, weight is the number of 1-elements), then we have all five symbols match.
<table>
  <tr><th>Sequence 1</th><th><tt>TTGGAGATCG</tt></th></tr>
  <tr><th>Sequence 2</th><th><tt>TAGGTGCTCG</tt></th></tr>
  <tr><th>Seed</th>      <th><tt>1010010101</tt></th></tr>
  <tr><th>Match</th>     <th><tt>1_1__1_1_1</tt></th></tr>
</table>
So, if we form a library of records based on <tt>1010010101</tt>, then one record is (16, <tt>TGGTG</tt>) and a read containing the second sequence <tt>TAGGTGCTCG</tt> can be pre-aligned to <tt>ACGACAACCTTGTCGTTGGAGATCGGAAGAGCACACGTCTGAAC</tt>.

Ideally, we should find seeds of considerable weight since by increasing the weight by one, we reduce the number of candidate positions to be checked in 4 times (assuming that the chance to have any of the symbols <tt>A</tt>, <tt>C</tt>, <tt>G</tt> and <tt>T</tt> is the same). Codes to generate such seed can be found in <a href="https://github.com/vtman/PerFSeeB">https://github.com/vtman/PerFSeeB</a>.

Standard seeds are <i>binary</i> when there are only two states (<b><tt>1</tt></b> or <b><tt>#</tt></b> for <i>match</i> and <b><tt>0</tt></b> or <b><tt>&#95;</tt></b> for ``don't care symbol''). In genetics, the chance to have a <b>transition</b> mutation (<tt>A</tt> &harr; <tt>G</tt> or <tt>C</tt> &harr; <tt>T</tt>) is often twice higher than a <b>transversion</b> mutation (<tt>A</tt> &harr; <tt>C</tt>, <tt>A</tt> &harr; <tt>T</tt>, <tt>G</tt> &harr; <tt>C</tt>, <tt>G</tt> &harr; <tt>T</tt>). Transition-constrained seeds use ternary alphabet {<b><tt>#</tt></b>, <b><tt>@</tt></b>, <b><tt>&#95;</tt></b>} where <b><tt>@</tt></b> is for a match or a transition mismatch.

We construct ternary seeds as periodic seeds when there is a whole number of repetitions of a string followed by a remainder. 


<h2 id="link_binaryLevel">periodicBinaryBlockLevel: binary blocks of less than maximum weight</h2>

In the PerFSeeB project, we have generated periodic binary blocks of the maximum weight for a given number of mismatches. We may reuse results obtained in that project and generate binary blocks with less than the maximum weight.

<h3>Parameters</h3>

<ol>
  <li>Input folder (binary blocks obtained in PerFSeeB project, can be found in <a href="https://github.com/vtman/perlotSeeds/tree/main/binaryBlocksText">binaryBlocksText</a>)</li>
  <li>Output folder</li>
  <li>Number of mismatches (from 2 to 9)</li>
  <li>Size of blocks (from 3 to 50)</li>
  <li>Minimum level (0 when maximum weight)</li>
  <li>Maximum level</li>
</ol>

<tt>periodicBinaryBlockLevel.exe E:\Temp2\perlotSeeds\binaryBlocksText E:\Temp2\perlotSeeds\binaryBlocks 4 24 0 2</tt>

In most cases, it is enough to use level = 0 (92%); level = 2 is only required in a couple of cases, and the rest is for level = 1.

Output files can be downloaded from <a href="https://zenodo.org/record/10641412">Zenodo</a>, and example output files are <a href="https://github.com/vtman/perlotSeeds/tree/main/ExampleBinaryBlocks">ExampleBinaryBlocks</a>.

Output files are in binary format. For a given block size <b>B</b>, we find the smallest number <b>N</b> such that <b>B &#8804;8N</b>. So, for a block of length 30, we need 4 bytes to store. A hundred blocks require 100*4= 400 bytes.

<h2 id="link_ternaryBlock">periodicTernaryBlocks: ternary blocks of maximum weight</h2>


<h3>Parameters</h3>

<ol>
  <li>Input folder (files for binary blocks, like in <a href="https://github.com/vtman/perlotSeeds/tree/main/ExampleBinaryBlocks">ExampleBinaryBlocks</a>)</li>
  <li>Output folder</li>
  <li>Block size</li>
  <li>Number of mismatches (transition)</li>
  <li>Number of mismatches (transversion)</li>
<li>Level</li>
</ol>

<tt>periodicTernaryBlocks.exe E:\Temp2\perlotSeeds\binaryBlocks E:\Temp2\perlotSeeds\ternaryBlocks 30 2 3 0</tt>

Output files for ternary blocks are also in binary format. Each element of a block requires 2 bits (<tt>&#95;</tt> = <tt>0</tt> = <tt>00</tt>, <tt>#</tt> = <tt>1</tt> = <tt>01</tt>, <tt>@</tt> = <tt>2</tt> = <tt>10</tt>). So, a block of length 30 requires 8 bytes.

Output files can be downloaded from <a href="https://zenodo.org/record/8370909">Zenodo</a>, and example output files are <a href="https://github.com/vtman/perlotSeeds/tree/main/ExampleTernaryBlocks">ExampleTernaryBlocks</a>. Some files are huge (tens of gigabytes). 

<h2 id="link_bin2text">convertBin2Text: convert ternary blocks in binary format into a text</h2>

<h3>Parameters</h3>

<ol>
  <li>Input folder (ternary blocks in binary format)</li>
  <li>Output folder</li>
  <li>Block size</li>
  <li>Number of mismatches (transition)</li>
  <li>Number of mismatches (transversion)</li>
</ol>

<tt>convertBin2Text.exe E:\Temp2\perlotSeeds\ternaryBlocks E:\Temp2\perlotSeeds\ternaryBlocksText 30 2 3</tt>

<h2 id="link_maxWeight">ternarySeedMaxWeight: ternary seeds of maximum weight</h2>

<h3>Parameters</h3>

<ol>
  <li>Input folder (ternary blocks in binary format)</li>
  <li>Output folder</li>
  <li>Number of mismatches (transition)</li>
  <li>Number of mismatches (transversion)</li>
  <li>Block size (minimum)</li>
  <li>Block size (maximum)</li>
  <li>Read length (minimum)</li>
  <li>Read length (maximum)</li>
</ol>

<tt>ternarySeedMaxWeight.exe E:\Temp2\TSeeds\Zen\TernaryBlocks\T3V4 E:\Temp2\perlotSeeds\ternarySeeds 3 4 3 50 50 80</tt>

<h2 id="link_data">Linked data</h2>

Binary, ternary blocks and best ternary seeds for reads of length 30 to 500 can be found in Zenodo archives (<a href="https://zenodo.org/record/8395215">10.5281/zenodo.10641412</a> and <a href="https://zenodo.org/record/8395813">10.5281/zenodo.8395813</a>).

<hr>

<h2 id="link_fna2acgt">fna2acgt: convert a reference FNA file into a binary one</h2>

We used Human genome assembly <a href="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14">GRCh38.p14</a>

<h3>Parameters</h3>

<ol>
  <li>Input file (FNA file)</li>
  <li>Output folder (+ prefix)</li>
  <li>Gap size (> 0)</li>
</ol>

Two files are created: name.acgt (binary output file) and name.iacgt (information/text file).

*.acgt file constists of 16-byte blocks. Each block is four 32-bit componenets to store information about presence/absense of four nucleotides <tt>A</tt>, <tt>C</tt>, <tt>G</tt>, <tt>T</tt>.

<tt>*.iacgt</tt> file consists of 3 fields for each chromosome (each filed is on a new line): name of the chromosom as it is written in the FNA file, starting position in the *.acgt file (in bytes), number of nucleotides.

<tt>fna2acgt.exe C:\Data\chr4.fna D:\Genome\Library\chr4 200</tt>

There are two output files: <tt>D:\Genome\Library\chr4.acgt</tt> and <tt>D:\Genome\Library\chr4.iacgt</tt>


<h2 id="link_acgt2lib">acgt2lib: creating an unsorted library of records (signature, position)</h2>

The code is written for specific seeds. You should uncomment one of <tt>#define CASE_...</tt> (and comment out the other lines). You should use the same seed also for <b>alignReads</b>, so the same step for that code. You may also modify two lines

<tt>const int nLetters = 4;</tt>

This means there will be <tt>2<sup>8 nLetters</sup> = 2<sup>16</sup> = 65536</tt> output files created in <tt>original</tt> subfolder. Parameter <tt>nLetters = 4</tt> seems to be optimal for the Human genome (not so many small files). As each record in a file will have the same <tt>8 nLetters</tt> bits, we remove them to save space.

<tt>const int indexLevel = 14;</tt>

This parameter will be used when sorting operation is performed in <tt>sortLib</tt> (we create a list of starting positions of records).

<h3>Parameters</h3>

<ol>
  <li>Path to the ACGT (reference file, see <b>fna2acgt</b>)</li>
  <li>Output folder</li>
</ol>

<tt>acgt2lib.exe D:\Genome\Cram\ref38.acgt D:\Genome\library</tt>

The output folder MUST contain two subfolder <tt>original</tt> and <tt>sorted</tt>. The code will also produce <tt>info.txt</tt> file with short information about the parameters (this file will be used by <b>sortLib</b>).

<h2 id="link_sortLib">sortLib: sorting the library of records</h2>

For each binary file in <tt>original</tt> subfolder we sort th records and also create an index array of starting records (to simplify search in <b>alignReads</b>).

<h3>Parameters</h3>

<ol>
  <li>Input/output folder</li>
</ol>

<tt>sortLib.exe D:\Genome\library</tt>

The code use the file created by <tt>acgt2lib.exe</tt>, so both <tt>original</tt> and <tt>sorted</tt> subfolders should be present as well as <tt>info.txt</tt> file. An extra <tt>stat.txt</tt> will be created in the main folder to provide statistics how many signatures have a given number of same records. Note that for seeds of small weight you may need to increase the length of a storage array

<tt>const int NSIZE = 1000000;</tt>

<h2 id="link_fastq2bin">fastq2bin: convert FASTQ files (reads) into a binary ones</h2>

All tests were performed for <a href="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR016/ERR016118">ERR016118</a> data. The corresponding "exact" alignment can be found in <a href="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CHS/HG00513/alignment">HG00513</a> (Han Chinese South, 1000 genomes project). 

<h3>Parameters</h3>

<ol>
  <li>Input file (FASTQ format)</li>
  <li>Output file (FQB format)</li>
  <li>Length of reads</li>
</ol>

<tt>fastq2bin.exe D:\Genome\Cram\DRR346006_1.fastq D:\Genome\Cram\DRR346006_1.fqb 150</tt>

FQB files format: length of reads (32-bit integer), total number of reads (64-bit integer), reads in binary format. For example, if the length of reads is <tt>76</tt>, then we round it up to the nearest multiple of <tt>32</tt>, i.e. <tt>96</tt> (in case of reads' length <tt>150</tt>, we round it up to <tt>160</tt>). Then divide by <tt>32</tt>, get <tt>M</tt> (<tt>3</tt> or <tt>5</tt>, respectively). So, each read in binary format will require <tt>16*M</tt> bytes. Each 32 bits correspond to to presence/absense of <tt>A</tt>, <tt>C</tt>, <tt>G</tt>, <tt>T</tt>. 



<h2 id="link_alignReads">alignReads: simplistic alignment of paired-end reads</h2>

We load the reference library in memory. Then mutiple threads start to process chuncks of binary reads files (FQB files), the size of chunks is defined by <tt>MAX_READ_COUNT</tt> (set to <tt>10000</tt>). We read two sequences for each read, create their reversed counterparts and find <i>signatures</i> using the specified seeds. Then find the corresponding position in the reference sequence and merge them into four lists. Since we have paired-edn reads, we have to  check in the positions from the lists for two seqeucnes are within a specifeid range. For those successful pairs of positions we perform alignment and count the number of transitional and transersional mismataches. In the output folder we have files for each group of <tt>MAX_READ_COUNT</tt> reads with information related to the total number of signatures found in the merged lists and detailed information about sucessful alignments. There is also <tt>statInfo.txt</tt> where a shortened information is provided (number of candidate positions, numbe rof successful alignment and the best score).

<h3>Parameters</h3>

<ol>
  <li>Path to the first binary reads file (created by <b>fastq2bin</b>): <tt>D:\Genome\Data\ERR016118_1.fqb</tt></li>
  <li>Path to the second binary reads file: <tt>D:\Genome\Data\ERR016118_2.fqb</tt></li>
  <li>Path to the library (containing <tt>original</tt> and <tt>sorted</tt> folders): <tt>D:\Genome\DataB32</tt></li>
  <li>Reference ACGT file: <tt>D:\Genome\RefData\ref.acgt</tt></li>
  <li>Info reference file (IACGT): <tt>D:\Genome\RefData\ref.iacgt</tt></li>
  <li>Output folder: <tt>D:\Genome\output\C32</tt></li>
  <li>Distance between start position of read's seqeunces, minimum: <tt>200</tt></li>
  <li>Distance, maximum: <tt>800</tt></li>
</ol>

<tt>alignReads.exe ../reads/first.fqb ../reads/second.fqb ../data/library ../reference/ref.acgt ../reference/ref.iacgt ../output/C40 200 800</tt>

Examples of out files can be found in <a href="https://zenodo.org/records/10645042">https://zenodo.org/records/10645042</a>.
